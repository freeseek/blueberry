/* The MIT License

   Copyright (C) 2019 Giulio Genovese

   Author: Giulio Genovese <giulio.genovese@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "bcftools.h"
#include "rbuf.h"

#define VERSION     "2019-09-13"
#define SNP_DIST    "20"
#define SNP_STR     "SNP_DIST"
#define INDEL_DIST  "150"
#define INDEL_STR   "INDEL_DIST"

/******************************************
 * AUXILIARY BUFFER                       *
 ******************************************/

typedef struct
{
    int dist;
    int snp_dist;
    int indel_dist;
    const char *snp_str;
    const char *indel_str;
    bcf_hdr_t *hdr;
    bcf1_t **line;
    rbuf_t rbuf;

    int32_t *dst;
    int ndst;
}
auxbuf_t;

static auxbuf_t *auxbuf_init(bcf_hdr_t *hdr,
                             int snp_dist,
                             int indel_dist,
                             const char *snp_str,
                             const char *indel_str)
{
    auxbuf_t *buf = (auxbuf_t *)calloc(1,sizeof(auxbuf_t));
    buf->snp_dist = snp_dist;
    buf->indel_dist = indel_dist;
    buf->dist = indel_dist > snp_dist ? indel_dist : snp_dist;
    buf->snp_str = snp_str;
    buf->indel_str = indel_str;

    int snp_info_id = bcf_hdr_id2int(hdr, BCF_DT_ID, snp_str);
    if ( !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, snp_info_id) )
    {
        bcf_hdr_printf(hdr, "##INFO=<ID=%s,Number=1,Type=Integer,Description=\"Distance from closest SNP\">", snp_str);
    }

    int indel_info_id = bcf_hdr_id2int(hdr, BCF_DT_ID, indel_str);
    if ( !bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, indel_info_id) )
    {
        bcf_hdr_printf(hdr, "##INFO=<ID=%s,Number=1,Type=Integer,Description=\"Distance from closest indel\">", indel_str);
    }

    buf->hdr = hdr;
    rbuf_init(&buf->rbuf, 0);
    return buf;
}

static void auxbuf_destroy(auxbuf_t *buf)
{
    for (int i=0; i<buf->rbuf.m; i++)
    {
        if ( buf->line[i] ) bcf_destroy(buf->line[i]);
    }
    free(buf->line);
    free(buf->dst);
    free(buf);
}

// push a new record into the buffer
static void auxbuf_push(auxbuf_t *buf, bcf1_t *line)
{
    rbuf_expand0(&buf->rbuf, bcf1_t *, buf->rbuf.n+1, buf->line);
    int curr = rbuf_append(&buf->rbuf);

    // if a record is already present in the buffer, destroy it
    if ( buf->line[curr] ) bcf_destroy(buf->line[curr]);
    buf->line[curr] = bcf_dup(line);

    int i = curr;
    while ( rbuf_prev(&buf->rbuf, &i) )
    {
        int dist = buf->line[curr]->pos - buf->line[i]->pos;
        if ( dist <= buf->snp_dist )
        {
            if ( bcf_is_snp( buf->line[curr] ) )
            {
                if ( bcf_get_info_int32(buf->hdr, buf->line[i], buf->snp_str, &buf->dst, &buf->ndst) < 0 || dist < buf->dst[0] )
                {
                    bcf_update_info_int32(buf->hdr, buf->line[i], buf->snp_str, &dist, 1);
                }
            }
            if ( bcf_is_snp( buf->line[i] ) )
            {
                if ( bcf_get_info_int32(buf->hdr, buf->line[curr], buf->snp_str, &buf->dst, &buf->ndst) < 0 || dist < buf->dst[0] )
                {
                    bcf_update_info_int32(buf->hdr, buf->line[curr], buf->snp_str, &dist, 1);
                }
            }
        }
        if ( dist <= buf->indel_dist )
        {
            if ( bcf_get_variant_types( buf->line[curr] ) & VCF_INDEL )
            {
                if ( bcf_get_info_int32(buf->hdr, buf->line[i], buf->indel_str, &buf->dst, &buf->ndst) < 0 || dist < buf->dst[0] )
                {
                    bcf_update_info_int32(buf->hdr, buf->line[i], buf->indel_str, &dist, 1);
                }
            }
            if ( bcf_get_variant_types( buf->line[i] ) & VCF_INDEL )
            {
                if ( bcf_get_info_int32(buf->hdr, buf->line[curr], buf->indel_str, &buf->dst, &buf->ndst) < 0 || dist < buf->dst[0] )
                {
                    bcf_update_info_int32(buf->hdr, buf->line[curr], buf->indel_str, &dist, 1);
                }
            }
        }
    }
}

static bcf1_t *auxbuf_flush(auxbuf_t *buf, int flush_all)
{
    if ( buf->rbuf.n==0 ) return NULL;
    int first = rbuf_kth(&buf->rbuf, 0);
    int last = rbuf_last(&buf->rbuf);
    if ( !flush_all && buf->line[last]->pos - buf->line[first]->pos <= buf->dist ) return NULL;

    int i = rbuf_shift(&buf->rbuf);
    return buf->line[i];
}

/******************************************
 * PLUGIN                                 *
 ******************************************/

const char *about(void)
{
    return "Record variant distance from closest indel.\n";
}

static const char *usage_text(void)
{
    return
"\n"
"About: Record variant distance from closest indel. ("VERSION")\n"
"\n"
"Usage: bcftools +add-variant-dist [options] <in.vcf.gz>\n"
"\n"
"Plugin options:\n"
"    -p, --snp <int>                    maximum allowed distance to add distance to closest SNP ["SNP_DIST"]\n"
"    -P  --filter-snp <int>             ID of SNP INFO field ["SNP_STR"]\n"
"    -i, --indel <int>                  maximum allowed distance to add distance to closest indel ["INDEL_DIST"]\n"
"    -I  --filter-indel <int>           ID of indel INFO field ["INDEL_STR"]\n"
"        --no-version                   do not append version and command line to the header\n"
"    -o, --output <file>                write output to a file [standard output]\n"
"    -O, --output-type b|u|z|v          b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
"    -r, --regions <region>             restrict to comma-separated list of regions\n"
"    -R, --regions-file <file>          restrict to regions listed in a file\n"
"    -t, --targets [^]<region>          similar to -r but streams rather than index-jumps. Exclude regions with \"^\" prefix\n"
"    -T, --targets-file [^]<file>       similar to -R but streams rather than index-jumps. Exclude regions with \"^\" prefix\n"
"        --threads <int>                number of extra output compression threads [0]\n"
"    -s, --samples [^]<list>            comma separated list of samples to include (or exclude with \"^\" prefix)\n"
"    -S, --samples-file [^]<file>       file of samples to include (or exclude with \"^\" prefix)\n"
"        --force-samples                only warn about unknown subset samples\n"
"\n"
"Example:\n"
"    bcftools +add-variant-dist --snp 20 --indel 150 file.bcf\n"
"\n";
}

static void flush(auxbuf_t *buf,
                  htsFile *out_fh,
                  int flush_all)
{
    bcf1_t *line;
    while ( ( line = auxbuf_flush(buf, flush_all) ) )
    {
        if ( bcf_write(out_fh, buf->hdr, line) < 0 ) error("Unable to write to output VCF file\n");
    }
}

int run(int argc, char **argv)
{
    char *tmp;
    int snp_dist = (int)strtol(SNP_DIST, &tmp, 0);
    char *snp_filter_str = SNP_STR;
    int indel_dist = (int)strtol(INDEL_DIST, &tmp, 0);;
    char *indel_filter_str = INDEL_STR;
    char *output_fname = NULL;
    int output_type = FT_VCF;
    int n_threads = 0;
    int record_cmd_line = 1;
    char *targets_list = NULL;
    int targets_is_file = 0;
    char *regions_list = NULL;
    int regions_is_file = 0;
    char *sample_names = NULL;
    int sample_is_file = 0;
    int force_samples = 0;

    static struct option loptions[] =
    {
        {"snp", required_argument, NULL, 'p'},
        {"filter-snp", required_argument, NULL, 'P'},
        {"indel", required_argument, NULL, 'i'},
        {"filter-indel", required_argument, NULL, 'I'},
        {"samples", required_argument, NULL, 's'},
        {"samples-file", required_argument, NULL, 'S'},
        {"force-samples", no_argument, NULL, 1},
        {"output-type", required_argument, NULL, 'O'},
        {"output-file", required_argument, NULL, 'o'},
        {"threads", required_argument, NULL, 9},
        {"targets", required_argument, NULL, 't'},
        {"targets-file", required_argument, NULL, 'T'},
        {"regions", required_argument, NULL, 'r'},
        {"regions-file", required_argument, NULL, 'R'},
        {"no-version", no_argument, NULL, 8},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "p:P:i:I:s:S:t:T:r:R:h?o:O:89",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'p':
                snp_dist = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --snp %s\n", optarg);
                if ( snp_dist <= 0 ) error("Distance to extend calls needs to be positive: --snp %s\n", optarg);
                break;
            case 'P': snp_filter_str = optarg; break;
            case 'i':
                indel_dist = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --indel %s\n", optarg);
                if ( indel_dist <= 0 ) error("Distance to extend calls needs to be positive: --indel %s\n", optarg);
                break;
            case 'I': indel_filter_str = optarg; break;
            case 'o': output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': output_type = FT_BCF_GZ; break;
                          case 'u': output_type = FT_BCF; break;
                          case 'z': output_type = FT_VCF_GZ; break;
                          case 'v': output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case  9 : n_threads = strtol(optarg, 0, 0); break;
            case  8 : record_cmd_line = 0; break;

            case 't': targets_list = optarg; break;
            case 'T': targets_list = optarg; targets_is_file = 1; break;
            case 'r': regions_list = optarg; break;
            case 'R': regions_list = optarg; regions_is_file = 1; break;

            case 's': sample_names = optarg; break;
            case 'S': sample_names = optarg; sample_is_file = 1; break;
            case  1 : force_samples = 1; break;

            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }

    char *input_fname = NULL;
    if ( optind == argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) input_fname = "-"; // reading from stdin
        else { error("%s", usage_text()); }
    }
    else if ( optind+1 != argc ) error("%s", usage_text());
    else input_fname = argv[optind];

    bcf_srs_t *srs = bcf_sr_init();
    if ( regions_list )
    {
        if ( bcf_sr_set_regions(srs, regions_list, regions_is_file) < 0 )
            error("Failed to read the regions: %s\n", regions_list);
    }
    if ( targets_list )
    {
        if ( bcf_sr_set_targets(srs, targets_list, targets_is_file, 0) < 0 )
            error("Failed to read the targets: %s\n", targets_list);
    }
    if ( bcf_sr_set_threads(srs, n_threads) < 0 ) error("Failed to create threads\n");
    if ( !bcf_sr_add_reader(srs, input_fname) ) error("Failed to open %s: %s\n", input_fname, bcf_sr_strerror(srs->errnum));

    bcf_hdr_t *hdr = bcf_sr_get_header(srs, 0);
    if (record_cmd_line) bcf_hdr_append_version(hdr, argc, argv, "bcftools_+extendFMT");

    if (sample_names)
    {
        int ret = bcf_hdr_set_samples(hdr, sample_names, sample_is_file);
        if ( ret < 0 ) error("Error parsing the list of samples: %s\n", sample_names);
        else if ( force_samples && ret > 0 ) error("Sample name mismatch: sample #%d not found in the header\n", ret);
    }

    htsFile *out_fh  = hts_open(output_fname ? output_fname : "-", hts_bcf_wmode(output_type));
    if ( !out_fh ) error("Can't write to \"%s\": %s\n", output_fname, strerror(errno));
    if ( n_threads ) hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, srs->p);

    auxbuf_t *buf = auxbuf_init(hdr, snp_dist, indel_dist, snp_filter_str, indel_filter_str);

    if ( bcf_hdr_write(out_fh, hdr) < 0 ) error("Unable to write to output VCF file\n");

    int prev_rid = -1;
    while ( bcf_sr_next_line(srs) )
    {
        bcf1_t *line = bcf_sr_get_line(srs, 0);
        if (prev_rid != line->rid) flush(buf, out_fh, 1);
        auxbuf_push(buf, line);
        flush(buf, out_fh, 0);
        prev_rid = line->rid;
    }
    flush(buf, out_fh, 1);
    auxbuf_destroy(buf);

    hts_close(out_fh);
    bcf_sr_destroy(srs);

    return 0;
}
