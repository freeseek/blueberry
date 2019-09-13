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

#include <getopt.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/ksort.h>
#include "regidx.h"
#include "kmin.h"
#include "beta_binom.h"
#include "genome_rules.h"
#include "bcftools.h"
#include "filter.h"

#define VERSION                       "2019-09-13"
#define RHO_HOM_DFLT                  "0.001"
#define RHO_HET_DFLT                  "0.001"
#define MIN_DST_DFLT                  "100"
#define MATERNAL_PHASE_ERR_PRB_DFLT   "0.0"
#define PATERNAL_PHASE_ERR_PRB_DFLT   "0.0"
#define MIN_COV_DFLT                  "0"
#define MAX_COV_DFLT                  "0"
#define XY_PRB_DFLT                   "1e-08"
#define ERR_PRB_DFLT                  "1e-03"
#define CROSS_PRB_DFLT                "1e-05"
#define MATERNAL_SWITCH_ERR_PRB_DFLT  "1e-03"
#define PATERNAL_SWITCH_ERR_PRB_DFLT  "1e-03"
#define SEQ_ERR_PRB_DFLT              "0.003" // from "Lo, Y. M. D. et al. Maternal plasma DNA sequencing reveals the
                                              // genome-wide genetic and mutational profile of the fetus. Sci. Transl. Med. 2, 61ra91 (2010)"
#define TEL_PRB_DFLT                  "1e-05"
#define CEN_PRB_DFLT                  "1e-04"
#define SHORT_ARM_CHRS_DFLT           "13,14,15,21,22,chr13,chr14,chr15,chr21,chr22"

/****************************************
 * HMM MODEL DEFINITIONS                *
 ****************************************/

#define FLT_INCLUDE      (1<<0)
#define FLT_EXCLUDE      (1<<1)
#define NO_LOG           (1<<2)
#define FF_SNP_HOM       (1<<3)
#define FF_SNP_HET       (1<<4)
#define NO_INDELS        (1<<5)
#define NO_UNPHASED      (1<<6)
#define HETS_ONLY        (1<<7)
#define USE_SHORT_ARMS   (1<<8)
#define USE_CENTROMERES  (1<<9)

#define MOCHA_UNK 0
#define MOCHA_CNP_DEL 4
#define MOCHA_CNP_DUP 5
#define MOCHA_CNP_CNV 6

#define MOCHA_NAN 0
#define MOCHA_ARM 1
#define MOCHA_TEL 2

// HMM states
#define HMM_0  0
#define HMM_L  1
#define HMM_R  2
#define HMM_LL 3
#define HMM_LR 4
#define HMM_RR 5
#define N_HMM  6

static const int null_model[] = {HMM_L, HMM_R};
static const int full_model[] = {HMM_0, HMM_L, HMM_R, HMM_LL, HMM_LR, HMM_RR};
static const int monosomy_model[] = {HMM_0};
static const int trisomy_model[] = {HMM_LL, HMM_LR, HMM_RR};

// genotype states
#define HOM_REF      0
#define HET_REF_ALT  1
#define HET_ALT_REF  2
#define HOM_ALT      3
#define HET_UNPHASED 4
#define N_GT         5

// 0 - L - R - LL - LR - RR
static const int8_t gt2a[] = {0, 1, 0, 2, 1, 0};
static const int8_t gt2b[] = {0, 0, 1, 0, 1, 2};

// HMM transition matrix across states (without switch errors)
// HMM |   0 |   L |   R |  LL | LR |  RR
// ----+-----+-----+-----+-----+----+----
//   0 |   0 |   x |   x |  2x | 2x |  2x
//   L |   x |   0 |   c |   x |  x | c+x
//   R |   x |   c |   0 | c+x |  x |   x
//  LL |  2x |   x | c+x |   0 |  c |  2c
//  LR |  2x |   x |   x |   c |  0 |   c
//  RR |  2x | c+x |   x |  2c |  c |   0

// HMM transition matrix across states (without crossovers)
// HMM |   0 |   L |   R |  LL | LR |  RR
// ----+-----+-----+-----+-----+----+----
//   0 |   0 |   x |   x |  2x | 2x |  2x
//   L |   x |   0 |   f |   x |  x | f+x
//   R |   x |   f |   0 | f+x |  x |   x
//  LL |  2x |   x | f+x |   0 |  0 |   f
//  LR |  2x |   x |   x |   0 |  0 |   0
//  RR |  2x | f+x |   x |   f |  0 |   0

// copy number cost matrix
static const int hmm_x_trans[] = {0, 1, 1, 2, 2, 2,
                                  1, 0, 0, 1, 1, 1,
                                  1, 0, 0, 1, 1, 1,
                                  2, 1, 1, 0, 0, 0,
                                  2, 1, 1, 0, 0, 0,
                                  2, 1, 1, 0, 0, 0};

// crossover cost matrix
static const int hmm_c_trans[] = {0, 0, 0, 0, 0, 0,
                                  0, 0, 1, 0, 0, 1,
                                  0, 1, 0, 1, 0, 0,
                                  0, 0, 1, 0, 1, 2,
                                  0, 0, 0, 1, 0, 1,
                                  0, 1, 0, 2, 1, 0};

// phase error switch cost matrix
static const int hmm_f_trans[] = {0, 0, 0, 0, 0, 0,
                                  0, 0, 1, 0, 0, 1,
                                  0, 1, 0, 1, 0, 0,
                                  0, 0, 1, 0, 0, 1,
                                  0, 0, 0, 0, 0, 0,
                                  0, 1, 0, 1, 0, 0};

// array with the 15 beta binomial models required to run blueberry
static beta_binom_t *beta_binom_arr[15];

// look up table indicating what model to use for each genotype state combination
static int8_t idx_tbl[N_HMM * N_HMM * N_GT * N_GT * 4];

// look up table indicating initial state probabilities for the blueberry HMM model
static double hmm_init_tbl[N_HMM * N_HMM];

// look up table indicating transition probabilities for the blueberry HMM model
static double hmm_trans_tbl[N_HMM * N_HMM * N_HMM * N_HMM];

static const char *hmm_str[] = {" 0",
                                " L",
                                " R",
                                "LL",
                                "LR",
                                "RR"};

static const char *gt_str[] = {"0|0",
                               "0|1",
                               "1|0",
                               "1|1",
                               "0/1"};

static inline double log_mean_exp(double x, double y) { return x == y ? x : ( x > y ? x + log( 1 + expf(y-x) ) : y + log( 1 + expf(x-y) ) ) - M_LN2; }
// static inline double log_mean_exp2(double x, double y) { return log( ( exp(x) + exp(y) ) / 2.0 ); }

/*********************************
 * INITIALIZE LOOKUP TABLES      *
 *********************************/

// provided the mother is heterozygous, given the number of fetal reference and alternate alleles, returns the index of the model to use
static inline int8_t het_map(int8_t ref_cnt, int8_t alt_cnt)
{
    assert(ref_cnt>=0);
    assert(alt_cnt>=0);
    assert(ref_cnt+alt_cnt<=4);
    return ( ref_cnt == alt_cnt ) ? 0 : ( ref_cnt > alt_cnt ? ref_cnt + 3* alt_cnt : -( 3 * ref_cnt + alt_cnt ) );
}

// provided the mother is homozygous reference, given the number of fetal reference and alternate alleles, returns the index of the model to use
static inline int8_t hom_ref_map(int8_t ref_cnt, int8_t alt_cnt)
{
    assert(ref_cnt>=0);
    assert(alt_cnt>=0);
    assert(alt_cnt<=2);
    assert(ref_cnt+alt_cnt<=4);
    return alt_cnt == 0 ? 7 : 4 + ref_cnt + 4 * alt_cnt;
}

// hmm_mother is the HMM state for the fetal maternal allele
// hmm_father is the HMM state for the fetal paternal allele
// gt_mother is the genotype state of the mother
// gt_father is the genotype state of the father
static void beta_binom_idx_tbl_init(int8_t *idx_tbl)
{
    int8_t *idx = idx_tbl;
    for(int hmm_mother=0; hmm_mother<N_HMM; hmm_mother++)
    {
        for(int hmm_father=0; hmm_father<N_HMM; hmm_father++)
        {
            for(int gt_mother=0; gt_mother<N_GT; gt_mother++)
            {
                for(int gt_father=0; gt_father<N_GT; gt_father++)
                {
                    int8_t ref_cnt[4] = {0, 0, 0, 0};
                    int8_t alt_cnt[4] = {0, 0, 0, 0};

                    switch (gt_mother)
                    {
                        case HOM_REF:
                            ref_cnt[0] += gt2a[hmm_mother] + gt2b[hmm_mother];
                            ref_cnt[1] += gt2a[hmm_mother] + gt2b[hmm_mother];
                            break;
                        case HET_REF_ALT:
                            ref_cnt[0] += gt2a[hmm_mother];
                            alt_cnt[0] += gt2b[hmm_mother];
                            ref_cnt[1] += gt2a[hmm_mother];
                            alt_cnt[1] += gt2b[hmm_mother];
                            break;
                        case HET_ALT_REF:
                            ref_cnt[0] += gt2b[hmm_mother];
                            alt_cnt[0] += gt2a[hmm_mother];
                            ref_cnt[1] += gt2b[hmm_mother];
                            alt_cnt[1] += gt2a[hmm_mother];
                            break;
                        case HOM_ALT:
                            alt_cnt[0] += gt2a[hmm_mother] + gt2b[hmm_mother];
                            alt_cnt[1] += gt2a[hmm_mother] + gt2b[hmm_mother];
                            break;
                        case HET_UNPHASED:
                            ref_cnt[0] += gt2a[hmm_mother];
                            alt_cnt[0] += gt2b[hmm_mother];
                            ref_cnt[1] += gt2a[hmm_mother];
                            alt_cnt[1] += gt2b[hmm_mother];
                            ref_cnt[2] += gt2b[hmm_mother];
                            alt_cnt[2] += gt2a[hmm_mother];
                            ref_cnt[3] += gt2b[hmm_mother];
                            alt_cnt[3] += gt2a[hmm_mother];
                            break;
                    }

                    switch (gt_father)
                    {
                        case HOM_REF:
                            ref_cnt[0] += gt2a[hmm_father] + gt2b[hmm_father];
                            ref_cnt[2] += gt2a[hmm_father] + gt2b[hmm_father];
                            break;
                        case HET_REF_ALT:
                            ref_cnt[0] += gt2a[hmm_father];
                            alt_cnt[0] += gt2b[hmm_father];
                            ref_cnt[2] += gt2a[hmm_father];
                            alt_cnt[2] += gt2b[hmm_father];
                            break;
                        case HET_ALT_REF:
                            ref_cnt[0] += gt2b[hmm_father];
                            alt_cnt[0] += gt2a[hmm_father];
                            ref_cnt[2] += gt2b[hmm_father];
                            alt_cnt[2] += gt2a[hmm_father];
                            break;
                        case HOM_ALT:
                            alt_cnt[0] += gt2a[hmm_father] + gt2b[hmm_father];
                            alt_cnt[2] += gt2a[hmm_father] + gt2b[hmm_father];
                            break;
                        case HET_UNPHASED:
                            ref_cnt[0] += gt2a[hmm_father];
                            alt_cnt[0] += gt2b[hmm_father];
                            ref_cnt[1] += gt2b[hmm_father];
                            alt_cnt[1] += gt2a[hmm_father];
                            ref_cnt[2] += gt2a[hmm_father];
                            alt_cnt[2] += gt2b[hmm_father];
                            ref_cnt[3] += gt2b[hmm_father];
                            alt_cnt[3] += gt2a[hmm_father];
                            break;
                    }

                    for (int i=0; i<4; i++)
                    {
                        switch (gt_mother)
                        {
                            case HOM_REF:
                                idx[i] = hom_ref_map(ref_cnt[i], alt_cnt[i]);
                                break;
                            case HOM_ALT:
                                idx[i] = -hom_ref_map(alt_cnt[i], ref_cnt[i]);
                                break;
                            case HET_REF_ALT:
                            case HET_ALT_REF:
                            case HET_UNPHASED:
                                idx[i] = het_map(ref_cnt[i], alt_cnt[i]);
                                break;
                        }
                    }

                    if ( gt_mother != HET_UNPHASED || hmm_mother == HMM_0 || hmm_mother == HMM_LR )
                    {
                        idx[2] = bcf_int8_missing;
                        idx[3] = bcf_int8_missing;
                    }

                    if ( gt_father != HET_UNPHASED || hmm_father == HMM_0 || hmm_father == HMM_LR )
                    {
                        idx[1] = idx[2];
                        idx[2] = bcf_int8_missing;
                        idx[3] = bcf_int8_missing;
                    }

                    idx += 4;
                }
            }
        }
    }
}

// initialize lookup table for HMM initial probabilities
static void init_hmm_init_tbl(double *hmm_init_tbl)
{
    for (int hmm_mother=0; hmm_mother<N_HMM; hmm_mother++)
    {
        for (int hmm_father=0; hmm_father<N_HMM; hmm_father++)
        {
            int8_t ploidy_mother = gt2a[hmm_mother] + gt2b[hmm_mother];
            int8_t ploidy_father = gt2a[hmm_father] + gt2b[hmm_father];
            if ( ploidy_mother == 0 && ploidy_father == 0 )
            {
                hmm_init_tbl[ N_HMM * hmm_mother + hmm_father ] = 2;
            }
            if ( ploidy_mother == 0 || ploidy_father == 0 )
            {
                hmm_init_tbl[ N_HMM * hmm_mother + hmm_father ] = 1;
            }
            else
            {
                int8_t ploidy = ploidy_mother + ploidy_father;
                hmm_init_tbl[ N_HMM * hmm_mother + hmm_father ] = ploidy > 2 ? ploidy - 2 : 2 - ploidy;
            }
        }
    }
}

// initialize lookup table for HMM transition probabilities
static void init_hmm_trans_tbl(double *hmm_trans_tbl,
                               double xy_log_prb,
                               double cross_log_prb,
                               double maternal_switch_err_log_prb,
                               double paternal_switch_err_log_prb)
{
    double maternal_tmp_tbl[N_HMM * N_HMM] = {0.0};
    double paternal_tmp_tbl[N_HMM * N_HMM] = {0.0};
    for (int hmm_prev=0; hmm_prev<N_HMM; hmm_prev++)
    {
        for (int hmm_next=0; hmm_next<N_HMM; hmm_next++)
        {
            double xy_trans = xy_log_prb * hmm_x_trans[ N_HMM * hmm_prev + hmm_next ];
            double cross_trans = cross_log_prb * hmm_c_trans[ N_HMM * hmm_prev + hmm_next ];
            double maternal_switch_trans = maternal_switch_err_log_prb * hmm_f_trans[ N_HMM * hmm_prev + hmm_next ];
            double paternal_switch_trans = paternal_switch_err_log_prb * hmm_f_trans[ N_HMM * hmm_prev + hmm_next ];
            // apply the transition cost that is smaller between switch errors and crossover
            maternal_tmp_tbl[ N_HMM * hmm_prev + hmm_next ] = xy_trans + ( ( 0 > cross_trans && cross_trans > maternal_switch_trans ) || maternal_switch_trans == 0 ? cross_trans : maternal_switch_trans );
            paternal_tmp_tbl[ N_HMM * hmm_prev + hmm_next ] = xy_trans + ( ( 0 > cross_trans && cross_trans > paternal_switch_trans ) || paternal_switch_trans == 0 ? cross_trans : paternal_switch_trans );
        }
    }
    for (int hmm_mother_prev = 0; hmm_mother_prev < N_HMM; hmm_mother_prev++)
    {
        for (int hmm_father_prev = 0; hmm_father_prev < N_HMM; hmm_father_prev++)
        {
            for (int hmm_mother_next = 0; hmm_mother_next < N_HMM; hmm_mother_next++)
            {
                for (int hmm_father_next = 0; hmm_father_next < N_HMM; hmm_father_next++)
                {
                    int state_prev = N_HMM * hmm_mother_prev + hmm_father_prev;
                    int state_next = N_HMM * hmm_mother_next + hmm_father_next;
                    hmm_trans_tbl[ N_HMM * N_HMM * state_prev + state_next ] = maternal_tmp_tbl[ N_HMM * hmm_mother_prev + hmm_mother_next ] + paternal_tmp_tbl[ N_HMM * hmm_father_prev + hmm_father_next ];
                }
            }
        }
    }
}

/*********************************
 * MODEL STRUCTURES              *
 *********************************/

typedef struct
{
    double err_log_prb;
    double xy_log_prb;
    double cross_log_prb;
    double maternal_switch_err_log_prb;
    double paternal_switch_err_log_prb;
    double seq_err_log_prb;
    double tel_log_prb;
    double cen_log_prb;
    int min_dst;
    double maternal_phase_err_prb;
    double paternal_phase_err_prb;
    int min_cov;
    int max_cov;
    int flags;
    genome_rules_t *genome_rules;

    int cellfree_id;
    int mother_id;
    int father_id;

    int rid;
    int n;
    int *pos_arr;
    int m_pos;
} model_t;

typedef struct
{
    int rid;
    int type_snps[5];
    int maternal_switch_count;
    int paternal_switch_count;
    double mean_cov;
    double ff_hom;
    double ff_het;
    double rho_hom;
    double rho_het;
    double mother_monosomy_lod;
    double father_monosomy_lod;
    double mother_trisomy_lod;
    double father_trisomy_lod;
} stats_t;

typedef struct
{
    double x_nonpar_cov_mean;
    double y_nonpar_cov_mean;
    double mt_cov_mean;
    stats_t stats;
    stats_t *stats_arr;
    int m_stats, n_stats;

    double ff;
    double rho_hom;
    double rho_het;

    int n;
    int *vcf_imap_arr;
    int m_vcf_imap;
    int16_t *ad_ref_arr;
    int16_t *ad_alt_arr;
    int m_ad_ref;
    int m_ad_alt;
    int8_t *gt_mother_arr;
    int8_t *gt_father_arr;
    int m_gt_mother;
    int m_gt_father;
} trio_t;

typedef struct
{
    int rid;
    int beg_pos;
    int end_pos;
    int length;
    int8_t p_arm;
    int8_t q_arm;
    int nsites;
    double lod;
    int8_t hmm_mother;
    int8_t hmm_father;
} blueberry_t;

typedef struct
{
    int n, m;
    blueberry_t *a;
} blueberry_table_t;

/*********************************
 * DEBUGGING                     *
 *********************************/

// print likelihood matrix for debugging purposes
void debug_log_lkl(FILE *restrict stream,
                   double *emis_log_lkl,
                   const int16_t *ad_ref,
                   const int16_t *ad_alt,
                   const int8_t *gt_mother,
                   const int8_t *gt_father,
                   const int *imap_arr,
                   int n,
                   const int *pos,
                   const int *vcf_imap_arr,
                   const int *mother_arr,
                   int n_mother,
                   const int *father_arr,
                   int n_father)
{
    for (int t=0; t<n; t++)
    {
        if (imap_arr) fprintf( stream, "%d\t%s\t%s\t%d\t%d", pos[ vcf_imap_arr[ imap_arr[t] ] ], gt_str[gt_mother[imap_arr[t]]], gt_str[gt_father[imap_arr[t]]], ad_ref[imap_arr[t]], ad_alt[imap_arr[t]] );
        else fprintf( stream, "%d\t%s\t%s\t%d\t%d", pos[ vcf_imap_arr[ imap_arr[t] ] ], gt_str[gt_mother[t]], gt_str[gt_father[t]], ad_ref[t], ad_alt[t] );
        for (int hmm_mother=0; hmm_mother<n_mother; hmm_mother++)
        {
            for (int hmm_father=0; hmm_father<n_father; hmm_father++)
            {
                fprintf( stream, "\tf(%s,%s)=%.4f", hmm_str[mother_arr[hmm_mother]], hmm_str[father_arr[hmm_father]],
                emis_log_lkl[ N_HMM * N_HMM * t + N_HMM * mother_arr[hmm_mother] + father_arr[hmm_father] ] * M_LOG10E );
            }
        }
        fprintf( stream, "\n" );
    }
}


// print Viterbi path for debugging purposes
void debug_path(FILE *restrict stream,
                int rid,
                const int8_t *path,
                const double *emis_log_lkl,
                const int16_t *ad_ref,
                const int16_t *ad_alt,
                const int8_t *gt_mother,
                const int8_t *gt_father,
                const int *imap_arr,
                int T,
                const int *pos,
                const int *vcf_imap_arr)
{
    for (int t=0; t<T; t++)
    {
        if ( imap_arr )
        {
            fprintf( stream, "%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n", rid+1, pos[ vcf_imap_arr[ imap_arr[t] ] ],
                gt_str[gt_mother[imap_arr[t]]], gt_str[gt_father[imap_arr[t]]], ad_ref[imap_arr[t]], ad_alt[imap_arr[t]],
                hmm_str[path[t] / N_HMM], hmm_str[path[t] % N_HMM]);
        }
        else
        {
            fprintf( stream, "%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n", rid+1, pos[ vcf_imap_arr[ t ] ],
                gt_str[gt_mother[t]], gt_str[gt_father[t]], ad_ref[t], ad_alt[t],
                hmm_str[path[t] / N_HMM], hmm_str[path[t] % N_HMM]);
        }
    }
}

void debug_idx_lookup_tbl(FILE *restrict stream,
                          const int8_t *idx_tbl,
                          const beta_binom_t * const *beta_binom_arr,
                          double ff)
{
    const int8_t *idx = idx_tbl;
    fprintf(stream, "Fetal fraction in use is: %.2f%%\n", ff);
    for(int hmm_mother=0; hmm_mother<6; hmm_mother++)
    {
        for(int hmm_father=0; hmm_father<6; hmm_father++)
        {
            for(int gt_mother=0; gt_mother<5; gt_mother++)
            {
                for(int gt_father=0; gt_father<5; gt_father++)
                {
                    fprintf(stream, "%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",
                    hmm_str[hmm_mother], hmm_str[hmm_father], gt_str[gt_mother], gt_str[gt_father],
                    idx[0], idx[1], idx[2], idx[3],
                    idx[0] == bcf_int8_missing ? NAN : ( idx[0]>=0 ? beta_binom_arr[idx[0]]->p : 1.0f - beta_binom_arr[-idx[0]]->p ),
                    idx[1] == bcf_int8_missing ? NAN : ( idx[1]>=0 ? beta_binom_arr[idx[1]]->p : 1.0f - beta_binom_arr[-idx[1]]->p ),
                    idx[2] == bcf_int8_missing ? NAN : ( idx[2]>=0 ? beta_binom_arr[idx[2]]->p : 1.0f - beta_binom_arr[-idx[2]]->p ),
                    idx[3] == bcf_int8_missing ? NAN : ( idx[3]>=0 ? beta_binom_arr[idx[3]]->p : 1.0f - beta_binom_arr[-idx[3]]->p ) );
                    idx += 4;
                }
            }
        }
    }
}

void debug_trio(FILE *restrict stream,
                trio_t *trio,
                model_t *model)
{
    for (int i=0; i<trio->n; i++)
    {
        fprintf(stream, "%d %d %d %s %s\n", model->pos_arr[trio->vcf_imap_arr[i]], trio->ad_ref_arr[i], trio->ad_alt_arr[i],
        gt_str[trio->gt_mother_arr[i]], gt_str[trio->gt_father_arr[i]]);
    }
}

/*********************************
 * MODEL AUXILIARY               *
 *********************************/

// initialize the array with the 15 beta binomial models required to run the blueberry model
static void init_log_lkl_mdls(beta_binom_t **beta_binom_arr, double ff, double rho_hom, double rho_het, int n1, int n2, double seq_err_rate)
{
    int ref_cnt, alt_cnt;
    int8_t idx;
    double p;

    alt_cnt = 0;
    for (ref_cnt=0; ref_cnt<=4; ref_cnt++)
    {
        idx = het_map(ref_cnt, alt_cnt);
        p = ( ( 1.0f - ff ) + (double)ref_cnt * ff ) / ( 2.0f * ( 1.0f - ff ) + (double)(ref_cnt + alt_cnt) * ff );
        beta_binom_update(beta_binom_arr[idx], p, rho_het, n1, n2);
    }

    alt_cnt = 1;
    for (ref_cnt=2; ref_cnt<=3; ref_cnt++)
    {
        idx = het_map(ref_cnt, alt_cnt);
        p = ( ( 1.0f - ff ) + (double)ref_cnt * ff ) / ( 2.0f * ( 1.0f - ff ) + (double)(ref_cnt + alt_cnt) * ff );
        beta_binom_update(beta_binom_arr[idx], p, rho_het, n1, n2);
    }

    ref_cnt = 0;
    alt_cnt = 0;
    p = 1.0f - seq_err_rate;
    idx = hom_ref_map(ref_cnt, alt_cnt);
    beta_binom_update(beta_binom_arr[idx], p, rho_hom, n1, n2);

    alt_cnt = 1;
    for (ref_cnt=0; ref_cnt<=3; ref_cnt++)
    {
        idx = hom_ref_map(ref_cnt, alt_cnt);
        p = ( 2.0f * ( 1.0f - ff ) + (double)ref_cnt * ff ) / ( 2.0f * ( 1.0f - ff ) + (double)(ref_cnt + alt_cnt) * ff );
        beta_binom_update(beta_binom_arr[idx], p, rho_hom, n1, n2);
    }

    alt_cnt = 2;
    for (ref_cnt=0; ref_cnt<=2; ref_cnt++)
    {
        idx = hom_ref_map(ref_cnt, alt_cnt);
        p = ( 2.0f * ( 1.0f - ff ) + (double)ref_cnt * ff ) / ( 2.0f * ( 1.0f - ff ) + (double)(ref_cnt + alt_cnt) * ff );
        beta_binom_update(beta_binom_arr[idx], p, rho_hom, n1, n2);
    }
}

static inline const int8_t *state2idx(const int8_t *idx_tbl, int hmm_mother, int hmm_father, int gt_mother, int gt_father)
{
    return &idx_tbl[4 * ( N_GT * ( N_GT * ( N_HMM * hmm_mother + hmm_father ) + gt_mother ) + gt_father )];
}

static double beta_binom_idx_log_lkl(beta_binom_t **beta_binom_arr, const int8_t *idx, int16_t ad0, int16_t ad1)
{
    if ( ad0 == bcf_int16_missing || ad1 == bcf_int16_missing ) return 0.0f;
    double log_lkl = idx[0] >= 0 ? beta_binom_log_unsafe(beta_binom_arr[idx[0]], ad0, ad1) :
                                  beta_binom_log_unsafe(beta_binom_arr[-idx[0]], ad1, ad0);
    if ( idx[1] != bcf_int8_missing )
    {
        double log_lkl1 = idx[1] >= 0 ? beta_binom_log_unsafe(beta_binom_arr[idx[1]], ad0, ad1) :
                                       beta_binom_log_unsafe(beta_binom_arr[-idx[1]], ad1, ad0);
        log_lkl = log_mean_exp(log_lkl, log_lkl1);
        if ( idx[2] != bcf_int8_missing )
        {
            double log_lkl2 = idx[2] >= 0 ? beta_binom_log_unsafe(beta_binom_arr[idx[2]], ad0, ad1) :
                                           beta_binom_log_unsafe(beta_binom_arr[-idx[2]], ad1, ad0);
            double log_lkl3 = idx[3] >= 0 ? beta_binom_log_unsafe(beta_binom_arr[idx[3]], ad0, ad1) :
                                           beta_binom_log_unsafe(beta_binom_arr[-idx[3]], ad1, ad0);
            log_lkl = log_mean_exp( log_lkl, log_mean_exp(log_lkl2, log_lkl3) );
        }
    }

    return log_lkl;
}

// adapted Giulio Genovese's implementation in bcftools/vcfmocha.c
void get_max_sum(const int16_t *ad0,
                 const int16_t *ad1,
                 int n,
                 const int *imap,
                 int *n1,
                 int *n2)
{
    *n1 = 0;
    *n2 = 0;
    for (int i=0; i<n; i++)
    {
        int a = imap ? ad0[ imap[i] ] : ad0[i];
        int b = imap ? ad1[ imap[i] ] : ad1[i];
        if ( a!=bcf_int16_missing && b!=bcf_int16_missing )
        {
            if (a > *n1) *n1 = a;
            if (b > *n1) *n1 = b;
            if (a + b > *n2) *n2 = a + b;
        }
    }
}

static void get_emis_log_lkl(const int16_t *ad_ref,
                             const int16_t *ad_alt,
                             const int8_t *gt_mother,
                             const int8_t *gt_father,
                             int T,
                             const int *imap,
                             double maternal_phase_err_prb,
                             double paternal_phase_err_prb,
                             double err_log_prb,
                             double seq_err_log_prb,
                             double ff,
                             double rho_hom,
                             double rho_het,
                             const int *mother_arr,
                             int n_mother,
                             const int *father_arr,
                             int n_father,
                             double *emis_log_lkl)
{
    double log_maternal_phase_err_prb = 0.0, log_1_maternal_phase_err_prb = 0.0;
    if ( maternal_phase_err_prb > 0 )
    {
        log_maternal_phase_err_prb = log( maternal_phase_err_prb ) + M_LN2;
        log_1_maternal_phase_err_prb = log( 1.0 - maternal_phase_err_prb ) + M_LN2;
    }

    double log_paternal_phase_err_prb = 0.0, log_1_paternal_phase_err_prb = 0.0;
    if ( paternal_phase_err_prb > 0 )
    {
        log_paternal_phase_err_prb = log( paternal_phase_err_prb ) + M_LN2;
        log_1_paternal_phase_err_prb = log( 1.0 - paternal_phase_err_prb ) + M_LN2;
    }

    // initialize beta binomial likelihoods look up tables
    int n1, n2;
    get_max_sum(ad_ref, ad_alt, T, imap, &n1, &n2);
    init_log_lkl_mdls(beta_binom_arr, ff, rho_hom, rho_het, n1, n2, expf(seq_err_log_prb));
    for (int t=0; t<T; t++)
    {
        double min_thr = -INFINITY;
        for (int hmm_mother=0; hmm_mother<n_mother; hmm_mother++)
        {
            for (int hmm_father=0; hmm_father<n_father; hmm_father++)
            {
                int16_t ad_r = imap ? ad_ref[ imap[t] ] : ad_ref[t];
                int16_t ad_a = imap ? ad_alt[ imap[t] ] : ad_alt[t];
                int8_t gt_m = imap ? gt_mother[ imap[t] ] : gt_mother[t];
                int8_t gt_f = imap ? gt_father[ imap[t] ] : gt_father[t];
                assert( mother_arr[hmm_mother]>=0 && mother_arr[hmm_mother]<N_HMM );
                assert( father_arr[hmm_father]>=0 && father_arr[hmm_father]<N_HMM );
                assert( gt_m>=0 && gt_m<N_GT );
                assert( gt_f>=0 && gt_f<N_GT );
                const int8_t *idx = state2idx( idx_tbl, mother_arr[hmm_mother], father_arr[hmm_father], gt_m, gt_f );
                double log_prb = beta_binom_idx_log_lkl(beta_binom_arr, idx, ad_r, ad_a );
                if ( maternal_phase_err_prb > 0 && ( gt_m == HET_REF_ALT || gt_m == HET_ALT_REF ) && paternal_phase_err_prb > 0 && ( gt_f == HET_REF_ALT || gt_f == HET_ALT_REF ) )
                {
                    idx = state2idx( idx_tbl, mother_arr[hmm_mother], father_arr[hmm_father], gt_m == HET_REF_ALT ? HET_ALT_REF : HET_REF_ALT, gt_f );
                    double log_prb1 = beta_binom_idx_log_lkl(beta_binom_arr, idx, ad_r, ad_a );
                    idx = state2idx( idx_tbl, mother_arr[hmm_mother], father_arr[hmm_father], gt_m, gt_f == HET_REF_ALT ? HET_ALT_REF : HET_REF_ALT );
                    double log_prb2 = beta_binom_idx_log_lkl(beta_binom_arr, idx, ad_r, ad_a );
                    idx = state2idx( idx_tbl, mother_arr[hmm_mother], father_arr[hmm_father], gt_m == HET_REF_ALT ? HET_ALT_REF : HET_REF_ALT, gt_f == HET_REF_ALT ? HET_ALT_REF : HET_REF_ALT );
                    double log_prb3 = beta_binom_idx_log_lkl(beta_binom_arr, idx, ad_r, ad_a );
                    log_prb = log_mean_exp( log_mean_exp(log_prb + log_1_maternal_phase_err_prb, log_prb1 + log_maternal_phase_err_prb) + log_1_paternal_phase_err_prb,
                                            log_mean_exp(log_prb2 + log_1_maternal_phase_err_prb, log_prb3 + log_maternal_phase_err_prb) + log_paternal_phase_err_prb );
                }
                else if ( maternal_phase_err_prb > 0 && ( gt_m == HET_REF_ALT || gt_m == HET_ALT_REF ) )
                {
                    idx = state2idx( idx_tbl, mother_arr[hmm_mother], father_arr[hmm_father], gt_m == HET_REF_ALT ? HET_ALT_REF : HET_REF_ALT, gt_f );
                    double log_prb1 = beta_binom_idx_log_lkl(beta_binom_arr, idx, ad_r, ad_a );
                    log_prb = log_mean_exp(log_prb + log_1_maternal_phase_err_prb, log_prb1 + log_maternal_phase_err_prb);
                }
                else if ( paternal_phase_err_prb > 0 && ( gt_f == HET_REF_ALT || gt_f == HET_ALT_REF ) )
                {
                    idx = state2idx( idx_tbl, mother_arr[hmm_mother], father_arr[hmm_father], gt_m, gt_f == HET_REF_ALT ? HET_ALT_REF : HET_REF_ALT );
                    double log_prb1 = beta_binom_idx_log_lkl(beta_binom_arr, idx, ad_r, ad_a );
                    log_prb = log_mean_exp(log_prb + log_1_paternal_phase_err_prb, log_prb1 + log_paternal_phase_err_prb);
                }
                emis_log_lkl[ N_HMM * N_HMM * t + N_HMM * mother_arr[hmm_mother] + father_arr[hmm_father] ] = log_prb;
                if ( min_thr < log_prb ) min_thr = log_prb;
            }
        }
        min_thr += err_log_prb;
        for (int hmm_mother=0; hmm_mother<n_mother; hmm_mother++)
        {
            for (int hmm_father=0; hmm_father<n_father; hmm_father++)
            {
                int state = N_HMM * N_HMM * t + N_HMM * mother_arr[hmm_mother] + father_arr[hmm_father];
                if ( emis_log_lkl[ state ] < min_thr ) emis_log_lkl[ state ] = min_thr;
            }
        }
    }
}

/*********************************
 * HMM AND OPTIMIZATION METHODS  *
 *********************************/

// this needs to cycle through 36 hidden states
static void log_viterbi_run(const double *emis_log_lkl,
                            int T,
                            double xy_log_prb,
                            double tel_log_prb,
                            double cen_log_prb,
                            int last_p,
                            int first_q,
                            const int *mother_arr,
                            int n_mother,
                            const int *father_arr,
                            int n_father,
                            int8_t *path)
{
    if (T<1) return;
    double log_prb[N_HMM * N_HMM];
    double new_log_prb[N_HMM * N_HMM];
    // allocate memory necessary for running the algorithm
    int8_t *ptr = (int8_t *)malloc(N_HMM * N_HMM * (T-1) * sizeof(int8_t));

    // initialize and rescale the first state
    for (int hmm_mother=0; hmm_mother<n_mother; hmm_mother++)
    {
        for (int hmm_father=0; hmm_father<n_father; hmm_father++)
        {
            int state = N_HMM * mother_arr[hmm_mother] + father_arr[hmm_father];
            log_prb[state] = hmm_init_tbl[state] * (xy_log_prb - tel_log_prb) + emis_log_lkl[state];
        }
    }

    // compute best probabilities at each position
    for (int t=1; t<T; t++)
    {
        for (int hmm_mother=0; hmm_mother<n_mother; hmm_mother++)
        {
            for (int hmm_father=0; hmm_father<n_father; hmm_father++)
            {
                int state = N_HMM * mother_arr[hmm_mother] + father_arr[hmm_father];
                new_log_prb[state] = log_prb[state];
                ptr[N_HMM * N_HMM * (t-1) + state] = (int8_t)state;
            }
        }

        // compute whether a hidden state switch would be used
        for (int hmm_mother_prev=0; hmm_mother_prev<n_mother; hmm_mother_prev++)
        {
            for (int hmm_father_prev=0; hmm_father_prev<n_father; hmm_father_prev++)
            {
                for (int hmm_mother_next=0; hmm_mother_next<n_mother; hmm_mother_next++)
                {
                    for (int hmm_father_next=0; hmm_father_next<n_father; hmm_father_next++)
                    {
                        int state_prev = N_HMM * mother_arr[hmm_mother_prev] + father_arr[hmm_father_prev];
                        int state_next = N_HMM * mother_arr[hmm_mother_next] + father_arr[hmm_father_next];
                        double trans_log_prb = hmm_trans_tbl[N_HMM * N_HMM * state_prev + state_next]; // using a lookup table speeds things up a little bit
                        if ( log_prb[ state_prev ] + trans_log_prb > new_log_prb[ state_next ] )
                        {
                            new_log_prb[ state_next ] = log_prb[ state_prev ] + trans_log_prb;
                            ptr[N_HMM * N_HMM * (t-1) + state_next] = (int8_t)state_prev;
                        }
                    }
                }
            }
        }

        // update and rescale the current state
        double max = -INFINITY;
        for (int hmm_mother=0; hmm_mother<n_mother; hmm_mother++)
        {
            for (int hmm_father=0; hmm_father<n_father; hmm_father++)
            {
                int state = N_HMM * mother_arr[hmm_mother] + father_arr[hmm_father];
                log_prb[state] = new_log_prb[state] + emis_log_lkl[N_HMM * N_HMM * t + state];
                if ( max < log_prb[state] ) max = log_prb[state];
            }
        }

        // rescale Viterbi log probabilities to avoid underflow issues
        for (int hmm_mother=0; hmm_mother<n_mother; hmm_mother++)
            for (int hmm_father=0; hmm_father<n_father; hmm_father++)
                log_prb[ N_HMM * mother_arr[hmm_mother] + father_arr[hmm_father] ] -= max;
    }

    // add closing cost to the last state
    for (int hmm_mother=0; hmm_mother<n_mother; hmm_mother++)
    {
        for (int hmm_father=0; hmm_father<n_father; hmm_father++)
        {
            int state = N_HMM * mother_arr[hmm_mother] + father_arr[hmm_father];
            log_prb[state] += hmm_init_tbl[state] * (xy_log_prb - tel_log_prb);
        }
    }

    // retrace Viterbi path from probabilities
    path[T-1] = (int8_t)(N_HMM * mother_arr[0] + father_arr[0]);
    for (int hmm_mother=0; hmm_mother<n_mother; hmm_mother++)
        for (int hmm_father=0; hmm_father<n_father; hmm_father++)
        {
            int state = N_HMM * mother_arr[hmm_mother] + father_arr[hmm_father];
            if (log_prb[(int)path[T-1]] < log_prb[state]) path[T-1] = (int8_t)state;
        }

    // compute best path by tracing back the Markov chain
    for (int t=T-1; t>0; t--)
        path[t-1] = ptr[(t-1) * N_HMM * N_HMM + (int)path[t]];

    // free memory
    free(ptr);
}

// this needs to cycle through 36 hidden states
static double get_lod(const double *emis_log_lkl,
                      const int8_t *path,
                      int a,
                      int b)
{
    double lod = 0.0;
    for (int t=a; t<b; t++)
    {
        lod += emis_log_lkl[N_HMM * N_HMM * t + path[t]];
        if ( t > a && path[t] != path[t-1] ) lod += hmm_trans_tbl[N_HMM * N_HMM * path[t-1] + path[t]];
    }
    return lod * M_LOG10E;
}

// this needs to cycle through 36 hidden states
static void trans_count(const int8_t *path,
                        int a,
                        int b,
                        int *m,
                        int *p)
{
    *m = 0;
    *p = 0;
    for (int t=a+1; t<b; t++)
    {
        if ( path[t-1] != path[t] )
        {
            if ( path[t-1] / N_HMM != path[t] / N_HMM ) (*m)++;
            if ( path[t-1] % N_HMM != path[t] % N_HMM ) (*p)++;
        }
    }
}

/*********************************
 * ANALYZE CONTIG                *
 *********************************/

// return segments called by the HMM or state with consecutive call
static int get_path_segs(const int8_t *path,
                         int n,
                         int **beg,
                         int **end)
{
    int beg_m = 0, end_m = 0, a = 0, b = 0;
    *beg = NULL;
    *end = NULL;
    int nseg = 0;
    for (b=0; b<n; b++)
    {
        // check whether it is the end of a segment
        if ( b != n-1 )
        {
            int curr = path[b];
            int next = path[b+1];
            if ( curr == next ) continue;
        }

        nseg++;
        hts_expand(int, nseg, beg_m, *beg);
        (*beg)[nseg-1] = a;
        hts_expand(int, nseg, end_m, *end);
        (*end)[nseg-1] = b;
        a = b + 1;
    }
    return nseg;
}

// process one contig for one sample
static void contig_run(trio_t *self,
                       const model_t *model,
                       stats_t *stats,
                       blueberry_table_t *blueberry_table)
{
    // do nothing if chromosome Y or MT are being tested
    if ( model->rid == model->genome_rules->y_rid || model->rid == model->genome_rules->mt_rid ) return;

    int cen_beg = model->genome_rules->cen_beg[model->rid];
    int cen_end = model->genome_rules->cen_end[model->rid];
    int length = model->genome_rules->length[model->rid];
    double tel_log_prb = model->rid == model->genome_rules->y_rid ? 0.0f : model->tel_log_prb; // because of PAR regions in males
    int n = self->n;
    int *imap_arr = (int *)malloc(n * sizeof(int));

    int last_p = 0, first_q = 0;
    int n_imap = 0;
    for (int i=0; i<n; i++)
    {
        if ( model->pos_arr[ self->vcf_imap_arr[i] ] < cen_beg ) last_p++;
        if ( model->pos_arr[ self->vcf_imap_arr[i] ] < cen_end ) first_q++;
        if ( !(model->flags & HETS_ONLY) || self->gt_mother_arr[i] == HET_ALT_REF || self->gt_mother_arr[i] == HET_REF_ALT || self->gt_mother_arr[i] == HET_UNPHASED )
        {
            n_imap++;
            imap_arr[n_imap - 1] = i;
        }
    }

    double *emis_log_lkl = (double *)malloc(N_HMM * N_HMM * n_imap * sizeof(double));
    get_emis_log_lkl( self->ad_ref_arr,
                      self->ad_alt_arr,
                      self->gt_mother_arr,
                      self->gt_father_arr,
                      n_imap,
                      imap_arr,
                      model->maternal_phase_err_prb,
                      model->paternal_phase_err_prb,
                      model->err_log_prb,
                      model->seq_err_log_prb,
                      self->ff,
                      self->rho_hom,
                      self->rho_het,
                      full_model,
                      6,
                      full_model,
                      6,
                      emis_log_lkl );

    // null model
    int8_t *path0 = (int8_t *)malloc(n_imap * sizeof(int8_t));
    log_viterbi_run(emis_log_lkl, n_imap, model->xy_log_prb, tel_log_prb, model->cen_log_prb, last_p, first_q, null_model, 2, null_model, 2, path0);
    double lod0 = get_lod(emis_log_lkl, path0, 0, n_imap);

    int8_t *path = (int8_t *)malloc(n_imap * sizeof(int8_t));

    // maternal trisomy model
    log_viterbi_run(emis_log_lkl, n_imap, model->xy_log_prb, tel_log_prb, model->cen_log_prb, last_p, first_q, trisomy_model, 3, null_model, 2, path);
    stats->mother_trisomy_lod = get_lod(emis_log_lkl, path, 0, n_imap) - lod0;

    // paternal trisomy model
    log_viterbi_run(emis_log_lkl, n_imap, model->xy_log_prb, tel_log_prb, model->cen_log_prb, last_p, first_q, null_model, 2, trisomy_model, 3, path);
    stats->father_trisomy_lod = get_lod(emis_log_lkl, path, 0, n_imap) - lod0;

    // maternal monosomy model
    log_viterbi_run(emis_log_lkl, n_imap, model->xy_log_prb, tel_log_prb, model->cen_log_prb, last_p, first_q, monosomy_model, 1, null_model, 2, path);
    stats->mother_monosomy_lod = get_lod(emis_log_lkl, path, 0, n_imap) - lod0;

    // paternal monosomy model
    log_viterbi_run(emis_log_lkl, n_imap, model->xy_log_prb, tel_log_prb, model->cen_log_prb, last_p, first_q, null_model, 2, monosomy_model, 1, path);
    stats->father_monosomy_lod = get_lod(emis_log_lkl, path, 0, n_imap) - lod0;

    // full model
    log_viterbi_run(emis_log_lkl, n_imap, model->xy_log_prb, tel_log_prb, model->cen_log_prb, last_p, first_q, full_model, 6, full_model, 6, path);

    int *beg, *end;
    int nseg = get_path_segs(path, n_imap, &beg, &end);
    hts_expand(blueberry_t, blueberry_table->n + nseg, blueberry_table->m, blueberry_table->a);
    for (int i=0; i<nseg; i++)
    {
        blueberry_t *blueberry = &blueberry_table->a[blueberry_table->n];
        blueberry->rid = model->rid;
        if ( beg[i] == 0 )
            if ( model->pos_arr[ self->vcf_imap_arr[ imap_arr[ beg[i] ] ] ] < cen_beg )
                blueberry->beg_pos = 0;
            else
                blueberry->beg_pos = cen_end;
        else blueberry->beg_pos = model->pos_arr[ self->vcf_imap_arr[ imap_arr[ beg[i] ] ] ];
        if ( end[i] == n_imap-1 )
            if ( model->pos_arr[ self->vcf_imap_arr[ imap_arr[ end[i] ] ] ] >= cen_end )
                blueberry->end_pos = length;
            else
                blueberry->end_pos = cen_beg;
        else blueberry->end_pos = model->pos_arr[ self->vcf_imap_arr[ imap_arr[ end[i] ] ] ];
        blueberry->length =  blueberry->end_pos - blueberry->beg_pos;
        if ( blueberry->beg_pos == 0 ) blueberry->p_arm = MOCHA_TEL;
        else if ( blueberry->beg_pos < cen_beg ) blueberry->p_arm = MOCHA_ARM;
        else blueberry->p_arm = MOCHA_NAN;
        if ( blueberry->end_pos == length ) blueberry->q_arm = MOCHA_TEL;
        else if ( blueberry->end_pos > cen_end ) blueberry->q_arm = MOCHA_ARM;
        else blueberry->q_arm = MOCHA_NAN;
        blueberry->nsites = end[i] - beg[i] + 1;
        blueberry->hmm_mother = path[ beg[i] ] / 6;
        blueberry->hmm_father = path[ beg[i] ] % 6;
        if ( ( blueberry->hmm_mother == HMM_L || blueberry->hmm_mother == HMM_R ) && ( blueberry->hmm_father == HMM_L || blueberry->hmm_father == HMM_R ) ) blueberry->lod = NAN;
        else blueberry->lod = get_lod(emis_log_lkl, path, beg[i], end[i]+1) - get_lod(emis_log_lkl, path0, beg[i], end[i]+1);
        blueberry_table->n++;
    }
    free(emis_log_lkl);
    free(beg);
    free(end);
    free(path);
    free(path0);
    free(imap_arr);

    return;
}
/*********************************
 * COMPUTE SAMPLE STATISTICS     *
 *********************************/

// compute mean of an array
static double get_mean(const double *v,
                       int n,
                       const int *imap)
{
    double mean = 0.0;
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        double tmp = (double)(imap ? v[imap[i]] : v[i]);
        if ( !isnan(tmp) )
        {
            mean += tmp;
            j++;
        }
    }
    if ( j <= 0 ) return NAN;
    return mean / (double)j;
}

// this macro from ksort.h defines the function
// double ks_ksmall_double(size_t n, double arr[], size_t kk);
KSORT_INIT_GENERIC(double)

// adapted Giulio Genovese's implementation in bcftools/vcfmocha.c
// compute the median of a vector using the ksort library (with iterator)
static double get_median(const double *v,
                         int n,
                         const int *imap)
{
    if ( n == 0 ) return NAN;
    double tmp, *w = (double *)malloc(n * sizeof(double));
    int j = 0;
    for (int i = 0; i < n; i++)
    {
        tmp = imap ? v[imap[i]] : v[i];
        if (!isnan(tmp)) w[j++] = tmp;
    }
    if ( j == 0 ) { free(w); return NAN; }
    double ret = ks_ksmall_double((size_t)j, w, (size_t)j/2);
    if (j%2==0) ret = (ret + w[j/2-1]) * 0.5f;
    free(w);
    return ret;
}

// this function computes the median of contig stats
static void contig_summary(trio_t *self, model_t *model)
{
    double *tmp_arr = (double *)malloc(self->n_stats * sizeof(double));

    for (int i=0; i<self->n_stats; i++)
    {
        for (int j=0; j<5; j++)
            self->stats.type_snps[j] += self->stats_arr[i].type_snps[j];
        self->stats.maternal_switch_count += self->stats_arr[i].maternal_switch_count;
        self->stats.paternal_switch_count += self->stats_arr[i].paternal_switch_count;
    }

    for (int i=0; i<self->n_stats; i++) tmp_arr[i] = self->stats_arr[i].ff_hom;
    self->stats.ff_hom = get_median( tmp_arr, self->n_stats, NULL );

    for (int i=0; i<self->n_stats; i++) tmp_arr[i] = self->stats_arr[i].ff_het;
    self->stats.ff_het = get_median( tmp_arr, self->n_stats, NULL );

    for (int i=0; i<self->n_stats; i++) tmp_arr[i] = self->stats_arr[i].mean_cov;
    self->stats.mean_cov = get_median( tmp_arr, self->n_stats, NULL );

    if ( isnan( self->ff ) )
    {
        if ( model->flags & FF_SNP_HOM ) self->ff = self->stats.ff_hom;
        else if ( model->flags & FF_SNP_HET ) self->ff = self->stats.ff_het;
    }

    free(tmp_arr);
}

// the terminology is from "Lo, Y. M. D. et al. Maternal plasma DNA sequencing reveals the
// genome-wide genetic and mutational profile of the fetus. Sci. Transl. Med. 2, 61ra91 (2010)"
static int get_snp_type(int8_t gt_mother, int8_t gt_father)
{
    if ( ( gt_mother == HOM_REF && gt_father == HOM_ALT ) || ( gt_mother == HOM_ALT && gt_father == HOM_REF ) ) return 1;
    else if ( ( gt_mother == HOM_REF && gt_father == HOM_REF ) || ( gt_mother == HOM_ALT && gt_father == HOM_ALT ) ) return 2;
    else if ( ( gt_mother == HOM_REF || gt_mother == HOM_ALT ) &&
              ( gt_father == HET_REF_ALT || gt_father == HET_ALT_REF || gt_father == HET_UNPHASED ) ) return 3;
    else if ( gt_mother == HET_REF_ALT || gt_mother == HET_ALT_REF || gt_mother == HET_UNPHASED )
    {
        if ( gt_father == HOM_REF || gt_father == HOM_ALT ) return 4;
        else if ( gt_father == HET_REF_ALT || gt_father == HET_ALT_REF || gt_father == HET_UNPHASED ) return 5;
    }
    return -1;
}

// function that returns number of alternate alleles of euploid fetus
static inline int get_fetal_genotype(int hmm_mother, int hmm_father, int gt_mother, int gt_father)
{
    if ( gt_mother == HET_UNPHASED || gt_father == HET_UNPHASED ||
       ( hmm_mother != HMM_L && hmm_mother != HMM_R ) ||
       ( hmm_father != HMM_L && hmm_father != HMM_R ) ) return -1;
    int alt_count = 0;
    if ( gt_mother == HOM_ALT || ( hmm_mother == HMM_L && gt_mother == HET_ALT_REF ) || ( hmm_mother == HMM_R && gt_mother == HET_REF_ALT ) ) alt_count++;
    if ( gt_father == HOM_ALT || ( hmm_father == HMM_L && gt_father == HET_ALT_REF ) || ( hmm_father == HMM_R && gt_father == HET_REF_ALT ) ) alt_count++;
    return alt_count;
}

// this function computes several contig stats and then clears the contig data from the sample
static void contig_stats(trio_t *self, const model_t *model)
{
    int n = self->n;
    if (n == 0) return;

    double *cov = (double *)malloc(n * sizeof(double));
    int *imap_arr = (int *)malloc(n * sizeof(int));
    int n_imap = 0;
    for (int i=0; i<n; i++) cov[i] = (double)(self->ad_ref_arr[i] + self->ad_alt_arr[i]);

    if ( model->rid == model->genome_rules->x_rid )
    {
        for (int i=0; i<n; i++)
        {
            int pos = model->pos_arr[ self->vcf_imap_arr[i] ];
            if ( pos > model->genome_rules->x_nonpar_beg && pos < model->genome_rules->x_nonpar_end &&
                  ( pos < model->genome_rules->x_xtr_beg || pos > model->genome_rules->x_xtr_end ) )
            {
                n_imap++;
                imap_arr[n_imap - 1] = i;
            }
        }
        self->x_nonpar_cov_mean = get_mean( cov, n_imap, imap_arr );
    }
    else if ( model->rid == model->genome_rules->y_rid )
    {
        int n_imap = 0;
        for (int i=0; i<n; i++)
        {
            int pos = model->pos_arr[ self->vcf_imap_arr[i] ];
            if ( pos > model->genome_rules->y_nonpar_beg && pos < model->genome_rules->y_nonpar_end &&
                  ( pos < model->genome_rules->y_xtr_beg || pos > model->genome_rules->y_xtr_end ) )
            {
                n_imap++;
                imap_arr[n_imap - 1] = i;
            }
        }
        self->y_nonpar_cov_mean = get_mean( cov, n_imap, imap_arr );
    }
    else if ( model->rid == model->genome_rules->mt_rid )
    {
        self->mt_cov_mean = get_mean( cov, n, NULL );
    }

    self->n_stats++;
    hts_expand0(stats_t, self->n_stats, self->m_stats, self->stats_arr);
    self->stats_arr[self->n_stats - 1].rid = model->rid;
    self->stats_arr[self->n_stats - 1].mean_cov = get_mean( cov, n, NULL );
    free(imap_arr);
    free(cov);

    // compute number of type 1,3,4,5 SNPs, fetal fraction, and overdispersion
    int last_p = 0, first_q = 0;
    int fetal_specific_count = 0;
    int non_fetal_specific_count = 0;

    int ref_count_pat_ref = 0;
    int alt_count_pat_ref = 0;
    int ref_count_pat_alt = 0;
    int alt_count_pat_alt = 0;
    int cen_beg = model->genome_rules->cen_beg[model->rid];
    int cen_end = model->genome_rules->cen_end[model->rid];
    int *imap_hom_arr = (int *)malloc(n * sizeof(int));
    int *imap_het_arr = (int *)malloc(n * sizeof(int));
    int n_imap_hom = 0, n_imap_het = 0;
    for (int i=0; i<n; i++)
    {
        if ( model->pos_arr[ self->vcf_imap_arr[i] ] < cen_beg ) last_p++;
        if ( model->pos_arr[ self->vcf_imap_arr[i] ] < cen_end ) first_q++;
        if ( self->ad_ref_arr[i] == bcf_int16_missing || self->ad_alt_arr[i] == bcf_int16_missing ) continue;
        int8_t gt_mother = self->gt_mother_arr[i];
        int8_t gt_father = self->gt_father_arr[i];
        int snp_type = get_snp_type( gt_mother, gt_father );
        if ( snp_type > 0 ) self->stats_arr[self->n_stats - 1].type_snps[snp_type-1]++;

        // use snp_type 1 for fetal fraction estimation
        if ( ( gt_mother == HOM_ALT && gt_father == HOM_REF ) )
        {
            fetal_specific_count += self->ad_ref_arr[i];
            non_fetal_specific_count += self->ad_alt_arr[i];
        }
        else if ( ( gt_mother == HOM_REF && gt_father == HOM_ALT ) )
        {
            fetal_specific_count += self->ad_alt_arr[i];
            non_fetal_specific_count += self->ad_ref_arr[i];
        }
        else if ( ( gt_mother == HET_ALT_REF || gt_mother == HET_REF_ALT || gt_mother == HET_UNPHASED ) && gt_father == HOM_REF )
        {
            ref_count_pat_ref += self->ad_ref_arr[i];
            alt_count_pat_ref += self->ad_alt_arr[i];
        }
        else if ( ( gt_mother == HET_ALT_REF || gt_mother == HET_REF_ALT || gt_mother == HET_UNPHASED ) && gt_father == HOM_ALT )
        {
            ref_count_pat_alt += self->ad_ref_arr[i];
            alt_count_pat_alt += self->ad_alt_arr[i];
        }

        if ( snp_type == 1 || snp_type == 3 ) imap_hom_arr[n_imap_hom++] = i;
        if ( snp_type == 4 || snp_type == 5 ) imap_het_arr[n_imap_het++] = i;
    }

    double ff_hom = 2.0 * (double)fetal_specific_count / (double)( fetal_specific_count + non_fetal_specific_count );
    double ff_het = 2.0f * ( (double)alt_count_pat_alt / (double)( ref_count_pat_alt + alt_count_pat_alt ) - (double)alt_count_pat_ref / (double)( ref_count_pat_ref + alt_count_pat_ref ) );

    int n1, n2;
    get_max_sum(self->ad_ref_arr, self->ad_alt_arr, n, NULL, &n1, &n2);
    double *emis_log_lkl = (double *)malloc(N_HMM * N_HMM * n * sizeof(double));
    int8_t *path = (int8_t *)malloc(n * sizeof(int8_t));

    // count how many switch errors and crossovers are required for the path that best explains the data
    get_emis_log_lkl( self->ad_ref_arr,
                      self->ad_alt_arr,
                      self->gt_mother_arr,
                      self->gt_father_arr,
                      n,
                      NULL,
                      model->maternal_phase_err_prb,
                      model->paternal_phase_err_prb,
                      model->err_log_prb,
                      model->seq_err_log_prb,
                      ff_hom,
                      self->rho_hom,
                      self->rho_het,
                      null_model,
                      2,
                      null_model,
                      2,
                      emis_log_lkl );
    log_viterbi_run(emis_log_lkl, n, model->xy_log_prb, model->tel_log_prb, model->cen_log_prb, last_p, first_q, null_model, 2, null_model, 2, path);
    trans_count(path, 0, n, &self->stats_arr[self->n_stats - 1].maternal_switch_count, &self->stats_arr[self->n_stats - 1].paternal_switch_count);

    if ( model->rid == model->genome_rules->x_rid || model->rid == model->genome_rules->y_rid || model->rid == model->genome_rules->mt_rid ||
      ( self->stats_arr[self->n_stats - 1].type_snps[0] + self->stats_arr[self->n_stats - 1].type_snps[1] + self->stats_arr[self->n_stats - 1].type_snps[2] + self->stats_arr[self->n_stats - 1].type_snps[3] + self->stats_arr[self->n_stats - 1].type_snps[4] <= 1351 ) )
    {
        self->stats_arr[self->n_stats - 1].ff_hom = NAN;
        self->stats_arr[self->n_stats - 1].ff_het = NAN;
        self->stats_arr[self->n_stats - 1].rho_hom = NAN;
        self->stats_arr[self->n_stats - 1].rho_het = NAN;
    }
    else
    {
        self->stats_arr[self->n_stats - 1].ff_hom = ff_hom;
        self->stats_arr[self->n_stats - 1].ff_het = ff_het;
        self->stats_arr[self->n_stats - 1].rho_hom = NAN;
        self->stats_arr[self->n_stats - 1].rho_het = NAN;
    }

    free(emis_log_lkl);
    free(path);
    free(imap_het_arr);
    free(imap_hom_arr);
}

/*********************************
 * OUTPUT METHODS                *
 *********************************/

static void trio_print(FILE *restrict stream,
                       const trio_t *trio,
                       const bcf_hdr_t *hdr,
                       int verbose)
{
    if ( verbose )
    {
        fprintf(stream, "CHROM\tSNP_1\tSNP_2\tSNP_3\tSNP_4\tSNP_5\tFF_HOM\tFF_HET\tRHO_HOM\tRHO_HET\tM_FLIPS\tP_FLIPS\tMEAN_COV\tLOD_MAT3\tLOD_PAT3\tLOD_MAT1\tLOD_PAT1\n");
        for (int i=0; i<trio->n_stats; i++)
        {
            const char *seq_name = bcf_hdr_id2name(hdr, trio->stats_arr[i].rid);
            fprintf(stream, "%s\t%d\t%d\t%d\t%d\t%d\t%.2f%%\t%.2f%%\t%.4f\t%.4f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", seq_name,
            trio->stats_arr[i].type_snps[0], trio->stats_arr[i].type_snps[1], trio->stats_arr[i].type_snps[2],
            trio->stats_arr[i].type_snps[3], trio->stats_arr[i].type_snps[4], 100 * trio->stats_arr[i].ff_hom,
            100 * trio->stats_arr[i].ff_het, trio->stats_arr[i].rho_hom, trio->stats_arr[i].rho_het,
            trio->stats_arr[i].maternal_switch_count, trio->stats_arr[i].paternal_switch_count,
            trio->stats_arr[i].mean_cov, trio->stats_arr[i].mother_trisomy_lod, trio->stats_arr[i].father_trisomy_lod,
            trio->stats_arr[i].mother_monosomy_lod, trio->stats_arr[i].father_monosomy_lod);
        }
        fprintf(stream, "ALL\t%d\t%d\t%d\t%d\t%d\t%.2f%%\t%.2f%%\t%.4f\t%.4f\t%d\t%d\t%.2f\tNA\tNA\tNA\tNA\n",
        trio->stats.type_snps[0], trio->stats.type_snps[1], trio->stats.type_snps[2],
        trio->stats.type_snps[3], trio->stats.type_snps[4], 100 * trio->stats.ff_hom,
        100 * trio->stats.ff_het, trio->stats.rho_hom, trio->stats.rho_het,
        trio->stats.maternal_switch_count, trio->stats.paternal_switch_count, trio->stats.mean_cov);
    } else {
        fprintf(stream, "CHROM\tHETS\tFF\tM_FLIPS\tP_FLIPS\tMEAN_COV\tLOD_MAT3\n");
        for (int i=0; i<trio->n_stats; i++)
        {
            const char *seq_name = bcf_hdr_id2name(hdr, trio->stats_arr[i].rid);
            fprintf(stream, "%s\t%d\t%.2f%%\t%d\t%d\t%.2f\t%.2f\n", seq_name,
            trio->stats_arr[i].type_snps[3] + trio->stats_arr[i].type_snps[4], 100 * trio->stats_arr[i].ff_hom,
            trio->stats_arr[i].maternal_switch_count, trio->stats_arr[i].paternal_switch_count,
            trio->stats_arr[i].mean_cov, trio->stats_arr[i].mother_trisomy_lod);
        }
        fprintf(stream, "ALL\t%d\t%.2f%%\t%d\t%d\t%.2f\tNA\n",
        trio->stats.type_snps[3] + trio->stats.type_snps[4], 100 * trio->stats.ff_hom,
        trio->stats.maternal_switch_count, trio->stats.paternal_switch_count, trio->stats.mean_cov);
    }
    if ( stream != stdout && stream != stderr ) fclose(stream);
}

static void blueberry_print(FILE *restrict stream,
                            const blueberry_t *blueberry,
                            int n,
                            const bcf_hdr_t *hdr,
                            char *genome)
{
    if ( stream == NULL ) return;
    char arm_type[3];
    arm_type[MOCHA_NAN] = 'N';
    arm_type[MOCHA_ARM] = 'Y';
    arm_type[MOCHA_TEL] = 'T';
    fprintf(stream, "CHROM\tBEG_%s\tEND_%s\tLENGTH\tP_ARM\tQ_ARM\tNSITES\tLOD\tMOTHER\tFATHER\n", genome, genome);
    for (int i=0; i<n; i++)
    {
        const char *seq_name = bcf_hdr_id2name(hdr, blueberry->rid);
        fprintf(stream, "%s\t%d\t%d\t%d\t%c\t%c\t%d\t%.4f\t%s\t%s\n",
            seq_name, blueberry->beg_pos, blueberry->end_pos, blueberry->length,
            arm_type[blueberry->p_arm], arm_type[blueberry->q_arm], blueberry->nsites, blueberry->lod,
            hmm_str[blueberry->hmm_mother], hmm_str[blueberry->hmm_father]);
        blueberry++;
    }
    if ( stream != stdout && stream != stderr ) fclose(stream);
}

static void bed_print(FILE *restrict stream,
                      const blueberry_t *blueberry,
                      int n,
                      const bcf_hdr_t *hdr,
                      int paternal)
{
    if ( stream == NULL ) return;
    int prev_rid = blueberry->rid;
    int prev_beg_pos = blueberry->beg_pos;
    int prev_end_pos = blueberry->end_pos;
    int prev_state = paternal ? blueberry->hmm_father : blueberry->hmm_mother;;
    for (int i=1; i<=n; i++)
    {
        blueberry++;
        if ( i==n || prev_rid != blueberry->rid || prev_state != ( paternal ? blueberry->hmm_father : blueberry->hmm_mother ) )
        {
            const char *seq_name = bcf_hdr_id2name(hdr, prev_rid);
            fprintf(stream, "%s\t%d\t%d\t%s\n", seq_name, prev_beg_pos, prev_end_pos, hmm_str[prev_state]);
            if ( i == n ) break;
            prev_rid = blueberry->rid;
            prev_beg_pos = blueberry->beg_pos;
            prev_state = paternal ? blueberry->hmm_father : blueberry->hmm_mother;
        }
        prev_end_pos = blueberry->end_pos;
    }
    if ( stream != stdout && stream != stderr ) fclose(stream);
}

/*********************************
 * VCF READ METHODS              *
 *********************************/

// adapted Giulio Genovese's implementation in bcftools/vcfmocha.c
// retrieve phase information from BCF record
// bcf_int8_missing if phase does not apply
// 0 if phase is not available
// 1 if higher number allele received from the mother
// -1 if higher number allele received from the father
// assumes little endian architecture
int bcf_get_genotype_phase(bcf_fmt_t *fmt,
                           int8_t *gt_phase_arr,
                           int nsmpl)
{
    // bcf_fmt_t *fmt = bcf_get_fmt_id(line, id);
    if ( !fmt || fmt->n != 2 ) return 0;

    #define BRANCH(type_t, bcf_type_vector_end) { \
        type_t *p = (type_t *)fmt->p; \
        for (int i=0; i<nsmpl; i++, p+=2) \
        { \
            if ( p[0]==bcf_type_vector_end || bcf_gt_is_missing(p[0]) || \
                 p[1]==bcf_type_vector_end || bcf_gt_is_missing(p[1]) ) \
            { \
                gt_phase_arr[i] = bcf_int8_missing; \
            } \
            else \
            { \
                type_t gt0 = bcf_gt_allele(p[0]) > 0; \
                type_t gt1 = bcf_gt_allele(p[1]) > 0; \
                if ( gt0 == gt1 ) gt_phase_arr[i] = bcf_int8_missing; \
                else if ( !bcf_gt_is_phased(p[1]) ) gt_phase_arr[i] = 0; \
                else if ( gt1 > gt0 ) gt_phase_arr[i] = 1; \
                else gt_phase_arr[i] = -1; \
            } \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end); break;
        default: error("Unexpected type %d", fmt->type);
    }
    #undef BRANCH

    return 1;
}

// adapted Giulio Genovese's implementation in bcftools/vcfmocha.c
// retrieve genotype alleles information from BCF record
// assumes little endian architecture
int bcf_get_genotype_alleles(const bcf_fmt_t *fmt,
                             int16_t *gt0_arr,
                             int16_t *gt1_arr,
                             int nsmpl)
{
    if ( !fmt || fmt->n != 2 ) return 0;

    // temporarily store genotype alleles in AD array
    #define BRANCH(type_t, bcf_type_vector_end) { \
        type_t *p = (type_t *)fmt->p; \
        for (int i=0; i<nsmpl; i++, p+=2) \
        { \
            if ( p[0]==bcf_type_vector_end || bcf_gt_is_missing(p[0]) || \
                 p[1]==bcf_type_vector_end || bcf_gt_is_missing(p[1]) ) \
            { \
                gt0_arr[i] = bcf_int16_missing; \
                gt1_arr[i] = bcf_int16_missing; \
            } \
            else \
            { \
                gt0_arr[i] = (int16_t)bcf_gt_allele(p[0]); \
                gt1_arr[i] = (int16_t)bcf_gt_allele(p[1]); \
            } \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end); break;
        default: error("Unexpected type %d", fmt->type);
    }
    #undef BRANCH

    return 1;
}

// adapted Giulio Genovese's implementation in bcftools/vcfmocha.c
static inline int bcf_sr_next_line_reader0(bcf_srs_t *sr)
{
    int nret = bcf_sr_next_line(sr);
    while ( nret > 0 && !bcf_sr_has_line(sr, 0) ) nret = bcf_sr_next_line(sr);
    return nret;
}

// retrive allelic depth information from BCF record
// assumes little endian architecture
int bcf_get_ref_alt_depth(const bcf_fmt_t *fmt,
                          int16_t *ad_ref_arr,
                          int16_t *ad_alt_arr,
                          int nsmpl)
{
    if ( !fmt ) return 0;
    int nalleles = fmt->n;

    #define BRANCH(type_t, bcf_type_vector_end, bcf_type_missing) { \
        type_t *p = (type_t *)fmt->p; \
        for (int i=0; i<nsmpl; i++, p+=nalleles) \
        { \
            ad_ref_arr[i] = p[0] == bcf_type_missing ? bcf_int16_missing : (int16_t)p[0]; \
            ad_alt_arr[i] = (int16_t)0; \
            for (int j=1; j<nalleles; j++) if ( p[j] != bcf_type_vector_end && p[j] != bcf_type_missing ) ad_alt_arr[i] += (int16_t)p[j]; \
        } \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH(int8_t, bcf_int8_vector_end, bcf_int8_missing); break;
        case BCF_BT_INT16: BRANCH(int16_t, bcf_int16_vector_end, bcf_int16_missing); break;
        case BCF_BT_INT32: BRANCH(int32_t, bcf_int32_vector_end, bcf_int32_missing); break;
        default: error("Unexpected type %d", fmt->type);
    }
    #undef BRANCH

    return 1;
}

// read data from the VCF
static int get_contig(bcf_srs_t *sr,
                      trio_t *trio,
                      model_t *model,
                      filter_t *filter,
                      int filter_logic)
{
    int rid = model->rid;
    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    bcf_sr_seek(sr, bcf_hdr_id2name( hdr, rid ), 0);

    int nsmpl = bcf_hdr_nsamples(hdr);

    int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
    int ad_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "AD");

    // maybe this should be assert instead
    if ( gt_id < 0 ) error("Error: input VCF file has no GT format field\n");
    if ( ad_id < 0 ) error("Error: input VCF file has no AD format field\n");

    int8_t *phase_arr = (int8_t *)malloc(nsmpl * sizeof(int8_t));
    int16_t *gt0 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int16_t *gt1 = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int16_t *ad_ref = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int16_t *ad_alt = (int16_t *)malloc(nsmpl * sizeof(int16_t));
    int last_het_pos = 0;
    int last_pos = 0;
    trio->n = 0;

    for (model->n=0; bcf_sr_next_line_reader0(sr); model->n++)
    {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        if ( rid != line->rid ) break;
        int pos = line->pos + 1;

        hts_expand(int, model->n+1, model->m_pos, model->pos_arr);
        model->pos_arr[model->n] = pos;

        // if line fails filter clause
        if ( filter )
        {
            int ret = filter_test(filter, line, NULL);
            if ( filter_logic==FLT_INCLUDE ) { if ( !ret ) continue; }
            else if ( ret ) continue;
        }

        // if failing inclusion/exclusion requirement, skip line
        if ( ( model->flags & FLT_EXCLUDE ) && bcf_sr_get_line(sr, 1) ) continue;
        if ( ( model->flags & FLT_INCLUDE ) && !bcf_sr_get_line(sr, 1) ) continue;

        // if site falls in short arm or centromere regions skip line
        if ( !( model->flags & USE_SHORT_ARMS ) && model->genome_rules->is_short_arm[rid] && pos < model->genome_rules->cen_beg[rid] ) continue;
        if ( !( model->flags & USE_CENTROMERES ) && pos > model->genome_rules->cen_beg[rid] && pos < model->genome_rules->cen_end[rid] ) continue;
        if ( ( model->flags & NO_INDELS ) && bcf_get_variant_types(line) == VCF_INDEL ) continue;

        // if there are no genotypes, skip line
        bcf_fmt_t *gt_fmt = bcf_get_fmt_id(line, gt_id);
        bcf_fmt_t *ad_fmt = bcf_get_fmt_id(line, ad_id);

        // TODO change bcf_get_genotype_phase / bcf_get_genotype_alleles / bcf_get_ref_alt_depth into a single function
        if ( !bcf_get_genotype_phase(gt_fmt, phase_arr, nsmpl) ) continue;
        if ( !bcf_get_genotype_alleles(gt_fmt, gt0, gt1, nsmpl) ) continue;
        if ( !bcf_get_ref_alt_depth(ad_fmt, ad_ref, ad_alt, nsmpl) ) continue;

        // skip sites that have no cellfree DNA coverage
        int ref_cnt = ad_ref[model->cellfree_id];
        int alt_cnt = ad_alt[model->cellfree_id];
        if ( ref_cnt == bcf_int16_missing || alt_cnt == bcf_int16_missing || ref_cnt + alt_cnt == 0 ) continue;
        if ( model->min_cov > 0 && ref_cnt + alt_cnt < model->min_cov ) continue;
        if ( model->max_cov > 0 && ref_cnt + alt_cnt > model->max_cov ) continue;

        int8_t gt_mother;
        if ( gt0[model->mother_id] == bcf_int16_missing || gt1[model->mother_id] == bcf_int16_missing ) continue;
        if ( gt0[model->mother_id] == 0 && gt1[model->mother_id] == 0 ) gt_mother = (int8_t)HOM_REF;
        else if ( gt0[model->mother_id] > 0 && gt1[model->mother_id] > 0 ) gt_mother = (int8_t)HOM_ALT;
        else if ( gt0[model->mother_id] == 0 && gt1[model->mother_id] > 0 && phase_arr[model->mother_id] == 0 ) gt_mother = HET_UNPHASED;
        else if ( gt0[model->mother_id] > 0 && gt1[model->mother_id] == 0 && phase_arr[model->mother_id] == 0 ) gt_mother = HET_UNPHASED;
        else if ( gt0[model->mother_id] == 0 && gt1[model->mother_id] > 0 && phase_arr[model->mother_id] != bcf_int8_missing ) gt_mother = HET_REF_ALT;
        else if ( gt0[model->mother_id] > 0 && gt1[model->mother_id] == 0 && phase_arr[model->mother_id] != bcf_int8_missing ) gt_mother = HET_ALT_REF;
        else continue;

        int8_t gt_father;
        if ( gt0[model->father_id] == bcf_int16_missing || gt1[model->father_id] == bcf_int16_missing ) continue;
        if ( gt0[model->father_id] == 0 && gt1[model->father_id] == 0 ) gt_father = (int8_t)HOM_REF;
        else if ( gt0[model->father_id] > 0 && gt1[model->father_id] > 0 ) gt_father = (int8_t)HOM_ALT;
        else if ( gt0[model->father_id] == 0 && gt1[model->father_id] > 0 && phase_arr[model->father_id] == 0 ) gt_father = HET_UNPHASED;
        else if ( gt0[model->father_id] > 0 && gt1[model->father_id] == 0 && phase_arr[model->father_id] == 0 ) gt_father = HET_UNPHASED;
        else if ( gt0[model->father_id] == 0 && gt1[model->father_id] > 0 && phase_arr[model->father_id] != bcf_int8_missing ) gt_father = HET_REF_ALT;
        else if ( gt0[model->father_id] > 0 && gt1[model->father_id] == 0 && phase_arr[model->father_id] != bcf_int8_missing ) gt_father = HET_ALT_REF;
        else continue;

        // skip sites that are either unphased in the mother or in the father (limited use)
        if ( ( model->flags & NO_UNPHASED ) && ( gt_mother == HET_UNPHASED || gt_father == HET_UNPHASED ) ) continue;

        int is_het_pos = ( gt_mother == HET_ALT_REF || gt_mother == HET_REF_ALT );

        // site too close to last selected site
        if ( ( pos < last_het_pos + model->min_dst ) ||
           ( !is_het_pos && ( pos < last_pos + model->min_dst ) ) ) continue;

        // substitute the last site with the current site
        if ( pos < last_pos + model->min_dst ) trio->n--;

        if ( is_het_pos ) last_het_pos = pos;
        last_pos = pos;

        trio->n++;
        hts_expand(int, trio->n, trio->m_vcf_imap, trio->vcf_imap_arr);
        hts_expand(int16_t, trio->n, trio->m_ad_ref, trio->ad_ref_arr);
        hts_expand(int16_t, trio->n, trio->m_ad_alt, trio->ad_alt_arr);
        hts_expand(int8_t, trio->n, trio->m_gt_mother, trio->gt_mother_arr);
        hts_expand(int8_t, trio->n, trio->m_gt_father, trio->gt_father_arr);
        trio->vcf_imap_arr[trio->n - 1] = model->n;
        trio->ad_ref_arr[trio->n - 1] = ref_cnt;
        trio->ad_alt_arr[trio->n - 1] = alt_cnt;
        trio->gt_mother_arr[trio->n - 1] = gt_mother;
        trio->gt_father_arr[trio->n - 1] = gt_father;
    }
    free(phase_arr);
    free(gt0);
    free(gt1);
    free(ad_ref);
    free(ad_alt);

    return model->n;
}

/*********************************
 * PLUGIN CODE                   *
 *********************************/

const char *about(void)
{
    return "Runs SNP-based method for detection of fetal aneuploidies in maternal plasma.\n";
}

static const char *usage_text(void)
{
    return
"\n"
"About: Runs SNP-based method for detection of fetal aneuploidies in maternal plasma. ("VERSION")\n"
"\n"
"Usage: bcftools +blueberry [options] <in.vcf.gz> <cellfree> <mother> <father>\n"
"\n"
"Required options:\n"
"    -r, --rules <assembly>[?]          predefined genome reference rules, 'list' to print available settings, append '?' for details\n"
"    -R, --rules-file <file>            genome reference rules, space/tab-delimited CHROM:FROM-TO,TYPE\n"
"\n"
"General options:\n"
"    -t,   --targets [^]<region>        similar to -r but streams rather than index-jumps. Exclude regions with \"^\" prefix\n"
"    -T,   --targets-file [^]<file>     similar to -R but streams rather than index-jumps. Exclude regions with \"^\" prefix\n"
"    -f,   --apply-filters <list>       require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n"
"    -i/e, --include/--exclude <expr>   select/exclude sites for which the expression is true (see man page for details)\n"
"          --variants [^]<file>         tabix-indexed [compressed] VCF/BCF file containing variants\n"
"                                       to include (or exclude with \"^\" prefix) in the analysis\n"
"          --no-indels                  exclude indels from being used in the model\n"
"          --no-unphased                exclude unphased sites from being used in the model\n"
"          --threads <int>              number of extra output compression threads [0]\n"
"\n"
"Output options:\n"
//"    -o, --output <file>                write output to a file [standard output]\n"
//"    -O, --output-type b|u|z|v          b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
//"        --no-version                   do not append version and command line to the header\n"
"        --verbose <file>               whether output should be verbose\n"
"        --no-log                       suppress progress report on standard error\n"
"        --log <file>                   write log to file [standard error]\n"
"        --stats <file>                 write chromosome statistics to file [standard output]\n"
"        --viterbi <file>               write Viterbi path to file\n"
"        --maternal <file>              write maternal inheritance to file\n"
"        --paternal <file>              write paternal inheritance to file\n"
"\n"
"HMM Options:\n"
"        --ff [HOM,HET]/<float>         cellfree fetal fraction [estimated from HOM SNPs]\n"
"        --rho-hom <float>              beta-binomial overdispersion ["RHO_HOM_DFLT"]\n"
"        --rho-het <float>              beta-binomial overdispersion ["RHO_HET_DFLT"]\n"
"        --min-dist <int>               minimum base pair distance between consecutive sites for WGS data ["MIN_DST_DFLT"]\n"
"        --mat-phase-err-prob <float>   maternal phase switch error probability ["MATERNAL_PHASE_ERR_PRB_DFLT"]\n"
"        --pat-phase-err-prob <float>   paternal phase switch error probability ["PATERNAL_PHASE_ERR_PRB_DFLT"]\n"
"        --min-cov <int>                minimum read coverage for a site to be used (0 for no minimum) ["MIN_COV_DFLT"]\n"
"        --max-cov <int>                maximum read coverage for a site to be used (0 for no maximum) ["MAX_COV_DFLT"]\n"
"        --use-only-hets                use only sites where the mother is heterozygous after fetal fraction estimation\n"
"        --err-prob <float>             uniform error probability ["ERR_PRB_DFLT"]\n"
"        --xy-prob <float>              transition probability ["XY_PRB_DFLT"]\n"
"        --cross-prob <float>           crossover probability ["CROSS_PRB_DFLT"]\n"
"        --mat-switch-err-prob <float>  maternal switch error probability ["MATERNAL_SWITCH_ERR_PRB_DFLT"]\n"
"        --pat-switch-err-prob <float>  paternal switch error probability ["PATERNAL_SWITCH_ERR_PRB_DFLT"]\n"
"        --seq-err-prob <float>         sequencing error probability ["SEQ_ERR_PRB_DFLT"]\n"
"        --telomere-advantage <float>   telomere advantage ["TEL_PRB_DFLT"]\n"
"        --centromere-penalty <float>   centromere penalty ["CEN_PRB_DFLT"]\n"
"        --short_arm_chrs <list>        list of chromosomes with short arms ["SHORT_ARM_CHRS_DFLT"]\n"
"        --use_short_arms               use variants in short arms\n"
"        --use_centromeres              use variants in centromeres\n"
"\n"
"Example:\n"
"    bcftools +blueberry -r GRCh38 cfDNA mother father input.bcf\n"
"\n";
}

static FILE *get_file_handle(const char *str)
{
    FILE *ret;
    if ( strcmp(str, "-") == 0 )
        ret = stdout;
    else
    {
        ret = fopen(str, "w");
        if ( !ret ) error("Failed to open %s: %s\n", str, strerror(errno));
    }
    return ret;
}

int run(int argc, char *argv[])
{
    char *tmp = NULL;
    char *rules = NULL;
    int rules_is_file = 0;
    char *targets_list = NULL;
    int targets_is_file = 0;
    char *filter_fname = NULL;
    int n_threads = 0;
//    char *output_fname = NULL;
//    int output_type = FT_VCF;
//    int record_cmd_line = 1;
    int verbose = 0;
    const char *short_arm_chrs = SHORT_ARM_CHRS_DFLT;
    char *input_fname = NULL;
    const char *ff = "HOM";
    FILE *log_file = stderr;
    FILE *stats_file = stdout;
    FILE *viterbi_file = NULL;
    FILE *maternal_file = NULL;
    FILE *paternal_file = NULL;

    // model parameters
    model_t model;
    memset(&model, 0, sizeof(model_t));
    model.min_dst = (int)strtol(MIN_DST_DFLT, &tmp, 0);
    model.maternal_phase_err_prb = strtod(MATERNAL_PHASE_ERR_PRB_DFLT, &tmp);
    model.paternal_phase_err_prb = strtod(PATERNAL_PHASE_ERR_PRB_DFLT, &tmp);
    model.min_cov = (int)strtol(MIN_COV_DFLT, &tmp, 0);
    model.max_cov = (int)strtol(MAX_COV_DFLT, &tmp, 0);
    model.err_log_prb = log( strtod(ERR_PRB_DFLT, &tmp) );
    model.xy_log_prb = log( strtod(XY_PRB_DFLT, &tmp) );
    model.cross_log_prb = log( strtod(CROSS_PRB_DFLT, &tmp) );
    model.maternal_switch_err_log_prb = NAN;
    model.paternal_switch_err_log_prb = NAN;
    model.seq_err_log_prb = log( strtod(SEQ_ERR_PRB_DFLT, &tmp) );
    model.tel_log_prb = log( strtod(TEL_PRB_DFLT, &tmp) );
    model.cen_log_prb = log( strtod(CEN_PRB_DFLT, &tmp) );
    double maternal_switch_err_prb_dflt = strtod(MATERNAL_SWITCH_ERR_PRB_DFLT, &tmp);
    double paternal_switch_err_prb_dflt = strtod(PATERNAL_SWITCH_ERR_PRB_DFLT, &tmp);

    // trio_t structure
    trio_t trio;
    memset(&trio, 0, sizeof(trio_t));
    trio.stats.ff_hom = NAN;
    trio.stats.ff_het = NAN;
    trio.rho_hom = strtod(RHO_HOM_DFLT, &tmp);
    trio.rho_het = strtod(RHO_HET_DFLT, &tmp);
    trio.ff = NAN;

    // blueberry table
    blueberry_table_t blueberry_table = {0, 0, NULL};

    // initialize beta_binom objects
    for (int i=0; i<15; i++) beta_binom_arr[i] = beta_binom_init();

    // initialize look up table for beta binomials
    beta_binom_idx_tbl_init(idx_tbl);

    // create synced reader object
    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);

    // create filter object
    filter_t *filter = NULL;
    char *filter_str = NULL;
    int filter_logic = 0;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)

    static struct option loptions[] =
    {
        {"rules", required_argument, NULL, 'r'},
        {"rules-file", required_argument, NULL, 'R'},
        {"targets", required_argument, NULL, 't'},
        {"targets-file", required_argument, NULL, 'T'},
        {"apply-filters", required_argument, NULL, 'f'},
        {"exclude", required_argument, NULL, 'e'},
        {"include", required_argument, NULL, 'i'},
        {"variants", required_argument, NULL, 10},
        {"no-indels", no_argument, NULL, 11},
        {"no-unphased", no_argument, NULL, 12},
        {"threads", required_argument, NULL, 9},
        {"output", required_argument, NULL, 'o'},
        {"verbose", no_argument, NULL, 13},
//        {"output-type", required_argument, NULL, 'O'},
//        {"no-version", no_argument, NULL, 8},
        {"log", required_argument, NULL, 14},
        {"no-log", no_argument, NULL, 15},
        {"stats", required_argument, NULL, 16},
        {"viterbi", required_argument, NULL, 17},
        {"maternal", required_argument, NULL, 18},
        {"paternal", required_argument, NULL, 19},
        {"ff", required_argument, NULL, 20},
        {"rho-hom", required_argument, NULL, 21},
        {"rho-het", required_argument, NULL, 22},
        {"min-dist", required_argument, NULL, 23},
        {"mat-phase-err-prob", required_argument, NULL, 24},
        {"pat-phase-err-prob", required_argument, NULL, 25},
        {"min-cov", required_argument, NULL, 26},
        {"max-cov", required_argument, NULL, 27},
        {"use-only-hets", no_argument, NULL, 28},
        {"err-prob", required_argument, NULL, 29},
        {"xy-prob", required_argument, NULL, 30},
        {"cross-prob", required_argument, NULL, 31},
        {"mat-switch-err-prob", required_argument, NULL, 32},
        {"pat-switch-err-prob", required_argument, NULL, 33},
        {"seq-err-prob", required_argument, NULL, 34},
        {"telomere-advantage", required_argument, NULL, 35},
        {"centromere-penalty", required_argument, NULL, 36},
        {"short_arm_chrs", required_argument, NULL, 37},
        {"use_short_arms", no_argument, NULL, 38},
        {"use_centromeres", no_argument, NULL, 39},
        {NULL, 0, NULL, 0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "r:R:t:T:f:e:i:v:o:O:l:c:d:m:p:h?",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'r': rules = optarg; break;
            case 'R': rules = optarg; rules_is_file = 1; break;
            case 't': targets_list = optarg; break;
            case 'T': targets_list = optarg; targets_is_file = 1; break;
            case 'f': sr->apply_filters = optarg; break;
            case 'e': filter_str = optarg; filter_logic |= FLT_EXCLUDE; break;
            case 'i': filter_str = optarg; filter_logic |= FLT_INCLUDE; break;
            case 10 :
                if (optarg[0]=='^')
                {
                    filter_fname = optarg + 1;
                    model.flags |= FLT_EXCLUDE;
                }
                else
                {
                    filter_fname = optarg;
                    model.flags |= FLT_INCLUDE;
                }
                break;
            case  9 :
                n_threads = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --threads %s\n", optarg);
                break;
//            case 'o': output_fname = optarg; break;
//            case 'O':
//                      switch (optarg[0]) {
//                          case 'b': output_type = FT_BCF_GZ; break;
//                          case 'u': output_type = FT_BCF; break;
//                          case 'z': output_type = FT_VCF_GZ; break;
//                          case 'v': output_type = FT_VCF; break;
//                          default: error("The output type \"%s\" not recognised\n", optarg);
//                      }
//                      break;
//            case  8 : record_cmd_line = 0; break;
            case 11 : model.flags |= NO_INDELS; break;
            case 12 : model.flags |= NO_UNPHASED; break;
            case 13 : verbose = 1; break;
            case 14 : log_file = get_file_handle( optarg ); break;
            case 15 : model.flags |= NO_LOG; break;
            case 16 : stats_file = get_file_handle( optarg ); break;
            case 17 : viterbi_file = get_file_handle( optarg ); break;
            case 18 : maternal_file = get_file_handle( optarg ); break;
            case 19 : paternal_file = get_file_handle( optarg ); break;
            case 20: ff = optarg; break;
            case 21 :
                trio.rho_hom = strtod(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --rho %s\n", optarg);
                break;
            case 22 :
                trio.rho_het = strtod(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --rho %s\n", optarg);
                break;
            case 23:
                model.min_dst = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --min-dist %s\n", optarg);
                break;
            case 24:
                model.maternal_phase_err_prb = strtod(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --mat-phase-err %s\n", optarg);
                break;
            case 25:
                model.paternal_phase_err_prb = strtod(optarg, &tmp);
                if ( *tmp ) error("Could not parse: --pat-phase-err %s\n", optarg);
                break;
            case 26 :
                model.min_cov = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --min-cov %s\n", optarg);
                break;
            case 27 :
                model.max_cov = (int)strtol(optarg, &tmp, 0);
                if ( *tmp ) error("Could not parse: --max-cov %s\n", optarg);
                break;
            case 28 : model.flags |= HETS_ONLY; break;
            case 29 :
                model.err_log_prb = log( strtod(optarg, &tmp) );
                if ( *tmp ) error("Could not parse: --err-prob %s\n", optarg);
                break;
            case 30 :
                model.xy_log_prb = log( strtod(optarg, &tmp) );
                if ( *tmp ) error("Could not parse: --xy-prob %s\n", optarg);
                break;
            case 31 :
                model.cross_log_prb = log( strtod(optarg, &tmp) );
                if ( *tmp ) error("Could not parse: --cross-prob %s\n", optarg);
                break;
            case 32 :
                model.maternal_switch_err_log_prb = log( strtod(optarg, &tmp) );
                if ( *tmp ) error("Could not parse: --maternal-switch-prob %s\n", optarg);
                break;
            case 33 :
                model.paternal_switch_err_log_prb = log( strtod(optarg, &tmp) );
                if ( *tmp ) error("Could not parse: --paternal-switch-prob %s\n", optarg);
                break;
            case 34 :
                model.seq_err_log_prb = log( strtod(optarg, &tmp) );
                if ( *tmp ) error("Could not parse: --seq-err-prob %s\n", optarg);
                break;
            case 35 :
                model.tel_log_prb = log( strtod(optarg, &tmp) );
                if ( *tmp ) error("Could not parse: --telomere-advantage %s\n", optarg);
                break;
            case 36 :
                model.cen_log_prb = log( strtod(optarg, &tmp) );
                if ( *tmp ) error("Could not parse: --centromere-penalty %s\n", optarg);
                break;
            case 37 : short_arm_chrs = optarg; break;
            case 38 : model.flags |= USE_SHORT_ARMS; break;
            case 39 : model.flags |= USE_CENTROMERES; break;
            case 'h':
            case '?':
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( !rules )
    {
        fprintf(log_file, "Genome reference assembly was not specified with --rules or --rules-file\n");
        error("%s", usage_text());
    }
    int len = strlen(rules);
    if ( !rules_is_file && ( strncmp(rules, "list", 4) == 0 || rules[len-1]=='?' ) ) genome_init_alias(log_file, rules, NULL);
    if ( filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( argc-optind!=4 ) error("%s", usage_text());

    if ( strcmp(ff, "HOM") == 0 ) model.flags |= FF_SNP_HOM;
    else if ( strcmp(ff, "HET") == 0 ) model.flags |= FF_SNP_HET;
    else
    {
        trio.ff = strtod(ff, &tmp);
        if ( *tmp ) error("Could not parse: --ff %s\n", optarg);
    }

    // initialize HMM lookup tables
    init_hmm_init_tbl(hmm_init_tbl);
    init_hmm_trans_tbl(hmm_trans_tbl,
                       model.xy_log_prb,
                       model.cross_log_prb,
                       isnan(model.maternal_switch_err_log_prb) ? log(maternal_switch_err_prb_dflt) : model.maternal_switch_err_log_prb,
                       isnan(model.paternal_switch_err_log_prb) ? log(paternal_switch_err_prb_dflt) : model.paternal_switch_err_log_prb);

    if ( targets_list )
    {
        if ( bcf_sr_set_targets(sr, targets_list, targets_is_file, 0) < 0 )
            error("Failed to read the targets: %s\n", targets_list);
    }
    if ( bcf_sr_set_threads(sr, n_threads)<0 ) error("Failed to create threads\n");

    input_fname = argv[optind];
    if ( !bcf_sr_add_reader(sr, input_fname) ) error("Failed to open %s: %s\n", input_fname, bcf_sr_strerror(sr->errnum));
    optind++;
    if ( filter_fname && !bcf_sr_add_reader(sr, filter_fname) ) error("Failed to open %s: %s\n", filter_fname, bcf_sr_strerror(sr->errnum));

    // check whether the necessary information has been included in the VCF
    bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
    if ( filter_str ) filter = filter_init(hdr, filter_str);
    if ( bcf_hdr_nsamples(hdr) == 0 )  error("Error: input VCF file has no samples\n");
    if ( bcf_hdr_id2int( hdr, BCF_DT_ID, "GT" ) < 0 ) error("Error: input VCF file has no GT format field\n");
    if ( bcf_hdr_id2int( hdr, BCF_DT_ID, "AD" ) < 0 ) error("Error: input VCF file has no AD format field\n");

    model.cellfree_id = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, argv[optind]);
    if ( model.cellfree_id < 0 ) error("Cellfree sample %s not found in VCF %s\n", argv[optind], input_fname);
    optind++;
    model.mother_id = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, argv[optind]);
    if ( model.mother_id < 0 ) error("Mother sample %s not found in VCF %s\n", argv[optind], input_fname);
    optind++;
    model.father_id = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, argv[optind]);
    if ( model.father_id < 0 ) error("Father sample %s not found in VCF %s\n", argv[optind], input_fname);
    optind++;

    // initialize genome parameters
    if ( rules_is_file ) model.genome_rules = genome_init_file(rules, hdr);
    else model.genome_rules = genome_init_alias(log_file, rules, hdr);
    if ( !(model.flags & NO_LOG) ) fprintf(log_file, "Using genome assembly from %s\n", rules);
    readlist_short_arms(model.genome_rules, short_arm_chrs, hdr);

    for (int rid=0; rid < hdr->n[BCF_DT_CTG]; rid++)
    {
        model.rid = rid;
        int nret = get_contig(sr, &trio, &model, filter, filter_logic);
        if ( nret<=0 ) continue;
        if ( !(model.flags & NO_LOG) ) fprintf(log_file, "Read %d variants from contig %s while used %d\n", nret, bcf_hdr_id2name( hdr, rid ), trio.n);

        if( model.genome_rules->length[rid] < model.pos_arr[model.n-1] )
            model.genome_rules->length[rid] = model.pos_arr[model.n-1];
        contig_stats(&trio, &model);
    }

    contig_summary(&trio, &model);
    if ( !(model.flags & NO_LOG) ) fprintf(log_file, "Mean coverage: %.4f\n", trio.stats.mean_cov);
    if ( !(model.flags & NO_LOG) ) fprintf(log_file, "Cellfree fetal fraction: %.4f%%\n", 100.0f * trio.ff );
    if ( !(model.flags & NO_LOG) ) fprintf(log_file, "Rho at homozygous sites: %.4f\n", trio.rho_hom);
    if ( !(model.flags & NO_LOG) ) fprintf(log_file, "Rho at heterozygous sites: %.4f\n", trio.rho_het);

    // estimate maternal and paternal phase switch error
    if ( isnan(model.maternal_switch_err_log_prb) )
    {
        model.maternal_switch_err_log_prb = log( (double)( 1 + trio.stats.maternal_switch_count ) / (double)( 1 + trio.stats.type_snps[3-1] + trio.stats.type_snps[5-1] ) );
    }
    fprintf(log_file, "Maternal phase switch error probability: %.4f\n", exp(model.maternal_switch_err_log_prb) );

    if ( isnan(model.paternal_switch_err_log_prb) )
    {
        model.paternal_switch_err_log_prb = log( (double)( 1 + trio.stats.paternal_switch_count ) / (double)( 1 + trio.stats.type_snps[4-1] + trio.stats.type_snps[5-1] ) );
    }
    fprintf(log_file, "Paternal phase switch error probability: %.4f\n", exp(model.paternal_switch_err_log_prb) );

    init_hmm_trans_tbl(hmm_trans_tbl, model.xy_log_prb, model.cross_log_prb, model.maternal_switch_err_log_prb, model.paternal_switch_err_log_prb);

    for (int i=0; i<trio.n_stats; i++)
    {
        model.rid = trio.stats_arr[i].rid;
        int nret = get_contig(sr, &trio, &model, filter, filter_logic);
        if ( nret<=0 ) continue;
        if ( !(model.flags & NO_LOG) ) fprintf(log_file, "Read %d variants from contig %s while used %d\n", nret, bcf_hdr_id2name( hdr, model.rid ), trio.n);
        contig_run(&trio, &model, &trio.stats_arr[i], &blueberry_table);
    }

//     htsFile *out_fh = hts_open(output_fname ? output_fname : "-", hts_bcf_wmode(output_type));
//     if ( !out_fh ) error("Can't write to \"%s\": %s\n", output_fname, strerror(errno));
//     if ( n_threads ) hts_set_opt(out_fh, HTS_OPT_THREAD_POOL, sr->p);
//
//     if (record_cmd_line) bcf_hdr_append_version(out_hdr, argc, argv, "bcftools_+importFMT");
//     if ( bcf_hdr_write(out_fh, out_hdr) < 0 ) error("Unable to write to output VCF file\n");

    trio_print(stats_file, &trio, hdr, verbose);
    free(trio.stats_arr);
    free(trio.vcf_imap_arr);
    free(trio.ad_ref_arr);
    free(trio.ad_alt_arr);
    free(trio.gt_mother_arr);
    free(trio.gt_father_arr);

    blueberry_print(viterbi_file, blueberry_table.a, blueberry_table.n, hdr, rules);
    bed_print(maternal_file, blueberry_table.a, blueberry_table.n, hdr, 0);
    bed_print(paternal_file, blueberry_table.a, blueberry_table.n, hdr, 1);
    free(blueberry_table.a);

    if ( log_file != stdout && log_file != stderr ) fclose(log_file);

    // clear model data
    for (int i=0; i<15; i++) beta_binom_destroy(beta_binom_arr[i]);
    genome_destroy(model.genome_rules);
    free(model.pos_arr);

//     bcf_hdr_destroy(out_hdr);
//     hts_close(out_fh);
    bcf_sr_destroy(sr);

    return 0;
}
