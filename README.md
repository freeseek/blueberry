![](blueberry_logo.png)

This site contains the instructions necessary to reproduce the analyses that will be published in:
```
Genovese G., McCarroll S. et al. Chromosomal phasing improves aneuploidy
determination in non-invasive prenatal testing at low fetal fractions
```
For any feedback, send an email to giulio.genovese@gmail.com

Installation
============

Install basic tools (Debian/Ubuntu specific)
```
sudo apt install wget gzip unzip bcftools samtools bwa
```

Optionally, you can install these libraries to activate further bcftools features:
```
sudo apt install liblzma-dev libbz2-dev libgsl0-dev
```

Preparation steps
```
mkdir -p $HOME/bin && cd /tmp
```

Download latest version of <a href="https://github.com/samtools/htslib">HTSlib</a> and <a href="https://github.com/samtools/bcftools">BCFtools</a> (if not downloaded already)
```
git clone --branch=develop git://github.com/samtools/htslib.git
git clone --branch=develop git://github.com/samtools/bcftools.git
```

Add patch and code for plugins
```
/bin/rm -f bcftools/{vcfnorm.patch,{beta_binom,genome_rules}.{c,h}} bcftools/plugins/{trio-phase,add-variant-dist,blueberry}.c
wget -P bcftools https://raw.githubusercontent.com/freeseek/mocha/master/{vcfnorm.patch,{beta_binom,genome_rules}.{c,h}}
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/mocha/master/trio-phase.c
wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/blueberry/master/{add-variant-dist,blueberry}.c
cd bcftools && patch < vcfnorm.patch && cd ..
```
If for any reason the patch fails with an error message, contact the <a href="mailto:giulio.genovese@gmail.com">author</a> for a fix

Compile latest version of HTSlib (optionally disable bz2 and lzma) and BCFtools (make sure you are using gcc version 5 or newer or else include the -std=gnu99 compilation flag)
```
cd htslib && autoheader && (autoconf || autoconf) && ./configure --disable-bz2 --disable-lzma && make && cd ..
cd bcftools && make && cd ..
/bin/cp bcftools/{bcftools,plugins/{trio-phase,add-variant-dist,blueberry}.so} $HOME/bin/
```
Notice that you will need some functionalities missing from the base version of bcftools to run the pipeline

Make sure the directory with the plugins is available to bcftools
```
export BCFTOOLS_PLUGINS=$HOME/bin
```

Optionally you can also download the code to run the simulations described in the paper
```
sudo apt install libgsl-dev
wget https://raw.githubusercontent.com/freeseek/blueberry/master/simulate.c
gcc -g -O2 simulate.c -lm -lgsl -lgslcblas -o $HOME/bin/simulate
simulate -h
```

Download resources
==================

Download the GRCh37 human genome reference
```
wget -O- ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz | \
  gzip -d > human_g1k_v37.fasta
samtools faidx human_g1k_v37.fasta
```

Download the GRCh38 human genome reference and raw genotypes
```
wget -O- ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | \
  gzip -d > GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

Download software and resources necessary to run the analyses
```
wget https://github.com/broadinstitute/picard/releases/download/2.19.0/picard.jar

wget -O eagle https://data.broadinstitute.org/alkesgroup/Eagle/downloads/dev/eagle_v2.4.1
chmod a+x eagle
wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz

wget https://github.com/broadinstitute/gatk/releases/download/4.0.12.0/gatk-4.0.12.0.zip
unzip gatk-4.0.12.0.zip
https://github.com/broadinstitute/gatk/releases/download/4.1.3.0/gatk-4.1.3.0.zip
unzip gatk-4.1.3.0.zip
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/{hapmap_3.3,1000G_omni2.5,1000G_phase1.snps.high_confidence,Mills_and_1000G_gold_standard.indels,Axiom_Exome_Plus.genotypes.all_populations.poly}.hg38.vcf.gz{,.csi}
```

Download 1000 Genomes Project phase 3 (fixing contig names, removing duplicate variants, removing incomplete variants) for GRCh38
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr{{1..22},X,Y}_GRCh38.genotypes.20170504.vcf.gz{,.tbi}
for chr in chr{{1..22},X,Y}; do
  (bcftools view --no-version -h ALL.${chr}_GRCh38.genotypes.20170504.vcf.gz | \
    grep -v "^##contig=<ID=[GNh]" | sed 's/^##contig=<ID=MT/##contig=<ID=chrM/;s/^##contig=<ID=\([0-9XY]\)/##contig=<ID=chr\1/'; \
  bcftools view --no-version -H -c 2 ALL.${chr}_GRCh38.genotypes.20170504.vcf.gz | \
  grep -v "[0-9]|\.\|\.|[0-9]" | sed 's/^/chr/') | \
  bcftools norm --no-version -Ou -m -any | \
  bcftools norm --no-version -Ob -o ALL.${chr}_GRCh38.genotypes.20170504.bcf \
    -d none -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna && \
  bcftools index -f ALL.${chr}_GRCh38.genotypes.20170504.bcf
done
```

The following instructions were designed to work with the GRCh38 human genome reference but they should be easily adaptable to other human genome references if necessary. The computational steps to align and process the sequence data are very computationally intensive and should be parallelized to be run quickly

Microarray genotype data (pre-processing)
=========================================

These pre-processing steps were run by the owners of the microarray data to anonymize data. They are reported for full disclosure of the steps taken. Microarray data was manually downloaded from 23andMe, AncestryDNA, and FamilyTreeDNA on November 9th, 2018 and converted to VCF

23andMe genotypes
```
zcat genome_xxxxxx_xxxxxx_v#_Full_yyyymmddhhmmss.zip" | sed 's/\t\([ACGT]\)\r$/\t\1\1/' | \
  bcftools convert --no-version -Ou --tsv2vcf - -f human_g1k_v37.fasta -s xxxxxx | \
  bcftools sort -Ob -o xxxxxx.me.GRCh37.bcf -T ./bcftools-sort.XXXXXX
```

AncestryDNA genotypes
```
zcat xxxxxx_xxxxxx_dna-data-yyyy-mm-dd.zip | sed -e 's/\(.*\)\t/\1/' -e 's/00\r$/--/' | \
  bcftools convert --no-version -Ou --tsv2vcf - -f human_g1k_v37.fasta -s xxxxxx | \
  bcftools sort -Ob -o xxxxxx.ad.GRCh37.bcf -T ./bcftools-sort.XXXXXX
```

FamilyTreeDNA genotypes
```
zcat 37_x_xxxxxx_Chrom_Autoso_yyyymmdd.csv.gz | sed -e 's/"//g' -e 's/,/\t/g' | \
  bcftools convert --no-version -Ou --tsv2vcf - -f human_g1k_v37.fasta -s xxxxxx | \
  bcftools sort -Ob -o xxxxxx.ft.GRCh37.bcf -T ./bcftools-sort.XXXXXX
```

For each provider genotypes were merged into a single VCF
```
for vcf in *.{me,ad,ft}.GRCh37.bcf; do bcftools index -f $vcf; done
for pfx in me ad ft; do bcftools merge --no-version -Ob -o $pfx.GRCh37.bcf *.$pfx.GRCh37.bcf; done
/bin/rm *.{me,ad,ft}.GRCh37.bcf{,.csi}
```

The three VCFs above, together a pedigree file describing the trio relationships, are available for download in this repository

Microarray genotype data
========================

Download the microarray genotype data (~39MB)
```
wget https://raw.githubusercontent.com/freeseek/blueberry/master/{{me,ad,ft}.GRCh37.bcf,blueberry.ped}
for vcf in {me,ad,ft}.GRCh37.bcf; do bcftools index -f $vcf; done
```

Lift genotypes from GRCh37 to GRCh38
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod a+x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
for pfx in me ad ft; do
  bcftools query -f "%CHROM\t%POS\t%ID\t%REF,%ALT[\t%GT]\n" $pfx.GRCh37.bcf | sed 's/^MT/M/' | \
    awk '{split($4,a,","); printf "chr%s\t%d\t%d\t%s",$1,$2-1,$2,$3;
      for (i=5; i<=NF; i++) if ($i=="./.") printf "_--"; else printf "_"a[substr($i,1,1)+1]a[substr($i,3,1)+1];
      printf "\t0\t+\n"}' | \
    ./liftOver \
      /dev/stdin \
      hg19ToHg38.over.chain.gz \
      /dev/stdout \
      /dev/null | \
    awk 'BEGIN {rev["A"]="T"; rev["C"]="G"; rev["G"]="C"; rev["T"]="A"; rev["-"]="-"}
      {split($4,a,"_"); printf "%s %s %d",a[1],$1,$3}
      $6=="+" {for (i=2; i in a; i++) printf " %s",a[i]}
      $6=="-" {for (i=2; i in a; i++) printf " %s%s",rev[substr(a[i],1,1)],rev[substr(a[i],2,1)]}
      {printf "\n"}' | \
    bcftools convert --no-version -Ou --tsv2vcf - -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
      -s $(bcftools query -l $pfx.GRCh37.bcf | tr '\n' ',' | sed 's/,$//') | \
    bcftools sort -Ob -o $pfx.GRCh38.bcf -T ./bcftools-sort.XXXXXX && \
  bcftools index -f $pfx.GRCh38.bcf
done
```

Merge genotype data from the three providers
```
bcftools merge --no-version --force-samples {me,ad,ft}.GRCh38.bcf | sed 's/1\/0/0\/1/g' | \
  awk -F"\t" -v OFS="\t" -v n=$(bcftools query -l me.GRCh38.bcf | wc -l) '$0!~"^#" {
    for (me=10; me<10+n; me++) {ad=me+n; ft=ad+n;
    if ($me=="./." && $ad=="./.") gt=$ft;
    else if ($me=="./." && $ft=="./.") gt=$ad;
    else if ($ad=="./." && $ft=="./.") gt=$me;
    else if ($me==$ad) gt=$me;
    else if ($ad==$ft) gt=$ad;
    else if ($me==$ft) gt=$ft;
    else gt="./."; $me=gt} }
    {NF=9+n; print}' | \
  bcftools annotate --no-version -Ou -x ID | \
  bcftools norm --no-version -Ob -o blueberry.array.GRCh38.bcf -m -any && \
  bcftools index -f blueberry.array.GRCh38.bcf
```

Clean up
```
/bin/rm {me,ad,ft}.GRCh3[78].bcf{,.csi}
```

Process Personal Genome Project sequence data
=============================================

Download the sequence data available on the Personal Genome Project portal (~227GB)
```
wget https://my.pgp-hms.org/user_file/download/3728
unzip 56001801068754A.zip
wget https://my.pgp-hms.org/user_file/download/3729
unzip 56001801068814A.zip
wget https://my.pgp-hms.org/user_file/download/3733
unzip 2115760-20190109T222357Z-001.zip
wget https://my.pgp-hms.org/user_file/download/3734
unzip 20181127.zip
```

Align sequence data against the GRCh38 human genome reference
```
r1=( $(cut -f1 blueberry.tsv) )
r2=( $(cut -f2 blueberry.tsv) )
id=( $(cut -f3 blueberry.tsv) )
pu=( $(cut -f4 blueberry.tsv) )
lb=( $(cut -f5 blueberry.tsv) )
sm=( $(cut -f6 blueberry.tsv) )

for i in $(seq 1 $(tail -n+2 blueberry.tsv | wc -l)); do \
  str="@RG\tID:${id[$i]}\tPL:ILLUMINA\tPU:${pu[$i]}\tLB:${lb[$i]}\tSM:${sm[$i]}"; \
  out=${f1%.fastq.gz}; \
  out=${out%.fq.gz}; \
  out=${out%_R[12]_001}; \
  out=${out%_1}; \
  if [ ${r2[$i]} == "NA" ]; then \
    fastq="${r1[$i]}"; \
  else \
    fastq="${r1[$i]} ${r2[$i]}"; \
  fi; \
  bwa mem -t 4 -M -R "$str" GCA_000001405.15_GRCh38_no_alt_analysis_set.fna $fastq | samtools view -Sb - | \
    samtools sort - -o $out.bam && \
    samtools index $out.bam; \
done
```

Remove duplicate reads (skipped for SNP-based data)
```
for sm in 56001801068{75,81}4A 2115760 20181127; do \
  if [ $sm == "2115760" ]; then \
    samtools merge $sm.tmp.bam 45{45949_BC10,40024_BC264,54782_BC28}.bam; \
  else \
    java \
      -jar picard.jar \
      MarkDuplicates \
      $(grep $sm$ fastq.tsv | cut -f1 | sed 's/^/I=/;s/_R1_001.fastq.gz/.bam/;s/_1.fq.gz/.bam/;s/.fastq.gz/.bam/' | tr '\n' ' ') \
      O=$sm.tmp.bam \
      M=$sm.txt; \
  fi && \
  samtools index $sm.tmp.bam; \
done
```

Clean up
```
for i in $(seq 1 $(tail -n+2 blueberry.tsv | wc -l)); do \
  out=${f1%.fastq.gz}; \
  out=${out%.fq.gz}; \
  out=${out%_R[12]_001}; \
  out=${out%_1}; \
  /bin/rm $out.bam{,.bai}; \
done
```

Recalibrate base pairs
```
for sm in 56001801068{75,81}4A 2115760 20181127; do \
  gatk-4.1.3.0/gatk \
    BaseRecalibrator \
    -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    $sm.tmp.bam \
    --known-sites 1000G_phase1.snps.high_confidence.vcf.gz \
    -O $sm.grp && \
  gatk-4.1.3.0/gatk \
    ApplyBQSR \
    -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    $sm.tmp.bam \
    --bqsr-recal-file $sm.grp \
    -O $sm.bam && \
  samtools index $sm.bam; \
done
```

Clean up
```
/bin/rm {56001801068{75,81}4A,2115760,20181127}.tmp.bam{,.bai}
```

Genotype calling of Dante Labs MPSS data
========================================

```
for chr in chr{{1..22},X,Y,M}; do \
  gatk-4.1.3.0/gatk \
    HaplotypeCaller \
    --max-reads-per-alignment-start 0 \
    -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -I 56001801068814A.bam \
    -I 56001801068754A.bam \
    -O blueberry.hc.$chr.GRCh38.vcf.gz \
    -L $chr \
    --create-output-variant-index true; \
done
```

```
bcftools concat \
  --no-version \
  -Oz -o blueberry.hc.GRCh38.vcf.gz \
  blueberry.hc.chr{{1..22},X,Y,M}.GRCh38.vcf.gz && \
bcftools index -f -t blueberry.hc.GRCh38.vcf.gz
```

Clean up
```
/bin/rm blueberry.hc.chr{{1..22},X,Y,M}.GRCh38.vcf.gz{,.tbi}
```

Variant quality score recalibration
===================================

The following code was inspired by the VQSR <a href="https://software.broadinstitute.org/gatk/documentation/article.php?id=2805">howto</a>

Recalibrate SNP quality scores
```
gatk-4.1.3.0/gatk \
  VariantRecalibrator \
  -V blueberry.hc.GRCh38.vcf.gz \
  -O blueberry.hc.GRCh38.recalibrate_snp.recal \
  --tranches-file blueberry.hc.GRCh38.recalibrate_snp.tranches \
  --trust-all-polymorphic \
  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  -an DP \
  -an QD \
  -an FS \
  -an SOR \
  -an MQ \
  -an MQRankSum \
  -an ReadPosRankSum \
  -mode SNP \
  --max-gaussians 6 \
  -resource:hapmap,known=false,training=true,truth=true,prior=15 hapmap_3.3.hg38.vcf.gz \
  -resource:omni,known=false,training=true,truth=true,prior=12 1000G_omni2.5.hg38.vcf.gz \
  -resource:1000G,known=false,training=true,truth=false,prior=10 1000G_phase1.snps.high_confidence.hg38.vcf.gz
```

```
gatk-4.1.3.0/gatk \
  ApplyVQSR \
  -O blueberry.snp.recal.GRCh38.vcf.gz \
  -V blueberry.hc.GRCh38.vcf.gz \
  --recal-file blueberry.hc.GRCh38.recalibrate_snp.recal \
  --tranches-file blueberry.hc.GRCh38.recalibrate_snp.tranches \
  --truth-sensitivity-filter-level 99.0 \
  --create-output-variant-index true \
  -mode SNP
```

Recalibrate indel quality scores
```
gatk-4.1.3.0/gatk \
  VariantRecalibrator \
  -V blueberry.snp.recal.GRCh38.vcf.gz \
  -O blueberry.snp.recal.GRCh38.recalibrate_indel.recal \
  --tranches-file blueberry.snp.recal.GRCh38.recalibrate_indel.tranches \
  --trust-all-polymorphic \
  -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
  -an QD \
  -an DP \
  -an FS \
  -an SOR \
  -an MQRankSum \
  -an ReadPosRankSum \
  -mode INDEL \
  --max-gaussians 4 \
  -resource:mills,known=false,training=true,truth=true,prior=12 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  -resource:axiomPoly,known=false,training=true,truth=false,prior=10 Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
```

```
gatk-4.1.3.0/gatk \
  ApplyVQSR \
  -O blueberry.vqsr.GRCh38.vcf.gz \
  -V blueberry.snp.recal.GRCh38.vcf.gz \
  --recal-file blueberry.snp.recal.GRCh38.recalibrate_indel.recal \
  --tranches-file blueberry.snp.recal.GRCh38.recalibrate_indel.tranches \
  --truth-sensitivity-filter-level 99.0 \
  --create-output-variant-index true \
  -mode INDEL
```

Phase genotype data
===================

Normalize data
```
bcftools filter --no-version -Ou -e "FMT/DP<10 | FMT/GQ<30" --set-GT . blueberry.vqsr.GRCh38.vcf.gz | \
  bcftools annotate --no-version -Ou -x ID,QUAL,INFO,^FMT/GT,^FMT/AD | \
  bcftools norm --no-version -Ou -m -any -k | \
  bcftools norm --no-version -Ou -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | \
  bcftools view --no-version -Ou -e 'ALT=="*"' | \
  bcftools +add-variant-dist --no-version -Oz -o blueberry.unphased.GRCh38.vcf.gz && \
  bcftools index -f -t blueberry.unphased.GRCh38.vcf.gz
```

Generate list of variants for Eagle to exclude from the phasing process
```
bcftools view --no-version -Ob -o blueberry.exclude.GRCh38.bcf -G \
  -i 'FILTER!="PASS" || SNP_DIST<=20 || INDEL_DIST<=20' blueberry.unphased.GRCh38.vcf.gz && \
  bcftools index -f blueberry.exclude.GRCh38.bcf
```

Use Eagle to phase the autosomes and the X chromosome
```
for chr in chr{{1..22},X}; do \
  eagle \
    --geneticMapFile genetic_map_hg38_withX.txt.gz \
    --outPrefix blueberry.1000g_phased.$chr.GRCh38 \
    --vcfRef ALL.${chr}_GRCh38.genotypes.20170504.bcf \
    --vcfTarget blueberry.unphased.GRCh38.vcf.gz \
    --vcfOutFormat b \
    --noImpMissing \
    --outputUnphased \
    --vcfExclude blueberry.exclude.GRCh38.bcf \
    --chrom $chr \
    --pbwtIters 2 && \
  bcftools index -f blueberry.1000g_phased.$chr.GRCh38.bcf; \
done
```

```
bcftools view --no-version -Ob -o blueberry.1000g_phased.other.GRCh38.bcf blueberry.unphased.GRCh38.vcf.gz \
  -t ^$(seq -s, 1 22),X,$(seq -f chr%.0f -s, 1 22),chrX && \
  bcftools index -f blueberry.1000g_phased.other.GRCh38.bcf
```

```
bcftools concat --no-version -Ob -o blueberry.1000g_phased.GRCh38.bcf blueberry.1000g_phased.{chr{{1..22},X},other}.GRCh38.bcf && \
  bcftools index -f blueberry.1000g_phased.GRCh38.bcf
```

Clean up
```
/bin/rm blueberry.1000g_phased.{chr{{1..22},X},other}.GRCh38.bcf{,.csi}
```

Remove switch errors using family genotype data
```
bcftools merge --no-version -Ou -m none blueberry.{1000g_phased,array}.GRCh38.bcf | \
  bcftools +trio-phase --no-version -Ou -- -p blueberry.ped | \
  bcftools view --no-version -Ob -o blueberry.trio_phased.GRCh38.bcf -I \
  -S <(bcftools query -l blueberry.1000g_phased.GRCh38.bcf) && \
  bcftools index -f blueberry.trio_phased.GRCh38.bcf
```

Allelic count of cfDNA SNP-based and MPSS data
==============================================

Currently version 4.0.12.0 of Mutect2 is the only version that counts fragments rather than counting reads
```
for sm in 2115760 20181127; do \
  for chr in chr{{1..22},X,Y,M}; do \
    gatk-4.0.12.0/gatk \
      Mutect2 \
      --max-reads-per-alignment-start 0 \
      --max-mnp-distance 0 \
      --allow-non-unique-kmers-in-ref \
      -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
      -I $sm.bam \
      --tumor $sm \
      -O $sm.$chr.GRCh38.vcf.gz \
      -L $chr \
      --genotyping-mode GENOTYPE_GIVEN_ALLELES \
      --alleles blueberry.unphased.GRCh38.vcf.gz; \
  done; \
done
```
Due to bugs in the assembly graph determination in GATK, allelic counts at variants next to each other are not trusted and will be filtered out later

```
for sm in 2115760 20181127; do \
  bcftools concat --no-version -Ou $sm.chr{{1..22},X,Y,M}.GRCh38.vcf.gz | \
    bcftools annotate --no-version -Ou -x ID,QUAL,INFO,^FMT/GT,^FMT/AD | \
    bcftools norm --no-version -Ou -m -any -k | \
    bcftools norm --no-version -Ob -o $sm.GRCh38.bcf -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna && \
    bcftools index -f $sm.GRCh38.bcf; \
done
```

Clean up
```
/bin/rm blueberry.{2115760,20181127}.chr{{1..22},X,Y,M}.GRCh38.vcf.gz
```

Merge phased genotype data and allelic counts in single VCF
```
bcftools merge --no-version -Ob -o blueberry.GRCh38.bcf -m none --force-samples \
  {blueberry.{trio,1000g}_phased,{2115760,20181127}}.GRCh38.bcf && \
  bcftools index -f blueberry.GRCh38.bcf
```

Run Blueberry algorithm
=======================

Run cfDNA SNP-based data data without maternal phasing
```
bcftools +blueberry \
  --rules GRCh38 \
  --apply-filters "PASS,." \
  --min-cov 200 \
  --min-dist 50 \
  --cross-prob 1e-5 \
  --mat-phase-err-prob .5 \
  blueberry.GRCh38.bcf \
  2115760 56001801068{81,75}4A
```

Run cfDNA SNP-based data with 1000 Genomes Project phasing
```
bcftools +blueberry \
  --rules GRCh38 \
  --apply-filters "PASS,." \
  --min-cov 200 \
  --min-dist 50 \
  --mat-switch-err-prob 5e-2 \
  --cross-prob 1e-5 \
  blueberry.GRCh38.bcf \
  2115760 2:56001801068{81,75}4A
```

Run cfDNA SNP-based data with family based maternal phasing
```
bcftools +blueberry \
  --rules GRCh38 \
  --apply-filters "PASS,." \
  --min-cov 200 \
  --min-dist 50 \
  --cross-prob 1e-5 \
  blueberry.GRCh38.bcf \
  2115760 56001801068{81,75}4A
```

Run cfDNA MPSS data without maternal phasing
```
bcftools +blueberry \
  --rules GRCh38 \
  --apply-filters "PASS,." \
  --exclude "SNP_DIST<=20 || INDEL_DIST<=20" \
  --no-indels \
  --no-unphased \
  --min-dist 100 \
  --cross-prob 1e-5 \
  --mat-phase-err-prob .5 \
  blueberry.GRCh38.bcf \
  20181127 56001801068{81,75}4A
```

Run cfDNA MPSS data with 1000 Genomes Project phasing
```
bcftools +blueberry \
  --rules GRCh38 \
  --apply-filters "PASS,." \
  --exclude "SNP_DIST<=20 || INDEL_DIST<=20" \
  --no-indels \
  --no-unphased \
  --min-dist 100 \
  --mat-switch-err-prob 5e-2 \
  --cross-prob 1e-5 \
  blueberry.GRCh38.bcf \
  20181127 2:56001801068{81,75}4A
```

Run cfDNA MPSS data with 1000 Genomes Project phasing maternal phasing
```
bcftools +blueberry \
  --rules GRCh38 \
  --apply-filters "PASS,." \
  --exclude "SNP_DIST<=20 || INDEL_DIST<=20" \
  --no-indels \
  --no-unphased \
  --min-dist 100 \
  --cross-prob 1e-5 \
  blueberry.GRCh38.bcf \
  20181127 56001801068{81,75}4A
```

Acknowledgements
================

This work was supported by NIH grant <a href="https://projectreporter.nih.gov/project_info_description.cfm?aid=8852155">R01 HG006855</a> and the Stanley Center for Psychiatric Research and by US Department of Defense Breast Cancer Research Breakthrough Award W81XWH-16-1-0316 (project BC151244)
