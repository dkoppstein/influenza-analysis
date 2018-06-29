# In supplemental material of Djebali et al., 2012
# "The 5’ linker upper oligonucleotide (N6: 5' ‐  CCACCGACAGGTTCAGAGTTCTA CAGXXXCAGCAGNNNNNN 
# Phos  ‐ 3', 
# GN5: 5' ‐  CCACCGACAGGTTCAGAG  TTCTACAGXXXCAGCAGGNNNNN Phos  ‐ 3') 
# was mixed in a 1:1 ratio with lower oligonucleotide 
# (5' ‐ Phos CTGCTG XXXCTGTAGAACTCTGAACCTGTCGGTGG ‐ 3')
# RC: CCACCGACAGGTTCAGAGTTCTACAGXXXCAGCAG

# The 3 ́ linker ligation was performed using 10 to 100ng of 3 ́ linker (Upper:
# 5' ‐ Phos NNTCGTATGCCGTCTTCTGCTTG ‐ 3',
# Lower: 5' ‐ CAAGCAGAAGACGGCATACGA ‐ 3')

# ^ means the adaptor has to be anchored at the 5' end in the read
# FIVE_PRIME_CAGE_ADAPTOR = ACAGXXXCAGCAG
FIVE_PRIME_CAGE_ADAPTOR = CAGCAG
THREE_PRIME_CAGE_ADAPTOR = TCGTATGCCGTCTTCT

# Discard the five prime adaptor
# After this, the RNA-PET and CAGE data processing is identical -> go to five_prime_library.mk
$(ANNO_DIR)/$(CAGE_PREFIX)_%_before_g_strip.fastq: $(ANNO_DIR)/$(CAGE_PREFIX)_%_5prime_filtered.fastq
	$(BSUB) "cutadapt -a $(THREE_PRIME_CAGE_ADAPTOR) --minimum-length=14 --quality-base=33 -f fastq --trimmed-only $< > $@"

# Discard the three-prime adaptor
$(ANNO_DIR)/$(CAGE_PREFIX)_%_5prime_filtered.fastq: $(DATA_DIR)/$(CAGE_PREFIX)_%.fastq
	$(BSUB) "cutadapt -g $(FIVE_PRIME_CAGE_ADAPTOR) -e 0 -f fastq --trimmed-only $< > $@"

################## FIND START SITES #############################

# Find all the start sites from gtf file, then immediately 
# sort with bedtools
# % = ucsc, refseq, ensembl

# fastq_to_fasta_fast $< > $@
# bowtie-build $< > bowtie/whatever
# bowtie --all --norc --solexa1.3-quals --tryhard -v 0 -S $(CAGE_TSS) | samtools view -bS - > $@
# bowtie --all --norc --solexa1.3-quals --tryhard -v 0 -S bowtie/SRR530909_collapsed ../influenza/daybefore/intermediate/NS1_TS_par.fastq | samtools view -bS - > NS1_TS_to_SRR530909.bam
# cutadapt -g ^GATCAGCAG -f fasta --trimmed-only SRR530909_collapsed.fasta > SRR530909_5filtered.fasta
# python ~/science/software/count_transcripts/count_transcripts.py -n 0 -g --no-collapse -f fasta SRR530909_filtered.fasta > SRR530909_nogs.fasta
# python ~/science/software/count_transcripts/count_transcripts.py -n 0 -g --no-collapse -f fasta SRR530909_filtered.fasta > SRR530909_nogs.fasta

# GATCAGCAG on 5' end
# TCGTATGCCGTCTTCT on 3' end
# TCGTATGCCGTCTTC on 3' end
# cutadapt -a TCGTATGCCGTCTTC -f fasta --trimmed-only SRR530909_5filtered.fasta > SRR530909_53filtered.fasta

# $cmd = "cd $home/STAR_alignments; STAR --genomeDir /lab/solexa_bartel/agarwal/databases/STAR --readFilesIn $file1 $file2 ".
#       "--runThreadN 30 --outFilterType BySJout --outFilterMultimapNmax 20 --outSJfilterOverhangMin -1 12 12 12 ".
#       "--outSJfilterCountUniqueMin -1 5 10 10 --outSJfilterCountTotalMin -1 50 100 100 --outSJfilterIntronMaxVsReadN 5000 10000 15000 ".
#       "--scoreGapNoncan -100 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --readMatesLengthsIn Equal --readFilesCommand $readFilesCommand ". #
#       "--outFileNamePrefix $prefix --outStd SAM | samtools view -bS -@ 30 - | samtools sort -@ 30 - $prefix; samtools rmdup -S $prefix.bam $prefix.filtered.bam; ".
#       "mv $prefix.bam bamfiles/; mv $prefix.filtered.bam bamfiles_filtered/"; #--scoreDelOpen -100 --scoreInsOpen -100

# STAR --genomeDir /lab/solexa_bartel/agarwal/databases/STAR_hg19 --readFilesIn 5prime_sequences_trimmed_rc_nogscutadapt_collapsed.fasta --runThreadN `nproc` --outStd SAM | samtools view -bS -@ `nproc` - > SRR575143_mapped.bam

# --outSJfilterCountTotal 

# Then collect only the aligned reads
# samtools view -bh -F0x4 $< > $@
# and spit them out into a fasta file
# samtools view -bh -F0x4 SRR575143_mapped.bam > SRR575143_aligned.bam
# samtools view SRR575143_aligned.bam | awk '{OFS="\t"; print ">"$1"\n"$10}' - > SRR575143_reads.fasta
# Compress them:
# fastx_collapser -i SRR575143_reads.fasta -o SRR575143_reads_collapsed.fasta
# bsub "samtools view SRR575143_aligned.bam | awk '{OFS=\"\t\"; print \">\"\$1\"\n\"\$10}' - > SRR575143_reads.fasta"
# Make a bowtie index: 
# bowtie-build SRR575143_reads.fasta bowtie/SRR575143_bowtie
# bowtie --all --norc --solexa1.3-quals --tryhard -v 0 -S bowtie/SRR575143_bowtie /lab/solexa_bartel/koppstein/influenza/daybefore/intermediate/NS1_30min_par.fastq | samtools view -bS -@ `nproc` > NS1_30min_par_mapped_to_SRR575143.bam
# bowtie --all --norc --solexa1.3-quals --tryhard -v 0 -S -m 1 -p `nproc` bowtie/SRR575143_bowtie /lab/solexa_bartel/koppstein/influenza/daybefore/intermediate/NS1_30min_par.fastq | samtools view -bS -@ `nproc` - > NS1_30min_par_mapped_to_SRR575143_unique.bam

# bowtie --all --norc --solexa1.3-quals --tryhard -v 0 -S -p `nproc` bowtie/SRR575143_aligned_to_gencode_bowtie /lab/solexa_bartel/koppstein/influenza/daybefore/intermediate/NS1_30min_par.fastq | samtools view -bS -@ `nproc` - > NS1_30min_par_mapped_to_SRR575143_gencode_filtered.bam

# bowtie --all --norc --solexa1.3-quals --tryhard -v 0 -S -m 1 -p `nproc` bowtie/SRR575143_aligned_to_gencode_bowtie /lab/solexa_bartel/koppstein/influenza/daybefore/intermediate/NS1_30min_par.fastq | samtools view -bS -@ `nproc` - > NS1_30min_par_mapped_to_SRR575143_gencode_filtered_unique.bam
# awk ' {if ($4 == 1) print } ' test.bam
# R 
# read.table(infile, sep="\t", header=FALSE, comment.char="")
#  reheader 	samtools reheader <in.header.sam> <in.bam> 
# bsub "cat NS1_30min_par_mapped_to_SRR575143.sam_header > NS1_30min_par_mapped_to_SRR575143_first_only.sam && awk ' { if (\$4 == 1) print } ' NS1_30min_par_mapped_to_SRR575143.sam >> NS1_30min_par_mapped_to_SRR575143_first_only.sam"


# Unique hits:
# bsub "bowtie --all -m 1 -p `nproc` --norc --solexa1.3-quals --tryhard -v 0 -S /lab/bartel1_ata/koppstein/Tak_Only/5prime_annotations/bowtie_indices/gencode_tss_25 /lab/solexa_bartel/koppstein/influenza/daybefore/intermediate/PB2_TS_par.fastq | samtools view -bS -@ `nproc` - > PB2_TS_par_mapped_to_gencode_tss_25_unique.bam"
# bsub "samtools view -F0x4 PB2_TS_par_mapped_to_gencode_tss_25_unique.bam | cut -f 1-4 | awk ' { printf \"%s\t%s\t%s\t%s\n\", \$1, \$2, \$3, (\$4 - 26) } ' > PB2_TS_par_mapped_to_gencode_tss_25_unique.txt"
# bsub "samtools view -@ `nproc` -F0x4 NS1_30min_par_mapped_to_SRR575143_gencode_filtered.bam | cut -f 1-4 | awk ' { printf \"%s\t%s\t%s\t%s\n\", \$1, \$2, \$3, (\$4 - 26) } ' > NS1_30min_par_mapped_to_SRR575143_gencode_filtered_stats.txt"
# Rscript /lab/bartel1_ata/koppstein/science/software/dklib/influenza/throwaway.R 

# Bedtools: bedtools intersect -abam SRR575143_aligned.bam -b /lab/bartel1_ata/koppstein/Tak_Only/5prime_annotations/gencode_tss_50.gtf -wa -u -s > SRR575143_aligned_gencode_tss_50.bam

# /lab/bartel1_ata/koppstein/virtualenvs/my2.7/bin/pyfasta split -n 10 GIS_RnaPet_A549_longPolyA_cell_rep1_tss_intersect_offset.fasta

# outFilterIntronMotifs None outFilterMismatchNmax 0
