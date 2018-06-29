
# Get the data

GROSEQ_DIR = groseq
SMALL_RNA_THREEPRIME_ADAPTER_V1 = TCGTATGCCGTCTTCT
HG19_BOWTIE_INDEX = /nfs/genomes/human_gp_feb_09_no_random/bowtie/hg19

groseq_map: $(GROSEQ_DIR)/GROSEQ_A549_mapped.bam

# map
$(GROSEQ_DIR)/GROSEQ_A549_mapped.bam: $(GROSEQ_DIR)/GROSEQ_A549_more_trimmed.fastq
	$(BSUB) "python $(SPLIT_AND_MAP) -n 10000 --bsub-command '$(BSUB)' --infile $< -o $@ --merge-memory 120000 --mapping-command 'bowtie -S -p 24 --all --tryhard -v 2 $(HG19_BOWTIE_INDEX)'"

# further trim adapters
$(GROSEQ_DIR)/GROSEQ_A549_more_trimmed.fastq: $(GROSEQ_DIR)/GROSEQ_A549_trimmed.fastq
	$(BSUB) "fastx_trimmer -f 5 -Q 33 -m 15 -i $< -o tmp.fastq" # discard first 4 nucs
	$(BSUB) "fastx_trimmer -t 6 -Q 33 -m 15 -i tmp.fastq -m 15 -o $@" # discard last 6 nucs
	rm tmp.fastq

# trim adapters
$(GROSEQ_DIR)/GROSEQ_A549_trimmed.fastq: $(DATA_DIR)/GROSEQ_A549.fastq
	$(BSUB) "cutadapt --quality-base=64 -a $(SMALL_RNA_THREEPRIME_ADAPTER_V1) -m 15 $< > $@"

# Fastq-dump
$(DATA_DIR)/GROSEQ_A549.fastq: $(DATA_DIR)/GROSEQ_A549.sra
	$(BSUB) "$(FASTQ_DUMP) `pwd`/$<"

# Get GROSEQ data
$(eval $(call wget-target,\
$(DATA_DIR)/GROSEQ_A549.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX119%2FSRX119647/SRR408117/SRR408117.sra))
