# Process ribo-zero RNA-seq samples

RNASEQ_TARGETS = $(foreach bn,$(RNASEQ_BNS),$(INT_DIR)/$(bn)_influenza_gene_counts.sam)

INFLUENZA_RNASEQ_TARGETS = $(foreach bn,$(RNASEQ_BNS),$(INT_DIR)/$(bn)_influenza_gene_counts_stdout.txt $(INT_DIR)/$(bn)_hg19_gene_counts_stdout.txt)

RNASEQ_MAPPED_TARGETS = $(foreach bn,$(RNASEQ_BNS),$(INT_DIR)/$(bn)_mapped_to_hg19.bam)

rnaseq_mapped: $(RNASEQ_MAPPED_TARGETS)
clean_rnaseq_mapped:
	rm -f $(RNASEQ_MAPPED_TARGETS)

influenza_rnaseq: $(INFLUENZA_RNASEQ_TARGETS)
clean_influenza_rnaseq: 
	rm -f $(INFLUENZA_RNASEQ_TARGETS)

clean_tophat:
	rm -rf $(foreach bn,$(RNASEQ_BNS),$(INT_DIR)/$(bn)_tophat_tmp)

rnaseq: $(RNASEQ_TARGETS)
clean_rnaseq: 
	rm -f $(RNASEQ_TARGETS)

###### QUANTIFY GENOME MAPPING ########

CUFFLINKS_TARGETS = $(foreach bn,$(RNASEQ_BNS),$(INT_DIR)/$(bn)_cufflinks_output/genes.fpkm_tracking)

cufflinks: $(CUFFLINKS_TARGETS)
clean_cufflinks:
	rm -f $(CUFFLINKS_TARGETS)

QUANTIFY_MAPPING_TARGETS = $(foreach bn,$(RNASEQ_BNS),$(INT_DIR)/$(bn)_htseq_gene_counts.sam)
QUANTIFY_MAPPING_TARGETS_TEXTFILES = $(foreach bn,$(RNASEQ_BNS),$(INT_DIR)/$(bn)_htseq_gene_counts.txt)

quantify_mapping: $(QUANTIFY_MAPPING_TARGETS)
clean_quantify_mapping:
	rm -f $(QUANTIFY_MAPPING_TARGETS) $(QUANTIFY_MAPPING_TARGETS_TEXTFILES)

$(INT_DIR)/%_htseq_gene_counts.sam: | \
$(INT_DIR)/%_mapped_to_hg19.bam $(GENCODE_ANNO)
	$(BSUB) "samtools view $(word 1,$|) | htseq-count -m intersection-strict -o $@ - $(GENCODE_ANNO) > $(INT_DIR)/$*_htseq_gene_counts.txt"

# use for nucleotide distribution
$(INT_DIR)/%_top_genes_summary.txt.gz: $(INT_DIR)/%_top_genes.fasta.gz $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) --summarize -f fasta | gzip -9 > $@"


$(INT_DIR)/%_top_genes.fasta.gz: $(INT_DIR)/%_cufflinks_output/isoforms.fpkm_tracking $(ANNO_DIR)/gencode_offset_3p_25_nogs.fasta $(TOP_GENES)
	$(BSUB) "$(MY_PYTHON) $(TOP_GENES) -c $< -f $(word 2,$^) --num-genes 13000 | gzip -9 > $@"

# Alternatively, use cufflinks
$(INT_DIR)/%_cufflinks_output/isoforms.fpkm_tracking: $(INT_DIR)/%_mapped_to_hg19.bam $(GENCODE_ANNO)
	mkdir -p $(INT_DIR)/$*_cufflinks_output
	$(BSUB) "cufflinks -G $(GENCODE_ANNO) $< -o $(INT_DIR)/$*_cufflinks_output"

###### MAP TO GENOME #########

MAP_TO_GENOME_TARGETS = $(foreach bn,$(RNASEQ_BNS),$(INT_DIR)/$(bn)_mapped_to_hg19.bam)

map_to_genome: $(MAP_TO_GENOME_TARGETS)
clean_map_to_genome:
	rm -f $(MAP_TO_GENOME_TARGETS)

GENCODE_ANNO_PREBUILT_INDEX := anno/gencode/index

# 
$(GENCODE_ANNO_PREBUILT_INDEX): $(BOWTIE_GENOME_DIR) $(GENCODE_ANNO)
	$(BSUB) "tophat -G $(GENCODE_ANNO) --transcriptome-index=$@ $(BOWTIE_HG19_GENOME)"

# http://barcwiki.wi.mit.edu/wiki/SOPs/mapping
# uses tophat 2 (default in PATH)
$(INT_DIR)/%_mapped_to_hg19.bam: \
$(INT_DIR)/%_no_threeprime_adaptor.fastq.gz $(GENCODE_ANNO_PREBUILT_INDEX)
	$(BSUB) "zcat $< > $(basename $<) && tophat -o $(INT_DIR)/$*_tophat_tmp --phred64-quals --no-novel-juncs --transcriptome-index=$(GENCODE_ANNO_PREBUILT_INDEX) $(BOWTIE_HG19_GENOME) $(basename $<)"
	mv $(INT_DIR)/$*_tophat_tmp/accepted_hits.bam $@
	rm -rf $(INT_DIR)/$*_tophat_tmp
	rm $(basename $<)

TRUSEQ_THREEPRIME_ADAPTOR = AGATCGGAAGAGC

INFLUENZA_GTF = $(ANNO_DIR)/influenza_sequences.gtf

# Use HTSeq-Count to get the counts per gene
$(INT_DIR)/%_influenza_gene_counts_stdout.txt: $(INT_DIR)/%_mapped_to_influenza.bam $(INFLUENZA_GTF)
	$(BSUB) "samtools view -h $< | htseq-count - $(INFLUENZA_GTF) -i ID > $@"

# Make a dummy GTF that is the same length as the influenza sequences, for htseq-count
$(INFLUENZA_GTF): $(INFLUENZA_SEQUENCES) $(FASTX_TO_GTF)
	$(BSUB)	"python $(FASTX_TO_GTF) -i $< -f fasta > $@"

INFLUENZA_INDEX = $(BOWTIE_DIR)/influenza_sequences

# Map to influenza genome
$(INT_DIR)/%_mapped_to_influenza.bam: \
$(INT_DIR)/%_no_threeprime_adaptor.fastq $(call bowtie-target-pattern,influenza_sequences)
	$(BSUB) "bowtie --solexa1.3-quals -p `nproc` -v 1 --norc -S $(INFLUENZA_INDEX) $< | samtools view -bS -@ `nproc` - > $@"

$(call bowtie-target-pattern,influenza_sequences): $(INFLUENZA_SEQUENCES)
	$(BSUB) "bowtie-build $< $(BOWTIE_DIR)/influenza_sequences"

# Look at human gene expression
$(INT_DIR)/%_hg19_gene_counts_stdout.txt: $(INT_DIR)/%_mapped_to_hg19.bam
	$(BSUB) "samtools view -h $< | htseq-count - $(GENCODE_ANNO) -o $*_hg19_gene_counts.sam -i gene_id > $@"

# Map to human genome
#$(INT_DIR)/%_mapped_to_hg19.bam: $(INT_DIR)/%_no_threeprime_adaptor.fastq $(HG19_WHOLE_GENOME)
#	$(BSUB) "bowtie --solexa1.3-quals -p `nproc` -v 1 -m 100 --best --strata -S $(BOWTIE_HG19_GENOME) $< | samtools view -bS -@ `nproc` - > $@"

# Remove the five prime end of the five prime adaptor
$(INT_DIR)/%_no_threeprime_adaptor.fastq.gz: $(DATA_DIR)/%_data.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a $(TRUSEQ_THREEPRIME_ADAPTOR) -O 4 -f fastq --minimum-length=27 - | gzip -9 > $@"
