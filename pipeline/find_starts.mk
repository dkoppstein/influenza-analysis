################## FIND START SITES #############################

# Find all the start sites from gtf file, then immediately 
# sort with bedtools
# % = ucsc, refseq, ensembl

# make the targets of bowtie-build with pattern matching on the fly

TSS_WINDOW_LOWER = 5
TSS_WINDOW_UPPER = 25
TSS_WINDOW_SIZE = 25
GENCODE_OFFSET_NUCS = 25
GENCODE_PATTERN = gencode_tss
GENCODE_TSS = $(ANNO_DIR)/$(GENCODE_PATTERN).gtf
GENCODE_STARTS = $(ANNO_DIR)/$(GENCODE_PATTERN)_$(TSS_WINDOW_SIZE).gtf
GENCODE_STARTS_0 := $(ANNO_DIR)/$(GENCODE_PATTERN)_0.gtf

gencode_tss: $(GENCODE_STARTS) $(GENCODE_STARTS_0)

clean_gencode_tss: 
	rm -f $(GENCODE_STARTS) $(GENCODE_STARTS_0)

########## EXTRACT SEQUENCES FROM TSS WINDOW ANNOTATIONS #######
# % = ucsc, refseq, ensembl

# Need to convert to uppercase because fastx_collapser can't handle it
$(ANNO_DIR)/gencode_before_g_strip.fastq: $(ANNO_DIR)/gencode_offset.fasta $(UPPERCASE_FASTQ)
	$(BSUB) "fasta_to_fastq -Q 33 $< | python $(UPPERCASE_FASTQ) > $@"

# Build a bowtie index from the merged file
$(call bowtie-target-pattern,gencode_offset_25_nogs_merged): $(ANNO_DIR)/gencode_offset_25_nogs_merged.fasta
	mkdir -p $(BOWTIE_DIR)
	$(BSUB) "bowtie-build $< $(BOWTIE_DIR)/gencode_offset_25_nogs_merged"

# Get fasta from the merged file
$(ANNO_DIR)/gencode_offset_25_nogs_merged.fasta: $(ANNO_DIR)/gencode_offset_25_nogs_merged.bed
	$(BSUB) "bedtools getfasta -s -fi $(HG19_WHOLE_GENOME) -fo $@ -bed $<"

# another stab at the background distribution
$(ANNO_DIR)/gencode_offset_25_nogs_summary.txt.gz: $(ANNO_DIR)/gencode_offset_3p_25_nogs.fasta $(FASTQ_TO_TABLE)
	$(BSUB) "$(MY_PYTHON) $(FASTQ_TO_TABLE) -i $< --summarize -f fasta | gzip -9 > $@"

# Get fasta from unmerged file
$(ANNO_DIR)/gencode_offset_25_nogs.fasta: $(ANNO_DIR)/gencode_offset_25_nogs.bed
	$(BSUB) "bedtools getfasta -s -fi $(HG19_WHOLE_GENOME) -fo $@ -bed $<"

# Merge the files
$(ANNO_DIR)/gencode_offset_25_nogs_merged.bed: \
$(ANNO_DIR)/gencode_offset_25_nogs.bed
	$(BSUB) "sortBed -i $< > tmp.bed"
	$(BSUB) "bedtools merge -s -c 4,5,6 -o distinct -i tmp.bed > $@"
	rm tmp.bed

# or this
$(ANNO_DIR)/gencode_offset_3p_20_nogs_5ptrim_5_summary.txt.gz: $(ANNO_DIR)/gencode_offset_3p_25_nogs_summary.txt.gz $(TRIM_SEQUENCES)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(TRIM_SEQUENCES) -l 5 | gzip -9 > $@"

# or this
$(ANNO_DIR)/gencode_offset_3p_25_nogs_summary.txt.gz: $(ANNO_DIR)/gencode_offset_3p_25_nogs.fasta $(FASTQ_TO_TABLE)
	$(BSUB) "$(MY_PYTHON) $(FASTQ_TO_TABLE) -i $< -f fasta --summarize | gzip -9 > $@"

# use this for background distribution
$(ANNO_DIR)/gencode_offset_3p_25_nogs.fasta: $(ANNO_DIR)/gencode_offset_3p_25_nogs_renamed.bed
	$(BSUB) "bedtools getfasta -s -fi $(HG19_WHOLE_GENOME) -fo $@ -bed $< -name"

# convert name to transcript_id
$(ANNO_DIR)/gencode_offset_3p_25_nogs_renamed.bed: $(ANNO_DIR)/gencode_offset_3p_25_nogs.bed $(CONVERT_GTF_NAMES)
	$(BSUB) "$(MY_PYTHON) $(CONVERT_GTF_NAMES) -i $< -o $@ --to transcript_id"

$(ANNO_DIR)/gencode_offset_3p_25_nogs.bed: $(ANNO_DIR)/gencode_offset_25_nogs_intermediate.bed
	$(BSUB) "$(MY_PYTHON) $(OFFSET_BED) -i $< --offset-three-prime-wrt-five-prime 25 > $@"

# Recenter around the new start coordinate
$(ANNO_DIR)/gencode_offset_25_nogs.bed: $(ANNO_DIR)/gencode_offset_25_nogs_intermediate.bed
	$(BSUB) "python $(OFFSET_BED) -i $< --offset-three-prime-wrt-five-prime 25 > gencode_tmp.bed"
	$(BSUB) "python $(OFFSET_BED) -i gencode_tmp.bed --five-prime-offset -25 > $@"
	rm gencode_tmp.bed

$(ANNO_DIR)/gencode_offset_25_nogs_intermediate.bed: $(ANNO_DIR)/gencode_3p_offset_25.bed $(BED_WITH_NO_GS) $(HG19_WHOLE_GENOME)
	$(BSUB) "$(MY_PYTHON) $(BED_WITH_NO_GS) -i $< -f $(HG19_WHOLE_GENOME) > $@"

# Get 1024 base pairs to ensure the sequence is unique, so you're targetting the 
# correct TSS when you map back
$(ANNO_DIR)/gencode_3p_offset_25.bed: $(ANNO_DIR)/gencode_tss_0.bed $(OFFSET_BED)
	$(BSUB) "python $(OFFSET_BED) -i $< --offset-three-prime-wrt-five-prime 25 > $@"

$(ANNO_DIR)/gencode_tss_0.bed: $(ANNO_DIR)/gencode_tss_0.gtf
	$(BSUB) "gtf2bed < $< > $@"

define find-starts
$(ANNO_DIR)/%_tss_$(1).gtf: $(ANNO_DIR)/%.gtf $(TSS_SCRIPT)
	$(BSUB) "python $(TSS_SCRIPT) -n $(1) $$< | bedtools sort -i stdin > $$@"
# Ensembl does things a little differently, so we need to add "chr" to each
ifeq ($*,ensembl)
	$(BSUB) "sed -e 's/^/chr/' $$@ > $$@.tmp"
	$(BSUB) "mv -f $$@.tmp $$@"
endif
endef

# Use a window size of 0 or 25 for the GENCODE gtf
$(eval $(call find-starts,$(TSS_WINDOW_SIZE)))
$(eval $(call find-starts,0))
