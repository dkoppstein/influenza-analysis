FIVEPRIME_TARGETS = $(call bowtie-target-pattern,gencode_closest_to_tss)

####### SHUFFLE PRIME-AND-REALIGNED SEQUENCES ########
SHUFFLE_TARGETS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_shuffled_collapsed.fasta)

shuffle: $(SHUFFLE_TARGETS)

clean_shuffle:
	rm -f $(SHUFFLE_TARGETS)

$(INT_DIR)/%_before_mapping.fastq: $(INT_DIR)/%_before_mapping.fasta
	$(BSUB) "fasta_to_fastq $< > $@"

$(INT_DIR)/%_shuffled_before_mapping.fasta: $(INT_DIR)/%_collapsed_before_mapping.fasta $(SHUFFLE_SCRIPT)
	mkdir -p $(INT_DIR)/$*_shuffled_before_mapping
	ln -s `pwd`/$< $(INT_DIR)/$*_shuffled_before_mapping/
	cd $(INT_DIR)/$*_shuffled_before_mapping && $(BSUB) pyfasta split -n 100 $(notdir $<) && rm $(notdir $<) && cd ../..
	$(BSUB) parallel "$(BSUB) \"python $(SHUFFLE_SCRIPT) --preserve-cpgs --not-same --no-gs -f fasta -i {} > {.}_shuffled.fasta\"" ::: $(INT_DIR)/$*_shuffled_before_mapping/*.fasta
	$(BSUB) "cat $(INT_DIR)/$*_shuffled_before_mapping/*_shuffled.fasta > $@"
	rm -f $(INT_DIR)/$*_shuffled_before_mapping/* # && rmdir $(INT_DIR)/$*_shuffled_before_mapping

# Make a graph of the actual library vs. the shuffled library

define graph-of-five-prime
$(GRAPH_DIR)/%_mapped_to_$(1).pdf: \
$(INT_DIR)/%_mapped_to_$(1)_aggregated_mapping.bed $(HISTOGRAM_COLUMN)
	mkdir -p $$(dir $$@)
	$(BSUB) "cat $(INT_DIR)/$$*_split_mapped_to_$(1)/split*.txt > $$*_mapped_to_$(1)_aggregated.txt"
	$(BSUB) "Rscript $(HISTOGRAM_COLUMN) --column 3 -i $$*_mapped_to_$(1)_aggregated.txt -o $$@"
endef

# Map reads to bowtie indices

define get-closest-genes
$(INT_DIR)/%_mapped_to_$(1)_closest_genes.bed: \
$(INT_DIR)/%_mapped_to_$(1)_aggregated_mapping.bed $(GENCODE_STARTS_0)
	$(BSUB) "bedtools intersect -a $$< -b $(GENCODE_STARTS_0) -s -wo > $$@"
endef


# $(1) is the RNA-seq library to map influenza reads onto (can be an RNA-PET or CAGE dataset)
define map-to-five-prime
$(INT_DIR)/%_mapped_to_$(1)_aggregated_mapping.bed: \
$(INT_DIR)/%_before_mapping.fastq \
$(call bowtie-target-pattern,$(1)_closest_to_tss) $(MAPPING_TO_BED)
	if [ -d $(INT_DIR)/$$*_split_mapped_to_$(1) ]; then rm -f $(INT_DIR)/$$*_split_mapped_to_$(1)/* && rmdir $(INT_DIR)/$$*_split_mapped_to_$(1); fi;
	mkdir -p $(INT_DIR)/$$*_split_mapped_to_$(1)
	$(BSUB) "split -a 4 -l 1000 $$< $(INT_DIR)/$$*_split_mapped_to_$(1)/split"
	for file in `ls $(INT_DIR)/$$*_split_mapped_to_$(1)/split*`; do mv $$$$file $$$$file.fastq; done;
# Parallelize mapping onto the read dataset
# grep -v ^@ removes the header
# cut  -f 1,3,4 to get a subset of the data
# grep -P -v ...: Remove reads that don't map (end with \t0$)
	$(BSUB) parallel "$(BSUB) -R rusage[mem=110000] -M 110000 \"bowtie --all --norc -p `nproc` --solexa1.3-quals --tryhard -v 0 -S $(BOWTIE_DIR)/$(1)_closest_to_tss {} | grep -v ^@ | cut -f 1,3,4 | grep -P -v \\\"\t0\\$$$$\\\" > {.}.txt\"" ::: $(INT_DIR)/$$*_split_mapped_to_$(1)/split*.fastq
	$(BSUB) parallel "$(BSUB) \"python $(MAPPING_TO_BED) -i {} -o {.}.bed --filter-position 6\"" ::: $(INT_DIR)/$$*_split_mapped_to_$(1)/split*.txt
	$(BSUB) "cat $(INT_DIR)/$$*_split_mapped_to_$(1)/split*.bed >> $$@"
endef

FIVE_PRIME_ANNOTATIONS = gencode $(CAGE_PREFIX)_cell_rep1 $(RNAPET_PREFIX)_cell_rep1

# Call map-to-five-prime for gencode, CAGE, and RNA-PET
$(foreach five_prime_annotation,$(FIVE_PRIME_ANNOTATIONS),$(eval $(call map-to-five-prime,$(five_prime_annotation))))

$(foreach five_prime_annotation,$(FIVE_PRIME_ANNOTATIONS),$(eval $(call graph-of-five-prime,$(five_prime_annotation))))

$(foreach five_prime_annotation,$(FIVE_PRIME_ANNOTATIONS),$(eval $(call get-closest-genes,$(five_prime_annotation))))

# Influenza sequence mapped to five prime targets make rule
FIVE_PRIME_TARGETS = \
$(foreach five_prime_annotation,$(FIVE_PRIME_ANNOTATIONS),\
$(foreach bn,$(TARGET_BNS),\
$(INT_DIR)/$(bn)_mapped_to_$(five_prime_annotation)_aggregated_mapping.bed))

five_prime_targets: $(FIVE_PRIME_TARGETS)

clean_five_prime_targets:
	rm -f $(FIVE_PRIME_TARGETS)

# Create the bowtie indices for the RNA-PET targets
# Do rep1 for now; incorporate both reps later
RNAPET_TARGETS = $(call bowtie-target-pattern,$(RNAPET_PREFIX)_cell_rep1_closest_to_tss)

rnapet_targets: $(RNAPET_TARGETS)
clean_rnapet_targets:
	rm -f $(RNAPET_TARGETS)

$(call bowtie-target-pattern,%_closest_to_tss): $(ANNO_DIR)/%_closest_to_tss.fasta
	$(BSUB) "bowtie-build $< $(BOWTIE_DIR)/$*_closest_to_tss"

# Get fasta file corresponding to the offset BED file. 
$(ANNO_DIR)/%_closest_to_tss.fasta: \
$(ANNO_DIR)/%_offset_threeprime_wrt_fiveprime.bed $(HG19_WHOLE_GENOME)
	$(BSUB) "bedtools getfasta -s -fi $(HG19_WHOLE_GENOME) -fo $@ -bed $<"

# Offset three prime with respect to five prime
$(ANNO_DIR)/%_offset_threeprime_wrt_fiveprime.bed: \
$(ANNO_DIR)/%_offset_fiveprime.bed $(OFFSET_BED)
	$(BSUB) "python $(OFFSET_BED) --offset-three-prime-wrt-five-prime 50 -i $< > $@"

# Offset BED format
$(ANNO_DIR)/%_offset_fiveprime.bed: \
$(ANNO_DIR)/%_closest_to_tss.bed $(OFFSET_BED)
	$(BSUB) "python $(OFFSET_BED) --five-prime-offset -25 -i $< > $@"

# For each read, use only the one that maps closest to a TSS
$(ANNO_DIR)/%_closest_to_tss.bed: \
$(ANNO_DIR)/%_mapped_aligned.bam $(GENCODE_STARTS_0) $(FIND_CLOSEST_MAPPING)
	$(BSUB) "python $(FIND_CLOSEST_MAPPING) -b $< -g $(GENCODE_STARTS_0) > $@"

# Filter for the aligned reads
$(ANNO_DIR)/%_mapped_aligned.bam: \
$(ANNO_DIR)/%_mapped_all.bam
	$(BSUB) "samtools view -bh -F0x4 -@ `nproc` $< > $@"

# Map to the genome, allow no mismatches
$(ANNO_DIR)/%_mapped_all.bam: \
$(ANNO_DIR)/%_collapsed_before_map.fastq $(BOWTIE_GENOME_DIR)
	$(BSUB) "bowtie -p `nproc` --tryhard --all -v 0 -S $(BOWTIE_HG19_GENOME) $< | samtools view -bS -@ `nproc` - > $@"

# Convert to FASTQ to give fake bases so that HTSeq.BAM_Reader doesn't crash :/
$(ANNO_DIR)/%_collapsed_before_map.fastq: $(ANNO_DIR)/%_collapsed_before_map.fasta
	$(BSUB) "fasta_to_fastq $< > $@"

# Collapse identical reads
$(ANNO_DIR)/%_collapsed_before_map.fasta: \
$(ANNO_DIR)/%_no_gs.fastq
	$(BSUB) "fastx_collapser -Q 33 -i $< -o $@"

# Remove 5' Gs to match template-switching data
$(ANNO_DIR)/%_no_gs.fastq: \
$(ANNO_DIR)/%_before_g_strip.fastq
	$(BSUB) "cutadapt -g GGGGGGGGGGGGGG -e 0 -O 1 --minimum-length=15 $< > $@"
