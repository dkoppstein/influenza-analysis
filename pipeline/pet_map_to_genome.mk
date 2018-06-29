# OK, we're going to try mapping the influenza sequences directly to the genome
# in this iteration, and then find the reads that intersect with the 
# PET data. 

SHUFFLE_TARGETS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_shuffled.fastq)

shuffle: $(INT_DIR)/PB2_TS_shuffled.fastq

clean_shuffle:
	rm -f $(SHUFFLE_TARGETS)

####### SHUFFLE PRIME-AND-REALIGNED SEQUENCES #########
$(INT_DIR)/%_shuffled.fastq: $(INT_DIR)/%_par.fastq $(SHUFFLE_SCRIPT)
	$(BSUB) "python $(SHUFFLE_SCRIPT) --preserve-cpg -i $< > $@"

# For getting the stats of mapping targets, 
# 1. Find the aligned reads
# 2. Get the fourth column
# 3. Plot the distribution
$(INT_DIR)/%_distribution.txt: $(INT_DIR)/%.bam
	$(BSUB) "samtools view -F0x4 PB2_TS_par_mapped_to_gencode_tss_25_unique.bam | cut -f 1-4 | awk ' { printf \"%s\t%s\t%s\t%s\n\", \$1, \$2, \$3, (\$4 - 26) } ' > PB2_TS_par_mapped_to_gencode_tss_25_unique.txt"

RNAPET_TARGETS :=

# NUM_FASTA_FILES = 20
# FASTA_NUMS = $(shell seq 1 $(NUM_FASTA_FILES))

# Split de fasta file
# define fasta-split-pattern
# $(foreach num,$(FASTA_NUMS),$(INT_DIR)/$(1).$(num).fasta)
# endef

# Map prime-and-realign subtracted influenza sequencing to the new database. 

# First argument is the $(prefix).fastq in $(INT_DIR), 
# Second argument is the RNA PET file
# If the third argument is "unique", then only output uniquely mapping reads
define map-to-rnapet-unique

RNAPET_TARGETS += $(INT_DIR)/$(1)_mapped_to_$(2)_$(3)_unique_pet.sam

NUM_LINES_IN_SPLIT = 1000000
SPLIT_PREFIX = _split

# For example, intermediate/PB2_TS_par_mapped_to_GIS_RnaPet_A549_longPolyA_cell_rep1_3.bam
# We require samtools commit 44f25f37ceb6a053df14f53ac1369db7e707350d 
# because TAK version can't process huge headers
$(INT_DIR)/$(1)_mapped_to_$(2)_$(3)_unique_pet.sam: \
$(INT_DIR)/$(1).fastq \
$(call bowtie-target-pattern,$(2))
	$(BSUB) "split -l $(NUM_LINES_IN_SPLIT) $< $(INT_DIR)/$(1)$(SPLIT_PREFIX)"
$(foreach file,$(shell ls $(INT_DIR)/$(1)$(SPLIT_PREFIX)*),\
ifeq ($(3),y)
	$(BSUB) -R "rusage[mem=140000]" -M 140000 "bowtie -m 1 --all --norc -p `nproc` --solexa1.3-quals --tryhard -v 0 -S $(TSS_BOWTIE_DIR)/$(2) $(file) > $$@"
else
	$(BSUB) -R "rusage[mem=140000]" -M 140000 "bowtie --all --norc -p `nproc` --solexa1.3-quals --tryhard -v 0 -S $(TSS_BOWTIE_DIR)/$(2) $$(filter %.fastq,$$^) > $$@"
endif
endef

# For both PB2_TS and PB2_TS_shuffled, map uniquely and non-uniquely, 
# and to the TSS-restricted and unrestricted RNA PET reads. 
$(foreach sample,PB2_TS_par PB2_TS_shuffled,\
$(foreach rnapet,GIS_RnaPet_A549_longPolyA_cell_rep1 \
GIS_RnaPet_A549_longPolyA_cell_rep1_tss_intersect,\
$(foreach cond,y n,$(eval $(call map-to-rnapet-unique,$(sample),$(rnapet),$(cond))))))

rnapet_targets: $(RNAPET_TARGETS)
clean_rnapet_targets: 
	rm -f $(RNAPET_TARGETS)

# NUM_FASTA_FILES = 10
# FASTA_NUMS = $(shell seq 1 $(NUM_FASTA_FILES))

# # Split de fasta file
# define fasta-split-pattern
# $(foreach num,$(FASTA_NUMS),$(INT_DIR)/$(1).$(num).fasta)
# endef

# null  :=
# space := $(null) #
# comma := ,

# # Call bowtie-build on all de split fastas
# define bowtie-build-split-fasta
# $(call bowtie-target-pattern,%): $(call fasta-split-pattern,%_offset)
# 	mkdir -p $(TSS_BOWTIE_DIR)
# 	$(BSUB) "bowtie-build $(subst $(space),$(comma),$(strip $$^)) $(TSS_BOWTIE_DIR)/$$*"
# endef

# $(foreach num,$(FASTA_NUMS),$(eval $(call bowtie-build-split-fasta,$(num))))

# $(call fasta-split-pattern,%_offset): $(INT_DIR)/%_offset.fasta
# 	$(BSUB) "pyfasta split -n $(NUM_FASTA_FILES) $<"

$(call bowtie-target-pattern,%): $(INT_DIR)/%_offset.fasta
	$(BSUB) "bowtie-build $< $(TSS_BOWTIE_DIR)/$*"

# Get fasta file corresponding to the offset BED file. 
$(INT_DIR)/$(RNAPET_PREFIX)%_offset.fasta: \
$(INT_DIR)/$(RNAPET_PREFIX)%_offset.bed $(HG19_WHOLE_GENOME)
	$(BSUB) "bedtools getfasta -s -fi $(HG19_WHOLE_GENOME) -fo $@ -bed $<"

# Offset BED format, then sort and get only unique reads
# Sort trick from https://groups.google.com/forum/#!topic/bedtools-discuss/2o7oUgBwebw
$(INT_DIR)/$(RNAPET_PREFIX)%_offset.bed: \
$(INT_DIR)/$(RNAPET_PREFIX)%_mapped_aligned.bed $(OFFSET_BED)
	$(BSUB) "python $(OFFSET_BED) --five-prime-offset -20 --three-prime-offset 9 -i $< | sort -k1,1 -k2,2n -k3,3n -u > $@"

# Convert to bed file
$(INT_DIR)/$(RNAPET_PREFIX)%_mapped_aligned.bed: \
$(INT_DIR)/$(RNAPET_PREFIX)%_mapped_aligned.bam
	$(BSUB) "bamToBed -i $< > $@"

# In parallel, find those that intersect annotated transcription start sites
$(INT_DIR)/$(RNAPET_PREFIX)%_tss_intersect_mapped_aligned.bam: \
$(INT_DIR)/$(RNAPET_PREFIX)%_mapped_aligned.bam $(GENCODE_STARTS)
	$(BSUB) "bedtools intersect -abam $< -b $(GENCODE_STARTS) -wa -u -s > $@"

$(INT_DIR)/$(RNAPET_PREFIX)%_mapped_aligned.bam: \
$(INT_DIR)/$(RNAPET_PREFIX)%_mapped_all.bam
	$(BSUB) "samtools view -bh -F0x4 -@ `nproc` $< > $@"

$(INT_DIR)/$(RNAPET_PREFIX)%_mapped_all.bam: \
$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed_rc_nogs_collapsed.fastq $(BOWTIE_GENOME_DIR)
	$(BSUB) "bowtie -p `nproc` --tryhard --all -v 0 -S $(BOWTIE_HG19_GENOME) $< | samtools view -bS -@ `nproc` - > $@"

# Convert to FASTQ to give fake bases so that HTSeq.BAM_Reader doesn't crash
$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed_rc_nogs_collapsed.fastq: \
$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed_rc_nogs_collapsed.fasta
	$(BSUB) "fasta_to_fastq $< > $@"

$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed_rc_nogs_collapsed.fasta: \
$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed_rc_nogs.fastq
	$(BSUB) "fastx_collapser -i $< -o $@"

$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed_rc_nogs.fastq: \
$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed_rc.fastq
	$(BSUB) "cutadapt -g GGGGGGGGGG -e 0 -O 1 --minimum-length=20 $< > $@"

$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed_rc.fastq: \
$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed.fastq
	$(BSUB) "fastx_reverse_complement -i $< -o $@"

$(INT_DIR)/$(RNAPET_PREFIX)%_5prime_trimmed.fastq: \
$(INT_DIR)/$(RNAPET_PREFIX)%_5prime.fastq
	$(BSUB) "cutadapt --trimmed-only -a CTGCTGATGG --overlap=9 -f fastq --quality-base=64 --minimum-length=25 $< > $@"

# Find the reads corresponding to the 5' end of the transcript, as 
# described in Ruan and Ruan, Methods Mol. Biol., 2012. 
$(INT_DIR)/$(RNAPET_PREFIX)%_5prime.fastq: \
$(DATA_DIR)/$(RNAPET_PREFIX)%_1.fastq $(DATA_DIR)/$(RNAPET_PREFIX)%_2.fastq $(PROCESS_PET)
	$(BSUB) "python $(PROCESS_PET) -i $(filter %.fastq,$^) > $@"

# http://stackoverflow.com/questions/16273724/side-by-side-histograms-in-the-same-graph-in-r

# for file in `ls prefix_PB2*`; do bsub -R "rusage[mem=140000]" -M 140000 "bowtie --all --norc -p `nproc` --solexa1.3-quals --tryhard -v 0 -S ../anno/bowtie_indices/GIS_RnaPet_A549_longPolyA_cell_rep1 $file > $file_mapped_to_GIS_RnaPet_A549_longPolyA_cell_rep1_byhand.sam"; done;

# split -l 100000 PB2_TS_par.fastq PB2_TS_par_split_prefix.

# for file in `ls PB2_TS_par_split_prefix*`; do bsub -R "rusage[mem=200000]" -M 200000 "bowtie --all --norc -p `nproc` --solexa1.3-quals --tryhard -v 0 -S ../anno/bowtie_indices/GIS_RnaPet_A549_longPolyA_cell_rep1 $file > $file.mapped_to_GIS_RnaPet_A549_longPolyA_cell_rep1_byhand.sam"; done;
