# This file collapses the fastq data file using the barcode 

##### MAKE TARGETS #########

# Remember to declare targets *before* making the general rule! 
COLLAPSE_TARGETS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_collapsed.fastq)

# This is declared as .PHONY in common.mk
collapsed: $(COLLAPSE_TARGETS)

clean_collapsed:
	rm -f $(COLLAPSE_TARGETS)

##### MAKE RULES ###########



# Specialized rule for dummy files

define five-g-trimmed-rule
$(INT_DIR)/$(1)_5GTRIMMED_collapsed.fastq.gz: $(DATA_DIR)/$(1)_data.fastq.gz
	$(BSUB) "zcat $$< | fastx_trimmer -f $(2) | cutadapt --match-read-wildcards -g GGGGGGGGGGGGGGGGGG -e 0 -O 1 - | gzip -9 > $$@"
endef

$(foreach bn,$(GEN2_TS_BNS) $(GEN2_TIME_COURSE_BNS),$(eval $(call five-g-trimmed-rule,$(bn),13)))

$(foreach bn,$(GEN1_TS_BNS),$(eval $(call five-g-trimmed-rule,$(bn),9)))

COUNT_PARAMS := 

# General rule for collapsing transcripts by barcode
$(INT_DIR)/%_collapsed.fastq.gz: $(DATA_DIR)/%_data.fastq.gz $(COUNT_TRANSCRIPTS)
	mkdir -p $(INT_DIR)
	$(BSUB) "zcat $< > $(basename $<) && python $(COUNT_TRANSCRIPTS) $(COUNT_PARAMS) $(basename $<) | gzip -9 > $@"
	rm $(basename $<)

TIME_COURSE_BNS := 

define collapse-rule
$(INT_DIR)/$(1)_collapsed.fastq.gz: COUNT_PARAMS := $(2)
endef

# Set of parameters for the second generation template-switching experiments
# Since it's second generation, 12 nucleotides are used for the barcode
$(foreach bn,$(GEN2_TIME_COURSE_BNS) $(GEN2_TS_BNS),\
$(eval $(call collapse-rule,$(bn),\
--nucs 12 --remove-guanines --remove-ns --ns-are-gs --append-guanines-to-header)))

# Set of parameters for the first generation TS experiments
$(foreach bn,$(GEN1_TS_BNS),$(eval $(call collapse-rule,$(bn),\
--nucs 8 --remove-guanines --no-ns --append-guanines-to-header)))

# Set of parameters for the first generation TS experiments
$(foreach bn,$(GEN1_TS_BNS),$(eval $(call collapse-rule,$(bn),\
--nucs 8 --remove-guanines --no-ns --append-guanines-to-header)))

# Set of parameters for the first generation TS 
$(foreach bn,$(GEN1_5R_BNS),$(eval $(call collapse-rule,$(bn),\
--nucs 8 --no-ns --append-guanines-to-header)))

# Set of parameters for the G-trimmed 5' RACE
$(foreach bn,$(GEN1_5R_NOGS_BNS),$(eval $(call collapse-rule,$(bn),\
--nucs 8 --remove-guanines --no-ns --append-guanines-to-header)))

# short dummy
$(eval $(call collapse-rule,DUMMY_NS1,\
--nucs 12 --remove-guanines --remove-ns --ns-are-gs --append-guanines-to-header))
