####### MAKE TARGETS ############

# Remember to declare targets *before* making the general rule! 
VFILTERED_TARGETS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_vfiltered.fastq.gz)

# This is declared as .PHONY in common.mk
vfiltered: $(VFILTERED_TARGETS)

VFILTERED_SUMMARY = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_vfiltered_summary.txt.gz)

vfiltered_summary: $(VFILTERED_SUMMARY)

clean_vfiltered_summary:
	rm -f $(VFILTERED_SUMMARY)

clean_vfiltered: 
	rm -f $(VFILTERED_TARGETS)

####### CONSTANTS ##########

# First n nucleotides to extract from the fasta file, keyed by a given ID
# Also, how many nucleotides of overlap are required in order for cutadapt to 
# trim the 3' end of the read
EXTRACT_FASTA_NUCS = 10

# Use cutadapt to filter the viral sequences; use error rate=0 as of 4/14/2014
# Argument $(1) is filter sequence
define cut-adapt
zcat $< | cutadapt --trimmed-only --quality-base=64 -O \
$(EXTRACT_FASTA_NUCS) -f fastq --error-rate=0 --minimum-length=1 \
-a $(1) - | gzip -9 > $@
endef

$(INT_DIR)/%_vfiltered_summary.txt.gz: $(INT_DIR)/%_vfiltered.fastq.gz $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) --use-header --summarize | gzip -9 > $@"

$(INT_DIR)/%_vfiltered.fastq.gz: $(INT_DIR)/%_collapsed.fastq.gz $(INFLUENZA_TRANSCRIPTS)
	$(BSUB) "$(call cut-adapt,$(FILTER_SEQ))"

# CAACGGAATCCCAAAAGCAGCTGTCGTATGCC
# GGCATACGACAGC

define filter-rule
# if base name is in GCAAA... base names, use GCA..., and vice versa. 
ifneq ($(filter $(1),$(GCAAAAGCAG_BNS)),)
$(INT_DIR)/$(1)_vfiltered.fastq.gz: FILTER_SEQ := GCAAAAGCAG
else
$(INT_DIR)/$(1)_vfiltered.fastq.gz: FILTER_SEQ := GCGAAAGCAG
endif
endef

# Filter by gene. 
$(foreach gene,$(INFLUENZA_GENES),\
$(foreach bn,$($(gene)_BNS),\
$(eval $(call filter-rule,$(bn)))))
$(eval $(call filter-rule,DUMMY_NS1))

# for 5GTRIMMED files
define gcaaaagcag-filter-rule
$(INT_DIR)/$(1)_5GTRIMMED_vfiltered.fastq.gz: FILTER_SEQ := GCAAAAGCAG
endef

define gcgaaagcag-filter-rule
$(INT_DIR)/$(1)_5GTRIMMED_vfiltered.fastq.gz: FILTER_SEQ := GCGAAAGCAG
endef

$(foreach bn,NS1_TS HA_TS NP_TS NA_TS $(GEN2_TIME_COURSE_BNS),\
$(eval $(call gcaaaagcag-filter-rule,$(bn))))

$(foreach bn,PB2_TS PB1_TS PA_TS MP_TS,\
$(eval $(call gcgaaagcag-filter-rule,$(bn))))
