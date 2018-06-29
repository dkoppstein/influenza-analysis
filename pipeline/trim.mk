# Target files
TRIMMED_TARGETS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_trimmed.fastq)

TRIMMED_SUMMARY = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_trimmed_summary.txt)

# .PHONY targets
trimmed_summary: $(TRIMMED_SUMMARY)

trimmed: $(TRIMMED_TARGETS)

cleantrimmed: 
	rm -f $(TRIMMED_TARGETS)

# Minimum number of nucleotides in the FASTQ file before mapping
MIN_NUCS = 9

# Maximum number of nucleotides in the FASTQ file before mapping
MAX_NUCS = 15

# Make Rules

# Make a summary of the trimmed data, and sort it
$(INT_DIR)/%_trimmed_summary.txt: $(INT_DIR)/%_trimmed.fastq $(COUNT_TRANSCRIPTS)
	$(BSUB) "python $(COUNT_TRANSCRIPTS) --make-table $< | sort -r -n -k 2 > $@"

$(INT_DIR)/%_trimmed.fastq: $(INT_DIR)/%_qfiltered.fastq $(DISCARD_FASTQ)
	$(BSUB) "python $(DISCARD_FASTQ) -l $(MIN_NUCS) -i $(MAX_NUCS) $< > $@"
