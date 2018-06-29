QFILTERED_TARGETS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_qfiltered.fastq.gz)

QFILTERED_SUMMARY = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_qfiltered_summary.txt.gz)

qfiltered: $(QFILTERED_TARGETS)

NUC_DISTS = $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_qfiltered_length_distribution.eps)

nuc_dists: $(NUC_DISTS)

clean_nuc_dists:
	rm -f $(NUC_DISTS)

cleanqfiltered: 
	rm -f $(QFILTERED_TARGETS)

qfiltered_summary: $(QFILTERED_SUMMARY)

clean_qfiltered_summary:
	rm -f $(QFILTERED_SUMMARY)

QFILTERED_PSSM = $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_qfiltered_pssm.png)

qfiltered_pssm: $(QFILTERED_PSSM)

clean_qfiltered_pssm:
	rm -f $(QFILTERED_PSSM)

$(GRAPH_DIR)/%_qfiltered_length_distribution.eps: $(INT_DIR)/%_qfiltered_summary.txt $(PSSM_PLOT)
	$(BSUB) "$(MY_PYTHON) $(PSSM_PLOT) --table --no-plot-ns -i $< -o $@ -f eps --plot-nucleotide-dist --ncols-from-right 9 --ylabel Frequency --restrict-lengths 9 15 --negative-offsets --max-y 0.7"

$(INT_DIR)/%_qfiltered_summary.txt.gz: $(INT_DIR)/%_qfiltered.fastq.gz $(COUNT_TRANSCRIPTS) $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) --use-header --sort | gzip -9 > $@"

$(INT_DIR)/%_qfiltered.fastq.gz: $(INT_DIR)/%_vfiltered.fastq.gz $(FASTQ_QUALITY_FILTER)
	$(BSUB) "zcat $< | python $(FASTQ_QUALITY_FILTER) --no-ns | gzip -9 > $@"

$(GRAPH_DIR)/%_qfiltered_pssm.png: $(INT_DIR)/%_qfiltered_summary.txt $(PSSM_PLOT)
	$(BSUB) "$(MY_PYTHON) $(PSSM_PLOT) -i $< --table --ncols-from-right 9 --restrict-lengths 9 15 --log-likelihood --no-ic-scale --infer-background -o $@"
