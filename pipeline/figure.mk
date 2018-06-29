COMPARE_DATA_PREFIXES = $(CURDIR)/$(GRAPH_DIR)/compare_data_with_legend $(CURDIR)/$(GRAPH_DIR)/compare_data_no_legend $(CURDIR)/$(GRAPH_DIR)/compare_data_no_color 

COMPARE_DATA_FIGURES = $(foreach prefix,$(COMPARE_DATA_PREFIXES),$(prefix).eps $(prefix).png)

TREE_ALGORITHM_FIGURES = $(foreach bn,$(foreach gene,$(INFLUENZA_GENES),$(gene)_TS),$(GRAPH_DIR)/$(bn)_tree_output.eps)

tree_algorithm: $(TREE_ALGORITHM_FIGURES)
clean_tree_algorithm:
	rm -f $(TREE_ALGORITHM_FIGURES)


NOG_STATS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_nog_stats.txt)

nog_stats: $(NOG_STATS) $(INT_DIR)/combined_nog_stats.txt
clean_nog_stats:
	rm -f $(NOG_STATS)

FIGURES = $(COMPARE_DATA_FIGURES) $(CURDIR)/$(GRAPH_DIR)/all_length_distributions.eps $(CURDIR)/$(GRAPH_DIR)/each_length_distribution.eps $(TREE_ALGORITHM_FIGURES)

RESULTS = $(GRAPH_DIR)/NS1_TS_CAAAAGCAG_trimmed_nonG_fraction.txt

BOTH_UPSTREAM_AND_DOWNSTREAM := $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_both_upstream_and_downstream.eps $(GRAPH_DIR)/$(bn)_both_upstream_and_downstream_nucleotide_dist.eps $(GRAPH_DIR)/$(bn)_both_upstream_and_downstream_relaxed.eps $(GRAPH_DIR)/$(bn)_both_upstream_and_downstream_relaxed_top_genes.eps $(GRAPH_DIR)/$(bn)_both_upstream_and_downstream_relaxed_ic.eps $(GRAPH_DIR)/$(bn)_both_upstream_and_downstream_mapped_to_trimmed_tss.eps)

both_upstream_and_downstream: $(BOTH_UPSTREAM_AND_DOWNSTREAM)

clean_both_upstream_and_downstream:
	rm -f $(BOTH_UPSTREAM_AND_DOWNSTREAM)

UPSTREAM_OF_CLEAVAGE := $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_upstream_of_cleavage_site_before_trimming_nucleotide_dist.eps $(GRAPH_DIR)/$(bn)_upstream_of_cleavage_site_before_trimming.eps $(GRAPH_DIR)/$(bn)_upstream_of_cleavage_deterministic.eps $(GRAPH_DIR)/$(bn)_upstream_of_cleavage_nucleotide_dist.eps $(GRAPH_DIR)/$(bn)_upstream_of_cleavage_deterministic_tss_background.eps)

upstream_of_cleavage: $(UPSTREAM_OF_CLEAVAGE)
clean_upstream_of_cleavage:
	rm -f $(UPSTREAM_OF_CLEAVAGE)

GCFS = $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_upstream_of_cleavage_from_gcfs.eps)

gcfs: $(GCFS)
clean_gcfs:
	rm -f $(GCFS)

DOWNSTREAM_OF_CLEAVAGE := $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_downstream_of_cleavage_site_nucleotide_dist.eps $(GRAPH_DIR)/$(bn)_downstream_of_cleavage_site.eps)

downstream_of_cleavage: $(DOWNSTREAM_OF_CLEAVAGE)
clean_downstream_of_cleavage:
	rm -f $(DOWNSTREAM_OF_CLEAVAGE)

DINUCLEOTIDE_CONTENT = $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_dinucleotide_content.eps $(GRAPH_DIR)/$(bn)_dinucleotide_content_with_background.eps)

dinucleotide_content: $(DINUCLEOTIDE_CONTENT)
clean_dinucleotide_content:
	rm -f $(DINUCLEOTIDE_CONTENT)

figures: $(FIGURES)

U_SNRNAS = $(GRAPH_DIR)/U_snRNA_vs_literature.png

# fraction A at position 0 | G not at position 1 / [(fraction G at position 0 or 1) + fraction A at position 0 | G not at pos 1]

u_snrnas: $(U_SNRNAS)
clean_u_snrnas:
	rm -f $(U_SNRNAS)

clean_figures:
	rm -f $(FIGURES)

results: $(RESULTS)

clean_results:
	rm -f $(RESULTS)

# U1 vs U2

$(CURDIR)/$(GRAPH_DIR)/testy.png: $(foreach bn,$(TS_BNS) $(GEN1_5R_NOGS_BNS),$(INT_DIR)/$(bn)_vfiltered_trimmed_summary.txt.gz) $(U1_VS_U2)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "for FILE in $(basename $(filter-out $(U1_VS_U2),$^)); do zcat $$FILE.gz > $$FILE; done;"
	$(BSUB) "$(MY_PYTHON) $(U1_VS_U2) --ts-infiles $(foreach bn,$(TS_BNS),$(INT_DIR)/$(bn)_vfiltered_trimmed_summary.txt) --race-infiles $(foreach bn,$(GEN1_5R_NOGS_BNS),$(INT_DIR)/$(bn)_vfiltered_summary.txt)"
	$(BSUB) "rm $(basename $(filter-out $(U1_VS_U2),$^))"
# PSSMs

$(CURDIR)/$(GRAPH_DIR)/GCAAAAGCAG_pssms.eps: $(foreach bn,$(GCAAAAGCAG_BNS),$(INT_DIR)/$(bn)_par_min3_max9_trim4_summary.txt) $(PSSM_PLOTS)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "$(MY_PYTHON) $(PSSM_PLOTS) --infiles $(filter-out $(PSSM_PLOTS),$^) --ignore-ns --log-likelihood --no-ic-scale --ncols-from-right 9 --infer-background --table --negative-offsets --no-plot-ns --format eps -o $@"

$(CURDIR)/$(GRAPH_DIR)/GCGAAAGCAG_pssms.eps: $(foreach bn,$(GCGAAAGCAG_BNS),$(INT_DIR)/$(bn)_par_min2_max9_trim4_summary.txt) $(PSSM_PLOTS)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "$(MY_PYTHON) $(PSSM_PLOTS) --infiles $(filter-out $(PSSM_PLOTS),$^) --ignore-ns --log-likelihood --no-ic-scale --ncols-from-right 9 --infer-background --table --negative-offsets --no-plot-ns --format eps -o $@"

$(U_SNRNAS): $(PLOT_U_SNRNAS)
	$(BSUB) "$(MY_PYTHON) $(PLOT_U_SNRNAS)"

# Length Distributions

$(eval $(call link-file,paper/figure/all_length_distributions.eps,$(CURDIR)/$(GRAPH_DIR)/all_length_distributions.eps))

$(CURDIR)/$(GRAPH_DIR)/each_length_distribution.eps: $(foreach bn,$(TS_BNS),$(INT_DIR)/$(bn)_vfiltered_summary.txt.gz) \
$(PLOT_EACH_LENGTH_DISTRIBUTION)
#$(foreach bn,$(TS_BNS),$(INT_DIR)/$(bn)_vfiltered_trimmed_summary.txt.gz)
	$(BSUB) "$(MY_PYTHON) $(PLOT_EACH_LENGTH_DISTRIBUTION) --gzipped --infiles $(foreach gene,$(INFLUENZA_GENES),$(INT_DIR)/$(gene)_TS_vfiltered_summary.txt.gz) -o $@"

$(CURDIR)/$(GRAPH_DIR)/each_pelchat_length_distribution.eps: $(foreach bn,$(PELCHAT_BNS_REAL),$(INT_DIR)/$(bn)_vfiltered_summary.txt.gz) \
$(PLOT_EACH_LENGTH_DISTRIBUTION)
	$(BSUB) "$(MY_PYTHON) $(PLOT_EACH_LENGTH_DISTRIBUTION) --gzipped --infiles $(filter-out $(PLOT_EACH_LENGTH_DISTRIBUTION),$^) -o $@"


# how many don't start with a G?

define cut-adapt-nog
zcat $< | cutadapt --trimmed-only --quality-base=64 -O \
9 -f fastq --error-rate=0 --minimum-length=1 \
-a $(1) - | $(MY_PYTHON) $(DISCARD_FASTQ) -l 9 -g 15 | $(MY_PYTHON) $(NOG_SCRIPT) > $@
endef

define nog-filter-rule
ifneq ($(filter $(1),$(GCAAAAGCAG_BNS)),)
$(INT_DIR)/$(1)_nog_stats.txt: FILTER_SEQ := CAAAAGCAG
else
$(INT_DIR)/$(1)_nog_stats.txt: FILTER_SEQ := CGAAAGCAG
endif
endef

$(INT_DIR)/combined_nog_stats.txt: $(foreach gene,$(INFLUENZA_GENES),$(INT_DIR)/$(gene)_TS_nog_stats.txt) $(NOG_STATS_COMBINER)
	$(BSUB) "$(MY_PYTHON) $(NOG_STATS_COMBINER) -i $(foreach gene,$(INFLUENZA_GENES),$(INT_DIR)/$(gene)_TS_nog_stats.txt) -o $@"

$(INT_DIR)/%_nog_stats.txt: $(INT_DIR)/%_collapsed.fastq.gz $(NOG_SCRIPT)
	$(BSUB) "$(call cut-adapt-nog,$(FILTER_SEQ))"

$(foreach gene,$(INFLUENZA_GENES),\
$(foreach bn,$($(gene)_BNS),\
$(eval $(call nog-filter-rule,$(bn)))))

$(CURDIR)/$(GRAPH_DIR)/each_length_distribution_trimmed.eps: $(foreach gene,$(INFLUENZA_GENES),$(INT_DIR)/$(gene)_TS_vfiltered_trimmed_summary.txt.gz) \
$(PLOT_EACH_LENGTH_DISTRIBUTION)
	$(BSUB) "$(MY_PYTHON) $(PLOT_EACH_LENGTH_DISTRIBUTION) --gzipped --infiles $(foreach gene,$(INFLUENZA_GENES),$(INT_DIR)/$(gene)_TS_vfiltered_trimmed_summary.txt.gz) -o $@"

$(INT_DIR)/%_vfiltered_trimmed_summary.txt.gz: $(INT_DIR)/%_vfiltered_trimmed.fastq.gz $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) --summarize --use-header | gzip -9 > $@"

# In parallel, trim vfiltered
$(INT_DIR)/%_vfiltered_trimmed.fastq.gz: $(INT_DIR)/%_vfiltered.fastq.gz $(COLLAPSE_PAR_V2)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(COLLAPSE_PAR_V2) --header -s $(TRIM_SEQUENCE) --num-trims 3 --min-length 2 --max-length 9 | gzip -9 > $@"

TRIM_SEQUENCE :=

# Determines the sequence to trim from the 3' of each influenza sequence
define trim-rule
ifneq ($(filter $(1),$(GCAAAAGCAG_BNS)),)
$(INT_DIR)/$(1)_vfiltered_trimmed.fastq.gz: \
TRIM_SEQUENCE := GCAAAAGCAG
else
$(INT_DIR)/$(1)_vfiltered_trimmed.fastq.gz: \
TRIM_SEQUENCE := GCGAAAGCAG
endif
endef

$(foreach bn,$(TARGET_BNS),$(eval $(call trim-rule,$(bn))))

$(CURDIR)/$(GRAPH_DIR)/all_length_distributions.eps: $(foreach bn,$(TS_BNS),$(INT_DIR)/$(bn)_qfiltered_summary.txt) $(PLOT_LENGTH_DISTRIBUTION)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "$(MY_PYTHON) $(PLOT_LENGTH_DISTRIBUTION) -i $(filter-out $(PLOT_LENGTH_DISTRIBUTION),$^) -f table -o $@"

$(CURDIR)/$(GRAPH_DIR)/compare_data_5gtrimmed.png: $(INT_DIR)/NS1_TS_5GTRIMMED_noshuffle_before_mapping_summary.txt.gz $(INT_DIR)/NS1_240_5GTRIMMED_noshuffle_before_mapping_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(CURDIR)/$(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --cpm --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --x-col COUNT --y-col COUNT -o $@ -f png --max-val 100000 --spearman-r --xlabel 'Biological Replicate 1 (cpm)' --ylabel 'Biological Replicate 2 (cpm)' --restrict-lengths 8 50 --cutoff-min-val 1"

$(CURDIR)/$(GRAPH_DIR)/compare_NS1_PB2_5gtrimmed_vfiltered.png: $(INT_DIR)/NS1_TS_5GTRIMMED_vfiltered_summary.txt.gz $(INT_DIR)/PB2_TS_vfiltered_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(CURDIR)/$(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --cpm --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --x-col COUNT --y-col COUNT -o $@ -f png --max-val 100000 --spearman-r --xlabel 'NS1 (cpm)' --ylabel 'PB2 (cpm)' --restrict-lengths 8 50 --snRNAs all --cutoff-min-val 1"

$(CURDIR)/$(GRAPH_DIR)/compare_data_5gtrimmed_vfiltered.png: $(INT_DIR)/NS1_TS_5GTRIMMED_sized_summary.txt.gz $(INT_DIR)/NS1_240_5GTRIMMED_vfiltered_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(CURDIR)/$(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --cpm --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --x-col COUNT --y-col COUNT -o $@ -f png --max-val 100000 --spearman-r --xlabel 'Biological Replicate 1 (cpm)' --ylabel 'Biological Replicate 2 (cpm)' --restrict-lengths 8 50 --cutoff-min-val 1"

$(CURDIR)/$(GRAPH_DIR)/compare_data_no_color.png: $(INT_DIR)/NS1_TS_vfiltered_summary.txt.gz $(INT_DIR)/NS1_240_vfiltered_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --x-col COUNT --y-col COUNT -o $@ -f png --max-val 100000 --spearman-r --xlabel 'NS1, replicate 1 (log\$$_{10}\$$ counts)' --ylabel 'NS1, replicate 2 (log\$$_{10}\$$ counts)' --restrict-lengths 8 50 --n-stats --cutoff-reads 7 --dpi 1200 -f png --add-histograms "

$(CURDIR)/$(GRAPH_DIR)/compare_data_100_statistics.png: $(INT_DIR)/NS1_TS_vfiltered_summary.txt.gz $(INT_DIR)/NS1_240_vfiltered_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --x-col COUNT --y-col COUNT -o $@ -f png --max-val 100000 --spearman-r --xlabel 'NS1, replicate 1 (log\$$_{10}\$$ counts)' --ylabel 'NS1, replicate 2 (log\$$_{10}\$$ counts)' --restrict-lengths 8 50 --n-stats --cutoff-reads 100 --dpi 1200 -f png --add-histograms "

$(CURDIR)/$(GRAPH_DIR)/compare_data_with_legend.eps: $(INT_DIR)/NS1_TS_noshuffle_before_mapping_summary.txt.gz $(INT_DIR)/NS1_240_noshuffle_before_mapping_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --cpm --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --x-col COUNT --y-col COUNT --legend -o $@ -f eps --max-val 10000 --spearman-r --xlabel 'Biological replicate 1 (cpm)' --ylabel 'Biological replicate 2 (cpm)' --restrict-lengths 8 50 --snRNAs all --add-histograms"

$(GRAPH_DIR)/%_5R_nucleotide_dist.eps: $(INT_DIR)/%_5R_vfiltered_summary.txt.gz $(PSSM_PLOT)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --plot-nucleotide-dist --no-plot-ns --ncols-from-left 15 -o $@ -f eps --ylabel 'Nucleotide Frequency'"

$(CURDIR)/$(GRAPH_DIR)/compare_data_no_legend.png: $(INT_DIR)/NS1_TS_par_min3_max9_trim4_summary.txt.gz $(INT_DIR)/NS1_240_par_min3_max9_trim4_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --add-histograms --x-col COUNT --y-col COUNT -o $@ -f png --max-val 1000000 --spearman-r --xlabel 'NS1, replicate 1 (log\$$_{10}\$$ counts)' --ylabel 'NS1, replicate 2 (log\$$_{10}\$$ counts)' --restrict-lengths 8 50 --snRNAs all --n-stats --cutoff-reads 7 --dpi 1200"

$(CURDIR)/$(GRAPH_DIR)/compare_NS1_PB2_no_legend.png: $(INT_DIR)/NS1_TS_noshuffle_before_mapping_summary.txt.gz $(INT_DIR)/PB2_TS_noshuffle_before_mapping_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --add-histograms --x-col COUNT --y-col COUNT -o $@ -f png --max-val 1000000 --spearman-r --xlabel 'NS1 (log\$$_{10}\$$ counts)' --ylabel 'PB2 (log\$$_{10}\$$ counts)' --restrict-lengths 8 50 --snRNAs all --n-stats --cutoff-reads 7 --dpi 1200"

$(CURDIR)/$(GRAPH_DIR)/compare_NS1_PB2_no_legend_trimmed.png: $(INT_DIR)/NS1_TS_par_min3_max9_trim4_summary.txt.gz $(INT_DIR)/PB2_TS_par_min2_max9_trim4_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --add-histograms --x-col COUNT --y-col COUNT -o $@ -f png --max-val 1000000 --spearman-r --xlabel 'NS1 (log\$$_{10}\$$ counts)' --ylabel 'PB2 (log\$$_{10}\$$ counts)' --restrict-lengths 8 50 --snRNAs all --n-stats --cutoff-reads 7 --dpi 1200"

$(CURDIR)/$(GRAPH_DIR)/compare_NS1_PB2_with_legend_trimmed.eps: $(INT_DIR)/NS1_TS_par_min3_max9_trim4_summary.txt.gz $(INT_DIR)/PB2_TS_par_min2_max9_trim4_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --add-histograms --x-col COUNT --y-col COUNT -o $@ -f eps --legend --max-val 1000000 --spearman-r --xlabel 'NS1 (log\$$_{10}\$$ counts)' --ylabel 'PB2 (log\$$_{10}\$$ counts)' --restrict-lengths 8 50 --snRNAs all --n-stats --cutoff-reads 7 --dpi 1200"

$(CURDIR)/$(GRAPH_DIR)/compare_data_no_legend_untrimmed.png: $(INT_DIR)/NS1_TS_sized_summary.txt.gz $(INT_DIR)/NS1_240_sized_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --add-histograms --x-col COUNT --y-col COUNT -o $@ -f png --max-val 100000 --spearman-r --xlabel 'NS1, Replicate 1 (\$$log_{10}\$$ counts)' --ylabel 'NS1, Replicate 2 (\$$log_{10}\$$ counts)' --restrict-lengths 8 50 --snRNAs all --n-stats --cutoff-reads 7 --dpi 1200 --include-par GCAAAAGCAG"

$(CURDIR)/$(GRAPH_DIR)/compare_data_no_color.eps: $(INT_DIR)/NS1_TS_noshuffle_before_mapping_summary.txt.gz $(INT_DIR)/NS1_240_noshuffle_before_mapping_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --gzipped -a $(word 1,$^) -b $(word 2,$^) --on SEQUENCE --x-col COUNT --y-col COUNT -o $@ -f eps --max-val 10000 --spearman-r --xlabel 'Biological Replicate 1 (counts)' --ylabel 'Biological Replicate 2 (counts)' --restrict-lengths 8 50"

$(CURDIR)/$(GRAPH_DIR)/compare_data_with_legend.eps: $(INT_DIR)/NS1_TS_noshuffle_before_mapping_summary.txt.gz $(INT_DIR)/NS1_240_noshuffle_before_mapping_summary.txt.gz $(COMPARE_DATA)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) --gzipped -a $(word 1,$^) -b $(word 2,$^) --cpm --on SEQUENCE --x-col COUNT --y-col COUNT --legend -o $@ -f eps --max-val 10000 --spearman-r --xlabel 'Biological Replicate 1 (counts)' --ylabel 'Biological Replicate 2 (counts)' --restrict-lengths 8 50 --snRNAs all --add-histograms"

####### nucleotide composition near cleavage site

# upstream of cleavage site -- considering only none trimmed

$(GRAPH_DIR)/%_considering_only_nonetrimmed_nucleotide_dist.eps: $(INT_DIR)/%_noshuffle_before_mapping_with_header_nonetrimmed.txt.gz $(PSSM_PLOT)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --plot-nucleotide-dist --ignore-ns -o $@ -f eps --ncols-from-right 8 --negative-offset --ylabel 'Nucleotide Frequency'"

$(GRAPH_DIR)/%_considering_only_nonetrimmed.eps: $(INT_DIR)/%_noshuffle_before_mapping_with_header_nonetrimmed.txt.gz $(PSSM_PLOT)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --infer-background --ignore-ns --log-likelihood --no-ic-scale -o $@ -f eps --ncols-from-right 8 --negative-offset --ylabel 'log\$$_2\$$(enrichment)'"

$(INT_DIR)/%_noshuffle_before_mapping_with_header_nonetrimmed.txt.gz: $(INT_DIR)/%_noshuffle_before_mapping_with_header.txt.gz $(SUBSET_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(SUBSET_TABLE) --three-prime-trimmed '' | gzip -9 > $@"

$(INT_DIR)/%_noshuffle_before_mapping_with_header.txt.gz: $(INT_DIR)/%_noshuffle_before_mapping.fastq.gz $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) --use-header | gzip -9 > $@"

# upstream of cleavage site -- before trimming
$(GRAPH_DIR)/%_upstream_of_cleavage_site_before_trimming_nucleotide_dist.eps: $(INT_DIR)/%_vfiltered_summary.txt.gz
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --plot-nucleotide-dist --ignore-ns -o $@ -f eps --ncols-from-right 8 --negative-offset --ylabel 'Nucleotide Frequency' --max-y 0.7"

# upstream of cleavage site -- log2
$(GRAPH_DIR)/%_upstream_of_cleavage_site_before_trimming.eps: $(INT_DIR)/%_vfiltered_summary.txt.gz $(ANNO_DIR)/gencode_offset_3p_20_nogs_5ptrim_5_summary.txt.gz $(PSSM_PLOT)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns --log-likelihood --no-ic-scale -o $@ -f eps --ncols-from-right 8 --negative-offset --ylabel 'log\$$_2\$$(enrichment)' --max-y 2.0 --simple-background-parse --flat-background --background-file $(word 2,$^) --gzipped-background --offset-bars"

# upstream of cleavage sites -- deterministic
$(GRAPH_DIR)/%_upstream_of_cleavage_deterministic.eps: $(INT_DIR)/%_noshuffle_before_mapping_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --infer-background --ignore-ns --log-likelihood --no-ic-scale -o $@ -f eps --ncols-from-right 8 --negative-offset --ylabel 'log\$$_2\$$(enrichment)' --max-y 2.0"

# upstream of cleavage using tss background
$(GRAPH_DIR)/%_upstream_of_cleavage_deterministic_tss_background.eps: $(INT_DIR)/%_noshuffle_before_mapping_summary.txt.gz $(ANNO_DIR)/gencode_offset_3p_20_nogs_5ptrim_5_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns --log-likelihood --no-ic-scale -o $@ -f eps --ncols-from-right 8 --negative-offset --ylabel 'log\$$_2\$$(enrichment)' --max-y 2.0 --simple-background-parse --flat-background --background-file $(word 2,$^) --gzipped-background --offset-bars"

$(GRAPH_DIR)/%_upstream_of_cleavage_nucleotide_dist.eps: $(INT_DIR)/%_noshuffle_before_mapping_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns --plot-nucleotide-dist -o $@ -f eps --ncols-from-right 8 --negative-offset --ylabel 'Nucleotide Frequency' --max-y 0.7"

# dinucleotide content
$(GRAPH_DIR)/%_dinucleotide_content.eps: $(INT_DIR)/%_noshuffle_restrict_to_zero_length_10.txt.gz $(FIND_DINUCLEOTIDE_AT_JUNCTION)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FIND_DINUCLEOTIDE_AT_JUNCTION) -o $@ --stats-file $(INT_DIR)/$*_dinucleotide_stats.txt"

$(GRAPH_DIR)/%_dinucleotide_content_with_background.eps: $(INT_DIR)/%_noshuffle_restrict_to_zero_length_10.txt.gz $(ANNO_DIR)/gencode_offset_3p_20_nogs_5ptrim_5_summary.txt.gz $(FIND_DINUCLEOTIDE_AT_JUNCTION_WITH_BACKGROUND)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FIND_DINUCLEOTIDE_AT_JUNCTION_WITH_BACKGROUND) -o $@ --stats-file $(INT_DIR)/$*_dinucleotide_stats_with_background.txt --background-file $(word 2,$^) --gzipped-background-file"

# downstream of cleavage site
$(GRAPH_DIR)/%_downstream_of_cleavage_site.eps: $(INT_DIR)/%_noshuffle_restrict_to_zero_length_10.txt.gz $(ANNO_DIR)/gencode_offset_3p_20_nogs_5ptrim_5_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns --log-likelihood --no-ic-scale -o $@ -f eps --ncols-from-left 8 --sequence-column THREEPRIME_OF_CLEAVAGE --ylabel 'log\$$_2\$$(enrichment)' --max-y 0.6 --weight-column WEIGHT --simple-background-parse --flat-background --background-file $(word 2,$^) --gzipped-background --offset-bars"

$(GRAPH_DIR)/%_both_upstream_and_downstream_relaxed.eps: $(INT_DIR)/%_noshuffle_length_10_around_five.txt.gz $(ANNO_DIR)/gencode_offset_3p_25_nogs_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns --log-likelihood --no-ic-scale -o $@ -f eps --ncols-from-right 17 --sequence-column FINAL_SEQUENCE --ylabel 'log\$$_2\$$(enrichment)' --max-y 2.0 --label-xticks -8 9 --weight-column WEIGHT --simple-background-parse --flat-background --background-file $(word 2,$^) --gzipped-background"

$(GRAPH_DIR)/%_both_upstream_and_downstream_relaxed_top_genes.eps: $(INT_DIR)/%_noshuffle_length_10_around_five.txt.gz $(INT_DIR)/RNASEQ_ZERO_top_genes_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns --log-likelihood --no-ic-scale -o $@ -f eps --ncols-from-right 17 --sequence-column FINAL_SEQUENCE --ylabel 'log\$$_2\$$(enrichment)' --max-y 2.0 --label-xticks -8 9 --weight-column WEIGHT --simple-background-parse --flat-background --background-file $(word 2,$^) --gzipped-background"

$(GRAPH_DIR)/%_both_upstream_and_downstream_relaxed_ic.eps: $(INT_DIR)/%_noshuffle_length_10_around_five.txt.gz $(ANNO_DIR)/gencode_offset_3p_25_nogs_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns -o $@ -f eps --ncols-from-right 17 --sequence-column FINAL_SEQUENCE --label-xticks -8 9 --weight-column WEIGHT --simple-background-parse --flat-background --background-file $(word 2,$^) --gzipped-background"

$(GRAPH_DIR)/%_both_upstream_and_downstream.eps: $(INT_DIR)/%_noshuffle_restrict_to_zero_length_10.txt.gz $(ANNO_DIR)/gencode_offset_3p_25_nogs_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns --log-likelihood --no-ic-scale -o $@ -f eps --ncols-from-right 17 --sequence-column FINAL_SEQUENCE --ylabel 'log\$$_2\$$(enrichment)' --max-y 2.0 --label-xticks -8 9 --weight-column WEIGHT --simple-background-parse --flat-background --background-file $(word 2,$^) --gzipped-background"

$(GRAPH_DIR)/%_both_upstream_and_downstream_mapped_to_trimmed_tss.eps: $(INT_DIR)/%_noshuffle_restrict_to_zero_length_10.txt.gz $(ANNO_DIR)/gencode_offset_3p_20_nogs_5ptrim_5_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns --log-likelihood --no-ic-scale -o $@ -f eps --ncols-from-right 17 --sequence-column FINAL_SEQUENCE --ylabel 'log\$$_2\$$(enrichment)' --max-y 2.0 --label-xticks -8 9 --weight-column WEIGHT --simple-background-parse --flat-background --background-file $(word 2,$^) --gzipped-background"

$(GRAPH_DIR)/%_both_upstream_and_downstream_nucleotide_dist.eps: $(INT_DIR)/%_noshuffle_restrict_to_zero_length_10.txt.gz $(ANNO_DIR)/gencode_offset_25_nogs_summary.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --infer-background --ignore-ns --plot-nucleotide-dist -o $@ -f eps --ncols-from-right 17 --sequence-column FINAL_SEQUENCE --weight-column WEIGHT --label-xticks -8 9 --ylabel 'Frequency' --max-y 0.5"

# downstream of cleavage site
$(GRAPH_DIR)/%_downstream_of_cleavage_site_nucleotide_dist.eps: $(INT_DIR)/%_noshuffle_restrict_to_zero_length_10.txt.gz $(PSSM_PLOT)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --infer-background --ignore-ns --plot-nucleotide-dist -o $@ -f eps --ncols-from-left 8 --sequence-column THREEPRIME_OF_CLEAVAGE --ylabel 'Frequency' --max-y 0.5"

# make graph of upstream cleavage site
$(GRAPH_DIR)/%_upstream_of_cleavage_from_gcfs.eps: $(INT_DIR)/%_gcfs.txt.gz
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PSSM_PLOT) --table --ignore-ns --plot-nucleotide-dist --infer-background -o $@ -f eps --ncols-from-right 8 --negative-offset"

# make gcf table

define prime-and-realign
zcat $< | $(MY_PYTHON) $(PRIME_AND_REALIGN) --output-gcfs --as-table --has-children --with-end $(1) --gcf-not-end-with $(1) | gzip -9 > $@
endef

$(INT_DIR)/%_gcfs.txt.gz: $(INT_DIR)/%_noshuffle_full_table.txt.gz $(PRIME_AND_REALIGN)
	$(BSUB) "$(call prime-and-realign,$(GCF_SEQ))"

define gcf-rule
# if base name is in GCAAA... base names, use GCA..., and vice versa. 
ifneq ($(filter $(1),$(GCAAAAGCAG_BNS)),)
$(INT_DIR)/$(1)_vfiltered.fastq.gz: GCF_SEQ := GCA
else
$(INT_DIR)/$(1)_vfiltered.fastq.gz: GCF_SEQ := GC
endif
endef

$(foreach gene,$(INFLUENZA_GENES),\
$(foreach bn,$($(gene)_BNS),\
$(eval $(call gcf-rule,$(bn)))))

# supplemental table
$(CURDIR)/$(GRAPH_DIR)/supplemental_table_1.xls: $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_noshuffle_full_table_no_paralogs.txt) $(MAKE_SUPPLEMENTAL_TABLE)
	mkdir -p $(GRAPH_DIR)
	$(BSUB) "$(MY_PYTHON) $(MAKE_SUPPLEMENTAL_TABLE) --target-directory $(CURDIR) -o $@"

############# RESULTS ##############

# These scripts create intermediate files necessary for certain results embedded
# in the paper

# Is it possible that influenza additionally primes off of the antepenultimate nucleotide of the viral sequence, a guanosine? After all, it has been shown that GpC can prime influenza transcription {Plotch:1977vz, Eriksson:1977ut}. Our results indicate that this may occur, but only, which is in the range of sequencing errors from Illumina machines ().

$(GRAPH_DIR)/NS1_TS_CAAAAGCAG_trimmed_nonG_fraction.txt: $(INT_DIR)/NS1_TS_CAAAGCAG_trimmed.txt $(QUICK_COUNT_NON_CAAAAGCAG_TRIMMED)
	$(BSUB) "$(MY_PYTHON) $(QUICK_COUNT_NON_CAAAAGCAG_TRIMMED) -i $< > $@"

$(INT_DIR)/NS1_TS_CAAAGCAG_trimmed.txt: $(INT_DIR)/NS1_TS_CAAAAGCAG_trimmed.fastq $(COUNT_TRANSCRIPTS)
	$(BSUB) "$(MY_PYTHON) $(COUNT_TRANSCRIPTS) --make-table $< > $@"

$(INT_DIR)/NS1_TS_CAAAAGCAG_trimmed.fastq: $(INT_DIR)/NS1_TS_collapsed.fastq
	$(BSUB) "cutadapt --trimmed-only --quality-base=64 -O 9 -f fastq --error-rate=0 --minimum-length=1 -a CAAAAGCAG $< > $@"


$(GRAPH_DIR)/par_stats.eps: $(foreach bn,$(GEN2_TS_BNS) $(GEN1_TS_BNS),$(INT_DIR)/$(bn)_par.fastq.gz) $(PLOT_PAR_LENGTHS)
	$(BSUB) "$(MY_PYTHON) $(PLOT_PAR_LENGTHS) --GCAAAAGCAG-files $(foreach bn,HA_TS NA_TS NP_TS NS1_TS,$(INT_DIR)/$(bn)_par_stats.txt) --GCGAAAGCAG-files $(foreach bn,MP_TS PA_TS PB1_TS PB2_TS,$(INT_DIR)/$(bn)_par_stats.txt) -o $@ --min-num 2 --max-num 7 -f eps"

###########

# Compare unbarcoded to barcoded samples

$(GRAPH_DIR)/TS_vs_240_5gtrimmed.png: $(INT_DIR)/NS1_TS_5GTRIMMED_par_summary.txt.gz $(INT_DIR)/NS1_240_5GTRIMMED_par_summary.txt.gz $(COMPARE_DATA)
	$(BSUB) "zcat $(word 1,$^) > $(basename $(word 1,$^)) && zcat $(word 2,$^) > $(basename $(word 2,$^))"
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(COMPARE_DATA) -a $(basename $(word 1,$^)) -b $(basename $(word 2,$^)) --no-header --on SEQUENCE --x-col COUNT --y-col COUNT -o $@ -f png --max-val 100000 --spearman-r --xlabel 'Biological Replicate 1 (reads)' --ylabel 'Biological Replicate 2 (reads)' ---lengths 8 50"
	rm $(basename $(word 1,$^)) $(basename $(word 2,$^))

$(GRAPH_DIR)/u1_vs_u2.png: $(foreach bn,$(GEN2_TS_BNS) $(GEN1_TS_BNS),$(INT_DIR)/$(bn)_noshuffle_before_mapping_summary.txt.gz) $(U1_VS_U2)
	$(BSUB) -R rusage[mem=141000] "$(MY_PYTHON) $(U1_VS_U2) --gen-one-ts-infiles $(filter-out $(U1_VS_U2),$^) --gen-two-ts-infiles $(filter-out $(U1_VS_U2),$^) --gzipped -o $@"

$(GRAPH_DIR)/gen2_distributions.png: $(foreach bn,$(GEN2_TIME_COURSE_BNS),$(INT_DIR)/$(bn)_noshuffle_before_mapping_summary.txt.gz) $(PLOT_DISTS)
	$(BSUB) "$(MY_PYTHON) $(PLOT_DISTS) --gzipped --infiles $(filter-out $(PLOT_DISTS),$^) --xlabel 'log10(cpm)' --ylabel 'Frequency' --normalize --log-x --cutoff-reads 10 --dpi 600 -o $@"

$(GRAPH_DIR)/influenza_distributions.png: $(foreach bn,$(foreach gene,$(INFLUENZA_GENES),$(gene)_TS),$(INT_DIR)/$(bn)_noshuffle_before_mapping_summary.txt.gz) $(PLOT_DISTS)
	$(BSUB) "$(MY_PYTHON) $(PLOT_DISTS) --gzipped --infiles $(filter-out $(PLOT_DISTS),$^) --normalize --log-x --xlabel 'log10(cpm)' --ylabel 'Frequency' --cutoff-reads 10 --dpi 600 -o $@"

$(GRAPH_DIR)/%_tree_output.eps: $(INT_DIR)/%_noshuffle_full_table.txt.gz $(PRIME_AND_REALIGN)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(PRIME_AND_REALIGN) -f eps -o $@"

$(GRAPH_DIR)/pelchat_distributions.png: $(foreach bn,$(PELCHAT_BNS_REAL),$(INT_DIR)/$(bn)_noshuffle_before_mapping_summary.txt.gz) $(PLOT_DISTS)
	$(BSUB) "$(MY_PYTHON) $(PLOT_DISTS) --gzipped --infiles $(filter-out $(PLOT_DISTS),$^) --normalize --xlabel 'log10(cpm)' --ylabel 'Frequency' --log-x --cutoff-reads 10 --dpi 600 -o $@"

# Linear regression
$(GRAPH_DIR)/linear_regression.png: $(foreach bn,$(GEN2_TIME_COURSE_BNS),$(INT_DIR)/$(bn)_noshuffle_before_mapping_summary.txt.gz) $(LINEAR_REGRESSION)
	$(BSUB) "$(MY_PYTHON) $(LINEAR_REGRESSION) --gzipped --infile-prefix $(INT_DIR)/NS1_ --infile-suffix _noshuffle_before_mapping_summary.txt.gz --quantile-normalize -f png -o $@"

####
# weblogo -a AGTC < intermediate/NS1_TS_noshuffle_before_mapping_ncols_from_right_8.fasta > weblogo0.eps
