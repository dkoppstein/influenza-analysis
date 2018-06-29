########## FURTHER COLLAPSE PRIME-AND-REALIGNED SEQUENCES #####################

# biological replicates:
# /lab/bartel1_ata/koppstein/virtualenvs/my2.7/bin/python /lab/bartel1_ata/koppstein/science/software/dklib/influenza/compare_data_v2.py -a intermediate/NS1_TS_par_summary.txt -b intermediate/NS1_240_par_summary.txt --no-header --on Sequence --x-col Count --y-col Count -o test2.png -f png --max-val 100000 --spearman-r --xlabel 'Biological Replicate 1' --ylabel 'Biological Replicate 2' --title 'NS1 Sequence Abundances, 4 Hours Post-Infection'

PAR_TARGETS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_par.fastq)

PAR_COLLAPSED = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_par_collapsed_before_mapping.fasta)

PAR_SUMMARY = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_par_summary.txt)

PAR_STATS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_par_stats.txt)

PAR_GRAPHS = $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_par_stats.png $(GRAPH_DIR)/$(bn)_not_full_length_par_stats.png)

PAR_TRIMMED = $(foreach bn,$(GCAAAAGCAG_BNS),$(INT_DIR)/$(bn)_par_min3_max9_trim4.fastq.gz) $(foreach bn,$(GCGAAAGCAG_BNS),$(INT_DIR)/$(bn)_par_min2_max9_trim4.fastq.gz)

par_trimmed: $(PAR_TRIMMED)
clean_par_trimmed:
	rm -f $(PAR_TRIMMED)

BEFORE_MAPPING = $(foreach shuffle,noshuffle shuffled,\
$(foreach bn,$(TARGET_BNS),\
$(INT_DIR)/$(bn)_$(shuffle)_before_mapping.fasta.gz))

PELCHAT_BEFORE_MAPPING = $(foreach bn,$(PELCHAT_BNS),$(INT_DIR)/$(bn)_noshuffle_before_mapping_summary.txt.gz)

pelchat_before_mapping: $(PELCHAT_BEFORE_MAPPING)
clean_pelchat_before_mapping:
	rm -f $(PELCHAT_BEFORE_MAPPING)

# Summary: a table of the collapsed sequences
par_summary: $(PAR_SUMMARY)

par_graphs: $(PAR_GRAPHS)

clean_par_graphs: 
	rm -f $(PAR_GRAPHS)

SIZED_SUMMARY := $(foreach bn,$(TARGET_BNS) $(PELCHAT_BNS),$(INT_DIR)/$(bn)_sized_summary.txt.gz)

sized_summary: $(SIZED_SUMMARY)
clean_sized_summary:
	rm -f $(SIZED_SUMMARY)

# Stats: a table of the frequency with which we trimmed sequences
par_stats: $(PAR_STATS)

par: $(INFLUENZA_TRANSCRIPTS) $(PAR_TARGETS)

par_v3: $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_par_v3.fastq)

before_mapping: $(BEFORE_MAPPING)

clean_before_mapping: 
	rm -f $(BEFORE_MAPPING)

BEFORE_MAPPING_SUMMARY = $(foreach shuffle,noshuffle shuffled,\
$(foreach bn,$(TARGET_BNS) NS1_240_5GTRIMMED NS1_TS_5GTRIMMED,\
$(INT_DIR)/$(bn)_$(shuffle)_before_mapping_summary.txt.gz)) $(foreach gene,$(INFLUENZA_GENES),$(INT_DIR)/$(gene)_TS_5GTRIMMED_noshuffle_before_mapping_summary.txt.gz) $(foreach bn,$(GEN2_TIME_COURSE_BNS),$(INT_DIR)/$(bn)_5GTRIMMED_noshuffle_before_mapping_summary.txt.gz)

before_mapping_summary: $(BEFORE_MAPPING_SUMMARY)

clean_before_mapping_summary:
	rm -f $(BEFORE_MAPPING_SUMMARY)

clean_par: 
	rm -f $(PAR_TARGETS) $(PAR_SUMMARY) $(PAR_STATS) $(PAR_GRAPHS) $(PAR_COLLAPSED)

#####################################
######## PRIME-AND-REALIGN ##########
#####################################

min_three: $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_par_min3_max9_trim4.png)

min_two: $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_par_min2_max9_trim4.png)

$(GRAPH_DIR)/%_par_min2_max9_trim4.png: $(INT_DIR)/%_par_min2_max9_trim4_summary.txt $(PSSM_PLOT)
	$(BSUB) "$(MY_PYTHON) $(PSSM_PLOT) --ignore-ns --table --plot-nucleotide-dist --ncols-from-right 9 --ylabel Frequency --negative-offsets --max-y 0.7 -i $< -o $@"

$(INT_DIR)/%_par_min2_max9_trim4_summary.txt.gz: $(INT_DIR)/%_par_min2_max9_trim4.fastq.gz $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) -f fastq-illumina --sort --summarize --use-header | gzip -9 > $@"

$(GRAPH_DIR)/%_par_min3_max9_trim4.png: $(INT_DIR)/%_par_min3_max9_trim4_summary.txt $(PSSM_PLOT)
	$(BSUB) "$(MY_PYTHON) $(PSSM_PLOT) --ignore-ns --table --plot-nucleotide-dist --ncols-from-right 9 --ylabel Frequency --negative-offsets --max-y 0.7 -i $< -o $@"

$(INT_DIR)/%_par_min3_max9_trim4_summary.txt.gz: $(INT_DIR)/%_par_min3_max9_trim4.fastq.gz $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) --use-header -f fastq-illumina --sort --summarize | gzip -9 > $@"

# Subtract prime-and-realigned nucleotides
$(INT_DIR)/%_par_min3_max9_trim4.fastq.gz $(INT_DIR)/%_par_min3_max9_trim4_stats.txt: $(INT_DIR)/%_sized.fastq.gz \
$(COLLAPSE_PAR_V2) $(INFLUENZA_TRANSCRIPTS)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(COLLAPSE_PAR_V2) --header -s $(PAR_SEQUENCE) --stats-file $(INT_DIR)/$*_par_min3_max9_trim4_stats.txt --num-trims 4 --min-length 3 --max-length 9 | gzip -9 > $(INT_DIR)/$*_par_min3_max9_trim4.fastq.gz"

# Subtract prime-and-realigned nucleotides
$(INT_DIR)/%_par.fastq.gz $(INT_DIR)/%_par_stats.txt: $(INT_DIR)/%_sized.fastq.gz \
$(COLLAPSE_PAR_V2) $(INFLUENZA_TRANSCRIPTS)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(COLLAPSE_PAR_V2) --header -s $(PAR_SEQUENCE) --min-length 1 --max-length 9 --num-trims 1 --stats-file $(INT_DIR)/$*_par_stats.txt | gzip -9 > $@"

PAR_SEQUENCE :=

# Determines the sequence to trim from the 3' of each influenza sequence
define par-rule
ifneq ($(filter $(1),$(GCAAAAGCAG_BNS)),)
$(INT_DIR)/$(1)_par.fastq.gz $(INT_DIR)/$(1)_par_stats.txt: \
PAR_SEQUENCE := GCAAAAGCAG
else
$(INT_DIR)/$(1)_par.fastq.gz $(INT_DIR)/$(1)_par_stats.txt: \
PAR_SEQUENCE := GCGAAAGCAG
endif
endef

# Determines the sequence to trim from the 3' of each influenza sequence
define par-rule-v2
ifneq ($(filter $(1),$(GCAAAAGCAG_BNS)),)
$(INT_DIR)/$(1)_par_min3_max9_trim4.fastq.gz $(INT_DIR)/$(1)_par_min3_max9_trim4_stats.txt: \
PAR_SEQUENCE := GCAAAAGCAG
else
$(INT_DIR)/$(1)_par_min3_max9_trim4.fastq.gz $(INT_DIR)/$(1)_par_min3_max9_trim4_stats.txt: \
PAR_SEQUENCE := GCGAAAGCAG
endif
endef

# Determines the sequence to trim from the 3' of each influenza sequence
define gcgaaagcag-par-rule
$(INT_DIR)/$(1)_par_min2_max9_trim4.fastq.gz $(INT_DIR)/$(1)_par_min2_max9_trim4_stats.txt: $(INT_DIR)/$(1)_sized.fastq.gz \
$(COLLAPSE_PAR_V2) $(INFLUENZA_TRANSCRIPTS)
	$(BSUB) "zcat $$< | $(MY_PYTHON) $(COLLAPSE_PAR_V2) --header -s GCGAAAGCAG --stats-file $(INT_DIR)/$(1)_par_min2_max9_trim4_stats.txt --num-trims 4 --min-length 2 --max-length 9 | gzip -9 > $(INT_DIR)/$(1)_par_min2_max9_trim4.fastq.gz"
endef

define gcaaaagcag-par-rule
$(INT_DIR)/$(1)_par_min3_max9_trim4.fastq.gz $(INT_DIR)/$(1)_par_min3_max9_trim4_stats.txt: $(INT_DIR)/$(1)_sized.fastq.gz \
$(COLLAPSE_PAR_V2) $(INFLUENZA_TRANSCRIPTS)
	$(BSUB) "zcat $$< | $(MY_PYTHON) $(COLLAPSE_PAR_V2) --header -s GCAAAAGCAG --stats-file $(INT_DIR)/$(1)_par_min3_max9_trim4_stats.txt --num-trims 4 --min-length 3 --max-length 9 | gzip -9 > $(INT_DIR)/$(1)_par_min3_max9_trim4.fastq.gz"
endef

$(foreach bn,$(GCGAAAGCAG_BNS),$(eval $(call gcgaaagcag-par-rule,$(bn))))
$(foreach bn,$(GCAAAAGCAG_BNS),$(eval $(call gcaaaagcag-par-rule,$(bn))))

# Call rule
$(foreach bn,$(TARGET_BNS),$(eval $(call par-rule,$(bn))))
#$(foreach bn,$(TARGET_BNS),$(eval $(call par-rule-v2,$(bn))))
#$(foreach bn,$(TARGET_BNS),$(eval $(call par-rule-v3,$(bn))))

$(INT_DIR)/NS1_TS_5GTRIMMED_par_min3_max9_trim4.fastq.gz $(INT_DIR)/NS1_TS_5GTRIMMED_par_min3_max9_trim4_stats.txt: \
PAR_SEQUENCE := GCAAAAGCAG

$(INT_DIR)/NS1_240_5GTRIMMED_par_min3_max9_trim4.fastq.gz $(INT_DIR)/NS1_240_5GTRIMMED_par_min3_max9_trim4_stats.txt: \
PAR_SEQUENCE := GCAAAAGCAG


##################################################
######### SIZE, COLLAPSE, AND SHUFFLE ############
##################################################

MIN_NUCS = 9
MAX_NUCS = 15

$(INT_DIR)/%_sized_summary.txt.gz: $(INT_DIR)/%_sized.fastq.gz $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) --summarize --sort -f fastq-illumina | gzip -9 > $@"

# Three additions: noshuffle shuffled
$(INT_DIR)/%_sized.fastq.gz: $(INT_DIR)/%_qfiltered.fastq.gz $(DISCARD_FASTQ)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(DISCARD_FASTQ) -l $(MIN_NUCS) -g $(MAX_NUCS) | gzip -9 > $@"

define GCAAAAGCAG-rule
$(INT_DIR)/$(1)_noshuffle_before_mapping.fastq.gz: $(INT_DIR)/$(1)_par_min3_max9_trim4.fastq.gz $(CUSTOM_COLLAPSER) $(DISCARD_FASTQ)
	$(BSUB) "zcat $$< | $(MY_PYTHON) $(CUSTOM_COLLAPSER) | $(MY_PYTHON) $(DISCARD_FASTQ) -f fasta -l $(MIN_NUCS) | fasta_to_fastq | gzip -9 > $$@"
endef

define GCGAAAGCAG-rule
$(INT_DIR)/$(1)_noshuffle_before_mapping.fastq.gz: $(INT_DIR)/$(1)_par_min2_max9_trim4.fastq.gz $(CUSTOM_COLLAPSER) $(DISCARD_FASTQ)
	$(BSUB) "zcat $$< | $(MY_PYTHON) $(CUSTOM_COLLAPSER) | $(MY_PYTHON) $(DISCARD_FASTQ) -f fasta -l $(MIN_NUCS) | fasta_to_fastq | gzip -9 > $$@"
endef

$(foreach bn,$(GCAAAAGCAG_BNS) $(foreach bn,NS1_TS HA_TS NP_TS NA_TS $(GEN2_TIME_COURSE_BNS),$(bn)_5GTRIMMED),$(eval $(call GCAAAAGCAG-rule,$(bn))))
$(foreach bn,$(GCGAAAGCAG_BNS) $(foreach bn,PB2_TS PB1_TS PA_TS MP_TS,$(bn)_5GTRIMMED),$(eval $(call GCGAAAGCAG-rule,$(bn))))

$(INT_DIR)/%_shuffled_before_mapping.fastq.gz: $(INT_DIR)/%_noshuffle_before_mapping.fastq.gz $(SHUFFLE_SCRIPT)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(SHUFFLE_SCRIPT) -f fastq-illumina --preserve-cpgs --not-same --no-gs | gzip -9 > $@"

$(INT_DIR)/%_before_mapping_summary.txt.gz: $(INT_DIR)/%_before_mapping.fastq.gz $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) -f fastq-illumina --sort --use-header --summarize | gzip -9 > $@"

###### PLOTS #####

# Plots of the prime-and-realign
$(GRAPH_DIR)/%_par_stats.png: $(INT_DIR)/%_par.fastq \
$(PLOT_PAR_AND_LENGTH) $(LENGTH_VS_PAR)
# Normalize by the number of counts
	mkdir -p $(GRAPH_DIR)
	if [ -f $*_tmp.txt ]; then rm -f $*_tmp.txt; fi;
	$(BSUB) "python $(LENGTH_VS_PAR) --full-length -i $< > $*_tmp.txt"
	$(BSUB) "Rscript $(PLOT_PAR_AND_LENGTH) -i $*_tmp.txt -o $@ -y 400 --title $*_full_length_cellular_fragment -l $(MIN_NUCS) -x $(MAX_NUCS)"
	rm -f $*_tmp.txt

# Plots of the prime-and-realign
$(GRAPH_DIR)/%_par_stats_no_colors.png: $(INT_DIR)/%_par.fastq \
$(PLOT_PAR_AND_LENGTH) $(LENGTH_VS_PAR)
# Normalize by the number of counts
	mkdir -p $(GRAPH_DIR)
	if [ -f $*_tmp.txt ]; then rm -f $*_tmp.txt; fi;
	$(BSUB) "python $(LENGTH_VS_PAR) --full-length -i $< > $*_tmp.txt"
	$(BSUB) "Rscript $(PLOT_PAR_AND_LENGTH) -i $*_tmp.txt -o $@ -y 400 --title $*_full_length_cellular_fragment -l $(MIN_NUCS) -x $(MAX_NUCS)"
	rm -f $*_tmp.txt

$(GRAPH_DIR)/%_not_full_length_par_stats.png: $(INT_DIR)/%_par.fastq \
$(PLOT_PAR_AND_LENGTH) $(LENGTH_VS_PAR)
# Normalize by the number of counts
	mkdir -p $(GRAPH_DIR)
	if [ -f $*_tmp_not_full_length.txt ]; then rm -f $*_tmp_not_full_length.txt; fi;
	$(BSUB) "python $(LENGTH_VS_PAR) -i $< > $*_tmp_not_full_length.txt"
	$(BSUB) "Rscript $(PLOT_PAR_AND_LENGTH) -i $*_tmp_not_full_length.txt -o $@ -n $$((`wc -l < $<` / 4)) -y 400 --title $* -l 3 -x $(MAX_NUCS)"
	rm -f $*_tmp_not_full_length.txt

# Make a summary of the prime and realigned data, and sort it
$(INT_DIR)/%_par_summary.txt.gz: $(INT_DIR)/%_par.fastq.gz $(FASTQ_TO_TABLE)
	$(BSUB) "zcat $< | $(MY_PYTHON) $(FASTQ_TO_TABLE) --use-header --sort | gzip -9 > $@"
