########## MAKE PLOT OF LENGTH DISTRIBUTIONS ############################

distributions: $(LENGTH_DISTRIBUTIONS)

cleandistributions: 
	rm -f $(LENGTH_DISTRIBUTIONS)

LENGTH_DISTRIBUTIONS = $(foreach bn,$(TARGET_BNS),$(GRAPH_DIR)/$(bn)_length_distribution.png)

$(GRAPH_DIR)/%_length_distribution.png: $(INT_DIR)/%_qfiltered.fastq \
$(PLOT_LENGTH_DISTRIBUTION)
	$(BSUB) "python $(PLOT_LENGTH_DISTRIBUTION) -i $< -o $@"

# PLOT OF NUCLEOTIDE DISTRIBUTIONS, 3' and 5' NORMALIZED
PLOT_NUCLEOTIDE_DIST = $(DKLIB_SCRIPTS)/plot_nucleotide_dist.py

par_nucleotide_dist: $(PAR_NUCLEOTIDE_PLOTS)

PAR_NUCLEOTIDE_PLOTS = $(foreach bn,$(TARGET_BNS),$(INT_DIR)/$(bn)_par_nucleotide_dist.png)

$(INT_DIR)/%_nucleotide_dist.png: $(INT_DIR)/%_stats.txt
	$(BSUB) "python $(PLOT_NUCLEOTIDE_DIST)"

$(INT_DIR)/%_stats.txt: $(INT_DIR)/%_equalized.fastq
	$(BSUB) "fastx_quality_stats -N -i $< -o $@"

$(INT_DIR)/%_three_prime_equalized.fastq: $(INT_DIR)/%.fastq $(EQUALIZE_LENGTHS)
	$(BSUB) "python $(EQUALIZE_LENGTHS) --three-prime -i $< -o $@"

$(INT_DIR)/%_nucleotide_dist.png: $(INT_DIR) $(PLOT_NUCLEOTIDE_DIST)
	$(BSUB) "fastx_quality_stats -N -i $< -o $*_stats.txt"
	$(BSUB) "python $(PLOT_NUCLEOTIDE_DIST)"

# COMPARE DATA
