ANNO = gencode
ANNO_TSS = $(ANNO)$(TSS_SUFFIX)
GENCODE_LOCATION = /nfs/genomes/human_gp_feb_09_no_random/gtf/gencode.v17.annotation.gtf
SHUFFLED_OPTIONS = noshuffle shuffled

BEDTOOLS_INTERSECT_GRAPH_TARGETS = $(foreach bn,$(TARGET_BNS),\
$(GRAPH_DIR)/$(bn)_graph_of_tss_length_10.eps)

TABLES = $(foreach bn,$(TARGET_BNS),\
$(foreach shuffled,$(SHUFFLED_OPTIONS),\
$(INT_DIR)/$(bn)_$(shuffled)_full_table.txt.gz))

bedtools_graph_targets: $(BEDTOOLS_INTERSECT_GRAPH_TARGETS)
clean_bedtools_graph_targets: 
	rm -f $(BEDTOOLS_INTERSECT_GRAPH_TARGETS)

MAP_TARGETS = $(foreach bn,$(TARGET_BNS),\
$(foreach shuffled,$(SHUFFLED_OPTIONS),\
$(INT_DIR)/$(bn)_$(shuffled)_mapped_to_hg19_all.bam))

INTERSECT_WITH_GENCODE = $(foreach bn,$(TARGET_BNS),\
$(foreach shuffled,$(SHUFFLED_OPTIONS),\
$(INT_DIR)/$(bn)_$(shuffled)_mapped_to_hg19_intersect_with_gencode.bed \
$(INT_DIR)/$(bn)_$(shuffled)_mapped_to_hg19_intersect_with_gencode.fasta))

intersect_with_gencode: $(INTERSECT_WITH_GENCODE)
clean_intersect_with_gencode:
	rm -f $(INTERSECT_WITH_GENCODE)

TABLES_NO_PARALOGS = $(foreach bn,$(TARGET_BNS),\
$(foreach shuffled,$(SHUFFLED_OPTIONS),\
$(INT_DIR)/$(bn)_$(shuffled)_full_table_no_paralogs.txt.gz))

tables_no_paralogs: $(TABLES_NO_PARALOGS)
clean_tables_no_paralogs:
	rm -f $(TABLES_NO_PARALOGS)

MAPPED_ONLY_TABLES = $(foreach bn,$(TARGET_BNS),\
$(foreach shuffled,$(SHUFFLED_OPTIONS),\
$(INT_DIR)/$(bn)_$(shuffled)_full_table_mapped_only.txt.gz))

map: $(MAP_TARGETS)
clean_map: 
	rm -f $(MAP_TARGETS)

tables: | $(TABLES)
clean_tables:
	rm -f $(TABLES)

RESTRICTED_TABLES = $(foreach bn,$(TARGET_BNS),\
$(INT_DIR)/$(bn)_noshuffle_restricted_full_table.txt.gz)

restricted_tables: | $(RESTRICTED_TABLES)
clean_restricted_tables: 
	rm -f $(RESTRICTED_TABLES)

THREE_PRIME_CLEAVAGE_GRAPHS = $(foreach bn,$(TARGET_BNS),\
$(GRAPH_DIR)/$(bn)_threeprime_graph_length_9.png)
three_prime_cleavage_graphs: | $(THREE_PRIME_CLEAVAGE_GRAPHS)
clean_three_prime_cleavage_graphs:
	rm -f $(THREE_PRIME_CLEAVAGE_GRAPHS)

mapped_only_tables: | $(MAPPED_ONLY_TABLES)
clean_mapped_only_tables: 
	rm -f $(MAPPED_ONLY_TABLES)

RESTRICT_TO_ZERO := $(foreach bn,$(TARGET_BNS),$(foreach shuffle,noshuffle shuffled,$(INT_DIR)/$(bn)_$(shuffle)_restrict_to_zero_length_10.txt.gz))

restrict_to_zero: | $(RESTRICT_TO_ZERO)
clean_restrict_to_zero:
	rm -f $(RESTRICT_TO_ZERO)

########################## MAP TO CUSTOM BOWTIE INDEX ###########################
# % = PB2_5R
# output: PB2_5R_mapped_to_refseq_tss_50.bam

# Find all the alignments, and try extra-hard to find all of them
# But don't allow any mismatches; output to SAM (immediately convert to BAM)
# Don't map to the reverse complement of it; it's already in the right orientation

$(GRAPH_DIR)/%_graph_of_tss_length_10.eps: | \
$(INT_DIR)/%_noshuffle_restrict_to_zero_length_10.txt.gz $(INT_DIR)/%_shuffled_restrict_to_zero_length_10.txt.gz $(GRAPH_OF_BEDTOOLS_INTERSECT)
	$(BSUB) "$(MY_PYTHON) $(GRAPH_OF_BEDTOOLS_INTERSECT) --gzipped --use-counts -u $(word 1,$|) -s $(word 2,$|) -o $(GRAPH_DIR)/$*_graph_of_tss --lengths 10 15 --format eps"

# Intersect the noshuffled with the shuffled full table. 
$(INT_DIR)/%_noshuffle_restricted_full_table.txt.gz: | \
$(INT_DIR)/%_noshuffle_full_table.txt.gz $(INT_DIR)/%_shuffled_full_table.txt.gz \
$(INTERSECT_TABLES)
	$(BSUB) -R rusage[mem=$$((`du -m $(word 1,$|) | cut -f 1` + 100000))] "$(MY_PYTHON) $(INTERSECT_TABLES) -a $(word 1,$|) -b $(word 2,$|) --gzipped -ob /dev/null | gzip -9 > $@"

$(INT_DIR)/%_length_10_around_five.txt.gz: | $(INT_DIR)/%_full_table_no_paralogs.txt.gz $(SUBSET_TABLE)
	$(BSUB) -R rusage[mem=$$((`du -m $(word 1,$|) | cut -f 1` + 100000))] "zcat $(word 1,$|) | $(MY_PYTHON) $(SUBSET_TABLE) --restrict-to-zero --use-rank --min-length 10 --combine-sequences -8 9 --restrict-around -5 5 --three-prime-trimmed '' | gzip -9 > $@"

$(INT_DIR)/%_restrict_to_zero_length_10.txt.gz: | $(INT_DIR)/%_full_table_no_paralogs.txt.gz $(SUBSET_TABLE)
	$(BSUB) -R rusage[mem=$$((`du -m $(word 1,$|) | cut -f 1` + 100000))] "zcat $(word 1,$|) | $(MY_PYTHON) $(SUBSET_TABLE) --use-rank --min-length 10 --combine-sequences -8 9 --three-prime-trimmed '' -- | gzip -9 > $@"

$(INT_DIR)/%_full_table_no_paralogs.txt.gz: | $(INT_DIR)/%_full_table.txt.gz $(SUBSET_TABLE)
	$(BSUB) -R rusage[mem=$$((`du -m $(word 1,$|) | cut -f 1` + 100000))] "zcat $(word 1,$|) | $(MY_PYTHON) $(SUBSET_TABLE) --columns all --drop-paralogs | gzip -9 > $@"

# Combine the files
# Strip header from second file, then output first file, then append second file. 
$(INT_DIR)/%_full_table.txt.gz: \
$(INT_DIR)/%_mapped_only_table.txt.gz $(INT_DIR)/%_unmapped_only_table.txt.gz
	$(BSUB) "zcat $< > $(basename $<)"
	$(BSUB) "zcat $(word 2,$^) | tail -n +2 | cat $(basename $<) - | gzip -9 > $@"
	rm $(basename $<)

# Mapped reads to table
$(INT_DIR)/%_mapped_only_table.txt.gz: | \
$(INT_DIR)/%_mapped_to_hg19_intersect_with_gencode_preprocessed.txt.gz \
$(INT_DIR)/%_before_mapping.fastq.gz \
$(HG19_WHOLE_GENOME) \
$(INT_DIR)/%_mapped_to_hg19_all.bam \
$(DEPARSE_BEDTOOLS_OUTPUT)
	$(BSUB) "zcat $(word 2,$|) > $(basename $(word 2,$|))"
	$(BSUB) "-R rusage[mem=141000]" "zcat $(word 1,$|) | $(MY_PYTHON) $(DEPARSE_BEDTOOLS_OUTPUT) -r $(basename $(word 2,$|)) -f $(word 3,$|) | gzip -9 > $@"
	rm $(basename $(word 2,$|))

# Unmapped reads to table
# Use order-only prerequisites because it updates for some reason. 
$(INT_DIR)/%_unmapped_only_table.txt.gz: | \
$(INT_DIR)/%_mapped_to_hg19_all.bam $(BAM_TO_TABLE)
# request 120 GB of memory
	$(BSUB) "$(MY_PYTHON) $(BAM_TO_TABLE) -b $(word 1,$|) | gzip -9 > $@"

# First, preprocess the bedtools output
$(INT_DIR)/%_mapped_to_hg19_intersect_with_gencode_preprocessed.txt.gz: | \
$(INT_DIR)/%_mapped_to_hg19_intersect_with_gencode.bed.gz $(PREPROCESS_BEDTOOLS_OUTPUT)
	$(BSUB) "zcat $(word 1,$|) | $(MY_PYTHON) $(PREPROCESS_BEDTOOLS_OUTPUT) | gzip -9 > $@"

# Intersect with gencode TSS and keep only the reads that are closest
$(INT_DIR)/%_mapped_to_hg19_intersect_with_gencode.bed.gz: | \
$(INT_DIR)/%_mapped_to_hg19_all_human_coords.bed.gz $(ANNO_DIR)/gencode_offset_25_nogs.bed
	$(BSUB) "zcat $(word 1,$|) > $(basename $(word 1,$|))"
	$(BSUB) "bedtools intersect -loj -wo -s -a $(basename $(word 1,$|)) -b $(word 2,$|) | gzip -9 > $@ && rm $(basename $(word 1,$|))"

$(INT_DIR)/%_mapped_to_hg19_all_human_coords.bed.gz: | \
$(INT_DIR)/%_mapped_to_hg19_all_human_coords.bam
	$(BSUB) "bamToBed -i $(word 1,$|) | gzip -9 > $@"

$(INT_DIR)/%_mapped_to_hg19_all_human_coords.bam: | \
$(INT_DIR)/%_mapped_to_hg19_all.bam $(CUSTOM_BAM_TO_HUMAN_COORDS) $(BOWTIE_GENOME_DIR)
	$(BSUB) "$(MY_PYTHON) $(CUSTOM_BAM_TO_HUMAN_COORDS) -g $(BOWTIE_HG19_GENOME) $(word 1,$|) > $@"

# Map influenza reads to the genome in chunks
$(INT_DIR)/%_mapped_to_hg19_all.bam: | $(INT_DIR)/%_before_mapping.fastq.gz $(BOWTIE_DIR)/gencode_offset_25_nogs_merged.1.ebwt $(SPLIT_AND_MAP)
	$(BSUB) "zcat $(word 1,$|) > $(basename $(word 1,$|))"
	$(BSUB) "$(MY_PYTHON) $(SPLIT_AND_MAP) -n 10000 --bsub-command '$(BSUB)' --infile $(basename $(word 1,$|)) -o $@ --merge-memory 120000 --mapping-command 'bowtie -S --norc -p 24 --all --tryhard -v 0 $(BOWTIE_DIR)/gencode_offset_25_nogs_merged'"
	rm $(basename $(word 1,$|))
