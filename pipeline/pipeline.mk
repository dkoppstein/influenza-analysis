compressed: $(foreach target,$(BN_TARGETS),$(target)_compressed.fastq)

INT_DIR = intermediate
ANNO_DIR = anno

$(INT)/%_mapped_to_$(ANNO_TSS)_htseq_counted_aggregate.txt: \
$(INT)/%_mapped_to_$(ANNO_TSS)_htseq_counted.sam
	$(BSUB) "python $(AGGREGATE_HTSEQ_COUNTED_SAM) -i $< -o $@"

$(INT)/%_mapped_to_$(ANNO_TSS)_htseq_counted.sam: \
$(INT)/%_mapped_to_$(ANNO_TSS)_human_coords_sorted.bam \
$(TSS_DIR)/$(ANNO_TSS).gtf
	$(BSUB) "samtools view -h $< | htseq-count -m intersection-strict -i gene_name -o $@ - $(TSS_DIR)/$(ANNO_TSS).gtf > $@.txt"

############################## CONVERT TO BED FILE ##############################
# PB2_5R_mapped_to_refseq_tss_50_human_coords.bed

$(INTERMEDIATE)/%_human_coords.bed: $(INTERMEDIATE)/%_human_coords_sorted.bam
	$(BSUB) "bamToBed -i $< > $@"

################ CONVERT BAM FILE TO HUMAN COORDINATES AND SORT ################
# PB2_5R_mapped_to_refseq_tss_50_human_coords_sorted.bam

$(INTERMEDIATE)/%_human_coords_sorted.bam: $(INTERMEDIATE)/%_human_coords.bam
	$(BSUB) "samtools sort -o $< arbitraryprefix > $@"

$(INTERMEDIATE)/%_human_coords.bam: $(INTERMEDIATE)/%_aligned.bam
	$(BSUB) "python $(CUSTOM_BAM_TO_HUMAN_COORDS) --distances -g $(HG19_BOWTIE_INDEX) $< > $@"


########################## MAP TO CUSTOM BOWTIE INDEX ###########################
# % = PB2_5R
# output: PB2_5R_mapped_to_refseq_tss_50.bam

# Find all the alignments, and try extra-hard to find all of them
# But don't allow any mismatches; output to SAM (immediately convert to BAM)
# Don't map to the reverse complement of ; it's already in the right orientation

# -F0x4: only keep aligned reads

$(INT)/%_mapped_to_$(ANNO_TSS)_aligned.bam: $(INT)/%_mapped_to_$(ANNO_TSS)_all.bam
	$(BSUB) "samtools view -bh -F0x4 $< > $@"

$(INT)%_mapped_to_$(ANNO_TSS)_all.bam: $(call bowtie-target-pattern,$(ANNO_TSS)) \
%_compressed.fastq
	$(BSUB) "bowtie --all --norc --solexa1.3-quals --tryhard -v 0 -S $(BOWTIE_DIR)/$(ANNO_TSS) $(filter %.fastq,$^) | samtools view -bS - > $@"


%_reexpanded.fastq: %_compressed.fastq
	$(BSUB) "python $(EXPAND_FASTQ) $< > $@"

%_compressed.fastq: %_trimmed.fastq
	$(BSUB) "python $(COMPRESS_FASTQ_BY_COUNT) -i $< > $@"

########## FURTHER COLLAPSE PRIME-AND-REALIGNED SEQUENCES #####################
# % = PB2_5R
# filtered -> par
# filtered -> stats (OUTFILE)

PAR_PARAMS :=

# -c: use count= to guide statistics of prime-and-realignment
define par
python $(COLLAPSE_PAR) -c $(PAR_PARAMS) $< > $*_par.fastq 2> $*_par_stats.txt
endef

PAR_STATS = $(foreach bn,$(BASENAMES),$(bn)_par_stats.txt)

%_par.fastq %_par_stats.txt: %_filtered.fastq
	$(BSUB) "$(call par)"

$(PB2_TS_BN)_par.fastq $(PB2_TS_BN)_par_stats.txt: \
PAR_PARAMS := -s $(call extract-influenza-fasta,PB2)

$(PB2_5R_BN)_par.fastq $(PB2_5R_BN)_par_stats.txt: \
PAR_PARAMS := -s $(call extract-influenza-fasta,PB2)

$(NS1_TS_BN)_par.fastq $(NS1_TS_BN)_par_stats.txt: \
PAR_PARAMS := -s $(call extract-influenza-fasta,NS1)

$(NS1_5R_BN)_par.fastq $(NS1_5R_BN)_par_stats.txt: \
PAR_PARAMS := -s $(call extract-influenza-fasta,NS1)


######### TRIM OFF SEQUENCES THAT ARE TOO SHORT OR TOO LONG #############

%_trimmed.fastq: %_qfiltered.fastq
	$(BSUB) "python $(DISCARD_FASTQ) -l $(MIN_NUCS) -i $(MAX_NUCS) $< > $@"

# Get rid of sequences with Ns or less than QUALITY_CUTOFF

%_qfiltered.fastq: %_noillumina.fastq
	$(BSUB) "python $(FASTQ_QUALITY_FILTER) --no-ns --min-quality $(QUALITY_CUTOFF) -i $< > $@"

####### DISCARD READS THAT MATCH THE ILLUMINA SEQUENCE ##################

%_noillumina.fastq: %_vfiltered.fastq
	$(BSUB) "cutadapt --quality-base=64 -f fastq -e 0 --discard -O $(ILLUMINA_OVERLAP_REQD) -b $(ILLUMINA_SEQUENCE) $< > $@"

####### EXCISE INFLUENZA SEQUENCES AND ONLY KEEP READS THAT MATCH #######

# % = PB2_5R
# collapsed -> filtered
FILTER_PARAMS :=

# Function to use cutadapt to sort files based on whether they have the influenza
# sequence or not
define cut-adapt
cutadapt --trimmed-only --quality-base=64 -O \
$(CUTADAPT_OVERLAP_REQUIRED) -f fastq --suffix=$(CUTADAPT_SUFFIX)$(1) -a \
$(call extract-influenza-fasta,$(1)) $< > $@
endef

%_vfiltered.fastq: %_collapsed.fastq
	$(BSUB) "$(call cut-adapt,$(FILTER_PARAMS))"

$(PB2_TS_BN)_vfiltered.fastq: FILTER_PARAMS := PB2
$(NS1_TS_BN)_vfiltered.fastq: FILTER_PARAMS := NS1
$(PB2_5R_BN)_vfiltered.fastq: FILTER_PARAMS := PB2
$(NS1_5R_BN)_vfiltered.fastq: FILTER_PARAMS := NS1
$(PB2_TS_BN)_nc_vfiltered.fastq: FILTER_PARAMS := PB2
$(NS1_TS_BN)_nc_vfiltered.fastq: FILTER_PARAMS := NS1
$(PB2_5R_BN)_nc_vfiltered.fastq: FILTER_PARAMS := PB2
$(NS1_5R_BN)_nc_vfiltered.fastq: FILTER_PARAMS := NS1


# In parallel, just trim the 5' barcode, but don't collapse (to compare)
%_nc_collapsed.fastq: %_data.fastq
	$(BSUB) "python $(COUNT_TRANSCRIPTS) --no-ns --no-collapse $(COUNT_PARAMS) $< > $@"

$(PB2_TS_BN)_nc_collapsed.fastq: COUNT_PARAMS := -g
$(PB2_TS_BN)_nc_collapsed.fastq: COUNT_PARAMS := -g

%_barcodes_equalized.eps: %_barcodes_equalized_stats.txt
	$(BSUB) "python $(PLOT_NUCLEOTIDE_DIST) -f eps $< -o $@"

# covers both three prime and five prime
%_barcodes_equalized_stats.txt: %_barcodes_equalized.fasta
	$(BSUB) "fastx_quality_stats -N -i $< -o $@"

# Equalize barcode file
%_barcodes_equalized.fasta: %_barcodes.fasta
	$(BSUB) "python $(EQUALIZE_LENGTHS) -f fasta --three-prime -i $< -o $@"

# Collapse transcripts by barcode
%_collapsed.fastq %_barcodes.fasta: %_data.fastq
	$(BSUB) "python $(COUNT_TRANSCRIPTS) --no-ns --barcode-file $*_barcodes.fasta $(COUNT_PARAMS) $< > $*_collapsed.fastq"

# Additionally trim guanines if using template-switching strategy
$(PB2_TS_BN)_collapsed.fastq $(PB2_TS_BN)_barcodes.fasta: COUNT_PARAMS := -g
$(NS1_TS_BN)_collapsed.fastq $(NS1_TS_BN)_barcodes.fasta: COUNT_PARAMS := -g
