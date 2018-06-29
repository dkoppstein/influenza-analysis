
######## COMMON FLAGS FOR ALL MK FILES ##########

MAKEFLAGS = -j
BSUB = bsub -K -q bartel

# Don't delete intermediate files
.SECONDARY:

.PHONY: getdata cleandata collapsed cleancollapsed vfiltered cleanvfiltered trimmed cleantrimmed trimmed_summary qfiltered cleanqfiltered par cleanpar distributions cleandistributions cage_hmm clean_to_collapsed test_command

clean_to_collapsed:
	cd $(INT_DIR) && find . -type f ! -name "*_collapsed.fastq" -delete && cd ..

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  getdata   to get data from GEO, link files on TAK, and unzip FASTQ"
	@echo "  cleandata to remove the data directory"
	@echo "  collapsed to collapse the initial reads"
	@echo "  cleancollapsed to get rid of the collapsed files "

####### COMMON DIRECTORY PATHS ##########

SOLEXA_BARTEL = /lab/solexa_public/Bartel
DATA_DIR = data
INT_DIR = intermediate
ANNO_DIR = anno
GRAPH_DIR = graphs
TSS_BOWTIE_DIR = $(ANNO_DIR)/bowtie_indices

####### CONSTANTS ########

# Base names of the bowtie-build executable output files
BOWTIE_BUILD_BASE = .1.ebwt .2.ebwt .3.ebwt .4.ebwt .rev.1.ebwt .rev.2.ebwt

####### SEQUENCING DATA BASE NAMES #######

NS1_30_BN := NS1_30min
NS1_45_BN := NS1_45min
NS1_60_BN := NS1_60min
NS1_90_BN := NS1_90min
NS1_120_BN := NS1_120min
NS1_240_BN := NS1_240min
PB2_5R_BN := PB2_5R
PB2_TS_BN := PB2_TS
NS1_5R_BN := NS1_5R
NS1_TS_BN := NS1_TS

GEN2_TIME_COURSE_BNS = $(NS1_30_BN) $(NS1_45_BN) $(NS1_60_BN) $(NS1_90_BN) \
$(NS1_120_BN) $(NS1_240_BN)

GEN1_5R_BNS = $(PB2_5R_BN) $(NS1_5R_BN)

GEN1_TS_BNS = $(PB2_TS_BN) $(NS1_TS_BN)

TARGET_BNS = $(GEN2_TIME_COURSE_BNS) $(GEN1_5R_BNS) $(GEN1_TS_BNS)

######## USEFUL FUNCTIONS #######################

# Link a file from $(2) to $(1) -- useful for "pulling into the directory" from
# locations on the filesystem
define link-file
$(1): $(2)
	mkdir -p $$(dir $$@)
	ln -s $$< $$@
endef

# Extract the given sequence from the fasta file
define extract-influenza-fasta
$(shell $(EXTRACT_FASTA) -i $(1) -n $(EXTRACT_FASTA_NUCS) \
$(INFLUENZA_SEQUENCES))
endef

# make the targets of bowtie-build with pattern matching on the fly
define bowtie-target-pattern
$(foreach base,$(BOWTIE_BUILD_BASE),$(TSS_BOWTIE_DIR)/$(1)$(base))
endef

###### SYMLINK THE SCRIPTS ########

include scripts.mk

###### GET THE DATA #########

include get_data.mk

###### COLLAPSE FASTQ FILES ########

include collapse_fastq.mk

###### KEEP ONLY VIRAL SEQUENCES ########

include vfiltered.mk

###### REMOVE CAP SEQUENCES WITH NS #########

include qfiltered.mk

###### TRIM THE VIRAL SEQUENCES #######

include trim.mk

# FURTHER COLLAPSE PRIME-AND-REALIGNED SEQUENCES

include par.mk

# GRAPH EVERYTHING

include graphs.mk

# DO CAGE DATA

include cage.mk