######## COMMON FLAGS FOR ALL MAKE FILES ##########

# Run jobs in parallel
MAKEFLAGS = -j

# -K causes the bsub command to wait until the job is finished to return
# -o and -e redirect stdout and stderr
BSUB = bsub -K -q bartel -o /lab/bartel1_ata/koppstein/bjob_output -e /lab/bartel1_ata/koppstein/bjob_stderr

# For home vs. Tak compilation
UNAME := $(shell uname)

######## USEFUL FUNCTIONS #######################

# Link a file from $(2) to $(1) -- useful for "pulling into the directory" from
# locations on the filesystem
define link-file
$(1): $(2)
	mkdir -p $$(dir $$@)
# if the file isn't already there, then make a symlink
	if [ ! -L $(1) ]; then ln -s $$< $$@; fi;
endef

# Extract the given sequence from the influenza transcripts fasta file
define extract-influenza-fasta
$(shell $(EXTRACT_FASTA) -d $(1) -n $(EXTRACT_FASTA_NUCS) \
-i $(INFLUENZA_TRANSCRIPTS))
endef

# These are the files made by bowtie-build
define bowtie-target-pattern
$(foreach base,$(BOWTIE_BUILD_BASE),$(BOWTIE_DIR)/$(1)$(base))
endef

####### COMMON DIRECTORY PATHS ##########

SOLEXA_BARTEL := /lab/solexa_public/Bartel
DATA_DIR := data
INT_DIR := intermediate
ANNO_DIR := anno
GRAPH_DIR := graphs
BOWTIE_DIR := $(ANNO_DIR)/bowtie_indices

BJOB_OUTPUT = /lab/bartel1_ata/koppstein/bjob_output
BJOB_STDERR = /lab/bartel1_ata/koppstein/bjob_stderr

# Don't delete intermediate files
.SECONDARY:

# You can use these sub-commands to generate intermediates
.PHONY: get_data clean_to_collapsed clean_data collapsed clean_collapsed vfiltered clean_vfiltered trzdimmed cleantrimmed trimmed_summary qfiltered cleanqfiltered par clean_par distributions cleandistributions cage_hmm clean_to_collapsed test_command par_5R_nogs shuffle clean_shuffle rnapet_targets clean_rnapet_targets remove_empty five_prime_targets clean_five_prime_targets influenza_rnaseq clean_influenza_rnaseq clean_par_graphs bedtools_graph_targets clean_bedtools_graph_targets tables clean_tables mapped_only_tables clean_mapped_only_tables paper clean_paper

# Deletes everything in the intermediate dir that's not the collapsed fastq file 
# (as this takes a long time to generate)
# also deletes broken symlinks
clean_to_collapsed:
	cd $(INT_DIR) && find . -type f ! -name "*_collapsed.fastq.gz" -delete && find -L . -type l -delete && cd ..

# Remove files that are empty (usually Makefiles automatically delete these;
# unfortunately, calls to bsub complicate things, so you sometimes need to delete
# the empty filehandles by hand)
remove_empty:
	for DIRECTORY in $(INT_DIR) $(GRAPH_DIR) $(ANNO_DIR); do cd $$DIRECTORY && find . -maxdepth 1 -type f -empty -exec rm {} \; && cd ..; done;

remove_almost_empty:
	for DIRECTORY in $(INT_DIR) $(GRAPH_DIR) $(ANNO_DIR) $(PELCHAT_DIR); do cd $$DIRECTORY && find . -name "*.gz" -size -2 -delete && cd ..; done;

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  getdata   to get data from GEO, link files on TAK, and unzip FASTQ"
	@echo "  cleandata to remove the data directory"
	@echo "  collapsed to collapse the initial reads"
	@echo "  clean_collapsed to get rid of the collapsed files "
	@echo "  vfiltered"
	@echo " "
	@echo "  All of the intermediates can be preceded by clean_"
	@echo "  to remove those intermediate files specifically. "

###### SYMLINK THE SCRIPTS ########

include pipeline/scripts.mk

# has pandas 0.13.1 installed

ifeq ($(UNAME),Linux)
MY_PYTHON := /lab/bartel1_ata/koppstein/virtualenvs/my2.7/bin/python
endif
ifeq ($(UNAME),Darwin)
MY_PYTHON := python
endif

ifeq ($(UNAME),Linux)
MY_IPYTHON := /lab/bartel1_ata/koppstein/virtualenvs/my2.7/bin/ipython
endif
ifeq ($(UNAME),Darwin)
MY_IPYTHON := ipython
endif

MY_WEBLOGO := /lab/bartel1_ata/koppstein/virtualenvs/my2.7/bin/weblogo

####### CONSTANTS ########

# Base names of the bowtie-build executable output files
BOWTIE_BUILD_BASE = .1.ebwt .2.ebwt .3.ebwt .4.ebwt .rev.1.ebwt .rev.2.ebwt

####### SEQUENCING DATA BASE NAMES #######

# All second-generation time course datasets
GEN2_TIME_COURSE_BNS := NS1_30 NS1_45 NS1_60 NS1_90 NS1_120 NS1_240

# All second-generation template-switching datasets (strand-specific)
GEN2_TS_BNS := HA_TS MP_TS NA_TS NP_TS PA_TS PB1_TS

# All 5'-RACE-like datasets
GEN1_5R_BNS := PB2_5R NS1_5R

# All first-generation template-switching datasets (paired end 200x200)
GEN1_TS_BNS := PB2_TS NS1_TS

# All 5'-RACE-like datasets; will get guanosines trimmed in parallel 
# for fair comparison with the template-switching datasets
GEN1_5R_NOGS_BNS := PB2_RACE_NOGS NS1_RACE_NOGS

# All template-switching datasets
TS_BNS := $(GEN2_TIME_COURSE_BNS) $(GEN2_TS_BNS) $(GEN1_TS_BNS)

INFLUENZA_GENES := HA MP NA NP NS1 PA PB1 PB2

# Datasets sorted by gene name
HA_BNS := HA_TS 
MP_BNS := MP_TS 
NA_BNS := NA_TS 
NP_BNS := NP_TS 
NS1_BNS := NS1_TS $(GEN2_TIME_COURSE_BNS) NS1_5R NS1_RACE_NOGS
PA_BNS := PA_TS 
PB1_BNS := PB1_TS 
PB2_BNS := PB2_TS PB2_5R PB2_RACE_NOGS

GCAAAAGCAG_BNS := $(NS1_BNS) $(HA_BNS) $(NP_BNS) $(NA_BNS) NS1_TS_5GTRIMMED HA_TS_5GTRIMMED NP_TS_5GTRIMMED NA_TS_5GTRIMMED $(foreach bn,$(GEN2_TIME_COURSE_BNS),$(bn)_5GTRIMMED) $(foreach rep,1 2,PELCHAT_HA_rep_$(rep) PELCHAT_MP_rep_$(rep) PELCHAT_NA_rep_$(rep) PELCHAT_NP_rep_$(rep) PELCHAT_NS1_rep_$(rep) PELCHAT_PA_rep_$(rep)) DUMMY_NS1
GCGAAAGCAG_BNS := $(PB2_BNS) $(PB1_BNS) $(PA_BNS) $(MP_BNS) PB2_TS_5GTRIMMED PB1_TS_5GTRIMMED PA_TS_5GTRIMMED MP_TS_5GTRIMMED $(foreach rep,1 2,PELCHAT_PB1_rep_$(rep) PELCHAT_PB2_rep_$(rep))

PELCHAT_BNS := $(foreach gene,$(INFLUENZA_GENES),\
$(foreach rep,1 2,\
PELCHAT_$(gene)_rep_$(rep)))

PELCHAT_BNS_REAL := PELCHAT_HA_rep_2 PELCHAT_MP_rep_2 PELCHAT_NA_rep_1 PELCHAT_NP_rep_1 PELCHAT_NS1_rep_2 PELCHAT_PA_rep_1 PELCHAT_PB1_rep_2 PELCHAT_PB2_rep_1

RNASEQ_BNS = RNASEQ_MOCK RNASEQ_30 RNASEQ_45 RNASEQ_60 RNASEQ_90 RNASEQ_120 RNASEQ_240 RNASEQ_ZERO

# All real datasets
TARGET_BNS = $(foreach gene,$(INFLUENZA_GENES),$($(gene)_BNS))

# Make the paper

paper: paper.docx

clean_paper:
	rm -rf paper/*.aux paper/*.depytx paper/*.log paper/*.pdf paper/*.pytxcode paper/pythontex-files-paper paper/temp.tex paper.docx paper/temp_nohead.tex paper/figure/all_length_distributions.eps

paper.docx: paper/paper.tex
# $(MY_IPYTHON) nbconvert $< --to latex --SphinxTransformer.author='David Koppstein'
#	cd paper && pdflatex $(notdir $<) && $(MY_PYTHON) ../$(PYTHONTEX) $(notdir $<) && pdflatex $(notdir $<) && $(MY_PYTHON) ../$(DEPYTHONTEX) $(notdir $<) -o temp.tex && tail -n +2 temp.tex > temp_nohead.tex && cd ..
#	
	rm -f paper/final_papers.bib
	cat paper/papers.bib paper/extras.bib >> paper/final_papers.bib
	cd paper && pandoc --bibliography=final_papers.bib --csl=nar_revised.csl -t docx -o paper.docx paper.tex && cd ..
	mv paper/paper.docx .
	cp -f paper.docx ~/science/software/dklib/influenza/compiled_paper.docx

###### GET THE DATA #########

include pipeline/get_data.mk

###### COLLAPSE FASTQ FILES ########

include pipeline/collapse_fastq.mk

###### KEEP ONLY VIRAL SEQUENCES ########

include pipeline/vfiltered.mk

###### READ QUALITY FILTERING #########

include pipeline/qfiltered.mk

###### TRIM THE VIRAL SEQUENCES #######

# for now, don't trim; we want to see 
# include pipeline/trim.mk

# FURTHER COLLAPSE PRIME-AND-REALIGNED SEQUENCES

include pipeline/par.mk

# GRAPH EVERYTHING

# include pipeline/graphs.mk

# FIND ALL TRANSCRIPTION STARTS SITES

include pipeline/find_starts.mk

# DO CAGE DATA

include pipeline/cage.mk

# MAP READS TO PET DATA

include pipeline/pet.mk

# MAP READS TO TRANSCRIPTION START SITES

include pipeline/five_prime_library.mk

# MAP RNA-SEQ DATA TO THE GENOME
include pipeline/rnaseq.mk

# MAP INFLUENZA DIRECTLY TO THE GENOME
include pipeline/map.mk

# MAKE FIGURES
include pipeline/figure.mk

# GROSEQ ANALYSIS
include pipeline/groseq.mk

# PELCHAT ANALYSIS
include pipeline/pelchat.mk
