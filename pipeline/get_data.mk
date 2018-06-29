###############################################
############### ALL TARGETS ###################
###############################################

# PREFIXES AND DIRECTORIES
RNAPET_PREFIX := GIS_RnaPet_A549_longPolyA
CAGE_PREFIX := RIKEN_Cage_A549_longPolyA
COMPARTMENTS = cytosol nucleus cell
REPS = rep1 rep2

# TARGETS
CAGE_TARGETS = \
$(foreach compartment,$(COMPARTMENTS),\
$(foreach rep,$(REPS),\
$(DATA_DIR)/$(CAGE_PREFIX)_$(compartment)_$(rep).fastq))

RNAPET_TARGETS = \
$(foreach rep,$(REPS),\
$(foreach compartment,$(COMPARTMENTS),\
$(DATA_DIR)/$(RNAPET_PREFIX)_$(compartment)_$(rep)_1.fastq))

LINK_TARGETS = $(foreach bn,$(TARGET_BNS) $(RNASEQ_BNS),\
$(DATA_DIR)/$(bn)_data.fastq.gz)

###############################################
############### MAKE RULES ####################
###############################################

get_data: $(CAGE_TARGETS) $(RNAPET_TARGETS) $(LINK_TARGETS)

clean_data: 
	rm -rf $(DATA_DIR)/*.fastq.gz

###############################################
######## GET AND PROCESS GEO DATA #############
###############################################

# Get data from the internet. 
# Example: $(eval $(call wget-target,$(DATA_DIR)/RNAPet/GIS_RnaPet_A549_cytosol_longPolyArep1.fastq,www.google.com
# fetches data from $(2) and puts it in $(1)
define wget-target
$(1):
	mkdir -p $$(dir $$@)
	$(BSUB) "wget -O $$@ $(2)"
endef

############## SRA -> FASTQ ###################

FASTQ_DUMP := fastq-dump --outdir $(DATA_DIR)

# Split the RNA-PET files but not the CAGE files
$(DATA_DIR)/$(RNAPET_PREFIX)%_1.fastq $(DATA_DIR)/$(RNAPET_PREFIX)%_2.fastq : $(DATA_DIR)/$(RNAPET_PREFIX)%.sra
	$(BSUB) "$(FASTQ_DUMP) --split-files `pwd`/$<"

$(DATA_DIR)/$(CAGE_PREFIX)_%.fastq: $(DATA_DIR)/$(CAGE_PREFIX)_%.sra
	$(BSUB) "$(FASTQ_DUMP) `pwd`/$<"

####### GET RNA PET DATA #######
$(eval $(call wget-target,\
$(DATA_DIR)/$(RNAPET_PREFIX)_cytosol_rep1.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX188%2FSRX188980/SRR575156/SRR575156.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(RNAPET_PREFIX)_cytosol_rep2.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX188%2FSRX188980/SRR575157/SRR575157.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(RNAPET_PREFIX)_cell_rep1.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX188%2FSRX188973/SRR575142/SRR575142.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(RNAPET_PREFIX)_cell_rep2.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX188%2FSRX188973/SRR575143/SRR575143.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(RNAPET_PREFIX)_nucleus_rep1.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX188%2FSRX188974/SRR575144/SRR575144.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(RNAPET_PREFIX)_nucleus_rep2.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX188%2FSRX188974/SRR575145/SRR575145.sra))

###### GET PELCHAT DATA ######
$(eval $(call wget-target,\
$(DATA_DIR)/pelchat_raw_data.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR130/SRR1301043/SRR1301043.sra))

####### GET CAGE DATA #########
$(eval $(call wget-target,\
$(DATA_DIR)/$(CAGE_PREFIX)_cell_rep1.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX111%2FSRX111957/SRR390526/SRR390526.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(CAGE_PREFIX)_cell_rep2.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX111%2FSRX111957/SRR390527/SRR390527.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(CAGE_PREFIX)_cytosol_rep1.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX172%2FSRX172609/SRR530906/SRR530906.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(CAGE_PREFIX)_cytosol_rep2.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX172%2FSRX172609/SRR530907/SRR530907.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(CAGE_PREFIX)_nucleus_rep1.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX172%2FSRX172610/SRR530908/SRR530908.sra))

$(eval $(call wget-target,\
$(DATA_DIR)/$(CAGE_PREFIX)_nucleus_rep2.sra,\
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX172%2FSRX172610/SRR530909/SRR530909.sra))

# Gzipped bedRNAElements file
$(ANNO_DIR)/%.bedRnaElements: $(ANNO_DIR)/%.bedRnaElements.gz
	$(BSUB) "zcat $< > $@; rm -f $<"

$(eval $(call wget-target,\
$(ANNO_DIR)/$(CAGE_PREFIX)$(CELL)$(HMM).bedRnaElements.gz,\
ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM849nnn/GSM849359/suppl/GSM849359%5Fhg19%5FwgEncodeRikenCageA549CellPapTssHmm%2EbedRnaElements%2Egz))

#########################################################
########## SYMLINK GZIPPED DATA AND UNZIP################
#########################################################

GZIP_EXT = _data.txt.tar.gz

# Unzip the fastq
$(DATA_DIR)/%_data.fastq.gz: $(DATA_DIR)/%$(GZIP_EXT) $(UNZIP_FASTQ)
	$(BSUB) "python $(UNZIP_FASTQ) -i $< | gzip -9 > $@"

# First generation data
PB2_5R_DATA := $(SOLEXA_BARTEL)/130724_WIGTC-HISEQA_D29RNACXX/QualityScore/ATCACG-s_8_1_sequence.txt.tar.gz
PB2_TS_DATA := $(SOLEXA_BARTEL)/130802_WIGTC-HISEQA_C1YLGACXX/QualityScore/ATCACG-s_7_1_sequence.txt.tar.gz
NS1_5R_DATA := $(SOLEXA_BARTEL)/130724_WIGTC-HISEQA_D29RNACXX/QualityScore/ACTTGA-s_8_1_sequence.txt.tar.gz
NS1_TS_DATA := $(SOLEXA_BARTEL)/130802_WIGTC-HISEQA_C1YLGACXX/QualityScore/ACTTGA-s_7_1_sequence.txt.tar.gz

# Time course data (second generation TS)
NS1_30_DATA := $(SOLEXA_BARTEL)/131231_WIGTC-HISEQA_C3JTPACXX/QualityScore/GTGAAA-s_6_1_sequence.txt.tar.gz
NS1_45_DATA := $(SOLEXA_BARTEL)/131231_WIGTC-HISEQA_C3JTPACXX/QualityScore/GTGGCC-s_6_1_sequence.txt.tar.gz
NS1_60_DATA := $(SOLEXA_BARTEL)/131231_WIGTC-HISEQA_C3JTPACXX/QualityScore/GTTTCG-s_6_1_sequence.txt.tar.gz
NS1_90_DATA := $(SOLEXA_BARTEL)/131231_WIGTC-HISEQA_C3JTPACXX/QualityScore/CGTACG-s_6_1_sequence.txt.tar.gz
NS1_120_DATA := $(SOLEXA_BARTEL)/131231_WIGTC-HISEQA_C3JTPACXX/QualityScore/GAGTGG-s_6_1_sequence.txt.tar.gz
NS1_240_DATA := $(SOLEXA_BARTEL)/131231_WIGTC-HISEQA_C3JTPACXX/QualityScore/GGTAGC-s_6_1_sequence.txt.tar.gz

# Second generation TS data
HA_TS_DATA := $(SOLEXA_BARTEL)/140203_WIGTC-HISEQB_C3E3MACXX/QualityScore/TGACCA-s_7_1_sequence.txt.tar.gz
MP_TS_DATA := $(SOLEXA_BARTEL)/140203_WIGTC-HISEQB_C3E3MACXX/QualityScore/CAGATC-s_7_1_sequence.txt.tar.gz
NA_TS_DATA := $(SOLEXA_BARTEL)/140203_WIGTC-HISEQB_C3E3MACXX/QualityScore/GCCAAT-s_7_1_sequence.txt.tar.gz
NP_TS_DATA := $(SOLEXA_BARTEL)/140203_WIGTC-HISEQB_C3E3MACXX/QualityScore/ACAGTG-s_7_1_sequence.txt.tar.gz
PA_TS_DATA := $(SOLEXA_BARTEL)/140203_WIGTC-HISEQB_C3E3MACXX/QualityScore/TTAGGC-s_7_1_sequence.txt.tar.gz
PB1_TS_DATA := $(SOLEXA_BARTEL)/140203_WIGTC-HISEQB_C3E3MACXX/QualityScore/CGATGT-s_7_1_sequence.txt.tar.gz

# RNA-seq data
RNASEQ_ADAPTOR = TCGTATGCCGTCTTCTGCTTG
RNASEQ_MOCK_DATA := $(SOLEXA_BARTEL)/140226_WIGTC-HISEQB_C44CAACXX/QualityScore/TGACCA-s_6_1_sequence.txt.tar.gz
RNASEQ_30_DATA := $(SOLEXA_BARTEL)/140226_WIGTC-HISEQB_C44CAACXX/QualityScore/ACAGTG-s_6_1_sequence.txt.tar.gz
RNASEQ_45_DATA := $(SOLEXA_BARTEL)/140226_WIGTC-HISEQB_C44CAACXX/QualityScore/GCCAAT-s_6_1_sequence.txt.tar.gz
RNASEQ_60_DATA := $(SOLEXA_BARTEL)/140226_WIGTC-HISEQB_C44CAACXX/QualityScore/CAGATC-s_6_1_sequence.txt.tar.gz
RNASEQ_90_DATA := $(SOLEXA_BARTEL)/140226_WIGTC-HISEQB_C44CAACXX/QualityScore/CTTGTA-s_6_1_sequence.txt.tar.gz
RNASEQ_120_DATA := $(SOLEXA_BARTEL)/140226_WIGTC-HISEQB_C44CAACXX/QualityScore/ATCACG-s_6_1_sequence.txt.tar.gz
RNASEQ_240_DATA := $(SOLEXA_BARTEL)/140226_WIGTC-HISEQB_C44CAACXX/QualityScore/TTAGGC-s_6_1_sequence.txt.tar.gz
RNASEQ_ZERO_DATA := $(SOLEXA_BARTEL)/140226_WIGTC-HISEQB_C44CAACXX/QualityScore/ACTTGA-s_6_1_sequence.txt.tar.gz

# No Gs 5' RACE datasets (symlink to the same place)
NS1_RACE_NOGS_DATA := $(NS1_5R_DATA)
PB2_RACE_NOGS_DATA := $(PB2_5R_DATA)

# Symlink all the files to DATA_DIR
$(foreach bn,$(TARGET_BNS) $(RNASEQ_BNS),\
$(eval $(call link-file,$(DATA_DIR)/$(bn)$(GZIP_EXT),$($(bn)_DATA))))

$(DATA_DIR)/DUMMY_NS1_data.fastq.gz: $(DATA_DIR)/NS1_30_data.fastq.gz
	$(BSUB) "zcat $< | head -100 | gzip -9 > $@"

# Get the FASTA file with all the 5' ends of the influenza mRNA sequences
INFLUENZA_TRANSCRIPTS = $(DATA_DIR)/influenza_transcripts.fa
INFLUENZA_TRANSCRIPTS_LOCATION = $(INFLUENZA_SCRIPTS_DIR)/influenza_transcripts.fa
$(eval $(call link-file,$(INFLUENZA_TRANSCRIPTS),$(INFLUENZA_TRANSCRIPTS_LOCATION)))

####### GET ANNOTATIONS #########

# GENCODE annotation
GENCODE_LOCATION = /nfs/genomes/human_gp_feb_09_no_random/gtf/gencode.v17.annotation.gtf
GENCODE_ANNO = $(ANNO_DIR)/gencode.gtf
$(eval $(call link-file,$(GENCODE_ANNO),$(GENCODE_LOCATION)))

# Use the no-random HG19 genome
HG19_WHOLE_GENOME_LOCATION = /nfs/genomes/human_gp_feb_09_no_random/fasta_whole_genome/hg19.fa
HG19_WHOLE_GENOME = $(ANNO_DIR)/hg19.fa
$(eval $(call link-file,$(HG19_WHOLE_GENOME),$(HG19_WHOLE_GENOME_LOCATION)))

# Symlink the HG19 bowtie index and files
BOWTIE_GENOME_LOCATION = /nfs/genomes/human_gp_feb_09_no_random/bowtie
BOWTIE_GENOME_DIR = $(ANNO_DIR)/bowtie_hg19
BOWTIE_HG19_GENOME = $(BOWTIE_GENOME_DIR)/hg19
$(eval $(call link-file,$(BOWTIE_GENOME_DIR),$(BOWTIE_GENOME_LOCATION)))

INFLUENZA_SEQUENCES = $(DATA_DIR)/influenza_sequences.fa
INFLUENZA_SEQUENCES_LOCATION = $(INFLUENZA_SCRIPTS_DIR)/influenza_sequences.fa
$(eval $(call link-file,$(INFLUENZA_SEQUENCES),$(INFLUENZA_SEQUENCES_LOCATION)))
