# This file deals specifically with processing RNA-PET data. 
# The RNA-PET data is used for providing a high-confidence subset of influenza reads that map to 
# bona fide five prime ends of transcripts. 
# 
# After extracting the 5' end reads, trimming off the adaptor, and reverse complementing,
# the suffix "before_g_strip.fastq" is used and the file is passed on to five_prime_library.mk, 
# where it is mapped to the genome and additionally processed. 

# Reverse complement
$(ANNO_DIR)/$(RNAPET_PREFIX)_%_before_g_strip.fastq: \
$(ANNO_DIR)/$(RNAPET_PREFIX)_%_trimmed.fastq
	$(BSUB) "fastx_reverse_complement -i $< -o $@"

# Remove the adaptor
$(ANNO_DIR)/$(RNAPET_PREFIX)_%_trimmed.fastq: \
$(ANNO_DIR)/$(RNAPET_PREFIX)_%_5prime.fastq
	$(BSUB) "cutadapt --trimmed-only -a CTGCTGATGG --overlap=8 -f fastq --minimum-length=15 $< > $@"

# Find the reads corresponding to the 5' end of the transcript, as 
# described in Ruan and Ruan, Methods Mol. Biol., 2012. 
$(ANNO_DIR)/$(RNAPET_PREFIX)_%_5prime.fastq: \
$(DATA_DIR)/$(RNAPET_PREFIX)_%_1.fastq $(DATA_DIR)/$(RNAPET_PREFIX)_%_2.fastq $(PROCESS_PET)
	$(BSUB) "python $(PROCESS_PET) -i $(filter %.fastq,$^) > $@"
