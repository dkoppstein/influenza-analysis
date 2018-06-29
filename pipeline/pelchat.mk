
PELCHAT_DIR = pelchat

PELCHAT_VFILTERED = \
$(foreach gene,$(INFLUENZA_GENES),\
$(foreach rep,1 2,\
$(INT_DIR)/PELCHAT_$(gene)_rep_$(rep)_vfiltered.fastq.gz))

pelchat_vfiltered: $(PELCHAT_VFILTERED)
clean_pelchat_vfiltered:
	rm -f $(PELCHAT_VFILTERED)

$(PELCHAT_DIR)/pelchat_raw_data.fastq.gz: $(DATA_DIR)/pelchat_raw_data.sra
	mkdir -p $(PELCHAT_DIR)
	$(BSUB) "fastq-dump `pwd`/$< --outdir $(PELCHAT_DIR)"
	$(BSUB) "gzip -c $(PELCHAT_DIR)/pelchat_raw_data.fastq > $@"
	rm -f $(PELCHAT_DIR)/pelchat_raw_data.fastq

$(PELCHAT_DIR)/A549_rep_1.fastq.gz: $(PELCHAT_DIR)/pelchat_raw_data.fastq.gz
	$(BSUB) "zcat $< | cutadapt -g ^GCAT --trimmed-only --quality-base=64 - | gzip -9 > $@"

$(PELCHAT_DIR)/A549_rep_2.fastq.gz: $(PELCHAT_DIR)/pelchat_raw_data.fastq.gz
	$(BSUB) "zcat $< | cutadapt -g ^CATG --trimmed-only --quality-base=64 - | gzip -9 > $@"

$(PELCHAT_DIR)/A549_rep_%_5ptrimmed_rc.fastq.gz: $(PELCHAT_DIR)/A549_rep_%.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a TTTCCTGCAGGCGGCCGCTCCA -O 10 --trimmed-only --minimum-length=9 --quality-base=64 - | fastx_reverse_complement -Q33 | gzip -9 > $@"

$(PELCHAT_DIR)/HA_rep_%.fastq.gz: $(PELCHAT_DIR)/A549_rep_%_5ptrimmed_rc.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a GGATAATTCTA --trimmed-only --quality-base=64 - | gzip -9 > $@"

$(PELCHAT_DIR)/NA_rep_%.fastq.gz: $(PELCHAT_DIR)/A549_rep_%_5ptrimmed_rc.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a TGAGTGACATC --trimmed-only --quality-base=64 - | gzip -9 > $@"

$(PELCHAT_DIR)/MP_rep_%.fastq.gz: $(PELCHAT_DIR)/A549_rep_%_5ptrimmed_rc.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a AGTGAAAATGA --trimmed-only --quality-base=64 - | gzip -9 > $@"

$(PELCHAT_DIR)/NS1_rep_%.fastq.gz: $(PELCHAT_DIR)/A549_rep_%_5ptrimmed_rc.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a GTGACAAAGAC --trimmed-only --quality-base=64 - | gzip -9 > $@"

$(PELCHAT_DIR)/PB1_rep_%.fastq.gz: $(PELCHAT_DIR)/A549_rep_%_5ptrimmed_rc.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a CAAACCATTTG --trimmed-only --quality-base=64 - | gzip -9 > $@"

$(PELCHAT_DIR)/PB2_rep_%.fastq.gz: $(PELCHAT_DIR)/A549_rep_%_5ptrimmed_rc.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a TCAATTATATT --trimmed-only --quality-base=64 - | gzip -9 > $@"

$(PELCHAT_DIR)/NP_rep_%.fastq.gz: $(PELCHAT_DIR)/A549_rep_%_5ptrimmed_rc.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a GTAGATAATCA --trimmed-only --quality-base=64 - | gzip -9 > $@"

$(PELCHAT_DIR)/PA_rep_%.fastq.gz: $(PELCHAT_DIR)/A549_rep_%_5ptrimmed_rc.fastq.gz
	$(BSUB) "zcat $< | cutadapt -a TACTGATTCGA --trimmed-only --quality-base=64 - | gzip -9 > $@"

################
define pelchat-vfiltered-rule
$(INT_DIR)/PELCHAT_PA_rep_$(1)_vfiltered.fastq.gz: $(PELCHAT_DIR)/PA_rep_$(1).fastq.gz $(CONVERT_FORMATS)
	$(BSUB) "zcat $$< | cutadapt --trimmed-only -f fastq -a GCAAAAGCAGG --minimum-length=9 --quality-base=64 - | cutadapt -g GGGGGGGGGGGG -e 0 -O 1 --quality-base=64 - | $(MY_PYTHON) $(CONVERT_FORMATS) |  gzip -9 > $$@"

$(INT_DIR)/PELCHAT_HA_rep_$(1)_vfiltered.fastq.gz: $(PELCHAT_DIR)/HA_rep_$(1).fastq.gz
	$(BSUB) "zcat $$< | cutadapt --trimmed-only -f fastq -a GCAAAAGCAGG --minimum-length=9 --quality-base=64 - | $(MY_PYTHON) $(CONVERT_FORMATS) | gzip -9 > $$@"

$(INT_DIR)/PELCHAT_NA_rep_$(1)_vfiltered.fastq.gz: $(PELCHAT_DIR)/NA_rep_$(1).fastq.gz
	$(BSUB) "zcat $$< | cutadapt --trimmed-only -f fastq -a GCAAAAGCAGG --minimum-length=9 --quality-base=64 - | $(MY_PYTHON) $(CONVERT_FORMATS) | gzip -9 > $$@"

$(INT_DIR)/PELCHAT_NS1_rep_$(1)_vfiltered.fastq.gz: $(PELCHAT_DIR)/NS1_rep_$(1).fastq.gz
	$(BSUB) "zcat $$< | cutadapt --trimmed-only -f fastq -a GCAAAAGCAGG --minimum-length=9 --quality-base=64 - | $(MY_PYTHON) $(CONVERT_FORMATS) | gzip -9 > $$@"

$(INT_DIR)/PELCHAT_MP_rep_$(1)_vfiltered.fastq.gz: $(PELCHAT_DIR)/MP_rep_$(1).fastq.gz
	$(BSUB) "zcat $$< | cutadapt --trimmed-only -f fastq -a GCAAAAGCAGG --minimum-length=9 --quality-base=64 - | $(MY_PYTHON) $(CONVERT_FORMATS) | gzip -9 > $$@"

$(INT_DIR)/PELCHAT_NP_rep_$(1)_vfiltered.fastq.gz: $(PELCHAT_DIR)/NP_rep_$(1).fastq.gz
	$(BSUB) "zcat $$< | cutadapt --trimmed-only -f fastq -a GCAAAAGCAGG --minimum-length=9 --quality-base=64 - | $(MY_PYTHON) $(CONVERT_FORMATS) | gzip -9 > $$@"

$(INT_DIR)/PELCHAT_PB1_rep_$(1)_vfiltered.fastq.gz: $(PELCHAT_DIR)/PB1_rep_$(1).fastq.gz
	$(BSUB) "zcat $$< | cutadapt --trimmed-only -f fastq -a GCGAAAGCAGG --minimum-length=9 --quality-base=64 - | $(MY_PYTHON) $(CONVERT_FORMATS) | gzip -9 > $$@"

$(INT_DIR)/PELCHAT_PB2_rep_$(1)_vfiltered.fastq.gz: $(PELCHAT_DIR)/PB2_rep_$(1).fastq.gz
	$(BSUB) "zcat $$< | cutadapt --trimmed-only -f fastq -a GCGAAAGCAGG --minimum-length=9 --quality-base=64 - | $(MY_PYTHON) $(CONVERT_FORMATS) | gzip -9 > $$@"
endef

$(foreach rep,1 2,$(eval $(call pelchat-vfiltered-rule,$(rep))))
