STAR_GENOME_DIR = $(DATA_DIR)/STAR_hg19

$(STAR_GENOME_DIR)/Genome: $(HG19_WHOLE_GENOME)
	$(BSUB) mkdir -p $(STAR_GENOME_DIR)
	$(BSUB) "STAR --runMode genomeGenerate --genomeDir $(STAR_GENOME_DIR) --genomeFastaFiles $< --runThreadN `nproc`"

/pathToStarDir/STAR --runMode genomeGenerate --genomeDir
/path/to/GenomeDir --genomeFastaFiles /path/to/genome/fasta1
/path/to/genome/fasta2 --runThreadN <n> ...