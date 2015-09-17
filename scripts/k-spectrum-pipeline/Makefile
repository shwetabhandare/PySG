DATA_DIR=../../data/
PROJECTS_DIR=/projects/bhandare
KSPECTRUM_PIPELINE_DIR=$(PROJECTS_DIR)/k-spectrum-pipeline
RESULTS_DIR=`mkdate`

POS_FILE=$(DATA_DIR)/HuR_TTP/Only_to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative.fa
NEG_FILE=$(DATA_DIR)/HuR_TTP/NONOnly_to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative.fa

HuR_NEGATIVE_SET=$(PROJECTS_DIR)/workspace/scripts/NegFileCreator/NonHuRTranscripts.txt
COMBINED_FILE=$(DATA_DIR)/HuR_TTP/Only_to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_NONOnly_to_ucsc.fp_0.05.filtered.hg19.elavl1_conservative_3206_3206_CompleteSet.txt

FEATURE_FILE=$(DATA_DIR)/HuR_TTP/OnlyHuRFeatures.dat

BASIC_KSPECTRUM_SCRIPT=$(PROJECTS_DIR)/workspace/kernels/k-spectrum/scripts/basic_kspectrum.py
TOP_KMERS_SCRIPT=$(PROJECTS_DIR)/workspace/scripts/extractTopFeatures.py

ONLY_HUR_YML=$(DATA_DIR)/HuR_TTP/OnlyHuR.yml

$(FEATURE_FILE): $(COMBINED_FILE) $(BASIC_KSPECTRUM_SCRIPT) $(ONLY_HUR_YML)
	python $(BASIC_KSPECTRUM_SCRIPT) $(ONLY_HUR_YML)
	mv OnlyHuRFeatures.dat $(DATA_DIR)/HuR_TTP/
	python $(TOP_KMERS_SCRIPT) $(FEATURE_FILE) 25 $(DATA_DIR)/HuR_TTP/$(FEATURE_FILE).topfeatures.txt
	git add $(FEATURE_FILE)
	git add $(DATA_DIR)/HuR_TTP/$(FEATURE_FILE).topfeatures.txt


$(COMBINED_FILE): $(POS_FILE) $(NEG_FILE) seq_counter.py 
	python seq_counter.py $(POS_FILE) $(NEG_FILE)
	mv $(COMBINED_FILE) $(DATA_DIR)/HuR_TTP/

$(NEG_FILE): $(POS_FILE) LenMatchNegatives.py fasta.py
	python LenMatchNegatives.py $(POS_FILE) $(HuR_NEGATIVE_SET)
	mv $(NEG_FILE) $(DATA_DIR)/HuR_TTP/
	git add $(DATA_DIR)/HuR_TTP/$(NEG_FILE)
	git commit -av
	git push origin master



