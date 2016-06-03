import numpy as np;
import shuffle_utils;
from scipy import stats;
import SeqGenUtils;
import pyngram;
import generateGraphs;


def ComputeIndividualDistributions(posDistList):

	aa_values = [x['AA'] for x in posDistList]
	ac_values = [x['AC'] for x in posDistList]
	gt_values = [x['GT'] for x in posDistList]
	ag_values = [x['AG'] for x in posDistList]
	cc_values = [x['CC'] for x in posDistList]
	ca_values = [x['CA'] for x in posDistList]
	cg_values = [x['CG'] for x in posDistList]
	tt_values = [x['TT'] for x in posDistList]
	gg_values = [x['GG'] for x in posDistList]
	gc_values = [x['GC'] for x in posDistList]
	at_values = [x['AT'] for x in posDistList]
	ga_values = [x['GA'] for x in posDistList]
	gt_values = [x['GT'] for x in posDistList]
	ta_values = [x['TA'] for x in posDistList]
	tc_values = [x['TC'] for x in posDistList]
	ct_values = [x['CT'] for x in posDistList]

	diNucDist = dict();

	diNucDist['AA'] = [round(np.mean(aa_values), 4), round(stats.sem(aa_values), 4)]
	diNucDist['AC'] = [round(np.mean(ac_values), 4),  round(stats.sem(ac_values), 4)]
	diNucDist['GT'] = [round(np.mean(gt_values), 4),  round(stats.sem(gt_values), 4)]
	diNucDist['AG'] = [round(np.mean(ag_values), 4), round(stats.sem(ag_values), 4)]
	diNucDist['CC'] = [round(np.mean(cc_values), 4),  round(stats.sem(cc_values), 4)]
	diNucDist['CA'] = [round(np.mean(ca_values), 4),  round(stats.sem(ca_values), 4)]
	diNucDist['CG'] = [round(np.mean(cg_values), 4),  round(stats.sem(cg_values), 4)]
	diNucDist['TT'] = [round(np.mean(tt_values), 4),  round(stats.sem(tt_values), 4)]
	diNucDist['GG'] = [round(np.mean(gg_values), 4),  round(stats.sem(gg_values), 4)]
	diNucDist['GC'] = [round(np.mean(gc_values), 4),  round(stats.sem(gc_values), 4)]
	diNucDist['AT'] = [round(np.mean(at_values), 4),  round(stats.sem(at_values), 4)]
	diNucDist['GA'] = [round(np.mean(ga_values), 4),  round(stats.sem(ga_values), 4)]
	diNucDist['GT'] = [round(np.mean(gt_values), 4),  round(stats.sem(gt_values), 4)]
	diNucDist['TA'] = [round(np.mean(ta_values), 4),  round(stats.sem(ta_values), 4)]
	diNucDist['TC'] = [round(np.mean(tc_values), 4),  round(stats.sem(tc_values), 4)]
	diNucDist['CT'] = [round(np.mean(ct_values), 4),  round(stats.sem(ct_values), 4)]


	return diNucDist;

def GetDistribListForDirectory(directory, filesToFind):
	distribList = list();
	for seqFile in SeqGenUtils.findFiles(directory, filesToFind):
		seqs, gc_list, fg_lengths = shuffle_utils.get_seqs(seqFile)
		dinuc_distrib = shuffle_utils.compute_dinuc_distrib(seqs, True)			
		distribList.append(dinuc_distrib)		
	return distribList;

def Compute3UtrNgramDistibution(nLen=2):

	three_utr_file_dir = "/projects/bhandare/workspace/scripts/NegFileCreator/";
	three_utr_file = three_utr_file_dir + "/" + "3UTR_transcripts_Human.txt"

	threeUtrFreq, threeUtrProb = GetNgramDistributionForFile(three_utr_file, nLen);

	return threeUtrFreq, threeUtrProb;

def Compute3UtrDistibution():

	three_utr_file_dir = "/projects/bhandare/workspace/scripts/NegFileCreator/";
	three_utr_file = "3UTR_transcripts_Human.txt"

	threeUtrDist = GetDistributionForFiles(three_utr_file_dir, three_utr_file);
	threeUtrValues = [x[0] for x in threeUtrDist.values()]
	threeUtrErrors = [x[1] for x in threeUtrDist.values()]	

	return threeUtrValues, threeUtrErrors;

def GetNgramDistributionForFile(fastaFile, nLen):

	freq, prob = pyngram.ComputeNgramFrequencyAndProbability(fastaFile, nLen)
	return freq, prob;

def GetDistributionForFiles(resultDir, searchRe):

	posDistList = GetDistribListForDirectory(resultDir, searchRe)
	posDiNucDist = ComputeIndividualDistributions(posDistList)
	return posDiNucDist;

def GetDiNucleotideDistribution(resultDir):

	posDiNucDist = GetDistributionForFiles(resultDir, "Signal*.fa")
	negDiNucDist = GetDistributionForFiles(resultDir, "NoSignal*.fa")

	return posDiNucDist, negDiNucDist

def GetNGramDistribution(fastaFile, N=2):

	threeUtrFreq, threeUtrProb = Compute3UtrNgramDistibution();
	print "Three UTR Freq: ", threeUtrFreq;

	faFileFreq, faFileProb = GetNgramDistributionForFile(fastaFile, N);
	print "Fasta File Freq: ", faFileFreq;

	print "Three UTR FReq (2): ", threeUtrFreq;

	generateGraphs.PlotNGrams(faFileFreq, threeUtrFreq, "N-Gram Frequency", "ngram.png")


