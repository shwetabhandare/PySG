from __future__ import division
import sys
import re
import parseKspectrum
import SeqGenUtils

def createSharedAndSpecificKmers(combinedKmers):
	hurKmers = set();
	ttpKmers = set();
	sharedKmers = set();

	for kmer in combinedKmers:
		if "HuR_" in kmer:
			hurKmers.add(kmer[4:]);
		elif "TTP_"	 in kmer:
			ttpKmers.add(kmer[4:]);
		else:
			sharedKmers.add(kmer);

	return hurKmers, ttpKmers, sharedKmers;


def validateKmers(RBPString, specificKmers, hurInDomainKmers, ttpInDomainKmers, sharedKmers):
	specificInHuR = hurInDomainKmers.intersection(specificKmers)
	specificInTTP = ttpInDomainKmers.intersection(specificKmers)
	specificInShared = sharedKmers.intersection(specificKmers)

	print RBPString + "-specific in HuR In Domain: ", len(specificInHuR), specificInHuR
	print RBPString + "-specific in TTP In Domain: ", len(specificInTTP), specificInTTP
	print RBPString + "-specific in Shared Kmers: ", len(specificInShared), specificInShared

	return specificInHuR, specificInTTP, specificInShared;

def getTopXKmers(hurInDomainFile, ttpInDomainFile, combinedInDomainFile, numKmers):
	hurInDomainKmers = set(parseKspectrum.FindKspectrumKmers(hurInDomainFile, numKmers));
	ttpInDomainKmers = set(parseKspectrum.FindKspectrumKmers(ttpInDomainFile, numKmers));
	combinedKmers = set(parseKspectrum.FindKspectrumKmers(combinedInDomainFile, numKmers));

	return hurInDomainKmers, ttpInDomainKmers, combinedKmers;


def writeFastaFileFromKmers(kmerList, prefix, outFileName):
	fastaFileLines = SeqGenUtils.createFastaFileFromKmers(kmerList, prefix);
	SeqGenUtils.writeFastaLinesToFile(fastaFileLines, outFileName);

if __name__ == "__main__":
	import sys

	hurInDomainKmers, ttpInDomainKmers, combinedKmers = getTopXKmers(sys.argv[1], sys.argv[2], sys.argv[3], 100);
	hurKmers, ttpKmers, sharedKmers = createSharedAndSpecificKmers(combinedKmers)

	print "=== Number of HUR Specific KMERS == ", len(hurKmers)
	print "=== Number of TTP Specific KMERS == ", len(ttpKmers)
	print "=== Number of Shared KMERS == ", len(sharedKmers)

	hurSpecific, hurInTTP, hurInShared = validateKmers("HuR", hurKmers, hurInDomainKmers, ttpInDomainKmers, sharedKmers)
	ttpInHuR, ttpSpecific, ttpInShared = validateKmers("TTP", ttpKmers, hurInDomainKmers, ttpInDomainKmers, sharedKmers)
	sharedInHuR, sharedInTTP, sharedInShared = validateKmers("Shared", sharedKmers, hurInDomainKmers, ttpInDomainKmers, sharedKmers);

	writeFastaFileFromKmers(hurSpecific, "HuR_In_HuRInDomain_", "HuR_In_HuRInDomain1.fa");
	writeFastaFileFromKmers(hurInTTP, "HuR_In_TTPInDomain_", "HuR_In_TTPInDomain1.fa");
	writeFastaFileFromKmers(hurInShared, "HuR_In_Shared_", "HuR_In_Shared1.fa");

	writeFastaFileFromKmers(ttpSpecific, "TTP_In_TTPInDomain_", "TTP_In_TTPInDomain1.fa");
	writeFastaFileFromKmers(ttpInHuR, "TTP_In_HuRInDomain_", "TTP_In_HuRInDomain1.fa");
	writeFastaFileFromKmers(ttpInShared, "TTP_In_Shared_", "TTP_In_Shared1.fa");

	writeFastaFileFromKmers(sharedInHuR, "Shared_In_HuRInDomain_", "Shared_In_HuRInDomain.fa");
	writeFastaFileFromKmers(sharedInTTP, "Shared_In_TTPInDomain_", "Shared_In_TTPInDomain.fa");
