from Bio import SeqIO
from collections import defaultdict

from GenbankPredictions import GenbankPredictions
from GenbankFile import GenbankFile

#find all the features with a non-empty antismash prediction
CDS = [feature for feature in SeqIO.read("c00001_NODE_1_...cluster001.gbk", "genbank").features if "aSProdPred" in feature.qualifiers.keys()]
CDS = [feature for feature in CDS if feature.qualifiers["aSProdPred"][0] != ""]

#find all the antismash domains with individual predictions
aSDomains = [feature for feature in SeqIO.read("c00001_NODE_1_...cluster001.gbk", "genbank").features if "specificity" in feature.qualifiers.keys()]

#group aSDomains into a tuple of subfeatures for each CDS
#grouped_aSDomains = [list((aSDomain.qualifiers["specificity"] for aSDomain in aSDomains if (aSDomain.location.start in feature.location and aSDomain.location.end in feature.location))) for feature in CDS]
#print(list(grouped_aSDomains))

#print(aSDomains)

grouped_aSDomains = []
for feature in CDS:
	aSDomainGroup = []
	for aSDomain in aSDomains:
		if(aSDomain.location.start in feature.location and aSDomain.location.end in feature.location):
			aSDomainGroup.append(aSDomain)
	grouped_aSDomains.append(aSDomainGroup)
	
print()
print(CDS)
print()
print(aSDomains)
print()
print(grouped_aSDomains)

gbk_file = []
for zipped_CDS, zipped_aSDs in zip(CDS, grouped_aSDomains):
	zipped_aSDs = [aSD["specificity"] for aSD in zipped_aSDs]
	overall_prediction = zipped_CDS.qualifiers[aSProdPred]
	predictions = {}
	gbk_file.append(GenbankPredictions(overall_prediction, predictions))
GenbankFile([GenbankPredictions(zipped_CDS.qualifiers["aSProdPred"], zipped_aSD) ])

print(gbk_file)
	