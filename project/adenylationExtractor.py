from Bio import SeqIO
from collections import defaultdict

from CDSPrediction import CDSPrediction
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

#Gives a list of lists, where internal lists contain asDomain predictions, grouped by the CDS they are a subsection of
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
	overall_prediction = zipped_CDS.qualifiers["aSProdPred"]

	zipped_aSDs = [aSD.qualifiers["specificity"] for aSD in zipped_aSDs]
	predictions = [{(" ".join(prediction.split(' ')[:-1])).strip(':') : prediction.split(' ')[-1] for prediction in asD} for asD in zipped_aSDs]
	gbk_file.append(CDSPrediction(overall_prediction, predictions))
gbk_file = GenbankFile(gbk_file)

print(gbk_file)
	