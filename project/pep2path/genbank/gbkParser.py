from Bio import SeqIO
import os

from .CDSPrediction import CDSPrediction
from .GenbankFile import GenbankFile

'''Parses a single Genbank file in antiSMASH format into a GenbankFile object.'''
def parse_genbank(path, file):
	
	features = SeqIO.read(os.path.join(path, file), "genbank").features

	#find all the features with a non-empty antismash prediction
	CDS = [feature for feature in features if "aSProdPred" in feature.qualifiers.keys()]
	CDS = [feature for feature in CDS if feature.qualifiers["aSProdPred"][0] != ""]
	
	#find all the antismash domains with individual predictions
	aSDomains = [feature for feature in features if "specificity" in feature.qualifiers.keys()]
	
	#Gives a list of lists, where internal lists contain asDomain predictions, grouped by the CDS they are a subsection of
	grouped_aSDomains = []
	for feature in CDS:
		aSDomainGroup = []
		for aSDomain in aSDomains:
			if(aSDomain.location.start in feature.location and aSDomain.location.end in feature.location):
				aSDomainGroup.append(aSDomain)
		grouped_aSDomains.append(aSDomainGroup)
		
	gbk_file = []
	for zipped_CDS, zipped_aSDs in zip(CDS, grouped_aSDomains):
		overall_prediction = zipped_CDS.qualifiers["aSProdPred"]

		zipped_aSDs = [aSD.qualifiers["specificity"] for aSD in zipped_aSDs]
		predictions = [{(" ".join(prediction.split(' ')[:-1])).strip(':') : prediction.split(' ')[-1] for prediction in asD} for asD in zipped_aSDs]
		gbk_file.append(CDSPrediction(overall_prediction[0], predictions))
	
	return GenbankFile(gbk_file)