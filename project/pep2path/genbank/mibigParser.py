from Bio import SeqIO
import os

from .CDSPrediction import CDSPrediction
from .GenbankFile import GenbankFile

from .aanames import AA_names

def parse_genbank(path, filename):
				
	#collect all CDS predictions
	CDS = [feature.qualifiers["product"] for feature in SeqIO.read(os.path.join(path, filename), "genbank").features if "product" in feature.qualifiers.keys()]
	#filter the ones with amino acids
	CDS = [[AA_names[AA] for AA in AA_names if (AA in product)] for product in CDS]
	
	#filter to unique
	newCDS = []
	for productList in CDS:
		newProductList = []
		for product in productList:
			if(not product in newProductList):
				newProductList.append(product)
		newCDS.append(newProductList)
				
	#construct CDSPredictions
	newCDS = [CDSPrediction("-" + "|".join(productList) + "-", {}) for productList in CDS]
	return GenbankFile(newCDS)