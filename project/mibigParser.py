from Bio import SeqIO
from collections import defaultdict
import itertools
import copy
import random
import os

from CDSPrediction import CDSPrediction
from GenbankFile import GenbankFile

from spectra.MassSpectrum import MassSpectrum
from spectra.MassSpectraAggregate import MassSpectraAggregate
from spectra.spectraMain import find_longest_tag

def parse_genbank(path, filename):
	AA_names = {	
					"Ala" : "Ala",
					"Alanine" : "Ala",
					"Arg" : "Arg",
					"Arginine" : "Arg",
					"Asn" : "Asn",
					"Asparagine" : "Asn",
					"Asp" : "Asp",
					"Aspartic Acid" : "Asp",
					"Cys" : "Cys",
					"Cysteine" : "Cys",
					"Gln" : "Gln",
					"Glutamine" : "Gln",
					"Gly" : "Gly",
					"Glycine" : "Gly",
					"His" : "His",
					"Histidine" : "His",
					"Ile" : "Ile/Leu",
					"Isoleucine" : "Ile/Leu",
					"Leu" : "Ile/Leu",
					"Leucine" : "Leu",
					"Lys" : "Lys",
					"Lysine" : "Lys",
					"Met" : "Met",
					"Methionine" : "Met",
					"Phe" : "Phe",
					"Phenylalanine" : "Phe",
					"Pro" : "Pro",
					"Proline" : "Pro",
					"Ser" : "Ser",
					"Serine": "Ser",
					"Thr" : "Thr",
					"Threonine" : "Thr",
					"Sec" : "Sec",
					"Selenocysteine" : "Sec",
					"Trp" : "Trp",
					"Tryptophan" : "Trp",
					"Tyr" : "Tyr",
					"Tyrosine" : "Tyr",
					"Val" : "Val",
					"Valine" : "Val"
				}
				
	#collect all CDS predictions
	CDS = [feature.qualifiers["product"] for feature in SeqIO.read(os.path.join(path, filename), "genbank").features if "product" in feature.qualifiers.keys()]
	#filter the ones with amino acids
	print(CDS)
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
	
def random_component(gbk_files):
	if(len(gbk_files) - 1 <= 0): return -1, -1, -1
	file = random.randint(0, len(gbk_files) - 1) #generate an index for random file
	if(gbk_files[file] == []): return file, -1, -1
	cds = random.randint(0, len(gbk_files[file]) - 1) #generate an index for random cds in that file
	if(gbk_files[file][cds] == []): return file, cds, -1
	component = random.randint(0, len(gbk_files[file][cds]) - 1)
	return file, cds, component
	
def random_swap(gbk_files):
	file1, cds1, component1 = random_component(gbk_files)
	file2, cds2, component2 = random_component(gbk_files)
	if(component1 == -1 or component2 == -1): return gbk_files
	temp = gbk_files[file1][cds1][component1]
	gbk_files[file1][cds1][component1] = gbk_files[file2][cds2][component2]
	gbk_files[file2][cds2][component2] = temp
	return gbk_files
	
def compare_unique_components(spectra_components, gbk_components):
	correct = 0
	total_spectra = 0
	total_gbk = 0
	for gbk, spectrum in zip(gbk_components, spectra_components):
		total_spectra += len(spectrum)
		total_gbk += len(gbk)
		for component in spectrum:
			if(component in gbk):
				correct += 1
	return correct, total_spectra, total_gbk

def main():

	file = open("mibig_gnps_links_q3.csv", 'r')
	for line in file: 
		if(len(line.split(',')) < 5): print(line)
		connections = [line.split(',') for line in file] #read all genbank to spectra connections
	file.close()

	#filter to genbank and spectra filenames, and filter to correct compound class
	connections = [(connection[0], connection[1]) for connection in connections if ("NRP" in connection[4] or "RiPP" in connection[4])]
	
	#get all unique components in genbank files
	path = os.path.join(os.getcwd(), "mibigData")
	gbk_files = [parse_genbank(path, gbk + ".gbk") for gbk, spectrum in connections]
	print(gbk_files)
	
	mass_tolerance_modes = [MassSpectrum.STATIC_MASS_TOLERANCE, 
							MassSpectrum.STATIC_MASS_TOLERANCE, 
							MassSpectrum.MAX_PPM_MASS_TOLERANCE] * 2
							
	mass_thresholds = [0.01, 0.001, 0.00001] *2
	intensity_thresholds = [0.05, 0.005] *3
	outs = ["matching_0.01Mass_0.05Inten", "matching_0.001Mass_0.005Inten", "matching_0.00001ppmMass_0.05Inten", 
											"matching_0.01Mass_0.005Inten", "mibig_0.001Mass_0.05Inten"]
	zipped = zip(mass_tolerance_modes, mass_thresholds, intensity_thresholds, outs)
	
	path = os.path.join(os.path.join(os.getcwd(), "spectra"), "spectraData")
	tables = []
	for mass_tolerance_mode, mass_threshold, intensity_threshold, out in zipped:
		
		spectra_tags = [find_longest_tag(path=path, pattern=filename, intensity_thresholds=intensity_threshold, 
											mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)[1]
												for filename in [spectrum + ".ms" for gbk, spectrum in connections]]
												
		spectra_tags = itertools.chain.from_iterable(spectra_tags) #flatten list of spectra tags
		
		#convert spectra tags to their unique components
		spectra_tags = [tags.unique_components() for tags in spectra_tags]
		
		#compare unique components of spectra tags to unshuffled gbk
		gbk_components = [gbk.unique_components() for gbk in gbk_files]
		
		counts = []
		counts.append(compare_unique_components(spectra_tags, gbk_components))
		headers = ["Original:"]
		
		#a flattened version of gbk_files that flattens out objects into lists
		shuffled_gbks = [[cds.unique_components() for cds in gbk.predictions] for gbk in gbk_files] 
		for i in range(5):
			#shuffle gbks components
			for j in range(10000):
				random_swap(shuffled_gbks) #shuffle
			#append counts for shuffled version
			counts.append(compare_unique_components(spectra_tags, [list(set(itertools.chain.from_iterable(gbk))) for gbk in shuffled_gbks]))
			headers.append("Random %d:" % i)
			
		tables.append(counts)
		
	for counts in tables:
		print("---%s---" % out)
		for header, (correct, spectra_total, gbk_total) in zip(headers, counts):
			print("%s Correct: %d Spectra_Com: %d Gbk_Com: %d" % (header, correct, spectra_total, gbk_total))
	
if __name__ == "__main__":

    main()