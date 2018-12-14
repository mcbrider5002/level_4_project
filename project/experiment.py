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

from spectra.mgfParser import main as mgfparser

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

'''Parses a single Genbank file in antiSMASH format into a GenbankFile object.'''
def genbank_parser(path, file):
	
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
	
'''Given a twice-nested list, chooses a random element from any of the innermost lists (i.e. randomly selects an element with three indices required to access it).'''
def random_component(gbk_files):
	if(len(gbk_files) - 1 <= 0): return -1, -1, -1
	file = random.randint(0, len(gbk_files) - 1) #generate an index for random file
	if(gbk_files[file] == []): return file, -1, -1
	cds = random.randint(0, len(gbk_files[file]) - 1) #generate an index for random cds in that file
	if(gbk_files[file][cds] == []): return file, cds, -1
	component = random.randint(0, len(gbk_files[file][cds]) - 1)
	return file, cds, component
	
'''Given a twice nested list, swaps two of the innermost elements at random.'''
def random_swap(gbk_files):
	file1, cds1, component1 = random_component(gbk_files)
	file2, cds2, component2 = random_component(gbk_files)
	if(component1 == -1 or component2 == -1): return gbk_files
	temp = gbk_files[file1][cds1][component1]
	gbk_files[file1][cds1][component1] = gbk_files[file2][cds2][component2]
	gbk_files[file2][cds2][component2] = temp
	return gbk_files
	
'''Given a nested list where internal lists represent genbank files and their contents are their unique components, 
	rearranges all the inner lists (files) and then randomly selects an animo acid to replace each of those components
	before filtering them to be unique again.'''
def randomise_components(gbk_components):
	random.shuffle(gbk_components) #change order of files to be random
	return [list(set([AA_names[str(random.choice(list(AA_names.keys())))] for component in gbk])) for gbk in gbk_components] #randomise all components
	
'''Given a nested list where internal lists represent a mass spectrum and their contents are their unique components,
	and a nested list where internal lists represent genbank files and their contents are their unique components,
	and two lists of names where the names correspond to the files in each of these lists
	returns several statistics based on comparisons of the unique components of these nested lists.'''
def compare_unique_components(spectra_names, spectra_components, gbk_names, gbk_components):

	total_correct = 0
	max_correct = ("No spectrum!", "No gbk!", 0)
	total_comparisons = 0
	total_spectra = 0
	total_gbk_comp = 0
	total_gbk_files = len(gbk_components)
	
	for spectrum in spectra_components:
		total_spectra += len(spectrum)
	
	for gbk_name, gbk in zip(gbk_names, gbk_components):
	
		total_gbk_comp += len(gbk)
		if(len(gbk) == 0): continue
		
		for spectrum_name, spectrum in zip(spectra_names, spectra_components):
			correct = 0
			for component in spectrum:
				total_comparisons += 1 
				if(component.lower() in gbk or (component[0].upper() + component[:1].lower()) in gbk or component in gbk):
					correct += 1
			
			total_correct += correct
			max_correct = max(max_correct, (spectrum_name, gbk_name, correct), key=lambda t: t[2])		
					
	return total_correct, max_correct, total_comparisons, total_spectra, total_gbk_comp, total_gbk_files
	
'''Given a nested list of spectra then (ordered) components,
	and a twice-nested list of genbank files, then cds, then (ordered) components,
	attempts to score a match of each spectrum to each genbank file by comparing the components to each possible ordering of the cds,
	possibly with some cds reversed.'''
'''def compare_alignment(spectra_components, gbk_components):
	
	best_score = ("No spectrum!", "No gbk!", 0)
	for file in files:
		products = itertools.product([(cds, "-".join(cds.split("-").reverse())) for cds in file])
		for product in products:
			orderings = itertools.permutations(product, len(product))
			for ordering in orderings:
				for spectrum in spectra_components:
					new_score = spectrum.alignment_comparison(ordering)
					best_score = max(best_score, ("spectrum_name", "gbk_name", new_score))
			
	return best_score'''

msagg = mgfparser()

#get all unique components in genbank files
gbk_dataset = open("dataset.out", 'r')
genbank_names = [gbk_name.strip() for gbk_name in gbk_dataset]
path = os.path.join(os.getcwd(), "justin-20181022")
gbk_files = [genbank_parser(path, gbk_name) for gbk_name in genbank_names]
gbk_dataset.close()

mass_tolerance_modes = [MassSpectrum.STATIC_MASS_TOLERANCE, 
						MassSpectrum.STATIC_MASS_TOLERANCE, 
						MassSpectrum.MAX_PPM_MASS_TOLERANCE] * 2

							
mass_thresholds = [0.01, 0.001, 0.00001] *2
intensity_thresholds = [0.05, 0.005] *3
#intensity_thresholds = [0.05]
outs = ["matching_0.01Mass_0.05Inten", "matching_0.001Mass_0.005Inten", "matching_0.00001ppmMass_0.05Inten", 
											"matching_0.01Mass_0.005Inten", "matching_0.001Mass_0.05Inten"]
zipped = zip(mass_tolerance_modes, mass_thresholds, intensity_thresholds, outs)
	
path = os.path.join(os.path.join(os.getcwd(), "spectra"), "spectraData")
tables = []
for mass_tolerance_mode, mass_threshold, intensity_threshold, out in zipped:
	
	tables = []
	spectra_names = [str(spectrum.id) for spectrum in msagg.spectra]
	msagg.filter_intensity_local(intensity_thresholds=[intensity_threshold*mint for mint in msagg.max_intensity_local()])
	spectra_tags = msagg.find_sequence_tags(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
		
	#convert spectra tags to their unique components
	spectra_tags = [tags.unique_components() for tags in spectra_tags if tags.unique_components() != []]
		
	#compare unique components of spectra tags to unshuffled gbk
	gbk_components = [gbk.unique_components() for gbk in gbk_files]
	gbk_components = [[component for component in gbk if (component[0].upper() + component[1:].lower()) in AA_names] for gbk in gbk_components]
		
	counts = []
	counts.append(compare_unique_components(spectra_names, spectra_tags, genbank_names, gbk_components))
	headers = ["Original:"]
		
	#a flattened version of gbk_files that flattens out objects into lists
	shuffled_gbks = [[cds.unique_components() for cds in gbk.predictions] for gbk in gbk_files] 
	gbk_components = [[component for component in gbk if (component[0].upper() + component[1:].lower()) in AA_names] for gbk in gbk_components]
	for i in range(5):
		#shuffle gbks components
		for j in range(10000):
			random_swap(shuffled_gbks) #shuffle
		#append counts for shuffled version
		shuffled_components = [list(set(itertools.chain.from_iterable(gbk))) for gbk in shuffled_gbks]
		counts.append(compare_unique_components(spectra_names, spectra_tags, [str(component) for component in shuffled_components], shuffled_components))
		headers.append("Shuffled %d:" % i)
		
	for i in range(5):
		randomised_gbks = randomise_components(gbk_components)
		counts.append(compare_unique_components(spectra_names, spectra_tags, [str(gbk) for gbk in randomised_gbks], randomised_gbks))
		headers.append("Randomised %d:" % i)
			
	tables.append(counts)
		
	for counts in tables:
		print("---%s---\n" % out)
		for header, (total_correct, (max_spectrum, max_gbk, max_correct), comparisons, spectra_total, total_gbk_comps, total_gbk_files) in zip(headers, counts):
			print("%s Total Correct: %d\n" % (header, total_correct) +
					"Max Spectrum: %s Max Gbk: %s Max Correct: %d\n" % (max_spectrum, max_gbk, max_correct) +
					"Comparisons: %d Spectra_Components: %d Gbk_Components: %d Gbk_Files: %d\n" % (comparisons, spectra_total, total_gbk_comps, total_gbk_files))