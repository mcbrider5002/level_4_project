import itertools
import os
import random
import numpy as np
from statistics import mean

from .spectra.mgfParser import load_files_from_dir as mgfparser
from .spectra.MassSpectrum import MassSpectrum
from .spectra.Tag import Tag
from .spectra.SpectrumTags import SpectrumTags

from .genbank.gbkParser import parse_genbank
from .genbank.CDSPrediction import CDSPrediction
from .genbank.GenbankFile import GenbankFile

from .masstables import AA_names

from . import comparisons as compare

'''Given a list of strings filters out any that are not present in the AA table (case-insensitive).'''
def only_AAs(ls):
	return [component for component in ls if (component[0].upper() + component[1:].lower()) in AA_names]

'''Read according to antiSMASH format all genbank files contained within the specified directory matching a name in filenames_path.'''
def get_gbks(path=os.path.join(os.path.join(os.path.join(os.path.dirname(__file__), "genbank"), "justin-20181022")),
				filenames_path=os.path.join(os.path.dirname(__file__), "dataset.out")):
	gbk_dataset = open(filenames_path, 'r')
	gbk_names = [gbk_name.strip() for gbk_name in gbk_dataset]
	gbk_files = [parse_genbank(path, gbk_name) for gbk_name in gbk_names]
	gbk_dataset.close()
	return gbk_files, gbk_names	
	
'''Shuffles the contents about between Genbank files and CDSes then compares the unique components of these random Genbank files to spectra as a baseline
	for how well the actual comparison works. (We expect to lose a few components during this process for component comparisons because of filtering out non-unique components.)'''
def match_shuffled(spectra_names, spectra_tags, gbk_components, headers, counts, iterations, comparator):
	for i in range(iterations):
		shuffled_components = compare.shuffle_components(gbk_components)
		counts.append(comparator(spectra_names, spectra_tags, ["" for gbk in shuffled_components], shuffled_components))
		headers.append("Shuffled %d:" % i)

'''Randomises the Genbank files then compares the unique components of these random Genbank files to spectra as a baseline
	for how well the actual comparison works. (We expect to lose a few components during this process for component comparisons because of filtering out non-unique components.)'''
def match_randomised(spectra_names, spectra_tags, gbk_components, headers, counts, iterations, comparator):
	for i in range(iterations):
		randomised_gbks = compare.randomise_components(gbk_components)
		counts.append(comparator(spectra_names, spectra_tags, ["" for gbk in randomised_gbks], randomised_gbks))
		headers.append("Randomised %d:" % i)
			
'''Given a list of tables in the format [[(best_spectrum, best_gbk, best_score), total_correct, comparisons, spectra_total, total_gbk_comps, total_gbk_files)]]
	prints this information out.'''
def print_simple(headers, counts, out):	
	print("---%s---\n" % out)
	for header, ((best_spectrum, best_gbk, best_score), total_correct, comparisons, spectra_total, total_gbk_comps, total_gbk_files) in zip(headers, counts):
		print("%s\n Total Correct: %d\n" % (header, total_correct) +
				"Best Spectrum: %s Best Gbk: %s Best Correct: %d\n" % (best_spectrum, best_gbk, best_score) +
				"Comparisons: %d Spectra_Components: %d Gbk_Components: %d Gbk_Files: %d\n" % (comparisons, spectra_total, total_gbk_comps, total_gbk_files))

'''Given a list of tables in the format [[(best_spectrum, best_gbk, best_score)]], prints this information out.'''				
def print_alignment(headers, counts, out):
	print("---%s---\n" % out)
	for header, (best_spectrum, best_gbk, best_score) in zip(headers, counts):
		print("%s\n" % (header) + "Best Spectrum: %s Best Gbk: %s Best Correct: %d\n" % (best_spectrum, best_gbk, best_score))
		
def print_table(title, group_headers, column_headers, split_comparisons, label_width=6, column_width=14):

	format_length = lambda l: "%-" + str(l) + "s"
	column_no = len(column_headers)
	group_no = len(group_headers)
	most_samples = max([len(s_comp) for s_comp in split_comparisons])
	
	row_labels = format_length(label_width)
	column_format = "|" + format_length(column_width)
	group_format = "|" + format_length((1 + column_width) * column_no - 1)
	total_length = label_width + ((1 + column_width) * column_no) * group_no

	print("---%s---\n" % title)
	print(row_labels % "" + "".join([group_format % string for string in group_headers]))
	print(row_labels % "" + ("".join([(column_format % string) for string in column_headers]) * group_no))
	print("-".join(["" for i in range(total_length)]))

	stringFormat = lambda x: ["".join([column_format % column for column in row]) for row in x] \
									+ ["".join([column_format % "" for j in range(column_no)]) for i in range(most_samples - len(x))] #pad blank rows
									
	for row in zip([row_labels % str(i) for i in range(1, most_samples + 1)], *[stringFormat(s_comp) for s_comp in split_comparisons]):
		print("".join(row))
		
	print("-".join(["" for i in range(total_length)]))
	
	avgFormat = lambda x: "".join([column_format % mean([y[i] for y in x]) for i in range(column_no)])
	print((row_labels % "Avg.") + "".join(["%s" % avgFormat(group) for group in split_comparisons]))
	
	print()
		
def print_intersection_tables(shuffled_iters, random_iters, table_data, label_width=6, column_width=14):

	for headerList, countList, out in table_data:
	
		countList = [(best_score, total_correct, comparisons)
						for ((best_spectrum, best_gbk, best_score), total_correct, comparisons, spectra_total, total_gbk_comps, total_gbk_files) in countList]
		
		actualComps = [countList[0]]
		shuffledComps = countList[1:1+shuffled_iters]
		randomComps = countList[1+shuffled_iters:]
		
		print_table(out, ("Normal", "Shuffled", "Random"), 
						 ("Best Match", "Total Correct", "Comparisons"),
						 [actualComps, shuffledComps, randomComps], 
						 label_width=label_width, 
						 column_width=column_width)
		
def print_alignment_tables(shuffled_iters, random_iters, table_data, label_width=6, column_width=14):
	
	for headerList, countList, out in table_data:
	
		countList = [(best_score,) for (best_spectrum, best_gbk, best_score) in countList]
		
		actualComps = [countList[0]]
		shuffledComps = countList[1:1+shuffled_iters]
		randomComps = countList[1+shuffled_iters:]
		
		print_table(out, ("Normal", "Shuffled", "Random"), 
						 ("Best Score",),
						 [actualComps, shuffledComps, randomComps], 
						 label_width=label_width, 
						 column_width=column_width)

def experiment(shuffled_iters, random_iters, comparator, printer, table_printer):

	gbk_files, gbk_names = get_gbks()

	#gbks flattened out into twice-nested lists where inner lists represent files, in them are their CDS and in them their prediction's tag as a list of strings
	gbk_files = [gbk.decompose_tags() for gbk in gbk_files]
	gbk_files = [[only_AAs(cds) for cds in gbk if cds] for gbk in gbk_files if gbk]
	
	msagg = mgfparser()[0][1] #extract the (presumed) only mgf, ignore the filename

	mass_tolerance_modes = [MassSpectrum.STATIC_MASS_TOLERANCE, 
							MassSpectrum.STATIC_MASS_TOLERANCE, 
							MassSpectrum.MAX_PPM_MASS_TOLERANCE] * 2
							
	mass_thresholds = [0.01, 0.001, 0.00001] *2
	intensity_thresholds = [0.05, 0.005] *3

	outs = ["matching_0.01Mass_0.05Inten", "matching_0.001Mass_0.005Inten", "matching_0.00001ppmMass_0.05Inten", 
												"matching_0.01Mass_0.005Inten", "matching_0.001Mass_0.05Inten"]
	zipped = zip(mass_tolerance_modes, mass_thresholds, intensity_thresholds, outs)
		
	path = os.path.join("spectra", "spectraData")
	table_data = []
	for mass_tolerance_mode, mass_threshold, intensity_threshold, out in zipped:
		
		spectra_names = [str(spectrum.id) for spectrum in msagg.spectra]
		msagg.filter_intensity_local(intensity_thresholds=[intensity_threshold*mint for mint in msagg.max_intensity_local()])
		spectra = msagg.find_sequence_tags(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
		
		counts = []
		counts.append(comparator(spectra_names, spectra, gbk_names, gbk_files))
		headers = ["Original:"]
		
		match_shuffled(spectra_names, spectra, gbk_files, headers, counts, shuffled_iters, comparator)
		match_randomised(spectra_names, spectra, gbk_files, headers, counts, random_iters, comparator)
		
		printer(headers, counts, out)
		table_data.append((headers, counts, out))
	
	table_printer(shuffled_iters, random_iters, table_data)

def intersection_experiment(shuffled_iters=5, random_iters=5):
	experiment(shuffled_iters, random_iters, compare.compare_unique_components, print_simple, print_intersection_tables)

def simple_experiment(shuffled_iters=5, random_iters=5):
	def fixed_alignment(spectra_names, spectra, gbk_names, gbk_tags): 
		compare.compare_alignment(spectra_names, spectra, gbk_names, gbk_tags, scoring=compare.SIMPLE_SCORE)
	experiment(shuffled_iters, random_iters, compare.compare_alignment, print_alignment, print_alignment_tables)
	
def p2p_experiment(shuffled_iters=5, random_iters=5):
	def fixed_alignment(spectra_names, spectra, gbk_names, gbk_tags): 
		compare.compare_alignment(spectra_names, spectra, gbk_names, gbk_tags, scoring=compare.P2P_SCORE)
	experiment(shuffled_iters, random_iters, compare.compare_alignment, print_alignment, print_alignment_tables)
	
def main():
	s =	[SpectrumTags(0, 3, {3: [Tag("Ser-Val-Gly", [], [])]})]
	g = [[["Ser", "Thr"], ["Val"], ["Ile", "Ile", "Ile"], ["Thr", "Gly"]]]
	spectrum_name, genbank_name, score = compare.compare_alignment([str(s)], s, [str(g)], g)
	print("Spectrum: " + spectrum_name)
	print("Genbank: " + genbank_name)
	print("Score: " + str(score))

	#alignment_experiment()
	
if __name__ == "__main__":

    main()