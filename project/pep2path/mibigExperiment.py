import itertools
import os

from spectra.MassSpectrum import MassSpectrum
from spectra.spectraMain import find_longest_tag

from genbank.mibigParser import parse_genbank

import comparisons as compare

file = open("mibig_gnps_links_q3.csv", 'r')
for line in file: 
	if(len(line.split(',')) < 5): print(line)
	connections = [line.split(',') for line in file] #read all genbank to spectra connections
file.close()

#filter to genbank and spectra filenames, and filter to correct compound class
connections = [(connection[0], connection[1]) for connection in connections if ("NRP" in connection[4] or "RiPP" in connection[4])]
	
#get all unique components in genbank files
gbk_path = os.path.join(os.getcwd(), "genbank")
path = os.path.join(gbk_path, "mibigData")
gbk_files = [parse_genbank(path, gbk + ".gbk") for gbk, spectrum in connections]
	
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
	counts.append(compare.compare_unique_components(spectra_tags, gbk_components))
	headers = ["Original:"]
		
	#a flattened version of gbk_files that flattens out objects into lists
	shuffled_gbks = [[cds.unique_components() for cds in gbk.predictions] for gbk in gbk_files] 
	for i in range(5):
		#shuffle gbks components
		for j in range(10000):
			random_swap(shuffled_gbks) #shuffle
		#append counts for shuffled version
		counts.append(compare.compare_unique_components(spectra_tags, [list(set(itertools.chain.from_iterable(gbk))) for gbk in shuffled_gbks]))
		headers.append("Random %d:" % i)
			
	tables.append(counts)
		
for counts in tables:
	print("---%s---" % out)
	for header, (correct, spectra_total, gbk_total) in zip(headers, counts):
		print("%s Correct: %d Spectra_Com: %d Gbk_Com: %d" % (header, correct, spectra_total, gbk_total))