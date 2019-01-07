import itertools
import os

from spectra.mgfParser import main as mgfparser
from spectra.MassSpectrum import MassSpectrum

from genbank.gbkParser import parse_genbank
from genbank.aanames import AA_names

import comparisons as compare

msagg = mgfparser()

#get all unique components in genbank files
gbk_path = os.path.join(os.getcwd(), "genbank")
gbk_dataset = open("dataset.out", 'r')
genbank_names = [gbk_name.strip() for gbk_name in gbk_dataset]
path = os.path.join(os.path.join(gbk_path, "justin-20181022"))
gbk_files = [parse_genbank(path, gbk_name) for gbk_name in genbank_names]
gbk_dataset.close()

mass_tolerance_modes = [MassSpectrum.STATIC_MASS_TOLERANCE, 
						MassSpectrum.STATIC_MASS_TOLERANCE, 
						MassSpectrum.MAX_PPM_MASS_TOLERANCE] * 2
						
mass_thresholds = [0.01, 0.001, 0.00001] *2
intensity_thresholds = [0.05, 0.005] *3

outs = ["matching_0.01Mass_0.05Inten", "matching_0.001Mass_0.005Inten", "matching_0.00001ppmMass_0.05Inten", 
											"matching_0.01Mass_0.005Inten", "matching_0.001Mass_0.05Inten"]
zipped = zip(mass_tolerance_modes, mass_thresholds, intensity_thresholds, outs)
	
path = os.path.join("spectra", "spectraData")
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
	counts.append(compare.compare_unique_components(spectra_names, spectra_tags, genbank_names, gbk_components))
	headers = ["Original:"]
		
	#a flattened version of gbk_files that flattens out objects into lists
	shuffled_gbks = [[cds.unique_components() for cds in gbk.predictions] for gbk in gbk_files] 
	gbk_components = [[component for component in gbk if (component[0].upper() + component[1:].lower()) in AA_names] for gbk in gbk_components]
	for i in range(5):
		#shuffle gbks components
		for j in range(10000):
			compare.random_swap(shuffled_gbks) #shuffle
		#append counts for shuffled version
		shuffled_components = [list(set(itertools.chain.from_iterable(gbk))) for gbk in shuffled_gbks]
		counts.append(compare.compare_unique_components(spectra_names, spectra_tags, [str(component) for component in shuffled_components], shuffled_components))
		headers.append("Shuffled %d:" % i)
		
	for i in range(5):
		randomised_gbks = compare.randomise_components(gbk_components)
		counts.append(compare.compare_unique_components(spectra_names, spectra_tags, [str(gbk) for gbk in randomised_gbks], randomised_gbks))
		headers.append("Randomised %d:" % i)
			
	tables.append(counts)
		
	for counts in tables:
		print("---%s---\n" % out)
		for header, (total_correct, (max_spectrum, max_gbk, max_correct), comparisons, spectra_total, total_gbk_comps, total_gbk_files) in zip(headers, counts):
			print("%s Total Correct: %d\n" % (header, total_correct) +
					"Max Spectrum: %s Max Gbk: %s Max Correct: %d\n" % (max_spectrum, max_gbk, max_correct) +
					"Comparisons: %d Spectra_Components: %d Gbk_Components: %d Gbk_Files: %d\n" % (comparisons, spectra_total, total_gbk_comps, total_gbk_files))