import itertools
import os

from spectra.mgfParser import load_files_from_dir as mgfparser
from spectra.MassSpectrum import MassSpectrum

from genbank.gbkParser import parse_genbank
from genbank.aanames import AA_names

import comparisons as compare

'''Read according to antiSMASH format all genbank files contained within the specified directory matching a name in filenames_path.'''
def get_gbks(path=os.path.join(os.path.join(os.path.join(os.getcwd(), "genbank"), "justin-20181022")),
				filenames_path="dataset.out"):
	gbk_dataset = open(filenames_path, 'r')
	gbk_names = [gbk_name.strip() for gbk_name in gbk_dataset]
	gbk_files = [parse_genbank(path, gbk_name) for gbk_name in gbk_names]
	gbk_dataset.close()
	return gbk_files, gbk_names

'''Shuffles the contents about between Genbank files and CDSes then compares the unique components of these random Genbank files to spectra as a baseline
	for how well the actual comparison works. (We expect to lose a few components during this process because of filtering out non-unique components.)'''
def match_shuffled(spectra_names, spectra_tags, cmpts_with_cds, headers, counts, iterations):
	for i in range(iterations):
		#shuffle gbks components
		for j in range(10000):
			shuffled_gbks = compare.random_swap(cmpts_with_cds) #shuffle
		#append counts for shuffled version
		shuffled_components = [list(set(itertools.chain.from_iterable(gbk))) for gbk in shuffled_gbks]
		counts.append(compare.compare_unique_components(spectra_names, spectra_tags, [str(component) for component in shuffled_components], shuffled_components))
		headers.append("Shuffled %d:" % i)

'''Randomises the Genbank files then compares the unique components of these random Genbank files to spectra as a baseline
	for how well the actual comparison works. (We expect to lose a few components during this process because of filtering out non-unique components.)'''
def match_randomised(spectra_names, spectra_tags, gbk_components, headers, counts, iterations):
	for i in range(iterations):
		randomised_gbks = compare.randomise_components(gbk_components)
		counts.append(compare.compare_unique_components(spectra_names, spectra_tags, [str(gbk) for gbk in randomised_gbks], randomised_gbks))
		headers.append("Randomised %d:" % i)
			
'''Given a list of tables in the format [[(total_correct, (max_spectrum, max_gbk, max_correct), comparisons, spectra_total, total_gbk_comps, total_gbk_files)]]
	prints this information out.'''
def print_output(headers, counts, out):	
	print("---%s---\n" % out)
	for header, (total_correct, (max_spectrum, max_gbk, max_correct), comparisons, spectra_total, total_gbk_comps, total_gbk_files) in zip(headers, counts):
		print("%s Total Correct: %d\n" % (header, total_correct) +
				"Max Spectrum: %s Max Gbk: %s Max Correct: %d\n" % (max_spectrum, max_gbk, max_correct) +
				"Comparisons: %d Spectra_Components: %d Gbk_Components: %d Gbk_Files: %d\n" % (comparisons, spectra_total, total_gbk_comps, total_gbk_files))
	
def experiment():

	gbk_files, gbk_names = get_gbks()
	#gbks flattened out into lists of lists where inner lists represent files and in them is their unique components
	gbk_components = [gbk.unique_components() for gbk in gbk_files]
	gbk_components = [[component for component in gbk if (component[0].upper() + component[1:].lower()) in AA_names] for gbk in gbk_components]
	#gbks flattened out into lists of lists where inner lists represent files, in them and in them are their CDS and in them their unique components
	cmpts_with_cds = [[cds.unique_components() for cds in gbk.predictions] for gbk in gbk_files] 
	cmpts_with_cds = [[[component for component in cds if (component[0].upper() + component[1:].lower()) in AA_names]
																										for cds in gbk]
																										for gbk in cmpts_with_cds]
	
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
	for mass_tolerance_mode, mass_threshold, intensity_threshold, out in zipped:
		
		spectra_names = [str(spectrum.id) for spectrum in msagg.spectra]
		msagg.filter_intensity_local(intensity_thresholds=[intensity_threshold*mint for mint in msagg.max_intensity_local()])
		spectra_tags = msagg.find_sequence_tags(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
			
		#convert spectra tags to their unique components
		spectra_tags = [tags.unique_components() for tags in spectra_tags if tags.unique_components() != []]
			
		counts = []
		counts.append(compare.compare_unique_components(spectra_names, spectra_tags, gbk_names, gbk_components))
		headers = ["Original:"]
		
		match_shuffled(spectra_names, spectra_tags, cmpts_with_cds, headers, counts, 5)
		match_randomised(spectra_names, spectra_tags, gbk_components, headers, counts, 5)
		
		print_output(headers, counts, out)
		
experiment()