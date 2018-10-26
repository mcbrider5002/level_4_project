import os

import massSpectraParser as parser
from MassSpectrum import MassSpectrum
from MassSpectraAggregate import MassSpectraAggregate

import numpy as np

import random
import string

import itertools
from collections import defaultdict

import csv

##########
###Data###
##########

#uses avg mass instead of monoisotopic
AA_mass_table = {	
					"Ala" : 71.0779,
					"Arg" : 156.1857,
					"Asn" : 114.1026,
					"Asp" : 115.0874,
					"Cys" : 103.1429,
					"Gln" : 128.1292,
					"Gly" : 57.0513,
					"His" : 137.1393,
					"Ile/Leu" : 113.1576,
					"Lys" : 128.1723,
					"Met" : 131.1961,
					"Phe" : 147.1739,
					"Pro" : 97.1152,
					"Ser" : 87.0773,
					"Thr" : 101.1039,
					"Sec" : 150.0379,
					"Trp" : 186.2099,
					"Tyr" : 163.1733,
					"Val" : 99.1311
				}
				
###########################
###Convenience Functions###
###########################

def setup_mass_spectra(file):
	
	filename, dict = file
	
	return MassSpectrum(id=filename,
						compound=dict["compound"], 
						formula=dict["formula"], 
						parent_mass=dict["parentmass"], 
						ionization=dict["ionization"], 
						inchi=dict["InChI"], 
						inchi_key=dict["InChIKey"], 
						smiles=dict["smiles"], 
						ms2peaks=dict["ms2peaks"], 
						mass_table=AA_mass_table)
						
'''Convenience function to perform the multiple steps necessary to find tags and their ids.'''
def find_longest_tag(path=os.path.join(os.getcwd(), "spectraData"), pattern="*.ms", 
						intensity_thresholds=None,
						mass_tolerance_mode=None, mass_threshold=0.00001):
	
	import time
	
	print("Starting to parse for pattern \"" + pattern + "\"...")
	start = time.clock()
	records = parser.load_files_from_dir(path=path, pattern=pattern)
	print("Time taken: " + str(time.clock() - start))
	
	print("Converting to objects...")
	start = time.clock()
	records = [setup_mass_spectra(record) for record in records]
	spectra = MassSpectraAggregate(records)
	print("Time taken: " + str(time.clock() - start))
	
	print("Filtering intensity...")
	start = time.clock()
	intensity_thresholds = [intensity*0.05 for intensity in spectra.max_intensity_local()] if intensity_thresholds is None else [intensity*intensity_thresholds for intensity in spectra.max_intensity_local()]
	intensity_thresholds = [intensity*intensity_thresholds for intensity in spectra.max_intensity_local()] if isinstance(intensity_thresholds, float) else intensity_thresholds
	spectra.filter_intensity_local(intensity_thresholds=intensity_thresholds)
	print("Time taken: " + str(time.clock() - start))
	
	print("Sorting by mass...")
	start = time.clock()
	spectra.sort_by_mass()
	print("Time taken: " + str(time.clock() - start))
	
	print("Finding longest tag...")
	start = time.clock()
	longest_tag, tags = spectra.find_longest_tag_global(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
	print("Time taken: " + str(time.clock() - start))
	
	print("Longest tag found: " + str(longest_tag) + "\n")
	
	return longest_tag, spectra.attach_spectra_ids(tags)
	
'''Prints out all the tied-longest tags and their spectrum id.'''
def print_longest_tags(path=os.path.join(os.getcwd(), "spectraData"), pattern="*.ms", 
						intensity_thresholds=None,
						mass_tolerance_mode=None, mass_threshold=0.00001):
						
	longest_tag, spectra_tags = find_longest_tag(path=path, pattern=pattern, intensity_thresholds=intensity_thresholds, 
											mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
											
	spectra_tags = filter(lambda spectrum_tags : spectrum_tags[1][0] == longest_tag, spectra_tags) #first get all spectra with the maximum tag lengths
	for spectrum in spectra_tags: #print the id of each spectrum, its longest tag length, and all tags with the longest tag length (which should be identical to local longest)
		print(spectrum[0], spectrum[1][0])
		for tag in spectrum[1][1][longest_tag]:	
			print(tag)
	
'''Writes to "tags.out" on the specifed path all tags with the specified length and their spectrum id.'''	
def write_tags(path=os.path.join(os.getcwd(), "spectraData"), pattern="*.ms", 
						intensity_thresholds=None,
						mass_tolerance_mode=None, mass_threshold=0.00001,
						lengths_to_print=None, output_path=os.path.join(os.getcwd(), "out"),
						filename="tags.out"):
						
	longest_tag, spectra_tags = find_longest_tag(path=path, pattern=pattern, intensity_thresholds=intensity_thresholds, 
											mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
								
	lengths_to_print = [longest_tag] if lengths_to_print is None else lengths_to_print #default value
	lengths_to_print = [lengths_to_print] if isinstance(lengths_to_print, str) else lengths_to_print #take string input
	lengths_to_print = list(filter(lambda length : length <= longest_tag and length > 0, lengths_to_print)) #remove impossible indices
	
	os.chdir(output_path)
	file = open(filename, 'w')
	for spectrum in spectra_tags: #print the id of each spectrum, the current tag length, and all tags with that length
	
		for length in lengths_to_print:
		
			if(length in spectrum[1][1].keys()):
			
				file.write(str("---" + spectrum[0]) + " " + str(length) + "---\n")
				
				for tag in spectrum[1][1][length]:	
					file.write(str(tag) + "\n")
		
	file.close()

'''Helper to take the mibig file and automatically write the contents of corresponding files to tags.out.'''	
def mibig_parser():

	file = open("mibig_gnps_links_q3.csv", 'r')
	filenames = [line.split(',')[1] + ".ms" for line in file]
	file.close()
	
	mass_tolerance_modes = [MassSpectrum.STATIC_MASS_TOLERANCE, MassSpectrum.STATIC_MASS_TOLERANCE, MassSpectrum.PPM_MASS_TOLERANCE] * 2
	mass_thresholds = [0.01, 0.001, 0.00001] *2
	intensity_thresholds = [0.05, 0.005] *3
	outs = ["mibig_hMass_hInten", "mibi_lMass_lInten", "mibig_ppmMass_hInten", "mibig_hMass_lInten", "mibig_lMass_hInten"]
	zipped = zip(mass_tolerance_modes, mass_thresholds, intensity_thresholds, outs)
	
	for mass_tolerance_mode, mass_threshold, intensity_threshold, out in zipped:
		
		#(longest_tag, (id, (local_longest, {length: [tag, peaks]})))
		dataset = [find_longest_tag(pattern=filename, intensity_thresholds=intensity_threshold, 
											mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
												for filename in filenames]
			
		dataset = [spectra_data for longest_tag, spectra_data in dataset if spectra_data != []]
		
		os.chdir(os.path.join(os.getcwd(), os.path.join("..", "out")))
		file = open(out + ".out", 'w')
		if(len(dataset) > 0):
			for pattern_results in dataset:
				for id, (local_longest, tags) in pattern_results: #print the id of each spectrum, the current tag length, and all tags with that length
					for length in range(local_longest):
						if(length in tags.keys()):
							file.write("---" + str(id) + " " + str(length) + "---\n")
							for tag in tags[length]:	
								file.write(str(tag) + "\n")
		
		file.close()
	
	
###########
###Tests###
###########

'''Tests correctness of MassSpectraAggregate's get_spectra_ids, attach_spectra_ids and get_spectra_by_id.'''
def test_sequence_by_id_helpers():

	#generate a random number of random ids of random length
	ids = ["".join([random.choice(string.ascii_letters) for length in range(random.randint(1, 6))]) for number in range(random.randint(1, 6))] 
	
	spectra = MassSpectraAggregate([MassSpectrum(id=id) for id in ids]) #convert to spectra
	
	if(not (ids == spectra.get_spectra_ids())):
		print("get_spectra_ids fails test!")
	
	combined = [(id, index) for id, index in zip(ids, range(len(ids)))]
	if(not (combined == spectra.attach_spectra_ids(range(len(ids))))):
		print("attach_spectra_ids fails test!")
		
	if(not (ids == [spectrum.id for spectrum in spectra.get_spectra_by_id(ids)])):
		print("get_spectra_by_id fails test!")
	
'''Helper to make generating test spectra easier.'''		
def generate_random_readings(columns):
	
	noReadings = random.randint(10, 30) #choose length of data
	ms2peaks = np.zeros((noReadings, 2))
		
	lower_bound = random.randint(0, 1000)
	interval = (lower_bound, lower_bound + random.randint(0, 2000)) #interval to generate random numbers in
	ms2peaks[:, columns] = (interval[1]-interval[0]) * np.random.random(size=(noReadings,len(columns))) + interval[0] #fill in mass readings with random numbers
		
	return MassSpectrum(ms2peaks=ms2peaks)
	
'''Helper function to make creating spectra a bit neater.'''
def generate_random_spectra(columns):
	spectra = [generate_random_readings(columns) for element in range(random.randint(2, 5))]
	spectra_aggregate = MassSpectraAggregate(spectra)
	return spectra, spectra_aggregate
	
'''Helper to check all elements in some spectra are less-than-or-equal to some list of corresponding elements.'''
def check_all_are_leq(spectra, returned_maxes, column):
	return [np.all((spectrum.ms2peaks[:, column] <= maxm)) for spectrum, maxm in zip(spectra, returned_maxes)]
	
'''Helper to check an element is present in corresponding spectra.'''
def check_presence(spectra, returned_maxes, column):
	return [maxm in spectrum.ms2peaks[:, column] for spectrum, maxm in zip(spectra, returned_maxes)]
	
'''Helper to check values actually are max values.'''
def verify_max_local(spectra, returned_maxes, column):
	smaller_than_max = check_all_are_leq(spectra, returned_maxes, column)
	present = check_presence(spectra, returned_maxes, column)
	return not(all(smaller_than_max) and all(present))

'''Helper to check value actually is max value.'''
def verify_max_global(spectra, returned_maxes, column):
	smaller_than_max = check_all_are_leq(spectra, returned_maxes, column)
	present = check_presence(spectra, returned_maxes, column)
	return not(all(smaller_than_max) and any(present))
	
'''Tests MassSpectrum's max_mass and MassSpectraAggregate's max_mass_local and max_mass_global.

	Generates arrays of random masses with intensities zero,
	(and sets up new MassSpectraAggregate wrapped around these arrays,)
	and checks arrays are all <= result.'''
def test_max_mass():

	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.MASS])

	if(verify_max_local(spectra, [spectrum.max_mass() for spectrum in spectra], MassSpectrum.MASS)):
		print("max_mass fails test!")
	if(verify_max_local(spectra, spectra_aggregate.max_mass_local(), MassSpectrum.MASS)):
		print("max_mass_local fails test!")
	if(verify_max_global(spectra, [spectra_aggregate.max_mass_global() for spectrum in spectra], MassSpectrum.MASS)):
		print("max_mass_global fails test!")

'''Tests MassSpectrum's max_intensity and MassSpectraAggregate's max_intensity_local and max_intensity_global.

	Generates arrays of random intensities with masses zero,
	(and sets up new MassSpectraAggregate wrapped around these arrays,)
	and checks arrays are all <= result.'''
def test_max_intensity():
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	
	if(verify_max_local(spectra, [spectrum.max_intensity() for spectrum in spectra], MassSpectrum.INTENSITY)):
		print("max_mass fails test!")
	if(verify_max_local(spectra, spectra_aggregate.max_intensity_local(), MassSpectrum.INTENSITY)):
		print("max_mass_local fails test!")
	if(verify_max_global(spectra, [spectra_aggregate.max_intensity_global() for spectrum in spectra], MassSpectrum.INTENSITY)):
		print("max_mass_global fails test!")
		
'''Test MassSpectrum's sort_by_mass and MassSpectraAggregate's sort_by_mass.
	
	Generates arrays of random masses with intensities zero,
	(and sets up new MassSpectraAggregate wrapped around these arrays,)
	calls sort_by_mass,
	then checks each index is less than previous.'''
def test_sort_by_mass():
	
	def verify_mass_sorted(spectra):
		return all([np.all(spectrum.ms2peaks[0:-1, MassSpectrum.MASS] <= spectrum.ms2peaks[1:, MassSpectrum.MASS]) for spectrum in spectra])
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.MASS])
	for spectrum in spectra: spectrum.sort_by_mass()
	if(not(verify_mass_sorted(spectra))):
		print("MassSpectrum's sort_by_mass fails test!")
		
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.MASS])
	spectra_aggregate.sort_by_mass()
	if(not(verify_mass_sorted(spectra))):
		print("MassSpectraAggregate's sort_by_mass fails test!")

'''Tests MassSpectrum's filter_intensity and MassSpectraAggregate's filter_intensity_local and filter_intensity_global.

	Generates arrays of random intensities with masses zero,
	(and sets up new MassSpectraAggregate wrapped around these arrays,)
	finds how many values are below threshold
	calls filter intensity,
	then checks length and that all values are above threshold.'''		
def test_filter_intensity():

	'''Helper to count starting lengths and count how many elements should be filtered before calling filter.'''
	def setup(intensity_thresholds):
		starting_lengths = [spectrum.ms2peaks.shape[0] for spectrum in spectra]
		below_threshold = [np.count_nonzero(spectrum.ms2peaks[:, MassSpectrum.INTENSITY] <= intensity_threshold) for spectrum, intensity_threshold in zip(spectra, intensity_thresholds)]
		return starting_lengths, below_threshold

	'''Helper to ensure correct number of elements were filtered out.'''
	def check_length(starting_lengths, below_threshold, spectra):
		zipped = zip(starting_lengths, below_threshold, spectra)
		return [starting_length == below + spectrum.ms2peaks.shape[0] for starting_length, below, spectrum in zipped]

	'''Helper to check all elements in some spectra are greater than some list of corresponding elements.'''
	def check_all_are_greater(spectra, returned_mins, column):
		return [all((spectrum.ms2peaks[:, column] > minm)) for spectrum, minm in zip(spectra, returned_mins)]
		
	'''Helper function to make the if-condition a little neater.'''
	def verify(starting_lengths, below_threshold, spectra, intensity_thresholds):
		return not(all(check_length(starting_lengths, below_threshold, spectra)) and all(check_all_are_greater(spectra, intensity_thresholds, MassSpectrum.INTENSITY)))
	
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	intensity_thresholds = [0.5 * spectrum.max_intensity() for spectrum in spectra] 
	starting_lengths, below_threshold = setup(intensity_thresholds)
	for spectrum, intensity_threshold in zip(spectra, intensity_thresholds): spectrum.filter_intensity(intensity_threshold=intensity_threshold)
	if(verify(starting_lengths, below_threshold, spectra, intensity_thresholds)):
		print("filter_intensity fails test!")
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	intensity_thresholds = [0.5 * spectrum.max_intensity() for spectrum in spectra] 
	starting_lengths, below_threshold = setup(intensity_thresholds)
	spectra_aggregate.filter_intensity_local(intensity_thresholds=intensity_thresholds)
	if(verify(starting_lengths, below_threshold, spectra, intensity_thresholds)):
		print("filter_intensity_local fails test!")
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	intensity_thresholds = [0.5 * spectra_aggregate.max_intensity_global() for spectrum in spectra]
	starting_lengths, below_threshold = setup(intensity_thresholds)	
	spectra_aggregate.filter_intensity_global(intensity_threshold=intensity_thresholds[0])
	if(verify(starting_lengths, below_threshold, spectra, intensity_thresholds)):
		print("filter_intensity_global fails test!")
	
'''Tests MassSpectrum's normalise_intensity, and MassSpectraAggregate's normalise_intensity_local and normalise_intensity_global.
	
	Generates arrays of random intensities with masses zero,
	(and sets up new MassSpectraAggregate wrapped around these arrays,)
	calls normalise intensity,
	then membership of new max and that all values are below new max.'''
def test_normalise_intensity():
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	new_maxes = [random.randint(100, 10000) for spectrum in spectra]
	for spectrum, new_max in zip(spectra, new_maxes): spectrum.normalise_intensity(new_max=new_max)
	if(verify_max_local(spectra, new_maxes, MassSpectrum.INTENSITY)):
		print("normalise_intensity fails test!")
		
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	new_maxes = [random.randint(100, 10000) for spectrum in spectra]
	spectra_aggregate.normalise_intensity_local(new_maxes=new_maxes)
	if(verify_max_local(spectra, new_maxes, MassSpectrum.INTENSITY)):
		print("normalise_intensity_local fails test!")
		
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	new_max = random.randint(100, 10000)
	spectra_aggregate.normalise_intensity_global(new_max=new_max)
	if(verify_max_global(spectra, [new_max for spectrum in spectra], MassSpectrum.INTENSITY)):
		print("normalise_intensity_global fails test!")
	
'''Tests MassSpectrum's static_mass_tolerance, percentile_mass_tolerance and ppm_mass_tolerance.'''
def test_mass_tolerance_calculations():

	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.MASS])
	
	percentiles = [0.05 * np.array(spectrum.ms2peaks) for spectrum in spectra] #we'll just use percentages for the static values
	mass_tolerances = [(spectrum.ms2peaks - percentile, spectrum.ms2peaks + percentile) for spectrum, percentile in zip(spectra, percentiles)]
	if(not(all([np.allclose(mass_tolerance, MassSpectrum.static_mass_tolerance(spectrum.ms2peaks, percentile))
				for mass_tolerance, spectrum, percentile in zip(mass_tolerances, spectra, percentiles)]))):
		print("static_mass_tolerance fails test!")
	if(not(all([np.allclose(mass_tolerance, MassSpectrum.percentile_mass_tolerance(spectrum.ms2peaks, 0.05))
				for mass_tolerance, spectrum in zip(mass_tolerances, spectra)]))):
		print("percentile_mass_tolerance fails test!")
		
	percentiles = [0.05 * maxm for maxm in spectra_aggregate.max_mass_local()]
	mass_tolerances = [(spectrum.ms2peaks - percentile, spectrum.ms2peaks + percentile) for spectrum, percentile in zip(spectra, percentiles)]
	if(not(all([np.allclose(mass_tolerance, spectrum.ppm_mass_tolerance(spectrum.ms2peaks, 0.05))
				for mass_tolerance, spectrum in zip(mass_tolerances, spectra)]))):
		print("ppm_mass_tolerance fails test!")
	
'''Tests MassSpectrum's find_longest_tag and MassSpectraAggregate's longest_tag_local and longest_tag_global (and the tag-finding methods they call).'''
def test_find_longest_tag():

	mass_threshold = 0.001 #used to add and detect noise (static value)
	
	def generate_spectrum():
		#generate the masses of a few mass tags to become a sequence tag, generate several 'sequence tags' this way
		printable_tags = [[str(random.choice(list(AA_mass_table.keys()))) for length in range(random.randint(2, 10))] for number in range(random.randint(1, 5))]
		tags = [[AA_mass_table[element] for element in tag] for tag in printable_tags]
		
		#first element in generated tag elements doesn't actually represent a difference, just a starting mass, so we exclude it from tag
		printable_tags = [[element for element in tag][1:] for tag in printable_tags] 
		
		#turn tags into rolling sums to generate a 'walk' of mass peaks
		tags = [np.cumsum(np.array(tag)) for tag in tags]
		
		#add some readings which shouldn't correspond to anything
		upper_bound = max([np.max(tag) for tag in tags])
		tags.append((upper_bound) * np.random.random(size=(random.randint(10, 20))))
		
		#now collapse into a single array and filter duplicates
		tags = np.unique(np.concatenate(tags))
	
		#add noise in range [mass_tolerance/2, mass_tolerance/2)
		tags += (mass_threshold) * np.random.random(size=tags.shape) - mass_threshold
	
		ms2peaks = np.zeros((tags.shape[0], 2))
		ms2peaks[:, MassSpectrum.MASS] = tags
		
		printable_tags = ["-" + "-".join(tag) + "-" for tag in printable_tags]
		
		return MassSpectrum(id=printable_tags, ms2peaks=ms2peaks, mass_table=AA_mass_table)
		
	def generate_spectra():
		spectra = [generate_spectrum() for i in range(random.randint(2, 5))]
		spectra_aggregate = MassSpectraAggregate(spectra)
		return spectra, spectra_aggregate

	
	spectra, spectra_aggregate = generate_spectra()
	spectra_tags = [spectrum.find_longest_tag(mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=mass_threshold) for spectrum in spectra]
	returned_tags = []
	for (long_tag, tags), spectrum in zip(spectra_tags, spectra):
		for id in spectrum.id:
			id_in_tag = [[id in tag for (tag, peaks) in value] for key, value in tags.items()]
			returned_tags.append(any([any(element) for element in id_in_tag]))
	if(not(all([longest_tag == 0 or longest_tag == max(tags.keys()) for longest_tag, tags in spectra_tags]) and all(returned_tags))):
		print("find_longest_tag fails test!")
	
	spectra, spectra_aggregate = generate_spectra()
	spectra_tags = spectra_aggregate.find_longest_tag_local(mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=mass_threshold)
	returned_tags = []
	for (long_tag, tags), spectrum in zip(spectra_tags, spectra):
		for id in spectrum.id:
			id_in_tag = [[id in tag for (tag, peaks) in value] for key, value in tags.items()]
			returned_tags.append(any([any(element) for element in id_in_tag]))
	if(not(all([longest_tag == 0 or longest_tag == max(tags.keys()) for longest_tag, tags in spectra_tags]) and all(returned_tags))):
		print("find_longest_tag_local fails test!")
		
	spectra, spectra_aggregate = generate_spectra()
	longest_tag, spectra_tags = spectra_aggregate.find_longest_tag_global(mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=mass_threshold)
	returned_tags = []
	for (long_tag, tags), spectrum in zip(spectra_tags, spectra):
		for id in spectrum.id:
			id_in_tag = [[id in tag for (tag, peaks) in value] for key, value in tags.items()]
			returned_tags.append(any([any(element) for element in id_in_tag]))
	tag_lengths = [long_tags for long_tags,tags in spectra_tags]
	if(not((longest_tag == 0 or longest_tag == max(tag_lengths)) and all(returned_tags))):
		print("find_longest_tag_global fails test!")
	
'''Helper function to automate generating dicts full of random 'sequence tags'.'''
def generate_tag_dict():

	def generate_random_string(length):
		return [str([random.choice(string.ascii_letters) for j in range(random.randint(1, 4))]) for i in range(length)]

	tags = {}
	for length in range(1, random.randint(4, 11)):
		tags[length] = [("-".join(generate_random_string(length)), 0) for number in range(random.randint(1, 6))]
	
	return tags
	
'''Tests MassSpectrum's expand_tag_notation.'''
def test_expand_tag_notation():
	pass
	'''tags = generate_tag_dict()
	expanded_tags = defaultdict(list)
	for length, spectra in tags.items():
		for tag, peaks in spectra:
			print("meme")
			print(list(eval(split_tag)) for split_tag in tag.split('-'))
			print(list(itertools.product(eval(split_tag) for split_tag in tag.split('-'))))
			#expanded_tags[length] += list(itertools.product(eval(split_tag) for split_tag in tag.split('-')))
		
	expanded_tags = dict(expanded_tags)
	
	print("tags:")
	print(tags)
	print("expanded_tags:")
	print(expanded_tags)
	if(not(MassSpectrum.expand_tag_notation(tags) == expanded_tags)):
		print("expand_tag_notation fails test!")'''

'''Tests MassSpectrum's flatten_tags.'''
def test_flatten_tags():
	tags = generate_tag_dict()		
	if(not(MassSpectrum.flatten_tags(tags) == [element for element in (tags[key] for key in tags.keys())])):
		print("flatten_tags fails test!")

'''Run tests. They should be randomised, but they will give an indication of whether something is obviously wrong.'''
def tests():
	test_sequence_by_id_helpers()
	test_max_mass()
	test_max_intensity()
	test_sort_by_mass()
	test_filter_intensity()
	test_normalise_intensity()
	test_mass_tolerance_calculations()
	test_find_longest_tag()
	test_expand_tag_notation()
	test_flatten_tags()

##########
###Main###
##########
		
def main():

	path = os.path.join(os.getcwd(), "spectraData")
	
	tests()
	
	mibig_parser()
	
	#SecArgGlnPheArgIle
	#...223959.ms
	#prev settings: 
	'''print_longest_tags(path=path, pattern="*.ms", 
						intensity_thresholds=0.005,
						mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=0.001)'''
	#current settings:
	'''print_longest_tags(path=path, pattern="*.ms", 
						intensity_thresholds=None)'''
						
	write_tags(path=path, pattern="*.ms", 
						intensity_thresholds=0.005,
						mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=0.001,
						lengths_to_print=range(4,100), filename="tags_0.001Mass_0.005Int.out")
						
	write_tags(path=path, pattern="*.ms", 
						intensity_thresholds=0.05,
						mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=0.01,
						lengths_to_print=range(4,100), filename="tags_0.01Mass_0.05Int.out")
						
	write_tags(path=path, pattern="*.ms", 
						intensity_thresholds=0.05,
						mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=0.001,
						lengths_to_print=range(4,100), filename="tags_0.001Mass_0.05Int.out")
						
	#Too long
	'''write_tags(path=path, pattern="*.ms", 
						intensity_thresholds=0.005,
						mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=0.01,
						lengths_to_print=range(4,100), filename="tags_0.01Mass_0.005Int.out")'''
						
	write_tags(path=path, pattern="*.ms", 
						intensity_thresholds=0.05,
						mass_threshold = 0.00001,
						lengths_to_print=range(4,100), filename="tags_10ppmMass_0.05Int.out")
	
	'''print("Parsing, converting to object...")
	mibig_spectrum_1 = setup_mass_spectra(parser.load_files_from_dir(path=path, pattern="CCMSLIB00000223959.ms")[0])
	filename, mis1 = mibig_spectrum_1
	print("Filtering intensity...")
	mis1.filter_intensity()
	#print(mis1.ppm_mass_tolerance(0, 0.00001))
	#print(mis1.ms2peaks.shape)
	print("Finding longest tag...")
	longest, tags = mis1.find_longest_tag()
	print(filename, longest)
	if(longest > 0):
		for tag in tags[longest]:
			print(tag)'''

if __name__ == "__main__":

    main()