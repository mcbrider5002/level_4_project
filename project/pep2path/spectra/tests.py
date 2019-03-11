import os
import numpy as np
import random
import string
import itertools
from collections import Counter

from .MassSpectrum import MassSpectrum
from .MassSpectraAggregate import MassSpectraAggregate
from .Tag import Tag
from .SpectrumTags import SpectrumTags

from .masstables import AA_mass_table, AA_alphabet

###############
###Tag Tests###
###############

'''Helper to generate a list of lists of random components.'''
def generate_comp_lists(alphabet):
	return [[str(random.choice(list(alphabet))) for length in range(random.randint(20, 40))] for number in range(random.randint(10, 20))]

'''Helper to generate a list of random tags.'''	
def generate_tag(alphabet):
	generated = generate_comp_lists(alphabet)
	strings = ["-".join(g) for g in generated]
	tags = [Tag(s, [], []) for s in strings]
	return generated, strings, tags  
	
'''Helper to generate a dict of random tags, stored by length.'''	
def generate_tag_dict(alphabet):
	generated, strings, tags = generate_tag(alphabet)
	lengths = [len(g) for g in generated]
	longest = max(lengths)
	tag_dict = {length:[t for l, t in zip(lengths, tags) if l == length] for length in range(1, longest + 1) if length in lengths}
	return generated, strings, tags, lengths, longest, tag_dict
		   
def test_decompose_tag(alphabet=AA_alphabet):

	generated, strings, tags, lengths, longest, tag_dict = generate_tag_dict(alphabet)
	
	string = strings[0]
	tag = tags[0]
	assert tag.decompose_tag() == generated[0] , "decompose_tag doesn't match pregenerated tag!"
	assert Tag("-%s-" % string, [], []).decompose_tag() == generated[0], "decompose_tag (on string with trailing '-') doesn't match pregenerated tag!"
	
	s_tags = SpectrumTags("dummy_id", longest, tag_dict)
	
	tag_dict2 = {length:[Tag("-%s-" % tag.tag, [], []) for tag in tags] for length, tags in tag_dict.items()}
	s_tags2 = SpectrumTags("dummy_id", longest, tag_dict2)
	
	tdict = {length:[g for g in generated if len(g) == length] for length in range(1, longest + 1) if length in lengths}
	assert (s_tags.decompose_tags() == tdict), "decompose_tags doesn't match pregenerated tags!"
	assert (s_tags2.decompose_tags() == tdict), "decompose_tags (on string with trailing '-') doesn't match pregenerated tags!"
	
def test_component_counts(alphabet=AA_alphabet):
	
	generated, _, tags, _, longest, tag_dict = generate_tag_dict(alphabet)
	counts = Counter(itertools.chain.from_iterable(generated))
	
	c = Counter(generated[0])
	tag = tags[0]
	assert tag.component_counts() == c, "component_counts doesn't match component counts of pregenerated tag!"
	
	s_tags = SpectrumTags("dummy_id", longest, tag_dict)
	assert s_tags.multi_component_counts() == counts, "multi_component_counts doesn't match component counts of pregenerated tag!"
	
def test_unique_components(alphabet=AA_alphabet):

	generated, _, tags, lengths, longest, tag_dict = generate_tag_dict(alphabet)
	uniques = list(map(set, generated))
	
	unique = uniques[0]
	tag = tags[0]
	assert set(tag.unique_components()) == unique, "Tag's unique_components doesn't match components of pregenerated tag!"
	
	s_tags = SpectrumTags("dummy_id", longest, tag_dict)
	assert set(s_tags.unique_components()) == set(itertools.chain.from_iterable(uniques)), "SpectrumTags' unique_components doesn't match components of pregenerated tags!"

'''Helper to ensure ordering of tags doesn't matter.'''
def srt(s): sorted(s, key = lambda tag : tag.tag)
	
def test_filter_by_length(alphabet=AA_alphabet):
	_, strings, tags, lengths, longest, tag_dict = generate_tag_dict(alphabet)
	s_tags = SpectrumTags("dummy_id", longest, tag_dict)
	
	assert all([srt(s_tags.filter_by_length(length)) == srt([t for l, t in zip(lengths, tags) if l == length]) for length in range(longest)]), \
			"filter_by_length doesn't match pregenerated tags filtered by length!"
		   
def test_expand_tag_notation(alphabet=AA_alphabet):
	tags = [[random.choice(list(AA_alphabet)) for i in range(4)] for i in range(4)]
	true_tags = ["-".join(tag) for tag in tags]
	uncertain_tag = [str(list(comps)) for comps in zip(*tags)]
	
	e_tags = Tag("-".join(uncertain_tag), [], []).expand_tag_notation()
	assert all([tag in e_tags for tag in true_tags]), "expand_tag_notation did not have the true tags within the expanded uncertain tag!"
	
	e_tags = Tag("-%s-" % "-".join(uncertain_tag), [], []).expand_tag_notation()
	assert all([tag in e_tags for tag in true_tags]), "expand_tag_notation did not have the true tags within the expanded uncertain tag (with dashes at the ends)!"
		
def test_flatten_tags(alphabet=AA_alphabet):
	_, _, tags, _, longest, tag_dict = generate_tag_dict(alphabet)
	s_tags = SpectrumTags("dummy_id", longest, tag_dict)
	assert srt(s_tags.flatten_tags()) == srt(tags), "flatten_tags doesn't match pregenerated tags!"

###################
###Spectra Tests###
###################	

def test_mappings(trials=1000):
	spectrum = MassSpectrum()
	array = np.array(range(trials))
	for i in range(5):
		indices = random.sample(range(trials), trials)
		spectrum.update_mappings(indices)
		array = array[indices]
	assert list(array) == list(spectrum.original_indices(range(trials))), "MassSpectrum does not map indices correctly!"
	
'''Tests correctness of MassSpectraAggregate's get_spectra_ids, attach_spectra_ids and get_spectra_by_id.'''
def test_sequence_by_id_helpers():

	#generate a random number of random ids of random length
	ids = ["".join([random.choice(string.ascii_letters) for length in range(random.randint(1, 6))]) for number in range(random.randint(1, 6))] 
	
	spectra = MassSpectraAggregate([MassSpectrum(id=id) for id in ids]) #convert to spectra
	assert (ids == spectra.get_spectra_ids()), "get_spectra_ids fails test!"
	
	combined = [(id, index) for id, index in zip(ids, range(len(ids)))]
	assert (combined == spectra.attach_spectra_ids(range(len(ids)))), "attach_spectra_ids fails test!"
		
	assert (ids == [spectrum.id for spectrum in spectra.get_spectra_by_id(ids)]), "get_spectra_by_id fails test!"
	
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
	return all(smaller_than_max) and all(present)

'''Helper to check value actually is max value.'''
def verify_max_global(spectra, returned_maxes, column):
	smaller_than_max = check_all_are_leq(spectra, returned_maxes, column)
	present = check_presence(spectra, returned_maxes, column)
	return all(smaller_than_max) and any(present)
	
'''Tests MassSpectrum's max_mass and MassSpectraAggregate's max_mass_local and max_mass_global.

	Generates arrays of random masses with intensities zero,
	(and sets up new MassSpectraAggregate wrapped around these arrays,)
	and checks arrays are all <= result.'''
def test_max_mass():

	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.MASS])

	assert verify_max_local(spectra, [spectrum.max_mass() for spectrum in spectra], MassSpectrum.MASS), "max_mass fails test!"
	assert verify_max_local(spectra, spectra_aggregate.max_mass_local(), MassSpectrum.MASS), "max_mass_local fails test!"
	assert verify_max_global(spectra, [spectra_aggregate.max_mass_global() for spectrum in spectra], MassSpectrum.MASS), "max_mass_global fails test!"

'''Tests MassSpectrum's max_intensity and MassSpectraAggregate's max_intensity_local and max_intensity_global.

	Generates arrays of random intensities with masses zero,
	(and sets up new MassSpectraAggregate wrapped around these arrays,)
	and checks arrays are all <= result.'''
def test_max_intensity():
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	
	assert verify_max_local(spectra, [spectrum.max_intensity() for spectrum in spectra], MassSpectrum.INTENSITY), "max_mass fails test!"
	assert verify_max_local(spectra, spectra_aggregate.max_intensity_local(), MassSpectrum.INTENSITY), "max_mass_local fails test!"
	assert verify_max_global(spectra, [spectra_aggregate.max_intensity_global() for spectrum in spectra], MassSpectrum.INTENSITY), "max_mass_global fails test!"
		
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
	assert verify_mass_sorted(spectra), "MassSpectrum's sort_by_mass fails test!"
		
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.MASS])
	spectra_aggregate.sort_by_mass()
	assert verify_mass_sorted(spectra), "MassSpectraAggregate's sort_by_mass fails test!"

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
		return all(check_length(starting_lengths, below_threshold, spectra)) and all(check_all_are_greater(spectra, intensity_thresholds, MassSpectrum.INTENSITY))
	
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	intensity_thresholds = [0.5 * spectrum.max_intensity() for spectrum in spectra] 
	starting_lengths, below_threshold = setup(intensity_thresholds)
	for spectrum, intensity_threshold in zip(spectra, intensity_thresholds): spectrum.filter_intensity(intensity_threshold=intensity_threshold)
	assert verify(starting_lengths, below_threshold, spectra, intensity_thresholds), "filter_intensity fails test!"
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	intensity_thresholds = [0.5 * spectrum.max_intensity() for spectrum in spectra] 
	starting_lengths, below_threshold = setup(intensity_thresholds)
	spectra_aggregate.filter_intensity_local(intensity_thresholds=intensity_thresholds)
	assert verify(starting_lengths, below_threshold, spectra, intensity_thresholds), "filter_intensity_local fails test!"
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	intensity_thresholds = [0.5 * spectra_aggregate.max_intensity_global() for spectrum in spectra]
	starting_lengths, below_threshold = setup(intensity_thresholds)	
	spectra_aggregate.filter_intensity_global(intensity_threshold=intensity_thresholds[0])
	assert verify(starting_lengths, below_threshold, spectra, intensity_thresholds), "filter_intensity_global fails test!"
	
'''Tests MassSpectrum's normalise_intensity, and MassSpectraAggregate's normalise_intensity_local and normalise_intensity_global.
	
	Generates arrays of random intensities with masses zero,
	(and sets up new MassSpectraAggregate wrapped around these arrays,)
	calls normalise intensity,
	then membership of new max and that all values are below new max.'''
def test_normalise_intensity():
	
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	new_maxes = [random.randint(100, 10000) for spectrum in spectra]
	for spectrum, new_max in zip(spectra, new_maxes): spectrum.normalise_intensity(new_max=new_max)
	assert verify_max_local(spectra, new_maxes, MassSpectrum.INTENSITY), "normalise_intensity fails test!"
		
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	new_maxes = [random.randint(100, 10000) for spectrum in spectra]
	spectra_aggregate.normalise_intensity_local(new_maxes=new_maxes)
	assert verify_max_local(spectra, new_maxes, MassSpectrum.INTENSITY), "normalise_intensity_local fails test!"
		
	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.INTENSITY])
	new_max = random.randint(100, 10000)
	spectra_aggregate.normalise_intensity_global(new_max=new_max)
	assert verify_max_global(spectra, [new_max for spectrum in spectra], MassSpectrum.INTENSITY), "normalise_intensity_global fails test!"
	
'''Tests MassSpectrum's static_mass_tolerance, percentile_mass_tolerance and ppm_mass_tolerance.'''
def test_mass_tolerance_calculations():

	spectra, spectra_aggregate = generate_random_spectra([MassSpectrum.MASS])
	
	percentiles = [0.05 * np.array(spectrum.ms2peaks) for spectrum in spectra] #we'll just use percentages for the static values
	mass_tolerances = [(spectrum.ms2peaks - percentile, spectrum.ms2peaks + percentile) for spectrum, percentile in zip(spectra, percentiles)]
	
	assert (all([np.allclose(mass_tolerance, MassSpectrum.static_mass_tolerance(spectrum.ms2peaks, percentile)) 
						for mass_tolerance, spectrum, percentile in zip(mass_tolerances, spectra, percentiles)])), \
			"static_mass_tolerance fails test!"
			
	assert (all([np.allclose(mass_tolerance, MassSpectrum.rel_ppm_mass_tolerance(spectrum.ms2peaks, 0.05))
						for mass_tolerance, spectrum in zip(mass_tolerances, spectra)])), \
		   "percentile_mass_tolerance fails test!"
		
	percentiles = [0.05 * maxm for maxm in spectra_aggregate.max_mass_local()]
	mass_tolerances = [(spectrum.ms2peaks - percentile, spectrum.ms2peaks + percentile) for spectrum, percentile in zip(spectra, percentiles)]
	
	assert (all([np.allclose(mass_tolerance, spectrum.max_ppm_mass_tolerance(spectrum.ms2peaks, 0.05))
				for mass_tolerance, spectrum in zip(mass_tolerances, spectra)])), \
		   "ppm_mass_tolerance fails test!"
	
'''Tests MassSpectrum's find_sequence_tags and MassSpectraAggregate's longest_tag.'''
def test_find_longest_tag(mass_table=AA_mass_table):

	mass_threshold = 0.001 #used to add and detect noise (static value)
	
	'''Helper to randomly generate tags in a spectrum, with extra peaks corresponding to nothing, and noise.'''
	def generate_spectrum():
		#generate the masses of a few mass tags to become a sequence tag, generate several 'sequence tags' this way
		printable_tags = [[str(random.choice(list(mass_table.keys()))) for length in range(random.randint(2, 10))] for number in range(random.randint(1, 5))]
		tags = [[mass_table[element] for element in tag] for tag in printable_tags]
		
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
		
		return MassSpectrum(ms2peaks=ms2peaks, mass_table=mass_table, misc={"printable_tags":printable_tags})
		
	'''Helper to generate spectra and spectra aggregate.'''
	def generate_spectra():
		spectra = [generate_spectrum() for i in range(random.randint(2, 5))]
		spectra_aggregate = MassSpectraAggregate(spectra)
		return spectra, spectra_aggregate
		
	'''Helper to check whether all generated tags are contained within the output from the tag finder for their respective spectrum.'''
	def search_for_tags(spectra_tags, spectra):
		returned_tags = []
		for spectrum_tags, spectrum in zip(spectra_tags, spectra):
			for ptag in spectrum.misc["printable_tags"]:
				#check if one of tags generated are a subsequence of the tags returned from their respective spectrum
				ptag_in_tag = [[(ptag in tag.tag) for tag in tag_list] for length, tag_list in spectrum_tags.tags.items()]
				#reduce whether tag was found in tag dict to a single boolean
				returned_tags.append(any([any(element) for element in ptag_in_tag]))
		return returned_tags
	
	spectra, spectra_aggregate = generate_spectra()
	spectra_tags = [spectrum.find_sequence_tags(mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=mass_threshold) for spectrum in spectra]
	returned_tags = search_for_tags(spectra_tags, spectra)
	assert (all([spectrum_tags.longest_tag == 0 or spectrum_tags.longest_tag == max(spectrum_tags.tags.keys()) 
						for spectrum_tags in spectra_tags]) and all(returned_tags)), \
		   "find_sequence_tags fails test!"
		
	spectra, spectra_aggregate = generate_spectra()
	longest_tag, spectra_tags = spectra_aggregate.find_longest_tag(mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE, mass_threshold=mass_threshold)
	returned_tags = search_for_tags(spectra_tags, spectra)
	tag_lengths = [spectrum_tags.longest_tag for spectrum_tags in spectra_tags]
	assert ((longest_tag == 0 or longest_tag == max(tag_lengths)) and all(returned_tags)), \
		   "find_longest_tag fails test!"

'''Run tests. They should be randomised, but they will give an indication of whether something is obviously wrong.'''
def tests():
	test_mappings()
	test_decompose_tag()
	test_component_counts()
	test_unique_components()
	test_filter_by_length()
	test_expand_tag_notation()
	test_flatten_tags()
	
	test_sequence_by_id_helpers()
	test_max_mass()
	test_max_intensity()
	test_sort_by_mass()
	test_filter_intensity()
	test_normalise_intensity()
	test_mass_tolerance_calculations()
	test_find_longest_tag()