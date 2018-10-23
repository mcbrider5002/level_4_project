import os

import massSpectraParser as parser
from MassSpectrum import MassSpectrum
from MassSpectraAggregate import MassSpectraAggregate

import numpy as np

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
	intensity_thresholds = [intensity*0.05 for intensity in spectra.max_intensity_local()] if intensity_thresholds is None else intensity_thresholds
	intensity_thresholds = [intensity*intensity_thresholds for intensity in spectra.max_intensity_local()] if isinstance(intensity_thresholds, float) else intensity_thresholds
	spectra.filter_intensity_local(intensity_thresholds=intensity_thresholds)
	print("Time taken: " + str(time.clock() - start))
	
	print("Finding longest tag...")
	start = time.clock()
	longest_tag, tags = spectra.find_longest_tag_global(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
	print("Time taken: " + str(time.clock() - start))
	
	return longest_tag, spectra.attach_spectra_ids(tags)
	
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
			
def write_tags(path=os.path.join(os.getcwd(), "spectraData"), pattern="*.ms", 
						intensity_thresholds=None,
						mass_tolerance_mode=None, mass_threshold=0.00001,
						lengths_to_print=None, output_path=os.getcwd()):
						
	longest_tag, spectra_tags = find_longest_tag(path=path, pattern=pattern, intensity_thresholds=intensity_thresholds, 
											mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
								
	lengths_to_print = [longest_tag] if lengths_to_print is None else lengths_to_print
	lengths_to_print = [lengths_to_print] if isinstance(lengths_to_print, str) else lengths_to_print
	lengths_to_print = list(filter(lambda length : length <= longest_tag and length > 0, lengths_to_print))
	
	os.chdir(output_path)
	file = open("tags.out", 'w')
	for spectrum in spectra_tags: #print the id of each spectrum, the current tag length, and all tags with that length
		for length in lengths_to_print:
			if(length in spectrum[1][1].keys()):
				file.write(str("---" + spectrum[0]) + " " + str(length) + "---\n")
				for tag in spectrum[1][1][length]:	
					file.write(str(tag) + "\n")
		
	file.close()
		
def test_max_mass():
	pass

def test_max_intensity():
	pass
		
def test_sort_by_mass():
	#generate numpy array of masses and intensities
	#call sort_by_mass
	#check each index is less than previous
	pass
	
def test_filter_intensity():
	pass
	
def test_normalise_intensity():
	pass
	
def test_mass_tolerance_calculations():
	pass
	
def test_find_longest_tag():
	pass
	
def test_expand_tag_notation():
	pass

def test_flatten():
	pass

'''Run tests. They should be randomised, but they will give an indication of whether something is obviously wrong.'''
def tests():
	pass
						
def main():

	path = os.path.join(os.getcwd(), "spectraData")
	
	#SecArgGlnPheArgIle
	#223959.ms
	#prev settings: 
	'''print_longest_tags(path=path, pattern="*.ms", 
						intensity_thresholds=0.005,
						mass_tolerance_mode=MassSpectrum().STATIC_MASS_TOLERANCE, mass_threshold=0.001)'''
	#current settings:
	'''print_longest_tags(path=path, pattern="*.ms", 
						intensity_thresholds=None)'''
						
	write_tags(path=path, pattern="*.ms", 
						intensity_thresholds=0.005,
						mass_tolerance_mode=MassSpectrum().STATIC_MASS_TOLERANCE, mass_threshold=0.001,
						lengths_to_print=range(4,100))
						
	'''write_tags(path=path, pattern="*.ms", 
						intensity_thresholds=None,
						mass_threshold = 0.00001,
						lengths_to_print=range(4,100))'''
	
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