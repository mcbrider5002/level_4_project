import os
import numpy as np
import itertools

from .massSpectraParser import load_files_from_dir
from .MassSpectrum import MassSpectrum
from .MassSpectraAggregate import MassSpectraAggregate
from .Tag import Tag
from .SpectrumTags import SpectrumTags
from .tests import tests as tests

from .masstables import AA_mass_table
				
def setup_mass_spectra(file, mass_table=AA_mass_table):
	
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
						mass_table=mass_table)
						
'''Convenience function to perform the multiple steps necessary to find tags and their ids.'''
def find_longest_tag(path=os.path.join(os.path.dirname(__file__), "spectraData"), pattern="*.ms", 
						intensity_thresholds=None,
						mass_tolerance_mode=None, mass_threshold=0.00001):
	
	import time
	
	print("Starting to parse for pattern \"" + pattern + "\"...")
	start = time.clock()
	records = load_files_from_dir(path=path, pattern=pattern)
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
	longest_tag, tags = spectra.find_longest_tag(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
	print("Time taken: " + str(time.clock() - start))
	
	print("Longest tag found: " + str(longest_tag) + "\n")
	
	return longest_tag, tags
	
'''Prints out all the tied-longest tags and their spectrum id.'''
def print_longest_tags(path=os.path.join(os.path.dirname(__file__), "spectraData"), pattern="*.ms", 
						intensity_thresholds=None,
						mass_tolerance_mode=None, mass_threshold=0.00001):
						
	longest_tag, spectra_tags = find_longest_tag(path=path, pattern=pattern, intensity_thresholds=intensity_thresholds, 
											mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
				
	spectra_tags = filter((lambda spectrum_tags : spectrum_tags.longest_tag == longest_tag), spectra_tags) #first get all spectra with the maximum tag lengths
	for spectrum_tags in spectra_tags: #print the id of each spectrum, its longest tag length, and all tags with the longest tag length (which should be identical to local longest)
		print(spectrum_tags.id, spectrum_tags.longest_tag)
		for tag in spectrum_tags.tags[longest_tag]:	
			print(tag)
	
'''Writes to "tags.out" on the specifed path all tags with the specified length and their spectrum id.'''	
def write_tags(path=os.path.join(os.path.dirname(__file__), "spectraData"), pattern="*.ms", 
						intensity_thresholds=None,
						mass_tolerance_mode=None, mass_threshold=0.00001,
						lengths_to_print=None, output_path=os.path.join(os.path.dirname(__file__), "out"),
						filename="tags.out"):
						
	longest_tag, spectra_tags = find_longest_tag(path=path, pattern=pattern, intensity_thresholds=intensity_thresholds, 
											mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
								
	lengths_to_print = [longest_tag] if lengths_to_print is None else lengths_to_print #default value
	lengths_to_print = [lengths_to_print] if isinstance(lengths_to_print, str) else lengths_to_print #take string input
	lengths_to_print = list(filter(lambda length : length <= longest_tag and length > 0, lengths_to_print)) #remove impossible indices
	
	file = open(os.path.join(output_path, filename), 'w')
	for spectrum_tags in spectra_tags: #print the id of each spectrum, the current tag length, and all tags with that length
	
		for length in lengths_to_print:
		
			if(length in spectrum_tags.tags.keys()):
			
				file.write(str("---" + spectrum_tags.id) + " " + str(length) + "---\n")
				
				for tag in spectrum_tags.tags[length]:	
					file.write(str(tag) + "\n")
		
	file.close()
	
def main():

	path = os.path.join(os.path.dirname(__file__), "spectraData")
	
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
	mibig_spectrum_1 = setup_mass_spectra(load_files_from_dir(path=path, pattern="CCMSLIB00000223959.ms")[0])
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