import os

from .msparsers import load_files_from_dir
from .MassSpectrum import MassSpectrum
from .MassSpectraAggregate import MassSpectraAggregate
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
					  intensity_thresholds=None, mass_tolerance_mode=None, mass_threshold=0.00001):
	
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
						intensity_thresholds=None, mass_tolerance_mode=None, mass_threshold=0.00001):
						
	longest_tag, spectra_tags = find_longest_tag(path=path, pattern=pattern, intensity_thresholds=intensity_thresholds, 
												  mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
				
	#print the id of each spectrum, its longest tag length, and all tags with the longest tag length
	for spectrum_tags in (s_t for s_t in spectra_tags if s_t.longest_tag == longest_tag):
		print(spectrum_tags.id, spectrum_tags.longest_tag)
		for tag in spectrum_tags.tags[longest_tag]:	print(tag)
	
'''Writes to "tags.out" on the specifed path all tags with the specified length and their spectrum id.'''	
def write_tags(path=os.path.join(os.path.dirname(__file__), "spectraData"), pattern="*.ms", 
				intensity_thresholds=None, mass_tolerance_mode=None, mass_threshold=0.00001,
				lengths_to_print=None, output_path=os.path.join(os.path.dirname(__file__), "out"), 
				filename="tags.out"):
						
	longest_tag, spectra_tags = find_longest_tag(path=path, pattern=pattern, intensity_thresholds=intensity_thresholds, 
												  mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
								
	lengths_to_print = [longest_tag] if lengths_to_print is None else lengths_to_print #default value
	lengths_to_print = [lengths_to_print] if isinstance(lengths_to_print, str) else lengths_to_print #take string input
	lengths_to_print = [length for length in lengths_to_print if (length <= longest_tag and length > 0)] #remove impossible indices
	
	#print the id of each spectrum, the current tag length, and all tags with that length
	with open(os.path.join(output_path, filename), 'w') as file:
		for spectrum_tags in spectra_tags: 
			for length in lengths_to_print:
				if(length in spectrum_tags.tags):
					file.write("---{} {}---\n".format(spectrum_tags.id, length))
					for tag in spectrum_tags.tags[length]: file.write(str(tag) + "\n")
	
def main(path=os.path.join(os.path.dirname(__file__), "spectraData")):
	
	def tag_write(ithres, mthres):
		fname = "tags_{}Mass_{}Int.out".format(mthres, ithres)
		
		write_tags(path=path, pattern="*.ms", intensity_thresholds=ithres, 
					mass_tolerance_mode=MassSpectrum.STATIC_MASS_TOLERANCE,
					mass_threshold=mthres, lengths_to_print=range(4,100), filename=fname)
						
	tag_write(0.005, 0.001)					
	tag_write(0.05, 0.01)
	tag_write(0.05, 0.001)				
	#Runtime too long
	'''tag_write(0.005, 0.01)'''
						
	write_tags(path=path, pattern="*.ms", intensity_thresholds=0.05,
				mass_tolerance_mode=MassSpectrum.MAX_PPM_MASS_TOLERANCE, mass_threshold=0.00001, 
				lengths_to_print=range(4,100), filename="tags_10ppmMass_0.05Int.out")

if __name__ == "__main__":
    main(cdir=os.getcwd())