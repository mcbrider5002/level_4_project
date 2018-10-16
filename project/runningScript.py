import os

import massSpectraParser as parser
from MassSpectrum import MassSpectrum
from MassSpectraAggregate import MassSpectraAggregate

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
					"Ile" : 113.1576,
					"Leu" : 113.1576,
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
	
	return MassSpectrum(compound=dict["compound"], 
						formula=dict["formula"], 
						parent_mass=dict["parentmass"], 
						ionization=dict["ionization"], 
						inchi=dict["InChI"], 
						inchi_key=dict["InChIKey"], 
						smiles=dict["smiles"], 
						ms2peaks=dict["ms2peaks"], 
						mass_table=AA_mass_table)

def main():

	path = os.path.join(os.getcwd(), "spectraData")
	'''file_tup = (parser.load_files_from_dir(path=path, pattern="CCMSLIB00000579925.ms"))[0]
	spectrum = setup_mass_spectra(file_tup)
	print(spectrum)
	
	file_tup = (parser.load_files_from_dir(path=path, pattern="CCMSLIB00000574963.ms"))[0]
	spectrum2 = setup_mass_spectra(file_tup)
	print(spectrum2)
	
	print(spectrum.max_intensity())
	print(spectrum2.max_intensity())'''
	
	'''spectrum.filter_intensity(intensity_threshold=40000)
	print(spectrum)
	spectrum2.filter_intensity()
	print(spectrum2)'''
	
	'''spectrum.normalise_intensity(intensity_threshold=1000)
	print(spectrum)
	spectrum2.normalise_intensity(intensity_threshold=88187.0)
	print(spectrum2)'''
	
	'''spectrum.filter_intensity(intensity_threshold=1)
	print(spectrum)
	print(spectrum.find_sequence_tags(mass_threshold=5))
	spectrum2.filter_intensity()
	print(spectrum2.find_sequence_tags(mass_threshold=0.05))
	print(spectrum2.find_longest_tag(mass_threshold=0.05))'''
	
	file_tup = (parser.load_files_from_dir(path=path, pattern="CCMSLIB00000001548.ms"))[0]
	spectrum3 = setup_mass_spectra(file_tup)
	spectrum3.normalise_intensity()
	max_tag, max_tag_length = spectrum3.find_longest_tag(mass_threshold=0.05)
	for tag_tup in spectrum3.find_sequence_tags(mass_threshold=0.05):
		tag, peaks = tag_tup
		if(max_tag_length <= len(peaks) - 1):
			print(tag, peaks)
	
	'''print("Starting to parse...")
	records = parser.load_files_from_dir(path=path)
	print("Converting to objects...")
	for index in range(len(records)):
		records[index] = setup_mass_spectra(records[index])
	print("Filtering intensity...")
	for record in records:
		record.filter_intensity()
	print("Finding longest tag...")
	longest_tag = ""
	longest_tag_length = 0
	for record in records:
		tag, tag_length = record.find_longest_tag()
		if(tag_length > longest_tag_length):
			longest_tag = tag
			longest_tag_length = tag_length
	print(longest_tag, longest_tag_length)'''

if __name__ == "__main__":

    main()