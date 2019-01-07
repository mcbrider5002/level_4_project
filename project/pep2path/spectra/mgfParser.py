import os
import glob
import numpy as np

from .MassSpectrum import MassSpectrum
from .MassSpectraAggregate import MassSpectraAggregate

'''Takes a directory and a pattern to find files by, 
parses them line by line according to parsing modes in dictionary (reads them as a single string by default)
and loads them into a list of MassSpectraAggregate, one for each .mgf'''

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

##################
###Line Parsers###
##################
	
'''Parses as a single line of text'''
def default_parser(dict, key, values):
	dict[key] = values[0].strip()

'''Parses as a single input of type float.'''
def float_parser(dict, key, values):
	dict[key] = float(values[0])
	
'''Parses line as one in a list of arrays of size 2 of related data items separated by spaces'''
def ms2peaks_parser(dict, key, values):
	dict[key] = np.array([np.array((line.split()[0], line.split()[1])) for line in values[1:]], dtype=float)
	
'''Parses as a series of lines of text (i.e. unchanged)'''
def block_parser(dict, key, values):
	dict[key] = values
	
##############################	
###Parser Lookup Dictionary###
##############################

'''Lookup dictionary for how to parse a field denoted by '>'name;
	these will continue to be used to parse new lines until a new field name is encountered;
	to extend the parser just write a new method and put it in the dictionary under the field's name.'''
parsing_modes = {	
					"PEPMASS": float_parser,
					"SCANS": ms2peaks_parser
				}

################################				
###Whole File Parsing Methods###
################################

'''Given a file handle, parses the file line-by-line into a record,
	using '=' to delineate field names,
	which are then used to access parsing methods for the dictionary above'''
def parse_file(file):

	named_fields = ["TITLE", "PEPMASS", "CHARGE", "SCANS"]

	spectrum = {}
	spectra = []
	line_parser = default_parser
	
	field = "" #names of data type
	data = []
	in_spectrum = False
	
	for line in file:
		
		if("BEGIN IONS" in line):
			in_spectrum = True
			spectrum = {}
			continue
			
		if("END IONS" in line):
			in_spectrum = False
			if(field != ""):
				spectrum[field] = data
			field = ""
			data = []
			misc = {}
			
			if(len(spectrum["SCANS"]) > 1):
				for field, values in spectrum.items():
					if(field in named_fields):
						parsing_modes.get(field, default_parser)(spectrum, field, values)
					else:
						parsing_modes.get(field, default_parser)(misc, field, values)		
					
				spectra.append(MassSpectrum(id=spectrum.get("TITLE", ""),
										parent_mass=spectrum.get("PEPMASS", 0.0),
										ionization=spectrum.get("CHARGE", ""),
										ms2peaks=spectrum.get("SCANS", np.array([])),
										mass_table=AA_mass_table,
										misc=misc))
			
			continue
			
		if(in_spectrum):
			if('=' in line): #new field
				if(field != ""):
					spectrum[field] = data
				field = line.split('=')[0]
				line = line.split('=')[1] #cut off field name, but use the rest of the data
				data = []
			
			data.append(line)
	
	return MassSpectraAggregate(spectra)

'''Takes a directory under path, and a file name pattern under pattern,
	and feeds them all to the file parsing method for parsing.'''
def load_files_from_dir(path=os.path.dirname(__file__), pattern="*.ms"):
	
	spectra_aggregates = [] #each element is a filename and spectra aggregate of a single .mgf

	for filename in glob.glob(os.path.join(path, pattern)):
		file = open(os.path.join(path, filename), 'r')
		spectra_aggregates.append((os.path.basename(filename), parse_file(file)))
		file.close()
		
	return spectra_aggregates
	
##########
###Main###
##########

def main():
	return (load_files_from_dir(path=os.path.join(os.path.dirname(__file__), "spectraData"), pattern="*.mgf"))[0][1]

if __name__ == "__main__":

    main()
			