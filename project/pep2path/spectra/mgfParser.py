import os
import glob
import numpy as np

from .MassSpectrum import MassSpectrum
from .MassSpectraAggregate import MassSpectraAggregate

from .masstables import AA_mass_table
	
'''Parses as a single line of text'''
def default_parser(dict, key, values):
	dict[key] = values[0].strip()

'''Parses as a single input of type float.'''
def float_parser(dict, key, values):
	dict[key] = float(values[0])
	
'''Parses line as one in a list of arrays of size 2 of related data items separated by spaces'''
def ms2peaks_parser(dict, key, values):
	dict[key] = np.array([line.split[:2] for line in values[1:]], dtype=float)
	
'''Parses as a series of lines of text (i.e. unchanged)'''
def block_parser(dict, key, values):
	dict[key] = values

'''Lookup dictionary for how to parse a field denoted by '>'name;
	these will continue to be used to parse new lines until a new field name is encountered;
	to extend the parser just write a new method and put it in the dictionary under the field's name.'''
parsing_modes = {	
					"PEPMASS": float_parser,
					"SCANS": ms2peaks_parser
				}

'''Given a file handle, parses the file line-by-line into a record,
	using '=' to delineate field names,
	which are then used to access parsing methods for the dictionary above'''
def parse_file(file):

	named_fields = ["TITLE", "PEPMASS", "CHARGE", "SCANS"]

	spectrum = {}
	spectra = []
	
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
			if(field):
				spectrum[field] = data
			field = ""
			data = []
			misc = {}
			
			if(len(spectrum["SCANS"]) > 1):
				for field, values in spectrum.items():
					out_dict = spectrum if (field in named_fields) else misc
					(parsing_modes.get(field, default_parser))(out_dict, field, values)
					
				spectra.append(MassSpectrum(id=spectrum.get("TITLE", ""),
										parent_mass=spectrum.get("PEPMASS", 0.0),
										ionization=spectrum.get("CHARGE", ""),
										ms2peaks=spectrum.get("SCANS", np.array([])),
										mass_table=AA_mass_table,
										misc=misc))
			continue
			
		if(in_spectrum):
			if('=' in line): #new field
				if(field):
					spectrum[field] = data
				field = line.split('=')[0]
				line = line.split('=')[1] #cut off field name, but use the rest of the data
				data = []
			data.append(line)
	
	return MassSpectraAggregate(spectra)

'''Takes a directory under path, and a file name pattern under pattern,
	and feeds them all to the file parsing method for parsing.'''
def load_files_from_dir(path=os.path.join(os.path.dirname(__file__), "spectraData"), pattern="*.mgf"):
	
	spectra_aggregates = [] #each element is a filename and spectra aggregate of a single .mgf

	for filename in glob.glob(os.path.join(path, pattern)):
		with open(os.path.join(path, filename), 'r') as file:
			spectra_aggregates.append((os.path.basename(filename), parse_file(file)))
	return spectra_aggregates

def main(path=os.path.join(os.path.dirname(__file__), "spectraData")):
	print(load_files_from_dir(path=path, pattern="*.mgf")[0][1])
	
if __name__ == "__main__":
    main(path=os.getcwd())