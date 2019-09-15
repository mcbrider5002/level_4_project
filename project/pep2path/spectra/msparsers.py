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
	
'''Parses data as list of pairs into np.array (ignoring the first line)'''
def ms2peaks_parser(dict, key, values):
	dict[key] = np.array([line.split()[:2] for line in values[1:]], dtype=float)
	
'''Leaves list of lines unchanged'''
def block_parser(dict, key, values):
	dict[key] = values

'''Mapping between headers and corresponding parser for data stored under them.'''
parsing_modes = {	
					"parentmass": float_parser,
					"ms2peaks": ms2peaks_parser,
					"PEPMASS": float_parser,
					"SCANS": ms2peaks_parser
				}

'''Convert names to standard names used inside MassSpectrum so they can be passed as kwargs to constructor.'''
def format_names(spectrum, name_dict):
	s = {"misc": {}}
	for k, v in spectrum.items():
		k = k.lower()
		if(k in name_dict): s[name_dict[k]] = v
		else: s["misc"][k] = v
	return s

'''Given the file handle of an .ms file, parses the file into a MassSpectrum'''
def parse_ms(id, file, mass_table=AA_mass_table):
	name_dict = {k:k for k in ["id", "compound", "formula", "ionization", "inchi", "smiles", "ms2peaks"]}
	name_dict.update({"parentmass": "parent_mass", "inchikey": "inchi_key"})
	
	spectrum = {}
	field = "" #name of current data type
	values = []
	def parse():
		parser = parsing_modes.get(field, default_parser)
		parser(spectrum, field, values)
	
	for line in filter(lambda ln: ln.strip() != "", file):
		if(line[0] == '>'): #new field
			if(field): parse()
			field = line.split()[0][1:]
			line = " ".join(line.split()[1:]) #cut off field name, but use the rest of the data
			values = []
		values.append(line)
	if(field): parse()
		
	kwargs = format_names(spectrum, name_dict)
	return MassSpectrum(id=id, mass_table=mass_table, **kwargs)
	
'''Given the file handle of a .mgf, parses the file into a MassSpectraAggregate using '=' to delineate field names'''
def parse_mgf(fname, file, mass_table=AA_mass_table):

	name_dict = {"pepmass": "parent_mass", "charge": "ionization", "scans": "ms2peaks"}

	spectrum = {}
	spectra = []
	field = "" #name of current data type
	values = []
	in_spectrum = False
	
	for line in filter(lambda ln: ln.strip() != "", file):
		
		if("BEGIN IONS" in line):
			in_spectrum = True
			spectrum = {}
			field = ""
			values = []
			
		elif("END IONS" in line):
			in_spectrum = False
			if(field): spectrum[field] = values
			
			if(len(spectrum["SCANS"]) > 1):
				for f, v in spectrum.items():
					parser = parsing_modes.get(f, default_parser)
					parser(spectrum, f, v)

				id = "{}\n{}".format(fname, spectrum.get("TITLE", "Scan Number: UNK"))
				kwargs = format_names(spectrum, name_dict)
				spectra.append(MassSpectrum(id=id, mass_table=mass_table, **kwargs))
			
		elif(in_spectrum):
			if('=' in line): #new field
				if(field): spectrum[field] = values
				field = line.split('=')[0]
				line = line.split('=')[1] #cut off field name, but use the rest of the data
				values = []
			values.append(line)
	
	return MassSpectraAggregate(spectra)

'''Parses files matching glob pattern with file parser in a directory specified with path '''
def load_files(file_parser, path, pattern):
	def get_file(filename):
		with open(os.path.join(path, filename), 'r') as file:
			fname = os.path.basename(filename)
			return (fname, file_parser(fname, file))
			
	return [get_file(filename) for filename in glob.glob(os.path.join(path, pattern))]
	
'''Loads ms files from dir at specified location: returns a list of tuples of (filename, MassSpectrum)'''
def load_files_from_dir(path, pattern="*.ms"):	
	return load_files(parse_ms, path, pattern)
	
'''Loads mgfs from dir at specified location: returns a list of tuples of (filename, MassSpectraAggregate)'''
def load_mgfs_from_dir(path, pattern="*.mgf"):
	return load_files(parse_mgf, path, pattern)

def main(path=os.path.join(os.path.dirname(__file__), "spectraData")):
	print("\n".join("{}\n{}".format(n, s) for (n, s) in load_files_from_dir(path, pattern="CCMSLIB00000001548.ms")))
	print((load_mgfs_from_dir(path=path, pattern="*.mgf")[0][1]).spectra[0])