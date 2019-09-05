import os
import glob
import numpy as np
from .MassSpectrum import MassSpectrum
from .MassSpectraAggregate import MassSpectraAggregate
	
'''Parses as a single line of text'''
def default_parser(dict, key, split_line):
	dict[key] = " ".join(split_line)

'''Parses as a single input of type float.'''
def float_parser(dict, key, split_line):
	dict[key] = float(split_line[0])
	
'''Parses line as one in a list of arrays of size 2 of related data items separated by spaces'''
def ms2peaks_parser(dict, key, split_line):
	dict[key].append(split_line[:2])
	
'''Parses line as a single instance of multiple data items, delimited by spaces'''
def list_parser(dict, key, split_line):
	dict[key] = split_line

'''Mapping between names of headers '>name' and the functions used to parse all data under them.'''
parsing_modes = {	
					"parentmass": float_parser,
					"ms2peaks": ms2peaks_parser
				}

'''Casts ms2peaks list to array at the end, once it's done being read in'''
def ms2peaks_cleanup(dict, key):
	dict[key] = np.array(dict[key], dtype=float)

'''Lookup dictionary for how to clean up after parsing a field denoted by '>'name;
	this is called once at the end of parsing a field name, when a new field is entered;
	to extend the parser just write a new method and put it in the dictionary under the field's name.'''
cleanup_dict = {
					"ms2peaks": ms2peaks_cleanup
			   }

'''Given a file handle, parses the file line-by-line into a record,
	using '>' as an escape character to delineate field names,
	which are then used to access parsing methods for the dictionary above'''
def parse_file(id, file):

	record = {}
	line_parser = default_parser
	
	field = "" #names of data type

	for line in file:
		
		split_line = line.split()
		if(not split_line): continue #empty line - next
			
		if(split_line[0][0] == '>'): #new field name
			
			if(field in cleanup_dict): #perform any necessary cleanup on last field before switching
				cleanup_dict[field](record, field)
		
			field = split_line[0][1:] #take first item in line, remove preceding '>' denoting new field name
			
			record[field] = []
			line_parser = parsing_modes.get(field, default_parser)
				
			split_line = split_line[1:] #remove new field name from the line before passing off for parsing
			if(not split_line): continue #empty line - next
			
		line_parser(record, field, split_line) #pass current line off for parsing
		
	if(field in cleanup_dict): #perform any necessary cleanup on last field before finishing with file
		cleanup_dict[field](record, field)
	
	return record
	
'''Takes a directory under path, and a file name pattern under pattern,
	and feeds them all to the file parsing method for parsing.'''
def load_files_from_dir(path=os.path.dirname(__file__), pattern="*.ms"):
	
	records = [] #each element is a dictionary representing a parsed file

	for filename in glob.glob(os.path.join(path, pattern)):
		with open(os.path.join(path, filename), 'r') as file:
			fname = os.path.basename(filename)
			records.append((fname, parse_file(fname, file)))
	return records

def main(path=os.path.join(os.path.dirname(__file__), "spectraData")):
	print(load_files_from_dir(path=path, pattern="*.ms"))

if __name__ == "__main__":
    main(path=os.getcwd())