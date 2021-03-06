import os
import glob
import numpy as np

'''Takes a directory and a pattern to find files by, 
parses them line by line according to parsing modes in dictionary (reads them as a single string by default)
and loads them into:
a list of objects containing:
a tuple of a filename and a dictionary, containing:
keys representing record headers in the file and various structures for values'''

##################
###Line Parsers###
##################
	
'''Parses as a single line of text'''
def default_parser(dict, key, split_line):
	dict[key] = " ".join(split_line)

'''Parses as a single input of type float.'''
def float_parser(dict, key, split_line):
	dict[key] = float(split_line[0])
	
'''Parses line as one in a list of arrays of size 2 of related data items separated by spaces'''
def ms2peaks_parser(dict, key, split_line):
	dict[key] += [np.array((split_line[0], split_line[1]), dtype=float)]
	
'''Parses line as a single instance of multiple data items, delimited by spaces'''
def list_parser(dict, key, split_line):
	dict[key] = split_line
	
##############################	
###Parser Lookup Dictionary###
##############################

'''Lookup dictionary for how to parse a field denoted by '>'name;
	these will continue to be used to parse new lines until a new field name is encountered;
	to extend the parser just write a new method and put it in the dictionary under the field's name.'''
parsing_modes = {	
					"parentmass": float_parser,
					"ms2peaks": ms2peaks_parser
				}
				
#########################
###Line Parser Cleanup###
#########################

'''Casts ms2peaks list to array at the end, once it's done being read in'''
def ms2peaks_cleanup(dict, key):
	dict[key] = np.array(dict[key], dtype=float)

###############################
###Cleanup Lookup Dictionary###
###############################

'''Lookup dictionary for how to clean up after parsing a field denoted by '>'name;
	this is called once at the end of parsing a field name, when a new field is entered;
	to extend the parser just write a new method and put it in the dictionary under the field's name.'''
cleanup_dict = {
					"ms2peaks": ms2peaks_cleanup
				}

################################				
###Whole File Parsing Methods###
################################

'''Given a file handle, parses the file line-by-line into a record,
	using '>' as an escape character to delineate field names,
	which are then used to access parsing methods for the dictionary above'''
def parse_file(file):

	record = {}
	line_parser = default_parser
	
	field = "" #names of data type

	for line in file:
		
		split_line = line.split()
		if(len(split_line) < 1): #empty line - next
			continue
	
		if(split_line[0][0] == '>'): #new field name
			
			if(field in cleanup_dict): #perform any necessary cleanup on last field before switching
				cleanup_dict[field](record, field)
		
			field = split_line[0][1:] #take first item in line, remove preceding '>' denoting new field name
			
			#set new field to store data under
			if(field in parsing_modes):
				record[field] = []
				line_parser = parsing_modes[field]
			else:
				line_parser = default_parser
				
			split_line = split_line[1:] #remove new field name from the line before passing off for parsing
			if(len(split_line) < 1): #empty line - next
				continue
			
		line_parser(record, field, split_line) #pass current line off for parsing
		
	if(field in cleanup_dict): #perform any necessary cleanup on last field before finishing with file
		cleanup_dict[field](record, field)
	
	return record

'''Takes a directory under path, and a file name pattern under pattern,
	and feeds them all to the file parsing method for parsing.'''
def load_files_from_dir(path=os.path.dirname(__file__), pattern="*.ms"):
	
	records = [] #each element is a dictionary representing a parsed file

	for filename in glob.glob(os.path.join(path, pattern)):

		file = open(os.path.join(path, filename), 'r')
		records += [(os.path.basename(filename),parse_file(file))]
		file.close()
		
	return records
	
##########
###Main###
##########

def main(parentpath=os.path.dirname(__file__)):
	print(load_files_from_dir(path=os.path.join(parentpath, "spectraData"), pattern="*.ms"))

if __name__ == "__main__":

    main(parentpath=os.getcwd())