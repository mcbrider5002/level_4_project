import os
import glob

#data structure is list of objects containing:
#a tuple of a filename and a dictionary, containing:
#keys representing record headers in the file and various structures for values

###	
	
#parses as a single line of text
def default_parser(dict, key, split_line):
	dict[key] = " ".join(split_line)

#parses line as one in a list of pairs of related data items separated by spaces
def ms2peaks_parser(dict, key, split_line):
	
	dict[key] += [(split_line[0], split_line[1])]
	
#parses line as a single instance of multiple data items, delimited by spaces
def list_parser(dict, key, split_line):

	dict[key] = split_line

###
	
parsing_modes = {	
					"ms2peaks": ms2peaks_parser
				}
				
###

def parse_file(file):

	record = {}
	line_parser = default_parser
	
	field = "" #names of data type

	for line in file:
		
		split_line = line.split()
		if(len(split_line) < 1): #empty line - next
			continue
	
		if(split_line[0][0] == '>'): #new field name
		
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
			
		line_parser(record, field, split_line)
	
	return record

def load_files_from_dir(path=os.getcwd(), pattern="*.ms"):
	
	os.chdir(path) #changes pwd if path is set
	
	records = [] #each element is a dictionary representing a parsed file
	
	for filename in glob.glob(pattern):
	
		file = open(filename, 'r')
		records += [(filename,parse_file(file))]
		file.close()
		
	return records
	
###
		
def main():

	import time
	times = []
	
	start = time.clock()
	print(load_files_from_dir(pattern="CCMSLIB00000579925.ms"))
	times += ["Time for small file: " + str(time.clock() - start)]
	
	start = time.clock()
	load_files_from_dir(pattern="CCMSLIB00000574963.ms")
	times += ["Time for large file: " + str(time.clock() - start)]
	
	start = time.clock()
	records = load_files_from_dir()
	times += ["Time for all files: " + str(time.clock() - start)]
	
	for time in times:
		print(time)
		
	counts = {}
	compoundCounts = {}
	formulaCounts = {}
	for record in records:
		name, dict = record
		for key in dict.keys():
			if(key in counts):
				counts[key] += 1
			else:
				counts[key] = 1
	
		if(dict["compound"] in formulaCounts):
			compoundCounts[dict["formula"]] += 1
		else:
			compoundCounts[dict["formula"]] = 1
	
		if(dict["formula"] in formulaCounts):
			formulaCounts[dict["formula"]] += 1
		else:
			formulaCounts[dict["formula"]] = 1
		
	print("Number of files: " + str(len(list(glob.glob("*.ms")))))

	print()
	print("field counts")	
	for key in counts:
		print(str(key) + " : " + str(counts[key]))
	
	print()
	print("formulaCounts")
	for key in formulaCounts:
		if(formulaCounts[key] > 1):
			print(str(key) + " : " + str(formulaCounts[key]))
	
if __name__ == "__main__":

    main()
			