import numpy as np
import copy

'''A wrapper class for mass-spectrometry files.'''
class MassSpectrum():
	
	#'Constants' to make selection of columns in ms2peaks more human-readable
	MASS = 0
	INTENSITY = 1
	
	#############
	###Methods###
	#############
	
	'''Initialiser
		All parameters are optional to allow leaving some blank, so be sure you fill everything you need'''
	def __init__(self, compound="", formula="", parent_mass=0.0, ionization="", inchi="", inchi_key="", smiles="", ms2peaks=np.array([]), mass_table={}, misc={}):
	
		#Although Python is dynamically-typed, I'll suggest data types here - I can't guarantee methods will work if performed on the wrong data type
		
		#String containing compound name
		self.compound = compound
		#String
		self.formula = formula
		#Double
		self.parent_mass = parent_mass
		#String
		self.ionization = ionization
		#String
		self.inchi = inchi
		#String
		self.inchi_key = inchi_key
		#String
		self.smiles = smiles
		#NumPy array containing masses of readings in first column, and intensity in second
		self.ms2peaks = ms2peaks
		#Linked mass residue table as dictionary, with compound names as keys and masses as (double) values
		self.mass_table = mass_table
		#A miscellaneous dictionary to hold any other data you need with respect to a particular set of mass-spectrometry data
		#Just put the names of any additional data types (e.g. "Annotations"...) as keys and list them under values, for easy extraction of all additional fields with .keys()
		#If you intend to create additional data to be held across many instances of this class with standard operations performed on it, 
		#I recommend you extend the class itself instead
		self.misc = misc
		
	'''Returns string representation of object.
		Not very elegant and mostly for debugging.
		I wouldn't use this for anything particularly heavy.'''
	def __str__(self):
		base_string = (">compound " + self.compound
						+ "\n>formula " + self.formula
						+ "\n>parentmass " + str(self.parent_mass)
						+ "\n>ionization " + self.ionization
						+ "\n>InChI " + self.inchi
						+ "\n>InChIKey " + self.inchi_key
						+ "\n>smiles " + self.smiles
						+ "\n>ms2peaks ")
						
		for index in range(len(self.ms2peaks[:, 0])):
			base_string += "\n" + str(self.ms2peaks[index, 0]) + " " + str(self.ms2peaks[index, 1])
			
		base_string += "\n>masstable" 
		for compound in self.mass_table.keys():
			base_string += "\n" + str(compound) + " " + str(self.mass_table[compound])
			
		#Don't know in advance what format misc data will be in so will just append it naively
		for key in self.misc.keys():
			base_string += "\n>" + str(key) + "\n" + str(self.misc[key])
			
		return base_string
	
	'''Returns the maximum intensity reading in our mass peak readings.'''
	def max_intensity(self):
		return np.max(self.ms2peaks[:, self.INTENSITY])
	
	'''Filters mass peak readings so as to remove any low intensity readings.
		By default the minimum threshold is 5% of the max intensity reading.'''
	def filter_intensity(self, intensity_threshold=None):
		intensity_threshold = 0.05*self.max_intensity() if intensity_threshold is None else intensity_threshold #default value
		self.ms2peaks = self.ms2peaks[self.ms2peaks[:, self.INTENSITY] > intensity_threshold, :]
	
	'''Normalises mass intensity readings so the max reading is equal to some value.
		Default value is 100. Doesn't try to avoid machine precision numerical errors.'''
	def normalise_intensity(self, intensity_threshold=100):
		self.ms2peaks[:, self.INTENSITY] = self.ms2peaks[:, self.INTENSITY] * (intensity_threshold / self.max_intensity())
		
	'''Calculates an upper and lower bound for a mass tolerance given a mass and a static value to vary by.'''
	@staticmethod
	def static_mass_tolerance(mass, mass_tolerance):
		return (mass - mass_tolerance), (mass + mass_tolerance)
		
	'''Calculates an upper and lower bound for a mass tolerance given a mass and a percentage value to vary by.'''
	@staticmethod
	def percentile_mass_tolerance(mass, mass_tolerance):
		return (mass * (1 - mass_tolerance)), (mass * (1 + mass_tolerance))
		
	#######################################################
	#Pass one of these names to anything that accepts a mass tolerance mode with instance.NAME
	#(This is just a convenience to make this class easier to read - just the function name without the brackets will do, if you want)
	STATIC_MASS_TOLERANCE = static_mass_tolerance
	PERCENTILE_MASS_TOLERANCE = percentile_mass_tolerance
	#######################################################
	
	########################	
	###Tag Search Methods###
	########################
		
	'''A helper method to unpack our provisional data structure when we perform breadth-first search and store all the connections between peaks in it.'''
	@staticmethod
	def __unpack_tag_paths(tag_paths):
	
		to_check = list(range(len(tag_paths))) #List of all tag-beginning points we are yet to check (we ignore any that are visited during the construction of another tag)
		tags = [] #List of pairs of tags and their corresponding peaks extracted from our data structure
		
		'''Internal function to help us recursively traverse the data structure.
			The for loop below starts it off with a 'source node' (i.e. a peak with no previous peaks that have a tag to it),
			then it calls itself for connections found, building a 'stack' of the current sequence and its corresponding peaks,
			and if it can't find any it adds the current state of this 'stack' to the list of tags.'''
		def recursive_structure_search(current_tag, peaks, tags):
			pushed = True #true if last operation was a push (not pop) preventing subsequences of tags from being registered as a separate tag
			
			for compound in tag_paths[peaks[-1]].keys():
				for index in tag_paths[peaks[-1]][compound]:
				
					pushed = False #If there are ANY additions to the current tag, the current tag is a subsequence and we ignore it
					
					if(index in to_check): #Any tag we visit as part of another tag's path can be ignored as a root/source for a path as it must be a subsequence
						to_check.remove(index)
						
					recursive_structure_search(current_tag + compound, copy.copy(peaks) + [index], tags) #visit the next tag
			
			if(pushed and len(peaks) > 1): #This is a maximum-length tag (i.e. not a subsequence) and isn't zero-length
				tags += [(current_tag, peaks)] #Save current tag, and the sequence of corresponding peaks
		
		index = 0
		while(index < len(to_check)): 
		
			recursive_structure_search("", [to_check[index]], tags)
			index += 1 
			
		return tags
		
	'''def __find_sequence_tags_naive_depth_first(self, mass_tolerance_mode=None, mass_threshold=0.05):
	
		mass_tolerance_mode = self.STATIC_MASS_TOLERANCE if mass_tolerance_mode is None else mass_tolerance_mode #default value
	
		tags = []
		
		stack = []
		mass_differences = []
		
		def recursive_tag_search():
			for compound in mass_table.keys():
				
				mass = mass_table[compound]
				
				lower_threshold, upper_threshold = mass_tolerance_mode(mass, mass_threshold)
				
				for next_peak in range(index, len(ms2peaks[:][self.MASS])):
					if(ms2peaks[next_peak, self.MASS] > ms2peaks[:
			
		
		for index in range(len(ms2peaks[:, self.MASS])):
	
	def __find_sequence_tags_naive_breadth_first():
	
	def __find_sequence_tags_depth_first():	'''	
	
	'''Sequence tag search method using NumPy arrays and breadth-first search, storing the results of discovering potential compounds between gaps
		in a provisional data structure and then constructing sequence tags by walking across it.'''
	def __find_sequence_tags_breadth_first(self, mass_tolerance_mode=None, mass_threshold=0.05):
	
		mass_tolerance_mode = self.STATIC_MASS_TOLERANCE if mass_tolerance_mode is None else mass_tolerance_mode #default value
	
		tag_paths = [{} for item in self.ms2peaks] #our provisional data structure to hold potential next nodes in sequence tags under the candidate compound name
		
		for index in range(len(tag_paths)):
		
			#calculate a table of mass differences between each peak and all following peaks
			mass_differences = np.array(self.ms2peaks)[(index+1):, self.MASS] - np.array(self.ms2peaks)[index, self.MASS] 
		
			#for each compound check if any mass differences are approximately equal to the mass of a compound in the mass table and add it to tag_paths if so
			for compound in self.mass_table.keys():
			
				mass = self.mass_table[compound]
				lower_threshold, upper_threshold = mass_tolerance_mode(mass, mass_threshold)
				candidate_positions, = (np.logical_and(mass_differences <= upper_threshold, mass_differences >= lower_threshold)).nonzero()
				candidate_positions += index + 1
				
				if(len(candidate_positions) > 0):
					tag_paths[index][compound] = candidate_positions
					
		return self.__unpack_tag_paths(tag_paths)
		
	##############################################################	
	#Pass one of these names to anything that accepts a  with instance.NAME
	#(This is just a convenience to make this class easier to read - just the function name without the brackets will do, if you want)
	BREADTH_FIRST = __find_sequence_tags_breadth_first
	'''DEPTH_FIRST = __find_sequence_tags_depth_first
	BREADTH_FIRST_NAIVE = __find_sequence_tags_naive_breadth_first
	DEPTH_FIRST_NAIVE = __find_sequence_tags_naive_depth_first'''
	##############################################################
				
	'''ms2peaks contains an array of mass spectrometry readings, with a mass value and an intensity value.
		We want to find gaps in the masses that correspond to the masses in mass_table (with some tolerance, either percentile or absolute).
		Then we want to find all sequences of masses where one mass begins on a peak as another ends, with no subsequences - this is a sequence tag. 
		This function returns a list of pairs of all sequence tags and the corresponding peaks in between whose gaps one finds the constitutents of the sequence tag.
		
		This problem can be modelled as a DAG (directed acyclic graph), where each vertex is enumerated and vertices can only have edges to vertices with higher numbers,
		edges have a label (and a given vertex may have several edges with different labels to another, although this is unlikely in practise), and these edges are unknown,
		but discoverable. In this case, vertices are individual mass spectrometry readings, and edges are potential compounds from gaps between mass peaks.
		This is the reason for any reference to computational graph terminology.
		
		So-called 'breadth-first' search methods in this scheme differ not only by being breadth-first searches, but they also perform the computation in two phases,
		discovering the graph's edges (i.e. any potential candidate compounds between peaks), storing the complete graph in a list of dictionaries,
		then traversing it to extract all sequence tags. So-called 'depth-first' algorithms are indeed 'depth-first', but do the calculation all in one stage,
		potentially needing to recompute mass-differences but avoiding the need to traverse the data structure.
	
		mass_tolerance_mode allows you to specify one of this function's static methods to use to calculate the upper and lower bounds for the mass tolerance
		for a given compound. (Standard names for passing these are defined with the functions, for readability's sake.) These methods support both fixed mass tolerance and 
		percentile. mass_threshold specifies the value to be used for the mass tolerance. So in the case of a fixed mass tolerance of 0.05 it will give +-0.05 to the value
		and in the case of percentile it will give +-5%.
	
		This method is actually just an interface around the variant sequence tag methods to reduce the complexity of this class from the outside.
		Rather than having various different methods to be called, this method wraps around them and chooses a method to call with the search_mode
		argument, which is by default breadth-first non-naive search. (Standard names are defined just above this time.)'''
	def find_sequence_tags(self, mass_tolerance_mode=None, mass_threshold=0.05, search_mode=None):
		mass_tolerance_mode = self.STATIC_MASS_TOLERANCE if mass_tolerance_mode is None else mass_tolerance_mode #default values
		search_mode = self.BREADTH_FIRST if search_mode is None else search_mode
		return search_mode(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
		
	'''Finds the longest tag present using the find_sequence_tags method.
		Just wraps around it and then traverses its output with a for loop.'''
	def find_longest_tag(self, mass_tolerance_mode=None, mass_threshold=0.05, search_mode=None):
	
		mass_tolerance_mode = self.STATIC_MASS_TOLERANCE if mass_tolerance_mode is None else mass_tolerance_mode #default values
		search_mode = self.BREADTH_FIRST if search_mode is None else search_mode
		
		tags = self.find_sequence_tags(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold, search_mode=search_mode)
		
		longest_tag_length = 0 #stores the length of the current longest tag
		longest_tag = "" #stores the readout of compounds in current longest tag
		
		for tag in tags:
			tag_rep, peaks = tag
			if(len(peaks) - 1 > longest_tag_length):
				longest_tag_length = len(peaks) - 1
				longest_tag = tag_rep
				
		return (longest_tag, longest_tag_length)
		