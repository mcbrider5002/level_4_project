import numpy as np
import copy
from decimal import Decimal
from collections import defaultdict

from Tag import Tag
from SpectrumTags import SpectrumTags

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
	def __init__(self, id="", compound="", formula="", parent_mass=0.0, ionization="", inchi="", inchi_key="", smiles="", ms2peaks=np.array([]), mass_table={}, misc={}):
	
		#String containing some identifier, probably file name
		self.id = id
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
		
		#An empty data-structure that is populated as we change the order of readings in mass-peaks
		#It allows the new indices to be mapped back to their original positions in the data
		self.mappings = []
		
	'''Returns string representation of object.
		Not very elegant and mostly for debugging.
		I wouldn't use this for anything particularly heavy.'''
	def __str__(self):
		base_string = (str(self.id)
						+"\n>compound " + str(self.compound)
						+ "\n>formula " + str(self.formula)
						+ "\n>parentmass " + str(self.parent_mass)
						+ "\n>ionization " + str(self.ionization)
						+ "\n>InChI " + str(self.inchi)
						+ "\n>InChIKey " + str(self.inchi_key)
						+ "\n>smiles " + str(self.smiles)
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
		
	'''When indices in ms2peaks are changed, add a mapping to translate them back to the indices in the original layout.
		Takes what the changed indices would be in the original enumeration and adds a mapping based on this.'''
	def update_mappings(self, old_indices):
		self.mappings = [{new_index:old_index for new_index, old_index in enumerate(old_indices)}] + self.mappings
	
	'''Returns the current transformed indices to the original indices, using the mappings in mapping.'''
	def original_indices(self, indices):
		for mapping in self.mappings:
			indices = [mapping[index] for index in indices]
		return indices
			
	'''Returns the maximum mass reading in our mass peak readings.'''
	def max_mass(self):
		return np.max(self.ms2peaks[:, self.MASS])
	
	'''Returns the maximum intensity reading in our mass peak readings.'''
	def max_intensity(self):
		return np.max(self.ms2peaks[:, self.INTENSITY])
	
	'''Sorts readings in non-decreasing order of mass.'''
	def sort_by_mass(self):
	
		sorted_indices = np.argsort(self.ms2peaks[:, self.MASS], axis=0)
		self.update_mappings(sorted_indices)
		self.ms2peaks = self.ms2peaks[sorted_indices, :]
	
	'''Filters mass peak readings so as to remove any intensity readings less than the number given.
		By default the minimum threshold is 5% of the max intensity reading, but the function takes a fixed input, not a percentage, so be careful.'''
	def filter_intensity(self, intensity_threshold=None):
	
		intensity_threshold = 0.05*self.max_intensity() if intensity_threshold is None else intensity_threshold #default value
		above_threshold, = (self.ms2peaks[:, self.INTENSITY] > intensity_threshold).nonzero()
		self.update_mappings(above_threshold)
		self.ms2peaks = self.ms2peaks[above_threshold, :]
	
	'''Normalises mass intensity readings so the max reading is equal to some value.
		Default value is 100. Casts to Decimal to try to avoid machine precision numerical errors, then back to float, so this may slow it down a little..'''
	def normalise_intensity(self, new_max=100, old_max=None):
	
		old_max = self.max_intensity() if old_max is None else old_max
		new_intensities = np.array([Decimal(peak) for peak in list(self.ms2peaks[:, self.INTENSITY])]) #cast to decimal
		new_intensities *= ((Decimal(new_max) / Decimal(old_max))) #normalise
		self.ms2peaks[:, self.INTENSITY] = np.array(new_intensities, dtype=float) #cast back, return
		
	'''Calculates an upper and lower bound for a mass tolerance given a mass and a static value to vary by.'''
	@staticmethod
	def static_mass_tolerance(mass, mass_tolerance):
		return (mass - mass_tolerance), (mass + mass_tolerance)
		
	'''Calculates an upper and lower bound for a mass tolerance given a mass and a percentage value to vary by, with percentage relative to each mass reading.'''
	@staticmethod
	def rel_ppm_mass_tolerance(mass, mass_tolerance):
		return (mass * (1 - mass_tolerance)), (mass * (1 + mass_tolerance))
		
	'''Calculates an uppper and lower bound for a mass tolerance given a mass and a percentage value to vary by, with percentage relative to the largest mass reading.'''
	def max_ppm_mass_tolerance(self, mass, mass_tolerance):
		return mass - (self.max_mass() * mass_tolerance), mass + (self.max_mass() * mass_tolerance)
		
	#######################################################
	#Pass one of these names to anything that accepts a mass tolerance mode with MassSpectrum.NAME
	#(Lambdas allow us to fix length of parameter list for use as a variable while keeping static and percentile mass tolerance static)
	STATIC_MASS_TOLERANCE = lambda self, m, mt: self.static_mass_tolerance(m, mt)
	REL_PPM_MASS_TOLERANCE = lambda self, m, mt: self.rel_ppm_mass_tolerance(m, mt)
	MAX_PPM_MASS_TOLERANCE = max_ppm_mass_tolerance
	#######################################################
	
	########################	
	###Tag Search Methods###
	########################
		
	'''A helper method to unpack our provisional data structure when we perform breadth-first search and store all the connections between peaks in it.'''
	def __unpack_tag_paths(self, tag_paths):
	
		to_check = list(range(len(tag_paths))) #List of all tag-beginning points we are yet to check (we ignore any that are visited during the construction of another tag)
		tags = defaultdict(list) #List of pairs of tags and their corresponding peaks extracted from our data structure
		
		'''Internal function to help us recursively traverse the data structure.
			The for loop below starts it off with a 'source node' (i.e. a peak with no previous peaks that have a tag to it),
			then it calls itself for connections found, building a 'stack' of the current sequence and its corresponding peaks,
			and if it can't find any it adds the current state of this 'stack' to the list of tags.'''
		def recursive_structure_search(current_tag, peaks, tags):
			pushed = True #true if last operation was a push (not pop) preventing subsequences of tags from being registered as a separate tag
			
			for index in tag_paths[peaks[-1]].keys():
				compounds = tag_paths[peaks[-1]][index]
				
				pushed = False #If there are ANY additions to the current tag, the current tag is a subsequence and we ignore it
					
				if(index in to_check): #Any tag we visit as part of another tag's path can be ignored as a root/source for a path as it must be a subsequence
					to_check.remove(index)
					
				if(len(compounds) == 1): #if the list of compounds only has one element, just add that element rather than string representation of the list
					compounds = compounds[0]
				recursive_structure_search(current_tag + str(compounds) + "-", (copy.copy(peaks)) + [index], tags) #visit the next tag
			
			if(pushed and len(peaks) > 1): #This is a maximum-length tag (i.e. not a subsequence) and isn't zero-length
				masses = list(self.ms2peaks[peaks, MassSpectrum.MASS])
				(tags[len(peaks) - 1]).append(Tag(current_tag, self.original_indices(peaks), masses)) #Save current tag, and the sequence of corresponding peaks
	
	
		index = 0
		while(index < len(to_check)): 
		
			recursive_structure_search("-", [to_check[index]], tags)
			index += 1 
			
		longest_tag = (max(tags.keys()) if len(tags.keys()) > 0 else 0)
			
		return SpectrumTags(self.id, longest_tag, tags)
	
	'''ms2peaks contains an array of mass spectrometry readings, with a mass value and an intensity value.
		We want to find gaps in the masses that correspond to the masses in mass_table (with some tolerance, either percentile or absolute).
		Then we want to find all sequences of masses where one mass begins on a peak as another ends, with no subsequences - this is a sequence tag. 
		This function returns a dictionary of pairs of sequence tags and their corresponding peaks, where each key is the length of the sequence tag
		(which is equal to the number of peaks minus one).
		
		This problem can be modelled as a DAG (directed acyclic graph), where each vertex is enumerated and vertices can only have edges to vertices with higher numbers,
		edges have a label (and edges are labelled by a list of compounds they could possibly represent), and these edges are unknown,
		but discoverable. In this case, vertices are individual mass spectrometry readings, and edges are potential compounds from gaps between mass peaks.
		This is the reason for any reference to computational graph terminology.
	
		mass_tolerance_mode allows you to specify one of this function's methods to use to calculate the upper and lower bounds for the mass tolerance
		for a given compound. (Standard names for passing these are defined with the functions, for readability's sake.) These methods support fixed mass tolerance, 
		percentile and percentile relative to the largest value in the file. 
		mass_threshold specifies the value to be used for the mass tolerance. So in the case of a fixed mass tolerance of 0.05 it will give +-0.05 to the value
		and in the case of percentile it will give +-5%.
		
		Returns a dictionary of lists of sequence tags (with no subsequences because those can be reconstructed from longer tags), stored under their length.'''
	def find_sequence_tags(self, mass_tolerance_mode=None, mass_threshold=0.00001):
	
		mass_tolerance_mode = MassSpectrum.MAX_PPM_MASS_TOLERANCE if mass_tolerance_mode is None else mass_tolerance_mode #default value
	
		tag_paths = [defaultdict(list) for item in self.ms2peaks]  #our provisional data structure where each entry is a dictionary representing a peak,
																	#with keys indices of other (later) peaks and values a list of possible compounds which could be represented
																	#by the mass differences between these peaks
		
		for index in range(len(tag_paths)):
		
			#calculate a table of mass differences between each peak and all following peaks
			mass_differences = np.array(self.ms2peaks)[(index+1):, self.MASS] - np.array(self.ms2peaks)[index, self.MASS] 
		
			#for each compound check if any mass differences are approximately equal to the mass of a compound in the mass table and add it to tag_paths if so
			for compound, mass in self.mass_table.items():
			
				lower_threshold, upper_threshold = mass_tolerance_mode(self, mass, mass_threshold)
					
				candidate_positions, = (np.logical_and(mass_differences <= upper_threshold, mass_differences >= lower_threshold)).nonzero() #get indices where element between threshold 
				candidate_positions += index + 1 #adjust for the offset in only checking part of the array (since peaks are directional we only check succeeding ones)
				
				for next_pos in candidate_positions: tag_paths[index][next_pos] += [compound] #add to provisional data structure
					
		return self.__unpack_tag_paths(tag_paths)