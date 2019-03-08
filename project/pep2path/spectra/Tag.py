from collections import Counter

'''Class to contain a sequence tag.'''
class Tag():

	def __init__(self, tag, peaks, masses):
		#String containing components of the tag, separated by dashes
		self.tag = tag
		#List of indices of readings across which the tag is found
		self.peaks = peaks
		#List of masses of these readings
		self.masses = masses
		#Int containing the length of the tag
		self.length = len(peaks) - 1
		
	def __str__(self):
		return self.tag + " " + str(self.peaks) + " " + str(self.masses)
		
	'''Breaks up tag, returning a list of its components.'''
	def decompose_tag(self):
		return [split_item for split_item in self.tag.split('-') if split_item != ""]
		
	'''Returns a dictionary of counts of all unique components in this tag.'''
	def component_counts(self):
		return Counter(self.decompose_tag())
	
	'''Returns a list of all unique components in this tag.'''
	def unique_components(self):
		return list(self.component_counts().keys())
		
	'''A function that expands any instance of multiple tag possibility into several sequence tags with defininite possibilities.
		So an example tag -A-[A, B, C]-B- will expand into tags -A-A-B-, -A-B-B-, -A-C-B-.
		This will result in an exponentially large number of tags, so I don't recommend doing it on large tags with lots of uncertainty,
		as the computation will never finish.
		
		Takes own tag and outputs a list of new tags each with one of the expanded tags.'''
	def expand_tag_notation(self):
	
		'''Helper to expand tags.'''
		def recursive_expansion(new_tag, split_tag, index, new_tags):
			if(len(new_tag) >= len(split_tag)):
				output_list.append(Tag("-".join(new_tag), self.peaks, self.masses))
				return
		
			for element in split_tag[index]:
				recursive_expansion(new_tag.append(element), split_tag, index+1, new_tags)
			
			
		new_tags = []
		split_tag = self.tag.split('-')
		spilt_tag = [(list(char) if char[0] == '[' else [char]) for char in split_tag] #convert tags to lists
		recursive_expansion([], split_tag, 0, new_tags)
			
		return new_tags