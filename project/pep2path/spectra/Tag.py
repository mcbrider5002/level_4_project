import itertools
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
		
	'''Expands the multiple-possibility brackets (e.g. [A, B, C]) in tag strings that are used in place of storing every individual tag 
		(as this would cause them to be exponential in number).'''
	@staticmethod
	def expand_brackets(comp):
		return ["".join([chr for chr in comp if chr.isalnum() or chr == '/']) for comp in comp.split(',')]
		
	'''A function that expands any instance of multiple tag possibility into several sequence tags with definite possibilities.
		So an example tag -A-[A, B, C]-B- will expand into tags -A-A-B-, -A-B-B-, -A-C-B-.
		This will result in an exponentially large number of tags, so I don't recommend doing it on large tags with lots of uncertainty,
		as it will produce quite a lot of values.
		
		Takes own tag and outputs an iterator of new tags each with one of the expanded tags.'''
	def expand_tag_notation(self):
		split_tag = [self.expand_brackets(comp) for comp in self.tag.split('-') if comp != ""]
		return map(lambda t: "-".join(t), itertools.product(*split_tag))