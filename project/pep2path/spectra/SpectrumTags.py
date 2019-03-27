from collections import Counter
import itertools

'''Class to contain all sequence tags from a spectrum.'''
class SpectrumTags():

	def __init__(self, id, longest_tag, tags):
		#String containing id of the file this was generated from
		self.id = id
		#Int containing the length of the longest tag in the dictionary
		self.longest_tag = longest_tag
		#Dictionary containing tag lengths as keys and lists of tags as values i.e. {length:[tag]}
		self.tags = tags
		
	def __str__(self):
		return self.id + "\n" + "\n".join([("---Length: %d---" % length) + "\n" + "\n".join([str(tag) for tag in self.tags[length]]) for length in self.tags.keys()])
	
	'''Returns tags but where each is a list of components rather than a string.'''
	def decompose_tags(self):
		return {length:[tag.decompose_tag() for tag in tags] for length, tags in self.tags.items()}
		
	'''Returns counts of all unique components in all the tags in the spectrum combined.'''
	def multi_component_counts(self):
		return Counter(comp for tag in self.flatten_tags() for comp in tag.decompose_tag())
		
	'''Returns a list of all unique components in all the tags in the spectrum combined.'''
	def unique_components(self):
		return list(self.multi_component_counts().keys())
		
	'''Return all tags of the given length.'''
	def filter_by_length(self, length):
		return self.tags[length] if (length in self.tags) else []
		
	'''Expands all tag notation for all tags in dict, according to the specification in Tag, returning a new tag dict.'''
	def multi_expand_tag_notation(self):
		return {length:(itertools.chain.from_iterable([tag.expand_tag_notation() for tag in tag_list])) for length, tag_list in self.tags.items()}
		
	'''Convenience method to flatten the dictionary form of the tags (where they are stored by length) into a single list (in no guaranteed order).'''
	def flatten_tags(self): 
		return [tag for length, tags in self.tags.items() for tag in tags]