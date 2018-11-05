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
		
	'''Expands all tag notation for all tags in dict, according to the specification in Tag, returning a new tag dict.'''
	def expand_tag_notation(self):
		return {length:(itertools.chain.from_iterable([tag.expand_tag_notation() for tag in tag_list])) for length, tag_list in self.tags.items()}
		
	'''Convenience method to flatten the dictionary form of the tags (where they are stored by length) into a single list (in no guaranteed order).'''
	def flatten_tags(self):
		return list(itertools.chain.from_iterable([tags for length, tags in self.tags.items()]))
		
	'''Return all tags of the given length.'''
	def filter_by_length(self, length):
		return self.tag_dict[length] if (length in self.tags) else []
		