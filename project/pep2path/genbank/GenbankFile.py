import itertools

class GenbankFile():
	
	def __init__(self, predictions):
		#List of CDSPrediction objects contained in the file
		self.predictions = predictions
		
	def __str__(self):
		string = ""
		for prediction in self.predictions: string += str(prediction)
		return string

	'''Returns a list of lists, where each internal list represents a CDS in this file and within are the components of its tag.'''
	def decompose_tags(self):
		return [prediction.decompose_tag() for prediction in self.predictions]
		
	'''Returns a list of dictionaries of counts of all unique components possibly in each predicted tag.'''
	def multi_component_counts(self):
		return [prediction.component_counts() for prediction in self.predictions]
		
	'''Returns all unique components in the Genbank file.'''
	def unique_components(self):
		return list(set(itertools.chain.from_iterable([prediction.unique_components() for prediction in self.predictions])))