from collections import defaultdict

class CDSPrediction():
	
	def __init__(self, overall_prediction, predictions):
		#String in the form of a tag (components separated by | or -, surrounded on either side by '-')
		self.overall_prediction = overall_prediction
		#List of dictionaries with keys as name of predictor, and values the value of the prediction, one dictionary per domain in a coding region
		self.predictions = predictions
		
	def __str__(self):
		return "---" + str(self.overall_prediction) + "---\n" + str(self.predictions) + "\n\n"
		
	'''Breaks up predicted tag, returning a list of its components.'''
	def decompose_tag(self):
		return [split_item for split_item in self.overall_prediction.split('-') if split_item != ""]
		
	'''Returns a dictionary of counts of all unique components possibly in this predicted tag.'''
	def component_counts(self):
		dict = defaultdict(lambda: 0)
		for component in self.decompose_tag():
			for p_component in component.split('|'):
				dict[p_component] += 1
		return dict
	
	'''Gets all the unique components predicted to possibly be in the tag.'''
	def unique_components(self):
		return list(self.component_counts().keys())
		
	'''Simple comparison where we just check overlap between constituent parts of predicted tag and given tag.'''
	def compare_to_tag(self, tag):
		#self.overall_prediction.
		pass
		
	'''Compares to all tags in a SpectrumTags object returned by spectra tag finding methods.'''
	def compare_to_spectrum(self, tag_dict):
		pass