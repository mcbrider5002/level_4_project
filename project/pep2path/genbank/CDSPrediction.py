from collections import Counter

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
		return Counter([s_comp for comp in self.decompose_tag() for s_comp in comp.split('|')])
	
	'''Gets all the unique components predicted to possibly be in the tag.'''
	def unique_components(self):
		return list(self.component_counts().keys())