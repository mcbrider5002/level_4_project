class CDSPrediction():
	
	def __init__(self, overall_prediction, predictions):
		#String
		self.overall_prediction = overall_prediction
		#List of dictionaries with keys as name of predictor, and values the value of the prediction, one dictionary per domain in a coding region
		self.predictions = predictions
		
	def __str__(self):
		return "---" + str(self.overall_prediction) + "---\n" + str(self.predictions) + "\n\n"
		
	'''Simple comparison where we just check overlap between constituent parts of predicted tag and given tag.'''
	def compare_to_tag(self, tag):
		self.overall_prediction.
		
	'''Gets all the unique components predicted to possibly be in the tag.'''
	def unique_components(self):
		
	'''Compares to all tags in a SpectrumTags object returned by spectra tag finding methods.'''
	def compare_to_spectrum(self, tag_dict):
		