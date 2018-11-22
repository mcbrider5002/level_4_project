import itertools

class GenbankFile():
	
	def __init__(self, predictions):
		#List of CDSPrediction objects contained in the file
		self.predictions = predictions
		
	def __str__(self):
		string = ""
		for prediction in self.predictions: string += str(prediction)
		return string

	'''Returns all unique components in the Genbank file.'''
	def unique_components(self):
		return list(set(itertools.chain.from_iterable([prediction.unique_components() for prediction in self.predictions])))