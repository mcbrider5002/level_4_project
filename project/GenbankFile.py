class GenbankFile():
	
	def __init__(self, predictions):
		#List of CDSPrediction objects contained in the file
		self.predictions = predictions
		
	def __str__(self):
		string = ""
		for prediction in self.predictions: string += str(prediction)
		return string
