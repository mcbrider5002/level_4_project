class GenbankFile():
	
	def __init__(self, predictions):
		self.predictions = predictions
		
	def __str__(self):
		string = ""
		for prediction in self.predictions: string += str(prediction)
		return string
