class GenbankPredictions():
	
	def __init__(self, overall_prediction, predictions):
		self.overall_prediction = overall_prediction
		self.predictions = predictions
		
	def __str__(self):
		return "---" + self.overall_prediction + "---\n" + str(predictions) + "\n\n"