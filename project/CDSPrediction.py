class CDSPrediction():
	
	def __init__(self, overall_prediction, predictions):
		#String
		self.overall_prediction = overall_prediction
		#List of dictionaries with keys as name of predictor, and values the value of the prediction, one dictionary per domain in a coding region
		self.predictions = predictions
		
	def __str__(self):
		return "---" + str(self.overall_prediction) + "---\n" + str(self.predictions) + "\n\n"