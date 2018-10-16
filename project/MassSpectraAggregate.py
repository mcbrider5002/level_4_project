'''A wrapper class to make multiple operations on MassSpectra easier and with the use of less boilerplate.
	It is NOT a guarantee of more efficiency, it's just a convenience.'''
class MassSpectraAggregate():
	def __init__(spectra):
		this.spectra = spectra