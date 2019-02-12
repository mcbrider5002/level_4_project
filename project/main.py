from sys import argv

from pep2path.spectra.spectraMain import main as spectra_main
from pep2path.spectra.drawSpectra import main as graph_main
from pep2path.spectra.tests import tests as spectra_tests

def flag_parser(args):
	return args

def main():

	def spectra(args): spectra_main()
	def graph(args): graph_main()
	def tests(args):
		#to add other tests
		spectra_tests()
	def experiment(args):
		#tba
		pass
	def comparisons(args):
		#tba
		pass
	
	modes = {"spectra" : spectra,
			 "graph" : graph,
			 "tests" : tests,
			 "experiment" : experiment,
			 "comparisons" : comparisons}
			 
	invalid_mode_msg = "Not a valid mode, valid modes include: " + ", ".join(modes.keys())
	
	def output_func(args): print("Default functionality")
	args = flag_parser(argv)
	if(len(args) > 1):
		def invalid_mode(args): print(invalid_mode_msg)
		output_func = modes.get(args[1], invalid_mode)
		
	output_func(args[2:])
	
if __name__ == "__main__":

    main()