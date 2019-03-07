from sys import argv

from pep2path.spectra.spectraMain import main as spectra_main
from pep2path.spectra.drawSpectra import main as graph_main

from pep2path.spectra.tests import tests as spectra_tests
from pep2path.genbank.tests import tests as genbank_tests
from pep2path.tests import tests as comparison_tests

from pep2path.ripp2path import test_ripp2path

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
def testripp(args):
	test_ripp2path()

def flag_parser(args):
	return args

def main():
	
	modes = {"spectra" : spectra,
			 "graph" : graph,
			 "tests" : tests,
			 "experiment" : experiment,
			 "comparisons" : comparisons,
			 "testripp" : testripp}
			 
	invalid_mode_msg = "Not a valid mode, valid modes include: " + ", ".join(modes.keys())
	
	def output_func(args): print("Default functionality")
	args = flag_parser(argv)
	if(len(args) > 1):
		def invalid_mode(args): print(invalid_mode_msg)
		output_func = modes.get(args[1], invalid_mode)
		
	output_func(args[2:])
	
if __name__ == "__main__":

    main()