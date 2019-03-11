from sys import argv

from pep2path.spectra.spectraMain import main as spectra_main
from pep2path.spectra.drawSpectra import main as graph_main

from pep2path.experiment import intersection_experiment, simple_experiment, p2p_experiment
from pep2path.ripp2path import test_ripp2path

from pep2path.spectra.tests import tests as spectra_tests
from pep2path.genbank.tests import tests as genbank_tests
from pep2path.tests import nrp_tests
from pep2path.tests import ripp_tests

def spectra(args): spectra_main()
def graph(args): graph_main()

def tests(args):
	spectra_tests()
	genbank_tests()
	nrp_tests()
	ripp_tests()
	
def experiment(args):
	exp_modes = {"intersection" : intersection_experiment,
				 "simple" : simple_experiment,
				 "nrp2path" : p2p_experiment}
				 
	if(len(args) > 0):
		exp_modes.get(args[-1].lower(), lambda : print("Not a recognised type of experiment..."))()
	else:
		intersection_experiment()
	
def comparisons(args):
	pass
	
def testripp(args):
	test_ripp2path()
	
modes = {"spectra" : spectra,
		 "graph" : graph,
		 "tests" : tests,
		 "experiment" : experiment,
		 "comparisons" : comparisons,
		 "testripp" : testripp}
			 
def flag_parser(args):
	return args

def main():
			 
	invalid_mode_msg = "Not a valid mode, valid modes include: " + ", ".join(modes.keys())
	
	def output_func(args): print("Default functionality")
	args = flag_parser(argv)
	if(len(args) > 1):
		def invalid_mode(args): print(invalid_mode_msg)
		output_func = modes.get(args[1].lower(), invalid_mode)
		
	output_func(args[2:])
	
if __name__ == "__main__":

    main()