import random

from .masstables import AA_alphabet
from .comparisons import *
from .experiment import *
from .ripp2path import one_letter_AA, compare_ripp, ripp2path_scorer

'''Helper to generate a list of lists of random components.'''
def generate_comp_lists(alphabet):
	return [[str(random.choice(list(alphabet))) for length in range(random.randint(20, 40))] for number in range(random.randint(10, 20))]
	
'''Helper to generate a list of lists of lists of random components.'''	
def generate_multi_comps(alphabet):
	return [generate_comp_lists(alphabet) for number in range(random.randint(5, 50))]
	
def test_shuffle_components(alphabet=AA_alphabet):
	comps = generate_multi_comps(alphabet)
	shuffled = shuffle_components(comps)
	def srt(ls): sorted(map(sorted, map(lambda l : map(sorted, l), ls)))
	assert comps != shuffled and srt(comps) == srt(shuffled), "shuffle_components fails test!"
	
def test_randomise_components(alphabet=AA_alphabet):
	comps = generate_multi_comps(alphabet)
	randomised = randomise_components(comps, alphabet)
	assert comps != randomised and all([comp in alphabet for file in randomised for cds in file for comp in cds]), "randomise_components fails test!"
	
def test_compare_unique_components(alphabet=AA_alphabet):
	pass
	
def test_simple_score(alphabet=AA_alphabet):

	comps = list(itertools.chain.from_iterable(generate_comp_lists(alphabet)))
	length = len(comps)
	
	no_dummies = random.randint(0, length)
	dummies = random.sample(list(range(length)), no_dummies)
	
	r_padder = lambda r : ["dummy" for i in range(random.randint(0, r))]
	dummy_comps = r_padder(5) + [("dummy" if i in dummies else comp) for i, comp in enumerate(comps)] + r_padder(5)
	
	assert simple_score(comps, dummy_comps) == (length - no_dummies), "simple_score does not give the expected value!"
	
def test_exact_match_I(alphabet=AA_alphabet):
	aa = random.choice(list(AA_alphabet))
	M_lookup = {}
	assert exact_match_I(aa, aa, {}) == 1, "exact_match_I doesn't evaluate to 1 on a match!"
	assert exact_match_I(aa, None, {}) == 0, "exact_match_I doesn't evaluate to 0 on a mismatch!"
	
def test_SVM_I(alphabet=AA_alphabet):
	pass
	
def test_six_class_I(alphabet=AA_alphabet):
	pass
	
def test_alignment_score(alphabet=AA_alphabet):
	pass
	
def test_compare_alignment(alphabet=AA_alphabet):
	pass
	
def nrp_tests():
	test_shuffle_components()
	test_randomise_components()
	test_compare_unique_components()
	test_simple_score()
	test_exact_match_I()
	test_SVM_I()
	test_six_class_I()
	test_alignment_score()
	test_compare_alignment()
	
def	test_compare_ripp(alphabet=AA_alphabet):

	comps = list(itertools.chain.from_iterable(generate_comp_lists(alphabet)))
	length = len(comps)
	
	no_dummies = random.randint(0, length)
	dummies = random.sample(list(range(length)), no_dummies)
	dummy_comps = [("dummy" if i in dummies else comp) for i, comp in enumerate(comps)]
	
	assert compare_ripp(comps, dummy_comps) == (length - no_dummies) / length, "compare_ripp does not give the expected value!"
	
def test_ripp2path_scorer(alphabet=AA_alphabet):
	pass
	
def ripp_tests():
	test_compare_ripp()
	test_ripp2path_scorer()