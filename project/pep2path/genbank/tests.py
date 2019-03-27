import random
import itertools
from collections import Counter

from .aanames import AA_alphabet

from .CDSPrediction import CDSPrediction
from .GenbankFile import GenbankFile

'''Helper to generate a list of lists of random components.'''
def generate_comp_lists(alphabet):
	return [[str(random.choice(list(alphabet))) for length in range(random.randint(20, 40))] for number in range(random.randint(10, 20))]
	
'''Helper to generate a list of random gbks.'''	
def generate_gbk(alphabet):
	generated = generate_comp_lists(alphabet)
	strings = ["-".join(g) for g in generated]
	preds = [CDSPrediction(s, [{}]) for s in strings]
	return generated, strings, preds

def test_decompose_tag(alphabet=AA_alphabet):
	
	generated, strings, preds = generate_gbk(alphabet)
	
	string = strings[0]
	pred = preds[0]
	assert pred.decompose_tag() == generated[0], "decompose_tag doesn't match pregenerated prediction!"
	assert CDSPrediction("-%s-" % string, [{}]).decompose_tag() == generated[0], "decompose_tag (on string with trailing '-') doesn't match pregenerated prediction!"
	
	file = GenbankFile(preds)
	file2 = GenbankFile([CDSPrediction("-%s-" % s, [{}]) for s in strings])
	assert (file.decompose_tags() == generated), "decompose_tags doesn't match pregenerated predictions!"
	assert (file2.decompose_tags() == generated), "decompose_tags (on string with trailing '-') doesn't match pregenerated predictions!"
	
def test_component_counts(alphabet=AA_alphabet):
	
	generated, _, preds = generate_gbk(alphabet)
	counts = [Counter(comps) for comps in generated]
	
	c = Counter(generated[0])
	pred = preds[0]
	assert pred.component_counts() == c, "component_counts doesn't match component counts of pregenerated prediction!"
	
	file = GenbankFile(preds)
	assert file.multi_component_counts() == counts, "multi_component_counts doesn't match component counts of pregenerated predictions!"
	
def test_unique_components(alphabet=AA_alphabet):

	generated, _, preds = generate_gbk(alphabet)
	uniques = list(map(set, generated))
	
	unique = uniques[0]
	pred = preds[0]
	assert set(pred.unique_components()) == unique, "CDSPrediction's unique_components doesn't match components of pregenerated prediction!"
	
	file = GenbankFile(preds)
	assert set(file.unique_components()) == set(itertools.chain.from_iterable(uniques)), "GenbankFile's unique_components doesn't match components of pregenerated predictions!"

def tests():
	test_decompose_tag()
	test_component_counts()
	test_unique_components()