import itertools
import copy
import random
import os
import math
from statistics import mean
import numpy as np
from collections import Counter

from .masstables import kersten_names, kersten_alphabet, AA_names, AA_alphabet
	
'''Given a twice-nested list, returns a new list with the internal structure and the elements themselves shuffled.'''
def shuffle_components(gbk_files):

	reorder = lambda ls: random.sample(ls, len(ls))
		
	cds_numbers = reorder([len(cds) for gbk in gbk_files for cds in gbk])
	gbk_numbers = reorder([len(gbk) for gbk in gbk_files])
	
	reimpose_structure = lambda ls, indices: [ls[i:j] for i, j in zip([0] + list(np.cumsum(indices)[:-1]), np.cumsum(indices))]
	flattened = reorder(list(itertools.chain.from_iterable(itertools.chain.from_iterable(gbk_files))))
	return reimpose_structure(reimpose_structure(flattened, cds_numbers), gbk_numbers)
	
'''Given a nested list where internal lists represent genbank files containing CDSes containing their unique components, 
	rearranges all the inner lists (files) and then randomly selects an animo acid to replace each of those components
	before filtering them to be unique again.'''
def randomise_components(gbk_components, alphabet=AA_alphabet):
	gbk_components = random.sample(gbk_components, len(gbk_components)) #change order of files to be random
	return [[[str(random.choice(list(alphabet))) for component in cds] for cds in gbk] for gbk in gbk_components] #randomise all components
	
'''Given a tag in list form, normalise it for comparison. Also used for gbk where inner cds is removed.'''
def normalise_tag(tag):
	return [comp.lower() for comp in tag]
	
'''Given a gbk in list form, normalise it for comparison.'''
def normalise_gbk(gbk):
	return [[comp.lower() for comp in cds] for cds in gbk]

'''Jaccard Similarity Index of components. Apply unique_components method and normalise first. '''
def score_unique_components(spectrum, gbk):
	return len(set(spectrum) & set(gbk)) / len(set(spectrum) | set(gbk))
	
'''Given a nested list where internal lists represent a mass spectrum and their contents are their unique components,
	and a nested list where internal lists represent genbank files and their contents are their unique components,
	and two lists of names where the names correspond to the files in each of these lists
	returns several statistics based on comparisons of the unique components of these nested lists.'''
def compare_unique_components(spectra_names, spectra, gbk_names, gbk_components):

	#get unique components
	spectra_components = [normalise_tag(spectrum) for spectrum in spectra]
	gbk_components = [normalise_tag(list(set(itertools.chain.from_iterable(gbk)))) for gbk in gbk_components]

	best_score = ("No spectrum!", "No gbk!", float("-inf"))
	score_sum, total_comparisons = 0, 0
	
	for gbk_name, gbk in zip(gbk_names, gbk_components):
		for spectrum_name, spectrum in zip(spectra_names, spectra_components):
			score = score_unique_components(spectrum, gbk)
			total_comparisons += 1
			score_sum += score
			best_score = max(best_score, (spectrum_name, gbk_name, score), key=lambda t: t[2])		
					
	avg = score_sum / total_comparisons
					
	return best_score, score_sum, avg
		
svm_alphabet = ['ala', 'gly', 'val', 'leu', 'ile', 'abu', 'iva', 'ser', 'thr', 'hpg', 'dhpg', 'cys', 'pro', 'pip', 'arg', 'asp', 'glu', 'his', 'asn', 'lys', 'gln', 'orn', 'aad', 'phe', 'tyr', 'trp', 'dhb', 'phg', 'bht']
		
three_classes = {'hydrophobic-aliphatic':['ala', 'gly', 'val', 'leu', 'ile', 'abu', 'iva', 'ser', 'thr', 'hpg', 'dhpg', 'cys', 'pro', 'pip'],
                 'hydrophilic':['arg', 'asp', 'glu', 'his', 'asn', 'lys', 'gln', 'orn', 'aad'],
                 'hydrophobic-aromatic':['phe', 'tyr', 'trp', 'dhb', 'phg', 'bht']}

big_classes = ['gly,ala,val,leu,ile,abu,iva',
			   'ser,thr,dhpg,hpg', 
			   'asp,asn,glu,gln,aad', 
			   'dhb,sal', 
			   'cys', 
			   'phe,trp,phg,tyr,bht', 
			   'orn,lys,arg', 
			   'pro,pip']

small_classes = ['val,leu,ile,abu,iva', 
				 'ser', 
				 'gly,ala', 
				 'dhb,sal', 
				 'cys', 
				 'thr', 
				 'asp,asn', 
				 'phe,trp', 
				 'glu,gln', 
				 'pro', 
				 'tyr,bht', 
				 'orn,horn', 
				 'dhpg,hpg', 
				 'arg', 
				 'aad']

'''Creates a doubly-nested dictionary that maps a pair of amino acids to whether their product class matches, based on the classes used in the original Pep2Path.'''
def svm_class_matcher():
	dict = {comp1 : {comp2 : 0.25} for cls, cls_items in three_classes.items() for comp1, comp2 in itertools.permutations(cls_items, 2)}
	dict.update({comp1 : {comp2 : 0.5} for cls_items in big_classes for comp1, comp2 in itertools.permutations(cls_items.split(','), 2)})
	dict.update({comp1 : {comp2 : 0.75} for cls_items in small_classes for comp1, comp2 in itertools.permutations(cls_items.split(','), 2)})
	for aa, inner_dict in dict.items(): inner_dict[aa] = 1.0
	return dict
	
'''Dict to use for unknown components.'''
def svm_default_dict(component):
	return {component: 1.0}

'''Returns 1 on exact match, 0 otherwise. M_lookup isn't used - we just want to fix the length of parameter lists.'''
def exact_match_I(M, A, M_lookup):
	return 1.0 if M == A else 0.0

'''Returns an SVM I-value according to the rules in the original Pep2Path.'''
def svm_I(M, M_lookup=svm_class_matcher()):
	return M_lookup.get(M.lower(), 0.0)

six_classes = {
		"aliphatic": ["Ala", "Ile", "Leu", "Met", "Val"],
		"aromatic": ["Phe", "Trp", "Tyr"],
		"neutral_side_chains": ["Asn", "Cys", "Gln", "Ser", "Thr"],
		"acidic": ["Asp", "Glu"],
		"basic": ["Arg", "His", "Lys"],
		"unique": ["Gly", "Pro"]
	}
		   
'''Creates a doubly-nested dictionary that maps a pair of amino acids to whether their product class matches, 
	based on https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html'''
def six_class_matcher():
	
	class_list = [set(c) for c in six_classes]
	dict = {ele:{other_ele:0.5 for other_ele in (c - set(ele))} for c in class_list for ele in c}
	dict.update({ele:{ele:1.0} for ele in alphabet})
	return dict
	
'''Based on https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html
	Value is 1 for an exact match, 0.5 for a class match, 0 otherwise.'''
def six_class_I(M, A, M_lookup):
	return M_lookup.get(M, 0.0)
	
'''Dict to use for unknown components.'''
def six_class_default_dict(component):
	return {component: 1.0}
	
def create_pdict():
	with open("NRP_aa.txt", 'r') as file:
		counts = Counter([line.strip("\n").split("\t")[4] for line in file])
		return {AA.lower():(freq/(sum([freq for AA, freq in counts.items()]))) for AA, freq in counts.items()}
		
def create_svm_pdict():
	return {aa.lower():(1/len(svm_alphabet)) for aa in svm_alphabet}

def create_kersten_pdict():
	return {aa.lower():(1/len(kersten_alphabet)) for aa in kersten_alphabet}
			
def P(A, pdict):
	return pdict.get(A, 1/len(pdict.items()))
	
'''A naive matching scoring function that just counts the number of times seq1[i] == seq2[i] (if there's a choice '|', any of the possible options are acceptable).
	Calculates a normalised and inverted version of the Hamming Distance, (1-Hamming)'''
class simple_score():
		
	def setup(self, spectra, gbk):
		return None

	def score(self, seq1, seq2, spectrum_lookups=None):

		matches = mean([any([a in y for a in x.split('|')]) for x, y in zip(seq1, seq2)])
		
		if(len(seq1) == len(seq2)):
			return matches 
		elif(len(seq1) > len(seq2)):
			return max(matches, self.score(seq1[1:], seq2))
		else:
			return max(matches, self.score(seq1, seq2[1:]))
	
'''An implementation of the scoring mechanism described in the Pep2Path paper.
	In Is, give a list of functions meant to calculate I-values.
	The default settings might perform badly if you use things outside of the svm_alphabet.
	Be very careful with hashing values in dictionaries to make sure they've been normalised.'''
class p2p_score():

	def __init__(self, alphabet=AA_alphabet, Is=[svm_I], pdict=create_svm_pdict(), 
						I_dicts=[svm_class_matcher()], I_defaults=[svm_default_dict], 
						c=1, x=0.01, reg=2):
	
		self.alphabet = set([a.lower() for a in alphabet])
		self.Is = Is
		self.pdict = {k.lower() : v for k, v in pdict.items()}
		self.I_dicts = [{k.lower() : v for k, v in dict.items()} for dict in I_dicts]
		self.I_defaults = I_defaults
		self.c = c
		self.x = x
		self.reg = reg
		
		self.alphabet_dicts = [[dict.get(A, default(A)) for dict, default in zip(self.I_dicts, self.I_defaults)] for A in self.alphabet]
		self.spectrum_lookups = []
		
	def setup(self, spectrum, gbk):
		self.spectrum_lookups = [[dict.get(component.lower(), default(component.lower())) for component in spectrum] for dict, default in zip(self.I_dicts, self.I_defaults)]
		
	def score(self, spectrum_tag, genbank_ordering):
		
		I = lambda M, i: mean([I_func(M, lookups[i]) for I_func, lookups in zip(self.Is, self.spectrum_lookups)])
		denom_Is = lambda M: sum([mean([a_dict.get(M, 0) for a_dict in A]) for A in self.alphabet_dicts])
		
		s_ct = lambda M, a, i: math.log((P(a, self.pdict) + self.c * (I(M, i) ** self.reg + self.x * P(a, self.pdict)))
										/ (P(a, self.pdict) * (1 + self.c * (denom_Is(M) ** self.reg + self.x))))
		
		def best_alignment(genbank_ordering, spectrum_tag):
			score = sum([s_ct(M, a, i) for M, a, i in zip(genbank_ordering, spectrum_tag, range(len(genbank_ordering)))])
			if(len(genbank_ordering) == len(spectrum_tag)):
				return score 
			elif(len(genbank_ordering) > len(spectrum_tag)):
				return max(score, best_alignment(genbank_ordering[1:], spectrum_tag))
			else:
				return max(score, simple_score(genbank_ordering, spectrum_tag[1:]))

		return best_alignment(genbank_ordering, spectrum_tag)
	
'''Given a nested list of spectra then (ordered) components,
	and a twice-nested list of genbank files, then cds, then (ordered) components,
	attempts to score a match of each spectrum to each genbank file by comparing the components to each possible ordering of the cds,
	possibly with some cds reversed (as specified by the algorithm in the original Pep2Path paper).'''
def score_alignment(spectrum, gbk, scorer=simple_score()):
	
	spectrum = normalise_tag(spectrum)
	gbk = normalise_gbk(gbk)	
	
	scorer.setup(spectrum, gbk)
	
	best_score = float("-inf")
	
	mirror = lambda m : (m, list(reversed(m))) if (m != list(reversed(m))) else (m,)
	products = itertools.product(*[mirror(cds) for cds in gbk])
	orderings = [[comp for cds in ordering for comp in cds] for product in products for ordering in itertools.permutations(product, len(product))]
		
	for ordering in orderings:
		new_score = scorer.score(spectrum, ordering)
		best_score = max(best_score, new_score)
				
	return best_score
	
'''Used to score a group of spectral tags and BGCs against one another.'''
def score_alignments(spectrum_names, spectra, gbk_names, gbks, scorer=simple_score()):
	scores = [(s_name, g_name, score_alignment(s, g, scorer)) for s_name, s in zip(spectrum_names, spectra) for g_name, g in zip(gbk_names, gbks)]
	print(scores)
	best = max(scores, key=lambda s: s[2])
	avg = mean([score[2] for score in scores])
	return best, avg