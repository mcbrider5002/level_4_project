import itertools
import copy
import random
import os
import math
from statistics import mean
import numpy as np
from collections import Counter

from masstables import kersten_names, kersten_alphabet, AA_names, AA_alphabet
	
'''Given a twice nested list, returns a new list with the internal structure and the elements themselves shuffled.'''
def shuffle_components(gbk_files):

	reorder = lambda ls: random.sample(ls, len(ls))
		
	cds_numbers = reorder([len(cds) for gbk in gbk_files for cds in gbk])
	gbk_numbers = reorder([len(gbk) for gbk in gbk_files])
	
	reimpose_structure = lambda ls, indices: [ls[i:j] for i, j in zip([0] + list(np.cumsum(indices)[:-1]), np.cumsum(indices))]
	flattened = reorder(list(itertools.chain.from_iterable(itertools.chain.from_iterable(gbk_files))))
	return reimpose_structure(reimpose_structure(flattened, cds_numbers), gbk_numbers)
	
'''Given a nested list where internal lists represent genbank files and their contents are their unique components, 
	rearranges all the inner lists (files) and then randomly selects an animo acid to replace each of those components
	before filtering them to be unique again.'''
def randomise_components(gbk_components):
	gbk_components = random.sample(gbk_components, len(gbk_components)) #change order of files to be random
	return [[[AA_names[str(random.choice(list(AA_names.keys())))] for component in cds] for cds in gbk] for gbk in gbk_components] #randomise all components
	
'''Given a nested list where internal lists represent a mass spectrum and their contents are their unique components,
	and a nested list where internal lists represent genbank files and their contents are their unique components,
	and two lists of names where the names correspond to the files in each of these lists
	returns several statistics based on comparisons of the unique components of these nested lists.'''
def compare_unique_components(spectra_names, spectra, gbk_names, gbk_components):

	#convert spectra tags to their unique components
	spectra_components = [tags.unique_components() for tags in spectra if tags.unique_components() != []]

	#flatten out cds (not relevant here) get rid of any duplicates
	gbk_components = [list(set(itertools.chain.from_iterable(gbk))) for gbk in gbk_components] 

	best_score = ("No spectrum!", "No gbk!", 0)
	total_correct, total_comparisons, total_spectra, total_gbk_comp = 0, 0, 0, 0
	total_gbk_files = len(gbk_components)
	
	for spectrum in spectra_components:
		total_spectra += len(spectrum)
	
	for gbk_name, gbk in zip(gbk_names, gbk_components):
	
		total_gbk_comp += len(gbk)
		if(len(gbk) == 0): continue
		
		for spectrum_name, spectrum in zip(spectra_names, spectra_components):
			correct = 0
			for component in spectrum:
				total_comparisons += 1 
				if(component.lower() in gbk or (component[0].upper() + component[:1].lower()) in gbk or component in gbk):
					correct += 1
			
			total_correct += correct
			best_score = max(best_score, (spectrum_name + " " + str(spectrum), gbk_name + " " + str(gbk), correct), key=lambda t: t[2])		
					
	return best_score, total_correct, total_comparisons, total_spectra, total_gbk_comp, total_gbk_files

'''A naive matching scoring function that just counts the number of times seq1[i] == seq2[i] (if there's a choice '|', any of the possible options are acceptable).'''
def simple_score(seq1, seq2):

	matches = [any([a in y.lower() for a in (x.lower().split('|'))]) for x, y in zip(seq1, seq2)].count(True)
	
	if(len(seq1) == len(seq2)):
		return matches 
	elif(len(seq1) > len(seq2)):
		return max(matches, simple_score(seq1[1:], seq2))
	else:
		return max(matches, simple_score(seq1, seq2[1:]))
		
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
	get_alphabet = lambda classes : set(itertools.chain.from_iterable([c.split(',') for c in classes]))
	alphabet = set(itertools.chain.from_iterable([value for key, value in three_classes.items()])) & get_alphabet(big_classes) & get_alphabet(small_classes)
	class_normaliser = lambda classes : [set(prod_class.split(',')) for prod_class in classes]
	create_lookup = lambda classes, text : {ele:{other_ele:text for other_ele in (c - set(ele))} for c in classes for ele in c}
	
	dict = create_lookup([set([ele for ele in c]) for c_name, c in three_classes.items()], "three_class_match")
	big_dict = create_lookup(class_normaliser(big_classes), "big_class_match")
	small_dict = create_lookup(class_normaliser(small_classes), "small_class_match")
	
	big_alphabet = itertools.chain.from_iterable([c.split(',') for c in big_classes])
	small_alphabet = itertools.chain.from_iterable([c.split(',') for c in small_classes])
	
	for aa in big_alphabet:
		if(not aa in dict): dict[aa] = {}
		inner_dict = dict[aa]
		inner_dict.update(big_dict.get(aa, {}))
		
	for aa in small_alphabet:
		if(not aa in dict): dict[aa] = {}
		inner_dict = dict[aa]
		inner_dict.update(small_dict.get(aa, {}))
		
	for aa, inner_dict in dict.items(): inner_dict[aa] = "exact_match"
	
	return dict

'''Returns 1 on exact match, 0 otherwise.'''
def exact_match_I(M, A, M_lookup):
	return 1.0 if M == A else 0.0

'''Returns an SVM I-value according to the rules in the original Pep2Path.'''
def svm_I(M, M_lookup):
	match = M_lookup.get(M.lower(), "N/A")
	return 0.0 if match == "N/A" else (
		   0.25 if match == "three_class_match" else (
		   0.5 if match == "big_class_match" else (
		   0.75 if match == "small_class_match" else 
		   1.0 if match == "exact_match" else 0)))

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
	dict.update({ele:{ele:"exact_match"} for ele in alphabet})
	return dict
	
'''Based on https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html
	Value is 1 for an exact match, 0.5 for a class match, 0 otherwise.'''
def six_class_I(M, A, M_lookup):
	return M_lookup.get(M, 0.0)
	
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
	
'''An implementation of the scoring mechanism described in the Pep2Path paper.
	In Is, give a list of functions meant to calculate I-values.'''		
def alignment_score(spectrum_tag, spectrum_lookups, genbank_ordering, alphabet, Is, pdict, c=1, x=0.01, reg=2):

	I = lambda M, i: mean([I_func(M, lookups[i]) for I_func, lookups in zip(Is, spectrum_lookups)])
	
	s_ct = lambda M, a, i: math.log((P(a, pdict) + c * (I(M, i) ** reg + x * P(a, pdict)))
						/ (P(a, pdict) * (1 + c * (sum([svm_class_matcher().get(M, {M:1}).get(A, 0) for A in alphabet]) ** reg + x))))
	
	def best_alignment(genbank_ordering, spectrum_tag):
		score = sum([s_ct(M, a, i) for M, a, i in zip(genbank_ordering, spectrum_tag, range(len(genbank_ordering)))])
		if(len(genbank_ordering) == len(spectrum_tag)):
			return score 
		elif(len(genbank_ordering) > len(spectrum_tag)):
			return max(score, best_alignment(genbank_ordering[1:], spectrum_tag))
		else:
			return max(score, simple_score(genbank_ordering, spectrum_tag[1:]))

	return best_alignment(genbank_ordering, spectrum_tag)

############################
#Fix function arity
SIMPLE_SCORE = lambda s, s_lookups, gbk, alphabet, Is, pdict, c, x, reg: simple_score(s, gbk)
P2P_SCORE = lambda s, s_lookups, gbk, alphabet, Is, pdict, c, x, reg: alignment_score(s, s_lookups, gbk, alphabet, Is, pdict, c=c, x=x, reg=reg)
############################	

'''Given a nested list of spectra then (ordered) components,
	and a twice-nested list of genbank files, then cds, then (ordered) components,
	attempts to score a match of each spectrum to each genbank file by comparing the components to each possible ordering of the cds,
	possibly with some cds reversed (as specified by the algorithm in the original Pep2Path paper).'''
def compare_alignment(spectra_names, spectra, gbk_names, gbk_tags, alphabet=AA_alphabet, svm_dict=svm_class_matcher(),
														pdict=create_svm_pdict(), scoring_method=P2P_SCORE, c=1, x=0.01, reg=2):
	
	#convert spectra tags into a list of their tied-longest tags
	spectra_tags = [(tags.decompose_tags())[tags.longest_tag] for tags in spectra if tags.decompose_tags() != {}]
	
	best_score = ("No spectrum!", "No gbk!", 0)
	
	for (gbk_name, gbk) in zip(gbk_names, gbk_tags):
	
		mirror = lambda m : (m, list(reversed(m))) if (m != list(reversed(m))) else (m,)
		products = itertools.product(*[mirror(cds) for cds in gbk])
		orderings = [ordering for product in products for ordering in itertools.permutations(product, len(product))]
		
		for ordering in orderings:
			for (spectrum_name, spectrum) in zip(spectra_names, spectra_tags):
			
				spectrum_lookups = [[svm_dict.get(component.lower(), {component:"exact_match"}) for component in spectrum]]
				Is = [svm_I]
			
				new_score = scoring_method(spectrum, spectrum_lookups, list(itertools.chain.from_iterable(ordering)), alphabet, Is, pdict, c, x, reg)
				best_score = max(best_score, (spectrum_name + " " + str(spectrum), gbk_name + " " + str(ordering), new_score), key=lambda t: t[2])
				
	return best_score