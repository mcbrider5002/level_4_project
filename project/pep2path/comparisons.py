import itertools
import copy
import random
import os
import statistics
import numpy as np

from genbank.aanames import AA_names, AA_alphabet
	
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
	
############################
#Fix function arity
SIMPLE_SCORE = lambda s, gbk, alphabet, c, x, reg: simple_score(s, gbk)
#P2P_SCORE = lambda s, gbk, alphabet, c, x, reg: p2p_score(s, gbk, alphabet, c=c, x=x, reg=reg)
############################
	
'''Given a nested list of spectra then (ordered) components,
	and a twice-nested list of genbank files, then cds, then (ordered) components,
	attempts to score a match of each spectrum to each genbank file by comparing the components to each possible ordering of the cds,
	possibly with some cds reversed (as specified by the algorithm in the original Pep2Path paper).'''
def compare_alignment(spectra_names, spectra, gbk_names, gbk_tags, alphabet=AA_alphabet, scoring_method=SIMPLE_SCORE, c=1, x=0.01, reg=2):
	
	#convert spectra tags into a list of their tied-longest tags
	spectra_tags = [(tags.decompose_tags())[tags.longest_tag] for tags in spectra if tags.decompose_tags() != {}]
	
	best_score = ("No spectrum!", "No gbk!", 0)
	
	for (gbk_name, gbk) in zip(gbk_names, gbk_tags):
	
		mirror = lambda m : (m, list(reversed(m))) if (m != list(reversed(m))) else (m,)
		products = itertools.product(*[mirror(cds) for cds in gbk])
		orderings = [ordering for product in products for ordering in itertools.permutations(product, len(product))]
		
		for ordering in orderings:
			for (spectrum_name, spectrum) in zip(spectra_names, spectra_tags):
				new_score = scoring_method(spectrum, list(itertools.chain.from_iterable(ordering)), alphabet, c, x, reg)
				best_score = max(best_score, (spectrum_name + " " + str(spectrum), gbk_name + " " + str(ordering), new_score), key=lambda t: t[2])
			
	return best_score