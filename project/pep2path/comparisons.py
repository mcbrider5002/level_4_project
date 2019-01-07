import itertools
import random
import os

from genbank.aanames import AA_names

'''Given a twice-nested list, chooses a random element from any of the innermost lists (i.e. randomly selects an element with three indices required to access it).'''
def random_component(gbk_files):
	if(len(gbk_files) - 1 <= 0): return -1, -1, -1
	file = random.randint(0, len(gbk_files) - 1) #generate an index for random file
	if(gbk_files[file] == []): return file, -1, -1
	cds = random.randint(0, len(gbk_files[file]) - 1) #generate an index for random cds in that file
	if(gbk_files[file][cds] == []): return file, cds, -1
	component = random.randint(0, len(gbk_files[file][cds]) - 1)
	return file, cds, component
	
'''Given a twice nested list, swaps two of the innermost elements at random.'''
def random_swap(gbk_files):
	file1, cds1, component1 = random_component(gbk_files)
	file2, cds2, component2 = random_component(gbk_files)
	if(component1 == -1 or component2 == -1): return gbk_files
	temp = gbk_files[file1][cds1][component1]
	gbk_files[file1][cds1][component1] = gbk_files[file2][cds2][component2]
	gbk_files[file2][cds2][component2] = temp
	return gbk_files
	
'''Given a nested list where internal lists represent genbank files and their contents are their unique components, 
	rearranges all the inner lists (files) and then randomly selects an animo acid to replace each of those components
	before filtering them to be unique again.'''
def randomise_components(gbk_components):
	random.shuffle(gbk_components) #change order of files to be random
	return [list(set([AA_names[str(random.choice(list(AA_names.keys())))] for component in gbk])) for gbk in gbk_components] #randomise all components
	
'''Given a nested list where internal lists represent a mass spectrum and their contents are their unique components,
	and a nested list where internal lists represent genbank files and their contents are their unique components,
	and two lists of names where the names correspond to the files in each of these lists
	returns several statistics based on comparisons of the unique components of these nested lists.'''
def compare_unique_components(spectra_names, spectra_components, gbk_names, gbk_components):

	total_correct = 0
	max_correct = ("No spectrum!", "No gbk!", 0)
	total_comparisons = 0
	total_spectra = 0
	total_gbk_comp = 0
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
			max_correct = max(max_correct, (spectrum_name, gbk_name, correct), key=lambda t: t[2])		
					
	return total_correct, max_correct, total_comparisons, total_spectra, total_gbk_comp, total_gbk_files
	
'''Given a nested list of spectra then (ordered) components,
	and a twice-nested list of genbank files, then cds, then (ordered) components,
	attempts to score a match of each spectrum to each genbank file by comparing the components to each possible ordering of the cds,
	possibly with some cds reversed.'''
'''def compare_alignment(spectra_components, gbk_components):
	
	best_score = ("No spectrum!", "No gbk!", 0)
	for file in files:
		products = itertools.product([(cds, "-".join(cds.split("-").reverse())) for cds in file])
		for product in products:
			orderings = itertools.permutations(product, len(product))
			for ordering in orderings:
				for spectrum in spectra_components:
					new_score = spectrum.alignment_comparison(ordering)
					best_score = max(best_score, ("spectrum_name", "gbk_name", new_score))
			
	return best_score'''