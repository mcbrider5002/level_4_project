import itertools
import os
import random
import numpy as np
from statistics import mean

from masstables import AA_alphabet

from . import comparisons as compare

'''This experiment tests our scoring function's behaviour as we iteratively introduce random mutations into a perfect match.
	Scores should generally decrease as number of mutations increases and the first score (as it is a perfect match) should be the (tied-)maximum.'''
def scoring_experiment(scoring_function, alphabet=AA_alphabet, length=random.randint(2, 8), batch_size=6, mutations=10):

	alphabet = list(alphabet)

	'''Generate partitions for cds regions in gbk.'''
	def generate_cds(length): 
		cds_length = random.randint(1, length)
		return generate_cds(length - cds_length) + [cds_length] if length - cds_length > 0 else [cds_length]
		
	spectrum = [random.choice(alphabet) for i in range(length)]
	
	cds_indices = np.cumsum([0] + generate_cds(len(spectrum - 1)))
	gbk = [spectrum[start:end] for start, end in zip(cds_indices[:-1], cds_indices[1:])] #partition
	gbk = list(random.sample(gbk, len(gbk))) #shuffle
	gbk = [(reversed(cds) if r else cds) for r, cds in zip([random.randint(0, 1) for i in range(len(gbk))], gbk)] #reverse some
	
	def random_mutation(labelled_spectrum):
		mutations, spectrum = labelled_spectrum
		index = random.randint(0, len(spectrum - 1))
		return ((mutations + 1), (spectrum[:index] + random.choice(AA_alphabet) + spectrum[index+1:]))
	
	labelled_spectrum = (0, spectrum)
	batch = [random_mutation(labelled_spectrum) for batch in batch_size]
	scores = [(0, spectrum, gbk, scoring_function(spectrum, gbk))] + [(m, s, gbk, scoring_function(s, gbk)) for m, s in batch]
	for m in range(1, mutations):
		batch = [random_mutation(batch) for batch in batch_size]
		scores.extend([(m, s, gbk, scoring_function(s, gbk)) for m, s in batch])
		
	return sorted(scores, key=lambda t: t[3])