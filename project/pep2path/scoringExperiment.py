import itertools
import os
import random
import math
import numpy as np
from statistics import mean
import matplotlib
import matplotlib.pyplot as plt

from .masstables import AA_alphabet

from . import comparisons as compare

'''This experiment tests our scoring function's behaviour as we iteratively introduce random mutations into a perfect match.
	Scores should generally decrease as number of mutations increases and the first score (as it is a perfect match) should be the (tied-)maximum.'''
def scoring_experiment(scoring_function, alphabet=AA_alphabet, length=random.randint(2, 8), batch_size=6, mutations=10):

	alphabet = list(alphabet)

	'''Generate partitions for cds regions in gbk.'''
	def generate_cds(length): 
		cds_length = random.randint(math.floor(length/2)+1, length)
		return generate_cds(length - cds_length) + [cds_length] if length - cds_length > 0 else [cds_length]
		
	spectrum = [random.choice(alphabet) for i in range(length)]
	
	ies = generate_cds(len(spectrum) - 1)
	cds_indices = np.cumsum([0] + ies)
	gbk = [spectrum[start:end] for start, end in zip(cds_indices[:-1], cds_indices[1:])] + [[spectrum[-1]]]#partition
	gbk = list(random.sample(gbk, len(gbk))) #shuffle
	gbk = [(list(reversed(cds)) if r else cds) for r, cds in zip([random.randint(0, 1) for i in range(len(gbk))], gbk)] #reverse some
	
	def random_mutation(labelled_spectrum):
		mutations, spectrum = labelled_spectrum
		index = random.randint(0, len(spectrum) - 1)
		return ((mutations + 1), (spectrum[:index] + [random.choice(list(AA_alphabet))] + spectrum[index+1:]))
	
	labelled_spectrum = (0, spectrum)
	batch = [random_mutation(labelled_spectrum) for b in range(batch_size)]
	scores = [(0, spectrum, gbk, scoring_function(spectrum, gbk))] + [(m, s, gbk, scoring_function(s, gbk)) for m, s in batch]
	for m in range(mutations):
		batch = [random_mutation(b) for b in batch]
		scores.extend([(m, s, gbk, scoring_function(s, gbk)) for m, s in batch])
		
	return scores

def plot_results(scores, title, loc, cmap='plasma'):

	muts = [mutation for mutation, _, _, _ in scores]	
	xs = list(range(max(muts) + 1))
	
	ys = [mean([score for m, _, _, score in scores if m == x]) for x in xs]
	
	fig, ax = plt.subplots()
	ax.scatter(xs, ys, c=ys, cmap=cmap)
	ax.set(xlabel="No. Random Mutations", ylabel="Similarity Score", title=title)
	ax.set_ylim(ymin=0)
	plt.tight_layout()
	plt.savefig(loc)
	plt.show()
	
'''Need to flatten gbk first.'''
def score_unique_components(s, g):
	s = compare.normalise_tag(s)
	g = compare.normalise_tag(list(itertools.chain.from_iterable(g)))
	return compare.score_unique_components(s, g)[0]
		
def simple_alignment(s, g):
	s = compare.normalise_tag(s)
	g = compare.normalise_gbk(g)
	score = compare.score_alignment(s, g, compare.simple_score)	
	return score
	
def run_scoring_experiment():
	
	#intersection experiment
	def run_jaccard(b=8, length=4, mutations=100):
		scores = scoring_experiment(score_unique_components, batch_size=b, length=length, mutations=mutations)
		plot_results(scores, "Random Mutations Against Score, Jaccard Similarity, Length %d tag" % length, "jaccard%d.png" % length, cmap='plasma')
	
	run_jaccard(b=30, length=2, mutations=100)
	run_jaccard(b=30, length=4, mutations=100)
	run_jaccard(b=30, length=6, mutations=100)
	run_jaccard(b=30, length=8, mutations=100)
	run_jaccard(b=30, length=20, mutations=100)
	
	#simple alignment experiment
	scores = scoring_experiment(simple_alignment, batch_size=10, length=4, mutations=100)
	plot_results(scores, "Random Mutations Against Score, Simple Alignment Score", "simple_alignment.png", cmap='plasma')
	
	#p2p experiment