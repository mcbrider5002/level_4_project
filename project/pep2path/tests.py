import random
import itertools

from .masstables import AA_alphabet
from .comparisons import *
from .experiment import *
from .ripp2path import one_letter_AA, compare_ripp, ripp2path_scorer

from . import scoringExperiment as scoreexp

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
	
def test_score_unique_components(alphabet=AA_alphabet):
	scorer = scoreexp.score_unique_components
	scores = scoreexp.scoring_experiment(scorer, alphabet=AA_alphabet, length=random.randint(2, 8), batch_size=6, mutations=10)
	original = scores[0]
	assert original == sorted(scores, key=lambda t: t[3], reverse=True)[0], "score_unique_components doesn't give an exact match the tied-highest score!"
	
def test_compare_unique_components(alphabet=AA_alphabet):
	
	spectra = [normalise_tag(spectrum) for spectrum in generate_comp_lists(alphabet)]
	gbks = [normalise_gbk(gbk) for gbk in generate_multi_comps(alphabet)]
	
	flattened_gbks = [list(itertools.chain.from_iterable(gbk)) for gbk in gbks]
	scores = [(str(s), str(g), score_unique_components(s, g)) for s in spectra for g in flattened_gbks]
	
	best_score = max(scores, key=lambda s: s[2])
	score_sum = sum([s[2] for s in scores])
	avg = mean([s[2] for s in scores])
	
	result = compare_unique_components([str(s) for s in spectra], spectra, [str(g) for g in flattened_gbks], gbks)
	
	et = 0.001 #error tolerance - deal with float precision issues
	assert result[0][0] == best_score[0] and result[0][1] == best_score[1] and result[0][2] - best_score[2] < et, "compare_unique_components did not return the correct maximum score!"
	assert result[1] - score_sum < et, "compare_unique_components did not return the correct score sum!"
	assert result[2] - avg < et, "compare_unique_components did not return the correct average!" 
	
def test_simple_score(alphabet=AA_alphabet):

	comps = list(itertools.chain.from_iterable(generate_comp_lists(alphabet)))
	length = len(comps)
	
	no_dummies = random.randint(0, length)
	dummies = random.sample(list(range(length)), no_dummies)
	
	r_padder = lambda r : ["dummy" for i in range(random.randint(0, r))]
	dummy_comps = r_padder(5) + [("dummy" if i in dummies else comp) for i, comp in enumerate(comps)] + r_padder(5)
	
	assert simple_score().score(comps, dummy_comps) == (length - no_dummies) / length, "simple_score does not give the expected value!"
	
def test_exact_match_I(alphabet=AA_alphabet):

	aa = random.choice(list(AA_alphabet))
	M_lookup = {}
	assert exact_match_I(aa, aa, {}) == 1, "exact_match_I doesn't evaluate to 1 on a match!"
	assert exact_match_I(aa, None, {}) == 0, "exact_match_I doesn't evaluate to 0 on a mismatch!"
	
'''This is a test of both p2p_score and score_alignment; p2p_score is a complicated numerical mess and not amenable to individual testing,
	and alignment_score doesn't work without a component. So we test compare_alignment with the simple_score component which we've already tested, then the unverified
	p2p_score.'''
def test_score_alignment(alphabet=AA_alphabet):

	scorer = scoreexp.score_simple_alignment
	scores = scoreexp.scoring_experiment(scorer, alphabet=AA_alphabet, length=random.randint(2, 4), batch_size=6, mutations=10)
	original = scores[0]
	assert original == sorted(scores, key=lambda t: t[3], reverse=True)[0], "score_alignment (with simple score) doesn't give an exact match the tied-highest score!"

	scorer = scoreexp.score_p2p_alignment
	scores = scoreexp.scoring_experiment(scorer, alphabet=AA_alphabet, length=random.randint(2, 4), batch_size=6, mutations=10)
	original = scores[0]
	assert original == sorted(scores, key=lambda t: t[3], reverse=True)[0], "score_alignment (with p2p score) doesn't give an exact match the tied-highest score!"
	
def test_score_alignments(alphabet=AA_alphabet):
	
	def gen_comps(): return [[str(random.choice(list(alphabet))) for length in range(random.randint(2, 4))] for number in range(random.randint(2, 4))]
	
	spectra = gen_comps()
	gbks = [gen_comps() for number in range(random.randint(2, 4))]
	scorer = p2p_score()
		
	scores = [score_alignment(spectrum, gbk, scorer=scorer) for spectrum in spectra for gbk in gbks]
	best, avg = score_alignments(["" for s in spectra], spectra, ["" for g in gbks], gbks, scorer=scorer)
	
	assert best[2] == max(scores), "score_alignments doesn't give the same best score as score_alignment run multiple times!"
	assert avg == mean(scores), "score_alignments doesn't give the same mean as score_alignment run multiple times!"
	
def nrp_tests():
	test_shuffle_components()
	test_randomise_components()
	test_score_unique_components()
	test_compare_unique_components()
	test_simple_score()
	test_exact_match_I()
	test_score_alignment()
	test_score_alignments()
	
def	test_compare_ripp(alphabet=AA_alphabet):

	comps = list(itertools.chain.from_iterable(generate_comp_lists(alphabet)))
	length = len(comps)
	
	no_dummies = random.randint(0, length)
	dummies = random.sample(list(range(length)), no_dummies)
	dummy_comps = [("dummy" if i in dummies else comp) for i, comp in enumerate(comps)]
	
	assert compare_ripp(comps, dummy_comps) == (length - no_dummies) / length, "compare_ripp does not give the expected value!"
	
def test_ripp2path_scorer(alphabet=AA_alphabet):

	AAs = [v for k, v in one_letter_AA.items()]
	def gen_gene_seq(minm, maxm):
		return [str(random.choice(AAs)) for length in range(random.randint(minm, maxm))]
		
	loc = random.randint(1000, 2000)	
	true_match = gen_gene_seq(6, 10)
	genome = gen_gene_seq(loc, loc) + true_match + gen_gene_seq(1000, 2000)
	seq_len = len(genome) * 3
	
	frames = [genome] + [gen_gene_seq(2006, 4010) for i in range(5)] #we'll just spoof frames - those are BioPython's domain anyway
	scores = ripp2path_scorer(seq_len, [("id", true_match)], frames)
	scores = [score for score in scores if score[5] == 1]
	assert any([score[2] == true_match and score[3] == loc * 3 and score[4] == '+' for score in scores]), "ripp2path_scorer doesn't return a perfect match as the best score!"
	
def ripp_tests():
	test_compare_ripp()
	test_ripp2path_scorer()