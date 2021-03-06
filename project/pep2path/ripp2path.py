import itertools
import os
import math
from statistics import mean
from Bio import SeqIO

from .spectra.Tag import Tag
from .spectra.SpectrumTags import SpectrumTags

#Bio will return its translation as a string of single characters, so we use this to change the spectrum's format to match
one_letter_AA = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 
				 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 
				 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}

'''Load Genbank file for further processing.'''
def get_sequence_file(path, file, format="genbank"):
	return SeqIO.read(os.path.join(path, file), format)
	
'''Returns the six translation frames of a given Genbank file's sequence.'''
def six_frame_translation(gbk):
	return [gbk.seq[i:].translate() for i in range(3)] \
			+ [gbk.reverse_complement().seq[i:].translate() for i in range(3)]
	
'''Return the degree of matching between two tags.'''
def compare_ripp(spectrum, gbk):
	return mean([1.0 if s == g else 0.0 for s, g in zip(spectrum, gbk)])

'''Converts three-character AA names to single character codes used by BioPython.'''	
def convert_to_single_char(spectra_tags):

	'''Transforms tags with 'Ile/Leu' into separate tags.'''
	def filter_ile(tag): return itertools.product(*[component.split('/') for component in tag])
	
	s_tags = [(spectrum_tags.id, tag) for spectrum_tags in spectra_tags 
											for s_tag in spectrum_tags.flatten_tags()
												for tag in filter_ile(s_tag.decompose_tag())]

	return [(id, "".join([one_letter_AA.get(comp.upper(), 'X') for comp in tag])) for id, tag in s_tags]
	
'''Scoring mechanism for Ripp2Path algorithm (performs the actual work of the algorithm once all the data is formatted correctly).'''
def ripp2path_scorer(seq_len, single_char_tags, frames, no_results=100):

	#regular strands will be labelled 0 to 2 and reverse complement strands will be labelled from (seq_len) to (seq_len - 2)
	nucpos_func = lambda f_no, i: abs((-seq_len) * math.floor(f_no / 3) + (f_no % 3 + i))
	strands = (['+'] * 3) + (['-'] * 3)
	
	scores = [(spectra_name, tag, frame[i:i+len(tag)], nucpos_func(frame_no, i*3), strand, compare_ripp(tag, frame[i:i+len(tag)]))
					for (frame_no, frame, strand), (spectra_name, tag) in itertools.product(zip(range(6), frames, strands), single_char_tags) 
							for i in range(len(frame) - len(tag))]
	
	return sorted(scores, key=lambda t: t[5], reverse=True)[:no_results]
	
'''Given a list of SpectrumTags objects, compares them against a given Genbank file's sequence data using the RiPP2Path algorithm.'''	
def ripp2path(spectra_tags, path, file, no_results=100):
	gbk = get_sequence_file(path, file)
	seq_len = len(gbk.seq)
	frames = six_frame_translation(gbk)
	single_char_tags = convert_to_single_char(spectra_tags)
	return ripp2path_scorer(seq_len, single_char_tags, frames, no_results)

def ripp_printer(headers, scores, label_width=4, column_width=20, outpath=os.path.dirname(__file__), outfile="ripps.out"):
	no_entries = len(headers)
	row_labels = "%-" + str(label_width) + "s|"
	output_string = ("%-" + str(column_width) + "s|") * no_entries
	with open(os.path.join(outpath, outfile), 'w') as out:
		out.write(row_labels % "" + output_string % headers + "\n")
		for i, score in enumerate(scores): out.write(row_labels % i + output_string % score + "\n")
	
def test_ripp2path():
	tags = {6 : [Tag("Val-His-Phe-Val-Gly-Trp-Ile/Leu", [], [])]}
	spectra_tags = [SpectrumTags("Test Sequence", 6, tags)]
	ripp_printer(("Spectrum Location", "Search Tag", "Gene Tag", "Start Position", "Strand", "Score"), 
					ripp2path(spectra_tags, os.path.dirname(__file__), "S_coelicolor.gbk"))