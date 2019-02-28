import itertools
import os
import math
from statistics import mean
from Bio import SeqIO

from .spectra.Tag import Tag
from .spectra.SpectrumTags import SpectrumTags

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

'''Given a list of SpectrumTags objects, compares them against a given Genbank file's sequence data using the RiPP2Path algorithm.'''	
def ripp2path(spectra_tags, path, file, noResults=100):

	gbk = get_sequence_file(path, file)
	seq_len = len(gbk.seq)
	frames = six_frame_translation(gbk)
	
	'''Transforms tags with 'Ile/Leu' into separate tags.'''
	def filter_ile(tag):
		return itertools.product(*[component.split('/') for component in tag])
	
	#Bio will return its translation as a string of single characters, so let's change the spectrum's format to match
	one_letter_AA = {'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 'ASP':'D', 'ASN':'N', 
					 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 'SER':'S', 
					 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 'PRO':'P', 'CYS':'C'}
	
	s_tags = [(spectrum_tags.id, tag) for spectrum_tags in spectra_tags 
											for s_tag in spectrum_tags.flatten_tags()
												for tag in filter_ile(s_tag.decompose_tag())]
	single_char_tags = [(id, "".join([one_letter_AA.get(comp.upper, 'X') for comp in tag])) for id, tag in s_tags]
	
	#regular strands will be labelled 0 to 2 and reverse complement strands will be labelled from (seq_len) to (seq_len - 2)
	nucpos_f = lambda f_no, i: abs((-seq_len) * math.floor(f_no / 3) + (f_no % 3 + i))
	strands = (['+'] * 3) + (['-'] * 3)
	
	scores = [(spectra_name, tag, frame[i:i+len(tag)], nucpos_f(frame_no, i*3), strand, compare_ripp(tag, frame[i:i+len(tag)]))
					for (frame_no, frame, strand), (spectra_name, tag) in itertools.product(zip(range(6), frames, strands), single_char_tags) 
							for i in range(len(frame) - len(tag))]
	
	return sorted(scores, key=lambda t: t[2], reverse=True)[:noResults]
	
def ripp_printer(headers, scores, label_width=4, column_width=14, outpath=os.path.dirname(__file__), outfile="ripps.out"):
	no_entries = len(headers)
	row_labels = "%-" + str(label_width) + "s|"
	output_string = ("%-" + str(column_width) + "s|") * no_entries
	with open(os.path.join(outpath, outfile), 'w') as out:
		out.write(row_labels % "" + output_string % headers)
		for i, score in enumerate(scores): out.write(row_labels % i + output_string % score)
	
def test_ripp2path():
	tags = {6 : [Tag("Val-His-Phe-Val-Gly-Trp-Ile/Leu", [], [])]}
	spectra_tags = [SpectrumTags("Test Sequence", 6, tags)]
	ripp_printer(("Spectrum Location", "Search Tag", "Gene Tag", "Start Position", "Strand", "Score"), 
					ripp2path(spectra_tags, os.path.dirname(__file__), "S_coelicolor.gbk"))