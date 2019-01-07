from Bio import SeqIO
import os
from glob import glob

'''Provides methods to select antiSMASH-generated Genbank files based on certain properties.
	(Right now just generates list of all nrps files.)'''
	
AA_names = {	
					"Ala" : "Ala",
					"Alanine" : "Ala",
					"Arg" : "Arg",
					"Arginine" : "Arg",
					"Asn" : "Asn",
					"Asparagine" : "Asn",
					"Asp" : "Asp",
					"Aspartic Acid" : "Asp",
					"Cys" : "Cys",
					"Cysteine" : "Cys",
					"Gln" : "Gln",
					"Glutamine" : "Gln",
					"Gly" : "Gly",
					"Glycine" : "Gly",
					"His" : "His",
					"Histidine" : "His",
					"Ile" : "Ile/Leu",
					"Isoleucine" : "Ile/Leu",
					"Leu" : "Ile/Leu",
					"Leucine" : "Leu",
					"Lys" : "Lys",
					"Lysine" : "Lys",
					"Met" : "Met",
					"Methionine" : "Met",
					"Phe" : "Phe",
					"Phenylalanine" : "Phe",
					"Pro" : "Pro",
					"Proline" : "Pro",
					"Ser" : "Ser",
					"Serine": "Ser",
					"Thr" : "Thr",
					"Threonine" : "Thr",
					"Sec" : "Sec",
					"Selenocysteine" : "Sec",
					"Trp" : "Trp",
					"Tryptophan" : "Trp",
					"Tyr" : "Tyr",
					"Tyrosine" : "Tyr",
					"Val" : "Val",
					"Valine" : "Val"
				}
				
'''Given a prediction as string, checks it isn't empty and only contains components from the keys of the given table of component names.'''
def components_in_dict(pred, table):
	if(pred == ""): return false
	components = itertools.chain.from_iterable([c.split('|') for c in pred.split('-')])
	return all([c in table.keys() for c in components])

'''Given a SeqRecord object returns the combined product types of all clusters.'''
def get_product_class(seq_rec):	
	products = [feature.qualifiers["product"][0] for feature in seq_rec.features if feature.type == "cluster"]
	return "/".join(list(set(products)))
	
'''Checks all clusters have the product type nrps.'''
def is_nrps(seq_rec):
	return get_product_class(seq_rec) == "nrps"

'''Given a SeqRecord object and a length n returns whether it contains predictions of total length greater than or equal to n.'''
def has_long_predictions(seq_rec, n):
	predictions = [feature.qualifiers["aSProdPred"][0] for feature in seq_rec.features if "aSProdPred" in feature.qualifiers.keys()]
	return sum([(pred.count("-") + 1 if pred.strip() != "" else 0) for pred in predictions]) >= n

'''Given a SeqRecord object, returns whether it has exactly one feature called 'cluster'.'''
def has_one_cluster(seq_rec):	
	return len([feature for feature in seq_rec.features if feature.type == "cluster"]) == 1
	
'''Given a list of filenames and SeqRecords in the format (name, record), returns any that have a cluster feature without the 'product' qualifier.'''
def check_has_product(files):
	return [(name, file) for name, file in files if any([feature.type == "cluster" and (not "product" in feature.qualifiers.keys()) for feature in file.features])]
	
def main():
	path = os.path.join(os.path.join(os.getcwd(), "genbank"), "justin-20181022")
	outpath = os.getcwd()
	
	names = [name for name in glob(os.path.join(path, "*.gbk")) if (not "final" in os.path.basename(name))]
	names = [(os.path.basename(pattern), SeqIO.read(pattern, "genbank")) for pattern in names]

	#find files with one cluster with class nrps and has predictions > length one
	names = [(name, file) for name, file in names if (has_one_cluster(file) and is_nrps(file) and has_long_predictions(file, 3))]
	
	out = open(os.path.join(outpath, "dataset.out"), 'w')
	out.write("\n".join([name for name, file in names]))
	out.close()
	
if __name__ == "__main__":

    main()