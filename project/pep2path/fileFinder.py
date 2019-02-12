from Bio import SeqIO
import os
from glob import glob
from collections import Counter
from genbank.aanames import AA_names

'''Provides methods to select antiSMASH-generated Genbank files based on certain properties.'''
				
'''Given a prediction as string, checks it isn't empty and only contains components from the keys of the given table of component names.'''
def components_in_dict(pred, table):
	if(pred == ""): return false
	components = itertools.chain.from_iterable([c.split('|') for c in pred.split('-')])
	return all([c in table.keys() for c in components])

'''Given a SeqRecord object returns the combined product types of all clusters.'''
def get_product_class(seq_rec):	
	products = [feature.qualifiers["product"][0] for feature in seq_rec.features if feature.type == "cluster"]
	return "/".join(list(set(products)))
	
'''Checks all clusters have the product type NRPS.'''
def is_nrps(seq_rec):
	return get_product_class(seq_rec).lower() == "nrps"

'''Checks all clusters have the product type RiPP.'''	
def is_ripp(seq_rec):
	subclasses = ["ripp", "lantipeptide", "bottromycin", "cyanobactin", "glycocin", "lassopeptide", "linaridin", "linearazol", 
					"microcin", "sactipeptide", "thiopeptide", "auto_inducing_peptide", "comx", "bacterial_head_to_tail_cyclized"]
	return get_product_class(seq_rec).lower() in subclasses 

'''Given a SeqRecord object and a length n returns whether it contains predictions of total length greater than or equal to n.'''
def has_long_predictions(seq_rec, n):
	predictions = [feature.qualifiers["aSProdPred"][0] for feature in seq_rec.features if "aSProdPred" in feature.qualifiers.keys()]
	return sum([(pred.count("-") + 1 if pred.strip() != "" else 0) for pred in predictions]) >= n

'''Given a SeqRecord object, returns whether it has exactly one feature called 'cluster'.'''
def has_one_cluster(seq_rec):	
	return len([feature for feature in seq_rec.features if feature.type == "cluster"]) == 1
	
'''Given a list of filenames and SeqRecords in the format (name, record), returns any that have a cluster feature without the 'product' qualifier.'''
def has_no_product(files):
	return [(name, file) for name, file in files if any([feature.type == "cluster" and (not "product" in feature.qualifiers.keys()) for feature in file.features])]
	
'''Given a list of SeqRecord objects counts the unique combined product types for all clusters.'''
def product_class_freqs(seq_recs):
	return Counter([get_product_class(seq_rec) for seq_rec in seq_recs])	
	
'''Given a SeqRecord object returns whether it has a lanc_like pfam domain.'''	
def has_lanclike(seq_rec):
	return any(["lanc_like" in map(lambda x: x.lower(), feature.qualifiers.get("domain", [])) for feature in seq_rec.features])
	
def main():
	path = os.path.join(os.path.join(os.getcwd(), "genbank"), "justin-20181022")
	outpath = os.getcwd()
	
	names = [name for name in glob(os.path.join(path, "*.gbk")) if (not "final" in os.path.basename(name))]
	gbks = [(os.path.basename(pattern), SeqIO.read(pattern, "genbank")) for pattern in names]
	
	print("\n".join(sorted(["%s : %d" % (name, count) for name, count in product_class_freqs([gbk for name, gbk in gbks]).items()])))

	#find files with one cluster with class NRPS and has predictions > length two
	nrps_gbks = [(name, file) for name, file in gbks if (has_one_cluster(file) and is_nrps(file) and has_long_predictions(file, 3))]
	
	with open(os.path.join(outpath, "dataset.out"), 'w') as out:
		out.write("\n".join([name for name, file in nrps_gbks]))
	
	#find files with one cluster with class RiPP
	ripp_gbks = [(name, file) for name, file in gbks if (has_one_cluster(file) and is_ripp(file) and has_lanclike(file))]
	
	with open(os.path.join(outpath, "ripp_dataset.out"), 'w') as out:
		out.write("\n".join([name + " " + str(get_product_class(file)) for name, file in ripp_gbks]))
	
if __name__ == "__main__":

    main()
