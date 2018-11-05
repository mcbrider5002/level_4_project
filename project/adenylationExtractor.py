from Bio import SeqIO
from collections import defaultdict

from GenbankPredictions import GenbankPredictions

features = [feature for feature in SeqIO.read("c00001_NODE_1_...cluster001.gbk", "genbank").features if "aSProdPred" in feature.qualifiers.keys()]

domain_counts = defaultdict(lambda:0)
for feature in features: 
	for qualifier in feature.qualifiers.keys(): domain_counts[qualifier] += 1
	print(feature.qualifiers["aSProdPred"])
	for item in feature.qualifiers["sec_met"]:
		print(item)
		print()
	print("---")
#print(features[0].qualifiers)
print(domain_counts)

for feature in features: print("\n" + str(feature.qualifiers))
print([feature for feature in SeqIO.read("c00001_NODE_1_...cluster001.gbk", "genbank").features])