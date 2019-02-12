import csv

#amino acids taken from supplementary table provided by kersten et al.
#https://www.nature.com/articles/nchembio.684
kersten_masses, kersten_names = {}, {}
with open("nchembio.684-S3.csv", 'r') as file:
	lines = [line for line in csv.reader(file)][2:] #split lines, cutting out headers
	kersten_masses = {line[2]:line[0] for line in lines}
	kersten_names = {line[2]:line[2] for line in lines}
	kersten_names.update({line[1]:line[2] for line in lines})
kersten_alphabet = set([name for name in kersten_masses.keys()])

#proteinogenic AAs
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
				"Leucine" : "Ile/Leu",
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
			
#uses avg mass instead of monoisotopic
AA_mass_table = {	
					"Ala" : 71.0779,
					"Arg" : 156.1857,
					"Asn" : 114.1026,
					"Asp" : 115.0874,
					"Cys" : 103.1429,
					"Gln" : 128.1292,
					"Gly" : 57.0513,
					"His" : 137.1393,
					"Ile/Leu" : 113.1576,
					"Lys" : 128.1723,
					"Met" : 131.1961,
					"Phe" : 147.1739,
					"Pro" : 97.1152,
					"Ser" : 87.0773,
					"Thr" : 101.1039,
					"Sec" : 150.0379,
					"Trp" : 186.2099,
					"Tyr" : 163.1733,
					"Val" : 99.1311
				}
			
AA_alphabet = set([name for name in AA_mass_table.keys()])