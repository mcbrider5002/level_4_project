import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import copy

from .massSpectraParser import load_files_from_dir
from .MassSpectrum import MassSpectrum
from .Tag import Tag
from .spectraMain import setup_mass_spectra

#A dict specifying the default colours to use for each line
default_colours = {	
					"Ala" : (0.33, 0, 0),
					"Arg" : (0, 0.33, 0),
					"Asn" : (0, 0, 0.33),
					"Asp" : (0.66, 0, 0),
					"Cys" : (0, 0.66, 0),
					"Gln" : (0, 0, 0.66),
					"Gly" : (0.33, 0.33, 0),
					"His" : (0.33, 0, 0.33),
					"Ile/Leu" : (0, 0.33, 0.33),
					"Lys" : (0.15, 0.15, 0.15),
					"Met" : (0.33, 0.33, 0.33),
					"Phe" : (0.5, 0.5, 0.5),
					"Pro" : (0.7, 0.7, 0.7),
					"Ser" : (1, 0, 0),
					"Thr" : (0, 1, 0),
					"Sec" : (0, 0, 1),
					"Trp" : (0.5, 0.25, 0.25),
					"Tyr" : (0.25, 0.5, 0.25),
					"Val" : (0.25, 0.25, 0.5)
				}

'''Given a spectrum and a tag, plots the tag above the spectrum.
	If you want the spectrum in its original form without any filtered-out readings, make sure to give it in that form.'''
def drawSpectrum(spectrum, tag, colour_key=default_colours):
	fig, ax = plt.subplots(figsize=(15, 8))
	
	tag_height = (spectrum.max_intensity()) * 1.01
	text_height = tag_height * 1.01
	
	split_tag = tag.decompose_tag()
	for compound, mass1, mass2 in zip(split_tag, tag.masses[0:-1], tag.masses[1:]):
		intensity_min = np.min(spectrum.ms2peaks[:, MassSpectrum.INTENSITY])
		ax.plot(np.array([mass1, mass1]), np.array([tag_height, intensity_min]), color="r", linestyle="dashed")
		ax.plot(np.array([mass2, mass2]), np.array([tag_height, intensity_min]), color="r", linestyle="dashed")
		ax.annotate(s='', xy=(mass2,tag_height), xytext=(mass1,tag_height), arrowprops=dict(arrowstyle='<->', shrinkA=0, shrinkB=0, color=colour_key.get(compound, "b")))
		ax.text(mass1, text_height, compound)

	for mass, intensity in zip(spectrum.ms2peaks[:, MassSpectrum.MASS], spectrum.ms2peaks[:, MassSpectrum.INTENSITY]):
		ax.plot(np.array([mass, mass]), np.array([0, intensity]), color="black")
		
	p_mass = spectrum.parent_mass
	ax.plot(np.array([p_mass, p_mass]), np.array([0, spectrum.max_intensity()*1.05]), color="cyan", linestyle="dashed", label="Parent Mass")
	ax.legend()
	
	ax.set(xlabel="Mass (m/z)", ylabel="Intensity", title=spectrum.id)
	ax.set_ylim(ymin=0)
	
	plt.tight_layout()

	plt.savefig(str(spectrum.id) + ".png")
	plt.show()
	
'''Plot a test spectrum.'''
def main():
	record = load_files_from_dir(path=os.path.join(os.path.dirname(__file__), "spectraData"), pattern="CCMSLIB00000078177.ms")[0]
	spectrum = setup_mass_spectra(record)
	spectrum.filter_intensity(intensity_threshold=spectrum.max_intensity()*0.05)
	spectrum.sort_by_mass()
	spectrum_tags = spectrum.find_sequence_tags()
	tags = spectrum_tags.tags[spectrum_tags.longest_tag]
	if(len(tags) > 0):
		drawSpectrum(spectrum, tags[0]) #take the first tag of the longest length
	else:
		print("No tags found in this spectrum!")
	
if __name__ == "__main__":

    main()
	