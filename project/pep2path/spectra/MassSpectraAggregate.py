from .Tag import Tag
from .SpectrumTags import SpectrumTags

'''A wrapper class to make multiple operations on MassSpectra easier and with the use of less boilerplate.
	It is NOT a guarantee of more efficiency, it's just a convenience.
	Functions have both local and global variants, where 'local' performs the operation on each file individually whereas 'global' applies it to all
	using the data accessible from all of them.'''
class MassSpectraAggregate():

	def __init__(self, spectra):
		#List (or other iterable) containing spectra objects
		self.spectra = spectra
		
	def __repr__(self):
		return "\n\n".join(str(s) for s in self.spectra)
		
	'''Returns a list of spectra ids for the spectra in the order that they're held.
		You can use this to find out which spectrum a result was generated from.'''
	def get_spectra_ids(self):
		return [spectrum.id for spectrum in self.spectra]
	
	'''Given a list the same length as the internal list of spectra, returns a new list of tuples containing the spectra ids and the elements of the original list.
		You can use this to attach spectra ids to data generated from that spectrum.'''
	def attach_spectra_ids(self, attach_to):
		return list(zip(self.get_spectra_ids(), attach_to))
	
	'''Given a list of ids (or a single id) returns a new list of spectra from self.spectra with matching ids.'''
	def get_spectra_by_id(self, ids):
		ids = [ids] if isinstance(ids, str) else ids
		return list(filter(lambda spectrum: spectrum.id in ids, self.spectra))
	
	'''Returns a list of the maximum masses in each spectrum.'''
	def max_mass_local(self):
		return [spectrum.max_mass() for spectrum in self.spectra]
		
	'''Returns the maximum mass across all spectra.'''
	def max_mass_global(self):
		masses = self.max_mass_local()
		return max(masses) if len(masses) > 0 else 0
	
	'''Returns a list of the maximum intensities across all spectra.'''
	def max_intensity_local(self):
		return [spectrum.max_intensity() for spectrum in self.spectra]
		
	'''Returns the maximum intensity across all spectra.'''
	def max_intensity_global(self):
		intensities = self.max_intensity_local()
		return max(intensities) if len(intensities) > 0 else 0
		
	'''Sorts all spectra by mass.'''
	def sort_by_mass(self):
		for spectrum in self.spectra:
			spectrum.sort_by_mass()
		
	'''Takes a list of intensity thresholds, and filters the intensity of each item in spectra by the corresponding item in the list.
		Default is deferred to spectrum class.'''
	def filter_intensity_local(self, intensity_thresholds=None):
		intensity_thresholds = [None for spectrum in spectra] if intensity_thresholds is None else intensity_thresholds
		for spectrum, threshold in zip(self.spectra, intensity_thresholds):
			spectrum.filter_intensity(intensity_threshold=threshold)
		
	'''Takes an intensity threshold, and filters the intensity of each item in spectra by it.
		Default is 5% of the maximum intensity across all spectra.'''
	def filter_intensity_global(self, intensity_threshold=None):
		intensity_threshold = 0.05 * self.max_intensity_global() if intensity_threshold is None else intensity_threshold
		for spectrum in self.spectra:
			spectrum.filter_intensity(intensity_threshold=intensity_threshold)
	
	'''Takes a list of new maximum intensities, and normalises each spectrum to the corresponding item in the list.
		Default is deferred to the spectrum class.'''
	def normalise_intensity_local(self, new_maxes=None):
		new_maxes = [None for spectrum in spectra] if new_maxes is None else new_maxes
		for spectrum, maxm in zip(self.spectra, new_maxes):
			spectrum.normalise_intensity(new_max=maxm)
	
	'''Takes a new maximum intensity, and normalises all spectra according to it.
		Default is deferred to the spectrum class.'''
	def normalise_intensity_global(self, new_max=None):
		old_max = self.max_intensity_global()
		for spectrum in self.spectra:
			spectrum.normalise_intensity(new_max=new_max, old_max=old_max)
		
	'''Returns a list of SpectrumTags, containing the sequence tags for each spectrum.
		Defers mass_tolerance_mode default to the spectrum class, but the default value for the mass threshold is defined here as (10^(-5))
		(10ppm with ppm mode).'''
	def find_sequence_tags(self, mass_tolerance_mode=None, mass_threshold=0.00001):
		return [spectrum.find_sequence_tags(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold) for spectrum in self.spectra]
		
	'''Returns a pair of the longest sequence tag across all spectra, and a list of SpectrumTags.
		Defers mass_tolerance_mode default to the spectrum class, but the default value for the mass threshold is defined here as (10^(-5))
		(10ppm with ppm mode).'''
	def find_longest_tag(self, mass_tolerance_mode=None, mass_threshold=0.00001):
		spectra_tags = self.find_sequence_tags(mass_tolerance_mode=mass_tolerance_mode, mass_threshold=mass_threshold)
		return (max([spectrum_tags.longest_tag for spectrum_tags in spectra_tags]) if len(spectra_tags)>0 else 0), spectra_tags
		