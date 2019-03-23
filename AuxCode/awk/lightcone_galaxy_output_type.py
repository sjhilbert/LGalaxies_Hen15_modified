# numpy dtype for custom_lightcone_galaxy_output_type_
import numpy
struct_dtype = numpy.dtype([
('ending','i4',0) # padding for mem alignment (remove if not needed)
])
properties_used = {}
for el in struct_dtype.names:
	properties_used[el] = False
