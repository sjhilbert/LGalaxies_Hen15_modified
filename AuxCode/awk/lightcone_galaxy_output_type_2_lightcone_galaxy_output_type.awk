# from a cleaned list of variables defining the lightcone_galaxy_output_type C struct, create a .h file.
BEGIN{
print "typedef struct custom_lightcone_galaxy_output_type_ {"
}
{
	print $0 ";"
}
END{
print "} custom_lightcone_galaxy_output_type;"
}