# extract a block starting with a line containing lightcone_galaxy_output_type_ and ending with a }
BEGIN {
flag=0
line=""
d=""
}
/custom_lightcone_galaxy_output_type_/,/\}/ {   
line=$0
if(flag == 0 && line ~ /\{/) {
	ic=index(line,"{")
	line=substr(line,ic+1)
	flag=1
}
if(line ~ /\}/) {
	ic=index(line,"}")
	if(ic > 0)
	line=substr(line,0,ic-1)
}
if(flag==1 && line != "" && line !~ /^#/) 
{
    isc=index(line,";")
    if(isc == 0)
      d=d" "line
    else {
	    while(isc > 0)
	    {
	      d= d" "substr(line,0,isc-1)
	      # declaration d is now complete, do something with it
	      print d
	      # reset, see if multiple declarations on one line
	      d=""
	      line=substr(line,isc+1)
	      isc=index(line,";")
	    }
	    d=line
	}
}
}