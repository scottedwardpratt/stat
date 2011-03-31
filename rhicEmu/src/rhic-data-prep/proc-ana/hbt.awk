# match up the hbt stuff, sadly these don't have quite so easy a regex to match
#
# here we only want the zeroth and first moments
{
		# see if the line starts with R
		line = match($2, /^R\w*/)
		if(line){
				test = match($2, /\w*ELL[2-6]/)
				if(!test) {
						print $2 "\t" $3 "\t" $4
				}
		}
}