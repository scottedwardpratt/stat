# match up the spectra pars from a results.dat file
# now we only want the first few moments so we need to match on 
# the first column and reject ones which have ELL[3-6]
/SPECTRA/{
		test = match($2, /\w*ELL[3-6]/)
		if(!test) {
				print $2 "\t" $3 "\t" $4
		}
}