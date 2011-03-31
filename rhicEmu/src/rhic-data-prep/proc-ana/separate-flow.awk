# separate out the flow cpts
# we only want to keep the first three cpts again
/V2/{
		test = match($2, /\w*ELL[3-6]/) # see if the label strings match the cpts we don't need
		if(!test) { 
				print $2 "\t" $3 "\t" $4
		}
}
