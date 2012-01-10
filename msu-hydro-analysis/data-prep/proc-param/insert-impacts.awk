## we will treat the impact-parameters as 
## additional paramters for each data set, so we want to insert
## them into the stats.param-combined data
{
		test = match($1, /^\#id/)
		if(test){
				sub("#", "", $0) # the result of the sub is still in $0
				print $0 ", IMPACT"  ## label the impact column and remove the leading #
		} else { 
				test = match($1, /^\#/)
				if(!test){
						nimpacts = 5
						impacts[1] = 3.29
						impacts[2] = 5.70
						impacts[3] = 7.36
						impacts[4] = 8.70
						impacts[5] = 9.87

						for( x = 1; x <= nimpacts; x++){
								print $0 "\t" impacts[x]
						}
				} else {
						print $0
				}
		}
}
