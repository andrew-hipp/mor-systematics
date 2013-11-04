#This function was written on 10-07-2013 in order to pick&choose gbif species data already processed and present in the workspace, according to a master species list, so that they may be used to extract from WorldClim database & continue on in the processing chain of gbif functions
#specieslist is master species list containing new desired species for continued processing
#datafiles are the present objects in the workspace containing all gbif processed species & gbif data
collect_matches=function(specieslist, datafiles) {
new_specieslist=list()
temp=list()
for (i in datafiles) {
	data_working=get(i)
	temp=data_working[names(data_working) %in% specieslist]
	new_specieslist=c(new_specieslist,temp) }
return(new_specieslist)
}

#This function was written on 10-04-2013 to filter out NULL species matrices within the list, i.e. no specimen data was retrieved. 
#It will return a vector of the indices of all NULL matrices within the list
delete=function(gbifdata) {
collect=vector()
 for (i in 1:length(gbifdata)) {
	if (!is.logical(gbifdata[[i]]$precise_enough)) {
		collect=c(collect,i)
		} }
	return(gbifdata[-collect])
}
	
	#(empty matrices within my list 10-04-2013) 17,26,41,48,57,73,75,77,82, 104, 105, 107, 108, 122, 127, 128, 145, 146, 156, 186, 190, 195, 210, 216, 222, 224, 254, 314, 319, 367, 385, 395, 400