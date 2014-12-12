# generate adjacency matrices for input into igraph
# each node is a genotype
# connect nodes if they occur at the same site at the same time


library(igraph)
setwd("/Users/kimtsao/Desktop/Summer08/Mammal Data")
#setwd("/Users/Public/Github/ospC-communities")
#n08 <- read.delim("Nymphs08_ospC.txt", header=T, strip.white=T)
ticks <- read.delim("PeleTicksOspC.txt",header=T, colClasses="character") #colClasses avoids importing data as levels/factors
tissues <- read.csv("PeleTissueOspC.txt",header=T, colClasses="character")
#dates <- read.delim("cap-site-sess-date-tag.txt",header=T, colClasses="character")

# temp = unique(sort(c(ticks$OspCType, tissues$OspCType)))
# alltypes = temp[temp!=""]
# taglist = unique(sort(c(ticks$TagNumber, tissues$TagNumber)))

# tallies co-occurrences of types within inputted dataset (separated by timepoint)
makeAdjMatrix <- function(subdata,t){
	# make ntype x ntype matrix init with 0s
	ntypes = 17
	adjmat = matrix(0,nrow=ntypes,ncol=ntypes)
		rownames(adjmat) = c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","T","U")
		colnames(adjmat) = rownames(adjmat)
	
	# for specified time point t
	timesub = subset(subdata, subdata[,"SessionNumber"]==as.character(t))
		# for each individual
		mice = unique(timesub[,"TagNumber"])
		for (m in 1:length(mice)){
		# list all types in individual at timepoint
			mousesub = subset(timesub, timesub[,"TagNumber"]==mice[m])
			types = unique(mousesub[,"OspCType"])
				# remove blank "types"
				blanks = which(types==""); types = types[-blanks]
			if (length(types)==1){
				adjmat[types[1],types[1]] = adjmat[types[1],types[1]]+1
			}
			else if (length(types)>1){
			# all pairs of types found, add to adj matrix
				pairs = combn(types,2) # all possible unique pairs of genotypes in this indiv
				for (p in 1:ncol(pairs)){
					pair = pairs[,p]
					adjmat[pair[1],pair[2]] = adjmat[pair[1],pair[2]]+1
					adjmat[pair[2],pair[1]] = adjmat[pair[2],pair[1]]+1
				}
			}
		} # end "for each mouse"
 return(adjmat)
}

sites = as.character(unique(n08$Site))
	
# for each site
# table ospC types

for (s in sites){
	sitesub = subset(n08,n08$Site==s)
	dates = unique(sitesub[order(as.Date(sitesub$Date..Year.2008., format="%m/%d/%Y")),"Date..Year.2008."])
	for (d in dates){
		datesub = subset(sitesub,sitesub$Date..Year.2008.==d)
		d = gsub('/','.',d)
		d = gsub('2011','2008',d)
		mat = makeAdjMatrix(datesub)
		write(mat,paste("AM_Nymphs",s,"_",d,".txt",sep=""),ncolumns=17)
	}
}

# combine tick and tissue data
ticksandtissue = rbind(
	cbind(ticks$TagNumber,ticks$SiteName,ticks$SessionNumber,ticks$OspCType),
	cbind(tissues$TagNumber,tissues$SiteName,tissues$SessionNumber,tissues$OspCType))
colnames(ticksandtissue) = c("TagNumber","SiteName","SessionNumber","OspCType")

subdata = subset(ticksandtissue,ticksandtissue[,"SiteName"]==s)
mat = makeAdjMatrix(subdata,1)