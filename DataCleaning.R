
#working directory for the species data
setwd("C:/Users/admin/Documents/Felts/Puget Sound DATA CLEANING")


#load in the very raw data
rawDat  <- read.csv("rawdata.csv")
head(rawDat)
str(rawDat)
dim(rawDat)

install.packages("tidyr")
#load the tidyr package to clean data
library(tidyr)

#remove columns that definitely won't be used
dat2<- subset(rawDat, select=-c(ID,CollectionMethodCode,FieldReplicate,FieldSampleID, LifeStageCode, DistinctCode, Result, TaxonomicQualifier, QACode, ResQualCode, LabSampleID, BenthicResultsComments, Phylum.1, Class.1, Orders.1, Family.1))
dat2$ExcludedTaxa = NULL

#from here we can sort out for the things that I am interested in, comparing across sites (maybe based on watershed), but we 
#will definitely want to standardize for just a few months, so that we are comparing within the summer months
View(dat2)

#Only the rows I care about but with all Taxa
AllDat<- subset(dat2, select=c(CodeDate,Genus,BAResult))
View(AllDat)

#this is for collapsing sites with multiple of the same genus and summing the abundance numbers up
library(plyr)
SumAllDat<-ddply(AllDat,~CodeDate + Genus,summarise,BAResult=sum(BAResult)) #looking good, lets spread it

#Spreads into site x species matrix
AllSpread = spread(SumAllDat, Genus, BAResult, fill = 0)
head(AllSpread) #looks good!
dim(AllSpread)

#clean-up and turn CodeDate into a row.name
library(tibble)

row.names(AllSpread) = AllSpread$CodeDate
head(AllSpread)

## lets remove the Site column
AllSpread$CodeDate = NULL
head(AllSpread)

#load packages useful for community data
library(ade4)
library(vegan)
library(FD)

#now to make and load in the spatial data
spa <- subset(dat2, select=c(CodeDate,Latitude,Longitude)) #cool but we have a ton of duplicates we can get rid of 
spaSites <- unique(spa)
head(spaSites) #now we have the individual sites rather than a bunch of duplicates

#make the row names the site names
row.names(spaSites) = spaSites$CodeDate
spaSites$CodeDate = NULL

#we can plot our sites
plot (spaSites, asp=1, type="n",  main = "Sites",
      xlab = "X coordinate (km)", ylab = "Y coordinate (km)")
text(spaSites, row.names(spaSites), cex=0.8, col="red") #This is a ton of sites though 

plot(spaSites, asp=1, col="brown", cex=AllSpread$Procambarus, main="Procambarus") # this is relating our species sheet to the spatial data, here we show where an invasive species of crayfish is found in California!

#going back to the species data set, we can sum the captures of species across all sites
spe.pres <- apply(AllSpread > 0, 2, sum)
sort(spe.pres) #and sort them in increasing order
#for most studies of beta diversity we want to get rid of the extremely rare "singleton" species, we will do this later...

#we can easily calculate species relative abundance
spe.rel <- 100*spe.pres/nrow(AllSpread)
round(sort(spe.rel), 1) #this rounds the output to 1 decimal place

#plot histograms of rel.abundance
windows(title="Frequency Histograms",8,5)
par(mfrow=c(1,2))

hist(spe.pres, main="Species Occurrences", right=FALSE, las=1,
     xlab="Number of occurrences", ylab="Number of species",
     breaks=seq(0,2000,by=5), col="bisque")

hist(spe.rel, main="Species Relative Frequencies", right=FALSE, las=1,
     xlab="Frequency of occurrences (%)", ylab="Number of species",
     breaks=seq(0, 100, by=10), col="bisque")


#Describe the data with basic diversity metrices

N0 <- rowSums(AllSpread > 0) ## Species richness
H <- diversity(AllSpread) ## Shannon entropy
N1 <- exp(H) ## Shannon diversity (number of abundant species)
N2 <- diversity(AllSpread, "inv") ## Simpson diversity (number of dominant species)
J <- H/log(N0) ## Pielou evenness
E10 <- N1/N0 ## Shannon evenness (Hill's ratio)
E20 <- N2/N0 ## Simpson evenness (Hill's ratio)
div <- data.frame(N0, H, N1, N2, E10, E20, J) ## create a dataframe for the above measures
div

plot(spaSites, asp=1, col="blue", cex=H, main="Shannon Diversity") # plot of diversity in California, looks like the Sierras are pretty diverse

#Make a rank-abundance plot of the data
spe.tot <- apply(AllSpread, 2, sum)
sort(spe.tot)
plot(radfit(spe.tot)) #We get an error here which I believe is due to the singleton values or empty sites

sort(rowSums(AllSpread)) #there are a few sites with only one count in them
sort(colSums(AllSpread)) #and many species that only occur once

#First we will begin by removing singleton species 
removeSingle<-AllSpread[,colSums(AllSpread>0)>1,drop=FALSE]
removeSingle
colSums(removeSingle)
sort(colSums(removeSingle))
sort(rowSums(removeSingle))

spe.tot2 <- apply(removeSingle, 2, sum)
sort(spe.tot2)
plot(radfit(spe.tot2)) #still getting an error with the removed species so I am going to export to an Excel sheet and check for non-numeric data in my spreadsheet

write.table(removeSingle, "D:/EcoRSWAMP/checkforNA.txt", sep="\t") 
SpeMatrix<-as.matrix(removeSingle)
head(SpeMatrix)      

spe.tot3 <- apply(SpeMatrix, 2, sum)
sort(spe.tot3)
plot(radfit(spe.tot3))

#alright maybe the error doesn't really matter?
#let's read in the environmental variables
rawEnv <-read.csv("RawEnv.csv")
tEnv <- subset(rawEnv, select=c(CodeDate,Variable,Result))
#this looks better!

spreadEnv = spread(tEnv, Variable, Result, fill = NA)

row.names(spreadEnv) = spreadEnv$CodeDate
head(spreadEnv)

spreadEnv$CodeDate = NULL
head(spreadEnv)

pairs(spreadEnv)
#hard to tell if anything is correlated because nothing is standardized... will need to get to this later

#back to species matrices

#measure the dissimilarity in species community composition across sites
spe.db <- vegdist(removeSingle, 'bray')
head(spe.db)
str(spe.db)
  
#Bray-Curtis Dissimilarity on relative abundances
spe.rel <- decostand(removeSingle, "total")
spe.dbrel <- vegdist(spe.rel, 'bray')
str(spe.dbrel)

#Chord distance matrix
spe.norm <- decostand(removeSingle, "normalize")
spe.dc <- vegdist(spe.norm, "euclidean")
str(spe.dc)

#Hellinger distance matrix
spe.hel <- decostand(removeSingle, "hellinger")
spe.dh <- vegdist(spe.hel, "euclidean")
str(spe.dh)

#assign a fake environment to coerce betadisp to run (could use HUC watersheds to group in the future)
env1 = rep(x = "F1", times = 1950)
macroB = betadisper(spe.db, env1)
macroB
mean(spe.db)

#lets try to add the HUC10 here to give a grouping for our betadisper function
#going to use plyr for this!
library(data.table)
setDT(removeSingle, keep.rownames = TRUE)[]
colnames(removeSingle)[1] <- "CodeDate"
dfSingle<-as.data.frame(removeSingle)

HUCMatchDup<-read.csv("Huc10Match.csv")
HUCMatch<-unique(HUCMatchDup)

library(plyr)

#so this is pretty much the same matrix except I am adding in the HUC as a grouping factor
combinedData <- join(dfSingle, HUCMatchDup, by='CodeDate', type='left', match='first') 
rareSpe<-combinedData[,colSums(combinedData>0)>10,drop=FALSE]
SCcomp = rareSpe[c(2:284)]
sort(colSums(SCcomp)) #there are a lot of species with less than 10 that we could remove
#it worked! for some reason i was having trouble because the setDT step back there converts away from data.frame so I had to convert the file back using as.data.frame

#trying again with relative abundances
#for most studies of beta diversity we want to get rid of the extremely rare "singleton" species, we will do this later...
SIdis = vegdist(SCcomp, method = "bray")
str(SIdis)
HUC10 = combinedData$HUC10

SCbeta = betadisper(SCdis, HUC10) #my error here was that I was using sites as a group even though there is only one site per group, now it should run because i am grouping by watershed
SCbeta

plot(SCbeta)

#distances to centroid
str(SCbeta)
boxplot(SCbeta$distances) # this is just overall beta and doesn't tell us a whole lot

#here we will group by watershed
boxplot(SCbeta$distances ~factor(HUC10)) #Awesome! we definitely see some range (and also some outliers that should probably be removed)
boxplot(SCbeta) #same thing

#we can now look at how this relates to environment using the "spreadEnv" we created earlier
spreadEnv

#extract the mean for each watershed
DtC = tapply(SCbeta$distances, HUC10, mean) #okay to do this we'd need our data to look much different so we will skip for now.

#something like boxplot(SIenv$Area.km.2) could be used to rarefy to area of each watershed to see how this influences diversity

#on to Module 3 stuff... we will go back to module 2 later

#we will still use SCcomp and spreadEnv
SCcomp
spreadEnv

## lets relativize the data and focus on using Bray-Curtis dissimilarities:
SCcomp.rel <- decostand(SCcomp, 'total')
head(SCcomp.rel)

#Calculate B-C
SCcomp.rel.bray<-vegdist(SCcomp.rel, method = 'bray')

#First analysis we can do is ANOSIM
#Going to use HUC10 (loaded earlier) as my treatment
HUC10

SCcomp.anosim1 <- anosim(SCcomp.rel.bray, HUC10, permutations=999)
summary(SCcomp.anosim1)

SCcomp.anosim2 <- anosim(SCcomp.rel, HUC10, permutations = 999, distance = "bray")
summary(SCcomp.anosim2)

#we can use the mrpp function which is more similar to nmds
SCcomp.rel.mrpp <- mrpp(SCcomp.rel, HUC10) # this uses Euclidean distances and is as such uncomparable with the anosim tests
SCcomp.rel.mrpp

SCcomp.rel.bray.mrpp1<-mrpp(SCcomp.rel.bray, HUC10) # this one makes a dissimilarity matrix
SCcomp.rel.bray.mrpp1 #significant
SCcomp.rel.bray.mrpp2<-mrpp(SCcomp.rel, HUC10, distance = 'bray') #this one makes a data.frame
SCcomp.rel.bray.mrpp2

#calculate a mean of within-group dissimilarities
SCcomp.rel.bray.md1 <- meandist(vegdist(SCcomp.rel, method = "bray"), HUC10)
SCcomp.rel.bray.md1

SCcomp.rel.bray.md2 <- meandist(SCcomp.rel.bray, HUC10)
SCcomp.rel.bray.md2

#Permutational multivariate analysis of variance
#adonis is used to perform PERMANOVA 

adonis(SCcomp.rel ~ HUC10, data = rareSpe, permutations = 999, method = "bray")
adonis2(SCcomp.rel ~ HUC10, data = rareSpe, permutations = 999, method = "bray")

SCcomp.disper<-betadisper(SCcomp.rel.bray, HUC10)
SCcomp.disper #Having trouble here (SCcomp.disper does not exist)
#Also just having general problems with my environmental matrix so should go back and compare to the files that marko has
anova(SCcomp.disper)
permutest(SCcomp.disper)
TukeyHSD(SCcomp.)

PoffNexus<-read.nexus("PoffNexus.nex", tree.names = NULL, force.multi = FALSE)
plot(PoffNexus)
#need to prune this tree to something that contains just the species of interest 
str(PoffNexus)
#could export the tip.label as a csv, break to just the last piece after the underscore, read back in as tip.labels?
Genus.labels<-PoffNexus$tip.label
Genus.labels
write.csv(Genus.labels, "genus.labels.csv")
genusOnly<-read.csv("genusOnly.csv")
genusOnly2 = as.vector(genusOnly)
PoffNexus$tip.label = genusOnly2
  str(PoffNexus)                                            #genusOnly<-gsub("..")
plot(PoffNexus)
tiplabels(PoffNexus)
#str.split
#https://www.r-bloggers.com/changing-phylogeny-tip-labels-in-r/ (this should help)
plot(PoffNexus)




