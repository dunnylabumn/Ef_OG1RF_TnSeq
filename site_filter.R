
args = commandArgs(trailingOnly=TRUE)
# test if there is at least three argument: if not, return an error
if (length(args)<=2) {
  stop("At least three argument must be supplied (TA sites) (input file FW) (input file RW) [output file]", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  args[4] = "outsites.tab"
}

ta_sites= data.frame(pos=read.table(paste(args[1],"/TA_locations.tab",sep=""))[,1],flag=c(1)) #Read in known TA sites in the genome
FW_Pos = read.table(args[2], col.names = c("chr", "pos", "weight")) #Read in Forward Set
RW_Pos = read.table(args[3], col.names = c("chr", "pos", "weight")) #Read in Reverse Set
RW_Pos[,'pos'] = RW_Pos[,'pos'] -1 # Correct for single base shift on reverse strand

TA_FW = merge(ta_sites, FW_Pos[,-1], by.x='pos', by.y='pos',all.y=TRUE) #Merge against TA sites (retain all observed integration sites)
TA_RW = merge(ta_sites, RW_Pos[,-1], by.x='pos', by.y='pos',all.y=TRUE) #Merge against TA sites (retain all observed integration sites)
TA_all = merge(TA_FW,TA_RW[,-2], by='pos', all.x=TRUE) #Merge resulting lists together

TA_all[,'weight.x'] = TA_all[,'weight.x'] + TA_all[,'weight.y'] #Add both orientations together.
TA_all = TA_all[,-4]

total_pos= nrow(TA_all) #Compute total sites with integrations
TA_count = sum(TA_all[which(!is.na(TA_all[,2])),3],na.rm =TRUE) #Compute total TA sites with integrations
off_count = sum(TA_all[which(is.na(TA_all[,2])),3],na.rm =TRUE) #Compute total off target integrations

TA_all = TA_all[which(TA_all[,3]>0 | !is.na(TA_all[,2] )),] #Drop all rows that have no integrations and are NOT a TA site.
names(TA_all) = c("Position", "TA_site", "Integrations") #Set column names
write.table(TA_all, args[4], quote=FALSE, sep="\t", row.names=FALSE) #Write out tabular file

annots = read.csv(paste(args[1],"/E_faecalis_Full_Annot.csv",sep=""))  #Load annotations

library(sqldf)

annot_TAs =sqldf("select a.*, b.* FROM annots as a LEFT OUTER JOIN TA_all as b ON b.Position <= a.maxPos and b.Position >= a.minPos") #Merge in annotations

write.table(annot_TAs, paste(args[4],".annot",sep=""), quote=FALSE, sep="\t", row.names=FALSE) #Write out annotated sites

