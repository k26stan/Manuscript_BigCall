## Load HaplotypeCaller MEM Data
TAB.l <- read.table( "Dropbox/Schork/JNJ11/Big_Call_Paper/DATA/Revisions/Updates/HC_MEM_USED.txt",sep=" ",header=F, fill=T )

## Remove Incomplete Jobs
TAB <- TAB.l[ which(TAB.l[,2]!=""), ]

## Pull out Relevant Columns
TAB <- TAB[,c(1,4,7)]

## Calculate Mean MEM and VMEM Values
MEM <- as.numeric(gsub( "kb","", as.character(TAB$V4) )) / 1e6
VMEM <- as.numeric(gsub( "kb","", as.character(TAB$V7) )) / 1e6

## Calculate Mean MEM & VMEM Values
MEM.MN <- mean( MEM )
VMEM.MN <- mean( VMEM )
 # ...and Stanard Deviation
MEM.SD <- sd( MEM )
VMEM.SD <- sd( VMEM )


d <- 2
## Print for Table
HC.MEM.print <- paste( round(MEM.MN,d), "newline (", round(MEM.SD,d), ")", sep="" )
HC.VMEM.print <- paste( round(VMEM.MN,d), "newline (", round(VMEM.SD,d), ")", sep="" )
