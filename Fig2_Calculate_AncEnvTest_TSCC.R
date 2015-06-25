## COMPARE NA12878 CALLS UNDER DIFFERENT ANCESTRAL ENVIRONMENTS ##
## Calculate Sens, Spec, PPV, Accuracy For Revisions ##
## June 25, 2015 ##
## Kristopher Standish ##

## Game Plan ##
# This Script
  # Load VCF files for Hi-Conf Locations
  # Load BED files for Hi-Conf Locations
  # Compile Contingency Tables
  # Save Tables
# Next Script: Calculate Sens, Spec, PPV, Accuracy
  # Plot Metrics vs % European (Fig 2)
  # Plot ANOVA & Concordance Metrics

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("SNP","FIL")
SI <- LINE[1]
UF <- LINE[2]

#############################################################
## LOAD DATA ################################################
#############################################################
print("########## SPECIFYING UNIVERSAL VARIABLES ##########")

############################################
## Set Variables/Paths #####################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Path to vcf files for NA12878 calls (TSCC)
if ( UF=="UNF" ) { PathToFiles <- "/projects/janssen/Ancestry/ANC_ENV_TEST/VCFs/3-MERGE/" }
if ( UF=="FIL" ) { PathToFiles <- "/projects/janssen/Ancestry/ANC_ENV_TEST/VCFs/F3-MERGE/" }
PathToHiC <- "/projects/janssen/Ancestry/ANC_ENV_TEST/VCFs/0-HIC/"
PathToChroms <- "/home/kstandis/HandyStuff/Chrom_Lengths.txt"
PathToSave <- paste("/projects/janssen/Ancestry/ANC_ENV_TEST/PLOTS/",DATE,"/",sep="")
dir.create( PathToSave )

## File Name Identifiers for Diff Chroms
IDS <- 0:3 # paste("GR",0:3,sep="_")
CHRS <- c(17,14,22,20)
CHR_TAGS <- paste( "Chr",CHRS,sep="_" )
NUM <- length(IDS) # Number of different chromosomes I looked at

############################################
## Chrom Lengths ###########################

## Load/Specify Length of Chromosomes
Chroms <- read.table( PathToChroms, sep="\t",header=T )
LENS <- Chroms[ CHRS, "LEN" ]
names(LENS) <- CHRS

############################################
## Hi-Confidence Sites #####################

## Load Hi-Confidence BED Files
HI_BED <- list()
for ( i in 1:NUM ) { id <- IDS[i] ; chr_tag <- CHR_TAGS[i]
	HI_BED[[chr_tag]] <- read.table( paste(PathToHiC,"Hi_",id,".list",sep=""), sep=":",header=F )
}

############################################
## Load Hi-Confidence Variants #############

R <- list()
start <- proc.time()
for ( i in 1:NUM ) {
	PathToR <- paste(PathToFiles,"REF_ENV_MRGD_NA12878_",IDS[i],".",SI,".vcf", sep="")
	print(PathToR)
	R[[i]] <- read.table( PathToR, header=F, sep="\t", colClass=c("factor", "numeric", "character", "factor", "factor", "numeric", "factor", "factor", rep("character",5))) 
	colnames(R[[i]]) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "H", "N", "M", "A", "E")
	print(paste( "Done with:",IDS[i],"-",round((proc.time()-start)[3],2) ))
}
names(R) <- CHR_TAGS

#############################################################
## ORGANIZE DATA ############################################
#############################################################

## Split Bed Coordinates
HI_BED.split <- lapply( HI_BED, function(x) t(sapply(strsplit( as.character(x[,2]), "-" ),"[",1:2)) )
HI_BED.diffs <- lapply( HI_BED.split, function(x) as.numeric(x[,2])-as.numeric(x[,1]) )
 # Calculate Number of Bases in Hi-Confidence Regions
HI_LENS <- unlist(lapply( HI_BED.diffs, sum ))

## FCT: Check if groups don't have values in both lines
GET_RPTS <- function(X, RPT_X) {
X_RPTS <- X[0,] # array(,dim=c(length(RPT_R),14))
No_Good <- c()
row <- 0
for (rptd in RPT_X) {
	row <- row + 1
	TEMP <- X[which(X$POS==rptd),]
	which_dot <- list()
	for (temp_row in 1:nrow(TEMP)) {
		which_dot[[temp_row]] <- which(TEMP[temp_row,10:14]!="./.")
	}
	if (length(Reduce(intersect, which_dot))==0) {
		NEW_LINE <- TEMP[1,]
		NEW_LINE$REF <- paste(TEMP$REF, collapse=":")
		NEW_LINE$ALT <- paste(TEMP$ALT, collapse=":")
		for (temp_row in 1:nrow(TEMP)) {
			NEW_LINE[9+which(TEMP[temp_row,10:14]!="./.")] <- TEMP[temp_row,9+which(TEMP[temp_row,10:14]!="./.")]
		}
		X_RPTS <- rbind(X_RPTS, NEW_LINE)
	}
	else { No_Good <- c(No_Good, rptd) }
} # Closes loop that goes through repeated positions
COMPILE <- list(X_RPTS, No_Good)
names(COMPILE) <- c("Table", "No_Good")
return(COMPILE)
} # Closes Function "GET_RPTS"

## Fix Positions Listed Twice in VCF files
RPT_R <- RPT_Ri <- list()
start <- proc.time()
for ( i in 1:NUM ) {
	RPT_R[[i]] <- R[[i]]$POS[which(duplicated(R[[i]]$POS))]
	RPT_Ri[[i]] <- which( R[[i]]$POS %in% RPT_R[[i]] ) # Some positions are listed 3 times
	print(paste( "Done with:",IDS[i],"-",round((proc.time()-start)[3],2) ))
}

## Sets in the columen "INFO" refer to first REF/ALT value ; Sets NOT mentioned in "INFO" are AFTER the ":" in REF/ALT
FLAT_R <- list()
start <- proc.time()
for ( i in 1:NUM ) {
	if (length(RPT_Ri[[i]])>0) {
		FLAT_R[[i]] <- GET_RPTS(R[[i]], RPT_R[[i]])
		# Take out Rows in RPT_Xi, put in X_RPTs, sort by position
		R[[i]] <- rbind(R[[i]][-RPT_Ri[[i]],],FLAT_R[[i]]$Table)
		R[[i]] <- R[[i]][order(R[[i]]$POS),]
		print(paste( "Done with:",IDS[i],"-",round((proc.time()-start)[3],2) ))
	}
}

# ## SANITY CHECK: Check to make sure all the variants are in the regions specified in the BED file
# for ( i in 1:NUM ) {
# 	pos <- as.numeric( R[[i]][,"POS"] )
# 	breaks <- as.numeric(c(t( HI_BED.split[[i]] )))
# 	tab <- cbind(pos, findInterval(pos, breaks) )
# 	print( length(which(tab[,2]%%2==0)) )
# 	print( all( tab[ which(tab[,2]%%2==0), 1 ] %in% R[[i]][,"POS"] ) )
# }

## Pull out genotypes at each position (using R tables)
Genos <- list()
start <- proc.time()
for ( i in 1:NUM ) {
	Genos[[i]] <- list()	

	Genos[[i]]$H <- sapply( strsplit(as.character(R[[i]]$H), ":"), "[", 1)
	Genos[[i]]$H <- gsub("|", "/", Genos[[i]]$H, fixed=T)
	temp <- t(sapply( strsplit(Genos[[i]]$H,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[i]]$H[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[i]]$M <- sapply( strsplit(as.character(R[[i]]$M), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[i]]$M,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[i]]$M[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[i]]$N <- sapply( strsplit(as.character(R[[i]]$N), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[i]]$N,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[i]]$N[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[i]]$A <- sapply( strsplit(as.character(R[[i]]$A), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[i]]$A,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[i]]$A[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[i]]$E <- sapply( strsplit(as.character(R[[i]]$E), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[i]]$E,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[i]]$E[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	print(paste( "Done with:",CHR_TAGS[i],"-",round((proc.time()-start)[3],2) ))
}
names(Genos) <- CHR_TAGS

## Make array with easy-access genotypes
C_Rr <- SNP <- SNP_1 <- list()
start <- proc.time()
for ( i in 1:NUM ) {
	C_Rr[[i]] <- data.frame(R[[i]][,2:5], H=Genos[[i]]$H, M=Genos[[i]]$M, N=Genos[[i]]$N, A=Genos[[i]]$A, E=Genos[[i]]$E)	

	# Make table of alt allele count (SNP Array) for each position and group
	temp <- as.matrix(C_Rr[[i]][,5:9])
	SNP[[i]] <- array( 2, dim(C_Rr[[i]][,5:9]) )
	SNP[[i]][grep("./.",temp,fixed=T)] <- 0
	SNP[[i]][grep("0/",temp,)] <- 1
	SNP[[i]][grep("0/0",temp,)] <- 0
	SNP[[i]] <- data.frame(SNP[[i]])

	# Make table of variant positions
	SNP_1[[i]] <- data.frame(ceiling(SNP[[i]]/2))
	colnames(SNP[[i]]) <- colnames(SNP_1[[i]]) <- c("HI", "MIX", "NA", "ASN", "EUR")

	print(paste( "Done with:",IDS[i],"-",round((proc.time()-start)[3],2) ))
}
names(SNP) <- names(SNP_1) <- CHR_TAGS

#############################################################
## CALCULATE CONTINGENCY TABLE FOR EACH SET #################
#############################################################
## Compile, for each chromosome & group:
 # Total Positions
 # Actual Variants (Condition Positive)
 # Actual Non-Variants (Condition Negative)
 # Called Variants (Test Positive)
 # Called Non-Variants (Test Negative)
 # True Positives
 # False Negatives
 # False Positives
 # False Negatives
## Calculate, for each chromosome & group:
 # Sensitivity
 # Specificity
 # Accuracy
 # Positive Predictive Value

Total_Pos <- HI_LENS
Cond_Pos <- unlist(lapply( SNP, function(y) length(which(y[,1]!=0)) ))
Cond_Neg <- HI_LENS - Cond_Pos # unlist(lapply( SNP, function(y) length(which(y[,1]==0)) ))
Test_Pos <- lapply( SNP, function(y) apply(y[,2:5],2,function(x) length(which(x!=0)) ) )
Test_Neg <- # lapply( SNP, function(y) apply(y[,2:5],2,function(x) length(which(x==0)) ) )


TEST <- lapply( SNP, function(y) apply(y[,2:5],2,function(x) which(x!=y[,1]) ) )




#############################################################
## DETERMINE CONCORDANCE AMONGST GROUPS #####################
#############################################################

To_Ignore <- Groups <- C_Gr <- Conc_T <- P_Conc_T <- Conc_3 <- P_Conc_3 <- 
Conc_2 <- P_Conc_2 <- Conc_1 <- P_Conc_1 <- list()
Conc_Compile <- array(, dim=c(NUM,12))
start <- proc.time()
for ( chrom in 1:NUM ) {
	## Look only at the groups
	 # (i.e., remove positions where only Hi-Conf called a variant there)
	To_Ignore[[chrom]] <- which( SNP_1[[chrom]]$H==1 & apply(SNP_1[[chrom]][2:5],1,sum)==0 )
	Groups[[chrom]] <- list(H=Genos[[chrom]]$H[-To_Ignore[[chrom]]], M=Genos[[chrom]]$M[-To_Ignore[[chrom]]], N=Genos[[chrom]]$N[-To_Ignore[[chrom]]], A=Genos[[chrom]]$A[-To_Ignore[[chrom]]], E=Genos[[chrom]]$E[-To_Ignore[[chrom]]] )
	C_Gr[[chrom]] <- C_Vr[[chrom]][-To_Ignore[[chrom]],]

	# Total Concordance for Variant Positions (ignore Hi-Conf status)
	Conc_T[[chrom]] <- which( Groups[[chrom]]$M==Groups[[chrom]]$N & Groups[[chrom]]$M==Groups[[chrom]]$A & Groups[[chrom]]$M==Groups[[chrom]]$E )
	P_Conc_T[[chrom]] <- length(Conc_T[[chrom]]) / length(Groups[[chrom]]$M) # 85% of variants are called the same in all groups
	names(P_Conc_T[[chrom]]) <- "All"

	# 3-way Concordance for Variant Positions (ignore Hi-Conf status)
	Conc_3[[chrom]] <- list()
	Conc_3[[chrom]]$NoM <- which( Groups[[chrom]]$N==Groups[[chrom]]$A & Groups[[chrom]]$N==Groups[[chrom]]$E & Groups[[chrom]]$N!=Groups[[chrom]]$M )
	Conc_3[[chrom]]$NoN <- which( Groups[[chrom]]$M==Groups[[chrom]]$A & Groups[[chrom]]$M==Groups[[chrom]]$E & Groups[[chrom]]$M!=Groups[[chrom]]$N )
	Conc_3[[chrom]]$NoA <- which( Groups[[chrom]]$N==Groups[[chrom]]$M & Groups[[chrom]]$N==Groups[[chrom]]$E & Groups[[chrom]]$N!=Groups[[chrom]]$A )
	Conc_3[[chrom]]$NoE <- which( Groups[[chrom]]$N==Groups[[chrom]]$A & Groups[[chrom]]$N==Groups[[chrom]]$M & Groups[[chrom]]$N!=Groups[[chrom]]$E )
	P_Conc_3[[chrom]] <- c( length(Conc_3[[chrom]]$NoM), length(Conc_3[[chrom]]$NoN), length(Conc_3[[chrom]]$NoA), length(Conc_3[[chrom]]$NoE) ) / length(Groups[[chrom]]$M)
	names(P_Conc_3[[chrom]]) <- c("NoM", "NoN", "NoA", "NoE")
	Reduce(intersect, Conc_3[[chrom]])

	# Pairwise Concordance for Variant Positions (ignore Hi-Conf status)
	Conc_2[[chrom]] <- list()
	Conc_2[[chrom]]$M_N <- which( Groups[[chrom]]$M==Groups[[chrom]]$N & Groups[[chrom]]$M!=Groups[[chrom]]$A & Groups[[chrom]]$M!=Groups[[chrom]]$E )
	Conc_2[[chrom]]$M_A <- which( Groups[[chrom]]$M==Groups[[chrom]]$A & Groups[[chrom]]$M!=Groups[[chrom]]$N & Groups[[chrom]]$M!=Groups[[chrom]]$E )
	Conc_2[[chrom]]$M_E <- which( Groups[[chrom]]$M==Groups[[chrom]]$E & Groups[[chrom]]$M!=Groups[[chrom]]$A & Groups[[chrom]]$M!=Groups[[chrom]]$N )
	Conc_2[[chrom]]$N_A <- which( Groups[[chrom]]$N==Groups[[chrom]]$A & Groups[[chrom]]$N!=Groups[[chrom]]$M & Groups[[chrom]]$N!=Groups[[chrom]]$E )
	Conc_2[[chrom]]$N_E <- which( Groups[[chrom]]$N==Groups[[chrom]]$E & Groups[[chrom]]$N!=Groups[[chrom]]$M & Groups[[chrom]]$N!=Groups[[chrom]]$A )
	Conc_2[[chrom]]$A_E <- which( Groups[[chrom]]$A==Groups[[chrom]]$E & Groups[[chrom]]$A!=Groups[[chrom]]$M & Groups[[chrom]]$A!=Groups[[chrom]]$N )
	P_Conc_2[[chrom]] <- c( length(Conc_2[[chrom]]$M_N), length(Conc_2[[chrom]]$M_A), length(Conc_2[[chrom]]$M_E), length(Conc_2[[chrom]]$N_A), length(Conc_2[[chrom]]$N_E), length(Conc_2[[chrom]]$A_E) ) / length(Groups[[chrom]]$M)
	names(P_Conc_2[[chrom]]) <- c("MN", "MA", "ME", "NA", "NE", "AE")

	# Unique to a Single Group
	Conc_1[[chrom]] <- which( Groups[[chrom]]$M!=Groups[[chrom]]$N & Groups[[chrom]]$M!=Groups[[chrom]]$A & Groups[[chrom]]$M!=Groups[[chrom]]$E & Groups[[chrom]]$N!=Groups[[chrom]]$A & Groups[[chrom]]$N!=Groups[[chrom]]$E & Groups[[chrom]]$A!=Groups[[chrom]]$E )
	P_Conc_1[[chrom]] <- length(Conc_1[[chrom]]) / length(Groups[[chrom]]$M)
	names(P_Conc_1[[chrom]]) <- "UNIQ"

	# Compile into Array
	Conc_Compile[chrom,] <- c( P_Conc_T[[chrom]], P_Conc_3[[chrom]], P_Conc_2[[chrom]], P_Conc_1[[chrom]] )

	print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
}
rownames(Conc_Compile) <- IDS
colnames(Conc_Compile) <- c( names(P_Conc_T[[chrom]]), names(P_Conc_3[[chrom]]), names(P_Conc_2[[chrom]]), names(P_Conc_1[[chrom]]) )


# Barplot of how many variants fall into each category
COLS <- c("mediumpurple3", "dodgerblue4", "dodgerblue3", "dodgerblue2", "dodgerblue", colorRampPalette(c("cadetblue4","cadetblue1","white"))(8)[1:6], "chartreuse")
jpeg( paste(PathToSave,SI,"_",UF,"_1-Concord_Grps.jpeg",sep=""), height=2000, width=2000, pointsize=20)
par(mfrow=c(2,2))
for ( chrom in 1:NUM ) {
	mids <- barplot(Conc_Compile[chrom,], col=COLS, main=paste("Concordance of Calls from Groups -",IDS[chrom]), ylim=c(0,1), xlab="Concordance", ylab="Fraction of Variant Positions")
	text(  mids, .2, labels=round(Conc_Compile[chrom,],3), srt=90 )
	arrows( mids[c(2,6),1], .4, mids[c(5,11),1], .4, length=0, col=c("dodgerblue", "cadetblue"), lwd=5 )
	text( mean(mids[c(2,5),1]), c(.5,.45), labels=c("Called in 3 Groups",round(sum(Conc_Compile[chrom,2:5]),4)))
	text( mean(mids[c(6,11),1]), c(.5,.45), labels=c("Called in 2 Groups",round(sum(Conc_Compile[chrom,6:11]),4)))
}
dev.off()

# Barplot summarizing the 4 groups from the previous plot
COLS <- c("mediumpurple3", "dodgerblue4", "dodgerblue3", "dodgerblue2", "dodgerblue", colorRampPalette(c("cadetblue4","cadetblue1","white"))(8)[1:6], "chartreuse")
CONC_MNS <- apply(Conc_Compile,2,mean)
CONC_SEM <- apply(Conc_Compile,2,sd)/2
jpeg( paste(PathToSave,SI,"_",UF,"_1-Mn_Concord_Grps.jpeg",sep=""), height=1000, width=1000, pointsize=20)
mids <- barplot(CONC_MNS, col=COLS, main=paste("Concordance of Calls from Groups -",IDS[chrom]), ylim=c(0,1), xlab="Concordance", ylab="Fraction of Variant Positions")
text(  mids, .2, labels=round(CONC_MNS,3), srt=90 )
arrows( mids[,1], CONC_MNS-CONC_SEM, mids[,1], CONC_MNS+CONC_SEM, length=0, col="black", lwd=2 )
arrows( mids[c(2,6),1], .4, mids[c(5,11),1], .4, length=0, col=c("dodgerblue", "cadetblue"), lwd=5 )
text( mean(mids[c(2,5),1]), c(.5,.45), labels=c("Called in 3 Groups",round(sum(CONC_MNS[2:5]),4)))
text( mean(mids[c(6,11),1]), c(.5,.45), labels=c("Called in 2 Groups",round(sum(CONC_MNS[6:11]),4)))
dev.off()

# Of Variants omitted by one group, what any group consistently higher?
NORM_3 <- prop.table(Conc_Compile[,2:(1+NUM)],1)
MEAN_N3 <- apply(NORM_3, 2, mean)
SD_N3 <- apply(NORM_3, 2, sd) / sqrt(NUM)
NORM_3b <- data.frame( PERC=c(NORM_3),GRP=c(rep("M",4),rep("N",4),rep("A",4),rep("E",4)) )
MOD <- aov(PERC ~ GRP, data=NORM_3b)
P_VAL <- summary(MOD)[[1]]$`Pr(>F)`[1] # ; P_VAL <- drop1(MOD, ~., test="F")

# Barplot of how many variants are omitted by single group only
COLS <- c("dodgerblue4", "dodgerblue3", "dodgerblue2", "dodgerblue1")
jpeg( paste(PathToSave,SI,"_",UF,"_1-Concord_3.jpeg",sep=""), height=1000, width=1000, pointsize=20)
mids <- barplot(MEAN_N3, col=COLS, main="Variants Called by 3 Groups Only", ylim=c(0,1), xlab="Missing Group", ylab="Fraction of Vars Missed by Only 1 Group", names.arg=c("M","N","A","E"))
arrows( mids, MEAN_N3-SD_N3, mids, MEAN_N3+SD_N3, length=.2, angle=90, col="black", lwd=3, code=3 )
arrows( mids[1,1], max(MEAN_N3)+.1, mids[NUM,1], max(MEAN_N3)+.1, length=0, col="dodgerblue", lwd=5 )
text(  mean(mids[,1]), max(MEAN_N3)+.15, labels="*" )
text(  mean(mids[,1]), max(MEAN_N3)+.2, labels=paste("ANOVA: p=",round(P_VAL,4),sep="") )
dev.off()

#############################################################
## DETERMINE NUMBER OF MISSES FOR EACH GROUP ################
#############################################################

# Pull out HI-CONF variant positions only
Hi_Only <- G_Hi <- C_Hi <- list()
for ( chrom in 1:NUM ) {
	Hi_Only[[chrom]] <- which(Genos[[chrom]]$H!="./.")

	G_Hi[[chrom]] <- list()
	G_Hi[[chrom]]$H <- Genos[[chrom]]$H[Hi_Only[[chrom]]]
	G_Hi[[chrom]]$M <- Genos[[chrom]]$M[Hi_Only[[chrom]]]
	G_Hi[[chrom]]$N <- Genos[[chrom]]$N[Hi_Only[[chrom]]]
	G_Hi[[chrom]]$A <- Genos[[chrom]]$A[Hi_Only[[chrom]]]
	G_Hi[[chrom]]$E <- Genos[[chrom]]$E[Hi_Only[[chrom]]]

	C_Hi[[chrom]] <- data.frame(V[[chrom]][Hi_Only[[chrom]],2:5], H=Genos[[chrom]]$H[Hi_Only[[chrom]]], M=Genos[[chrom]]$M[Hi_Only[[chrom]]], N=Genos[[chrom]]$N[Hi_Only[[chrom]]], A=Genos[[chrom]]$A[Hi_Only[[chrom]]], E=Genos[[chrom]]$E[Hi_Only[[chrom]]])
}
names(Hi_Only) <- names(G_Hi) <- names(C_Hi) <- IDS

# Calculate the number of MISSES for each Group
HET_HI <- HOM_HI <- MISSES <- list()
for ( chrom in 1:NUM ) {
	HET_HI[[chrom]] <- grep("0/", G_Hi[[chrom]][[1]])
	HOM_HI[[chrom]] <- grep("0/", G_Hi[[chrom]][[1]], invert=T)
	MISSES[[chrom]] <- array( , dim=c(7,5))
	colnames(MISSES[[chrom]]) <- c("HIC", "MIX", "NA", "ASN", "EUR")
	rownames(MISSES[[chrom]]) <- c("Same", "Wrong", "HOM_Under", "HOM_Miss", "HET_Over", "HET_Miss", "Other")
	for ( i in 1:(ncol(C_Hi[[chrom]])-4) ) {
		MISSES[[chrom]][1,i] <- length(which( G_Hi[[chrom]][[1]]==G_Hi[[chrom]][[i]] ))
		MISSES[[chrom]][2,i] <- length(which( G_Hi[[chrom]][[1]]!=G_Hi[[chrom]][[i]] ))

		HOM_CL <- setdiff( 1:length(Hi_Only[[chrom]]), union( grep("0/",G_Hi[[chrom]][[i]]), grep("./.",G_Hi[[chrom]][[i]],fixed=T) ) )
		HET_CL <- setdiff( grep("0/",G_Hi[[chrom]][[i]]), grep("0/0",G_Hi[[chrom]][[i]]) )
		REF_CL <- union( grep("./.",G_Hi[[chrom]][[i]],fixed=T), grep("0/0",G_Hi[[chrom]][[i]]) )

		MISSES[[chrom]][3,i] <- length(intersect( HOM_HI[[chrom]], HET_CL ))
		MISSES[[chrom]][4,i] <- length(intersect( HOM_HI[[chrom]], REF_CL ))
		MISSES[[chrom]][5,i] <- length(intersect( HET_HI[[chrom]], HOM_CL ))
		MISSES[[chrom]][6,i] <- length(intersect( HET_HI[[chrom]], REF_CL ))
		MISSES[[chrom]][7,i] <- MISSES[[chrom]][2,i] - sum(MISSES[[chrom]][3:6,i])
	}
}
names(HOM_HI) <- names(HET_HI) <- names(MISSES) <- IDS

# Barplot of how many misses in each group (by type)
MAXES <- numeric(NUM) ; for (chrom in 1:NUM) { MAXES[chrom] <- max(MISSES[[chrom]][2,2:5])/nrow(C_Hi[[chrom]]) }
YMAX <- max(MAXES)
COLS <- c("darkorchid2", "chocolate1", "chartreuse2", "darkgoldenrod2", "deepskyblue2") # "brown2"
LEGEND <- c("HOM_Und", "HOM_Miss", "HET_Over", "HET_Miss", "Other")
jpeg( paste(PathToSave,SI,"_",UF,"_2-Misses.jpeg",sep=""), height=1000, width=1000, pointsize=20)
par(mfrow=c(2,2))
for ( chrom in 1:NUM ) {
	mids <- barplot(MISSES[[chrom]][3:7,2:5]/nrow(C_Hi[[chrom]]), ylim=c(0,YMAX), col=COLS, main=paste("Misses by Type and Group -",IDS[[chrom]]), xlab="Groups", ylab="Misses")
	legend(0.3*mids[1], YMAX, fill=COLS, legend=LEGEND, title="Type of Error" )
}
dev.off()

# Barplot of how many misses in each group (total)
CHR_HIC <- numeric(4) ; names(CHR_HIC) <- IDS
COLS <- c("darkorchid2")
LEGEND <- c("Total Misses")
jpeg( paste(PathToSave,SI,"_",UF,"_2-Misses_Total.jpeg",sep=""), height=1000, width=1000, pointsize=20)
par(mfrow=c(2,2))
for ( chrom in 1:NUM ) {
	CHR_HIC[chrom] <- nrow(C_Hi[[chrom]])
	mids <- barplot(MISSES[[chrom]][2,2:5]/nrow(C_Hi[[chrom]]), ylim=c(0,YMAX), col=COLS, main=paste("Total Misses by Group -",IDS[[chrom]]), xlab="Groups", ylab="Misses")
	legend(0.3*mids[1], YMAX, fill=COLS, legend=LEGEND, title="Type of Error" )
}
dev.off()

# Run 2-way ANOVA on Chrom and Group
NORM_MS <- MISSES[[1]]["Wrong",2:(NUM+1)]
for ( chrom in 2:NUM ) { NORM_MS <- cbind( NORM_MS, MISSES[[chrom]]["Wrong",2:(NUM+1)] ) }
colnames(NORM_MS) <- IDS
NORM_MSp <- t(t(NORM_MS) / CHR_HIC )
NORM_MSb <- data.frame( PERC=c(t(NORM_MSp)),GRP=c(rep("M",NUM),rep("N",NUM),rep("A",NUM),rep("E",NUM)), CHR=rep(c("0","1","2","3"),NUM) )
MOD <- lm( PERC ~ GRP + as.factor(CHR), data=NORM_MSb )
anova(MOD)
P_VAL2 <- anova(MOD)$`Pr(>F)`[1]

# Normalize Total Misses to Mean of Chrom & Compare
MEAN_MS <- apply(NORM_MS,2,mean)
SD_MS <- apply(NORM_MS,2,sd) / sqrt(NUM)
NORM_MSc <- t(t(NORM_MS) / MEAN_MS )
NORM_MSd <- data.frame( PERC=c(t(NORM_MSc)),GRP=c(rep("M",4),rep("N",4),rep("A",4),rep("E",4)), CHR=rep(c("0","1","2","3"),4) )
MOD <- aov(PERC ~ GRP, data=NORM_MSd)
MEAN_MS <- apply(NORM_MSc,1,mean)
SD_MS <- apply(NORM_MSc,1,sd) / sqrt(NUM)
P_VAL <- summary(MOD)[[1]]$`Pr(>F)`[1] # ; P_VAL <- drop1(MOD, ~., test="F")

# Barplot of how many variants are MISSED by each group
Y_MAX <- round(max(NORM_MSc)+.3,1)
COLS <- c("darkorchid2")
jpeg( paste(PathToSave,SI,"_",UF,"_2-Total_Misses.jpeg",sep=""), height=1000, width=1000, pointsize=20)
mids <- barplot(MEAN_MS, col=COLS, main="Proportion of Misses by Group Relative to Mean for Chrom", ylim=c(0,Y_MAX), xlab="Group", ylab="Proportion of Misses Relative to Mean", names.arg=c("M","N","A","E"))
arrows( mids, MEAN_MS-SD_MS, mids, MEAN_MS+SD_MS, length=.2, angle=90, col="black", lwd=3, code=3 )
arrows( mids[1,1], max(MEAN_MS)+.1, mids[NUM,1], max(MEAN_MS)+.1, length=0, col="darkorchid4", lwd=5 )
if (P_VAL<.05 | P_VAL2<.05) { text(  mean(mids[,1]), max(MEAN_MS)+.15, labels="*" ) }
text(  mean(mids[,1]), max(MEAN_MS)+.19, labels=paste("Norm ANOVA: p=",round(P_VAL,4),sep="") )
text(  mean(mids[,1]), max(MEAN_MS)+.23, labels=paste("2-way ANOVA: p=",round(P_VAL2,4),sep="") )
dev.off()

save(MISSES, file=paste(PathToSave,"Miss_Counts",SI,"_",UF,".Rdata",sep="")) 

##############################################################################################
##############################################################################################
##############################################################################################

#############################################################
## LOAD MERGED/VARIANT ONLY DATA ############################
#############################################################

###############################################
## Load MERGED vcfs with ONLY VARIANTS ########
R <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	PathToR <- paste(PathToFiles,"REF_ENV_MRGD_NA12878_",IDS[chrom],".",SI,".vcf", sep="")
	print(PathToR)
	R[[chrom]] <- read.table( PathToR, header=F, sep="\t", colClass=c("factor", "numeric", "character", "factor", "factor", "numeric", "factor", "factor", rep("character",5))) 
	colnames(R[[chrom]]) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "H", "N", "M", "A", "E")
	print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
}
names(R) <- IDS

###############################################
## Fix Positions Listed Twice #################
RPT_R <- RPT_Ri <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	RPT_R[[chrom]] <- R[[chrom]]$POS[which(duplicated(R[[chrom]]$POS))]
	RPT_Ri[[chrom]] <- which( R[[chrom]]$POS %in% RPT_R[[chrom]] ) # Some positions are listed 3 times
	print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
}

# Sets in the columen "INFO" refer to first REF/ALT value ; Sets NOT mentioned in "INFO" are AFTER the ":" in REF/ALT
FLAT_R <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	if (length(RPT_Ri[[chrom]])>0) {
		FLAT_R[[chrom]] <- GET_RPTS(R[[chrom]], RPT_R[[chrom]])
		# Take out Rows in RPT_Xi, put in X_RPTs, sort by position
		R[[chrom]] <- rbind(R[[chrom]][-RPT_Ri[[chrom]],],FLAT_R[[chrom]]$Table)
		R[[chrom]] <- R[[chrom]][order(R[[chrom]]$POS),]
		print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
	}
}

###############################################
## Organize into More Convenient Formats ######

# Pull out genotypes at each position (using R tables)
Genos <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	Genos[[chrom]] <- list()	

	Genos[[chrom]]$H <- sapply( strsplit(as.character(R[[chrom]]$H), ":"), "[", 1)
	Genos[[chrom]]$H <- gsub("|", "/", Genos[[chrom]]$H, fixed=T)
	temp <- t(sapply( strsplit(Genos[[chrom]]$H,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$H[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[chrom]]$M <- sapply( strsplit(as.character(R[[chrom]]$M), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[chrom]]$M,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$M[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[chrom]]$N <- sapply( strsplit(as.character(R[[chrom]]$N), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[chrom]]$N,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$N[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[chrom]]$A <- sapply( strsplit(as.character(R[[chrom]]$A), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[chrom]]$A,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$A[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[chrom]]$E <- sapply( strsplit(as.character(R[[chrom]]$E), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[chrom]]$E,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$E[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
}
names(Genos) <- IDS

# Make array with easy-access genotypes
C_Rr <- SNP <- SNP_1 <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	C_Rr[[chrom]] <- data.frame(R[[chrom]][,2:5], H=Genos[[chrom]]$H, M=Genos[[chrom]]$M, N=Genos[[chrom]]$N, A=Genos[[chrom]]$A, E=Genos[[chrom]]$E)	

	# Make table of alt allele count (SNP Array) for each position and group
	temp <- as.matrix(C_Rr[[chrom]][,5:9])
	SNP[[chrom]] <- array( 2, dim(C_Rr[[chrom]][,5:9]) )
	SNP[[chrom]][grep("./.",temp,fixed=T)] <- 0
	SNP[[chrom]][grep("0/",temp,)] <- 1
	SNP[[chrom]][grep("0/0",temp,)] <- 0
	SNP[[chrom]] <- data.frame(SNP[[chrom]])

	# Make table of variant positions
	SNP_1[[chrom]] <- data.frame(ceiling(SNP[[chrom]]/2))
	colnames(SNP[[chrom]]) <- colnames(SNP_1[[chrom]]) <- c("HI", "MIX", "NA", "ASN", "EUR")

	print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
}
names(SNP) <- names(SNP_1) <- IDS

#############################################################
## CALCULATE FALSE POSITIVES ################################
#############################################################

## Pull out positions where HIC Set is Reference (i.e., "./.")
H_REF <- list()
for ( chrom in 1:NUM ) {
	H_REF[[chrom]] <- which(Genos[[chrom]]$H=="./.")
}
names(H_REF) <- IDS

## Make Arrays for each CHR of only HIC Ref positions (e.g., potential false positives)
FP_CND <- list()
FP_CND_SNP <- list()
for ( chrom in 1:NUM ) {
	FP_CND[[chrom]] <- data.frame(M=Genos[[chrom]]$M[H_REF[[chrom]]],N=Genos[[chrom]]$N[H_REF[[chrom]]],A=Genos[[chrom]]$A[H_REF[[chrom]]],E=Genos[[chrom]]$E[H_REF[[chrom]]])
	FP_CND_SNP[[chrom]] <- SNP[[chrom]][H_REF[[chrom]],]
}
names(FP_CND) <- names(FP_CND_SNP) <- IDS

## Get Counts of how many False Positives there are (i.e., non- "0/0" or "./." in the dataframes)
FP_ARR <- list()
for ( chrom in 1:NUM ) {
	FP_ARR[[chrom]] <- array(,dim=c(4,4))
	colnames(FP_ARR[[chrom]]) <- c("MIX","NA","ASN","EUR")
	rownames(FP_ARR[[chrom]]) <- c("Ref","Var","Het","Hom")
	for ( grp in 1:4 ) {
		nREF <- length(which(FP_CND_SNP[[chrom]][,grp+1]==0))
		nVAR <- length(which(FP_CND_SNP[[chrom]][,grp+1]!=0))
		nHET <- length(which(FP_CND_SNP[[chrom]][,grp+1]==1))
		nHOM <- length(which(FP_CND_SNP[[chrom]][,grp+1]==2))
		FP_ARR[[chrom]][,grp] <- c(nREF,nVAR,nHET,nHOM)
	}
}
names(FP_ARR) <- IDS

#############################################################
## PLOT FALSE POSITIVES #####################################
#############################################################

# Barplot of how many false positives in each group (by type)
COLS <- c("chocolate1", "brown2") # "brown2"
LEGEND <- c("HET","HOM")
YMAX <- Reduce(max,FP_ARR)
jpeg( paste(PathToSave,SI,"_",UF,"_3-FPs_By_Type.jpeg",sep=""), height=1000, width=1000, pointsize=20)
par(mfrow=c(2,2))
for ( chrom in 1:NUM ) {
	mids <- barplot(FP_ARR[[chrom]][3:4,], ylim=c(0,YMAX), col=COLS, main=paste("False Positives by Type and Group -",IDS[[chrom]]), xlab="Groups", ylab="FPs")
	legend(0.3*mids[1], YMAX, fill=COLS, legend=LEGEND, title="Type of FP" )
}
dev.off()

# Barplot of how many FPs in each group (total)
COLS <- c("darkorchid2")
YMAX <- Reduce(max,FP_ARR)
jpeg( paste(PathToSave,SI,"_",UF,"_3-FPs_Total.jpeg",sep=""), height=1000, width=1000, pointsize=20)
par(mfrow=c(2,2))
for ( chrom in 1:NUM ) {
	mids <- barplot(FP_ARR[[chrom]][2,], ylim=c(0,YMAX), col=COLS, main=paste("Total False Positives by Group -",IDS[[chrom]]), xlab="Groups", ylab="FPs")
}
dev.off()

# Run 2-way ANOVA on Chrom and Group
NORM_FP <- FP_ARR[[1]]["Var",1:NUM]
for ( chrom in 2:NUM ) { NORM_FP <- cbind( NORM_FP, FP_ARR[[chrom]]["Var",1:NUM] ) }
colnames(NORM_FP) <- IDS
#NORM_FPp <- t(t(NORM_FP) / CHR_HIC )
NORM_FPb <- data.frame( PERC=c(t(NORM_FP)),GRP=c(rep("M",NUM),rep("N",NUM),rep("A",NUM),rep("E",NUM)), CHR=rep(c("0","1","2","3"),NUM) )
MOD <- lm( PERC ~ GRP + as.factor(CHR), data=NORM_FPb )
anova(MOD)
P_VAL2 <- anova(MOD)$`Pr(>F)`[1]

# Normalize Total Misses to Mean of Chrom & Compare
MEAN_FP <- apply(NORM_FP,2,mean)
SD_FP <- apply(NORM_FP,2,sd) / sqrt(NUM)
NORM_FPc <- t(t(NORM_FP) / MEAN_FP )
NORM_FPd <- data.frame( PERC=c(t(NORM_FPc)),GRP=c(rep("M",4),rep("N",4),rep("A",4),rep("E",4)), CHR=rep(c("0","1","2","3"),4) )
MOD <- aov(PERC ~ GRP, data=NORM_FPd)
MEAN_FP <- apply(NORM_FPc,1,mean)
SD_FP <- apply(NORM_FPc,1,sd) / sqrt(NUM)
P_VAL <- summary(MOD)[[1]]$`Pr(>F)`[1] # ; P_VAL <- drop1(MOD, ~., test="F")

# Barplot of how many variants are MISSED by each group
COLS <- c("darkorchid2")
jpeg( paste(PathToSave,SI,"_",UF,"_3-Total_FPs.jpeg",sep=""), height=1000, width=1000, pointsize=20)
mids <- barplot(MEAN_FP, col=COLS, main="Proportion of False Positives by Group Relative to Mean for Chrom", ylim=c(0,Y_MAX), xlab="Group", ylab="Proportion of False Positives Relative to Chromosome Mean", names.arg=c("M","N","A","E"))
arrows( mids, MEAN_FP-SD_FP, mids, MEAN_FP+SD_FP, length=.2, angle=90, col="black", lwd=3, code=3 )
arrows( mids[1,1], max(MEAN_FP)+.1, mids[NUM,1], max(MEAN_FP)+.1, length=0, col="darkorchid4", lwd=5 )
if (P_VAL<.05 | P_VAL2<.05) { text(  mean(mids[,1]), max(MEAN_FP)+.15, labels="*" ) }
text(  mean(mids[,1]), max(MEAN_FP)+.19, labels=paste("Norm ANOVA: p=",round(P_VAL,4),sep="") )
text(  mean(mids[,1]), max(MEAN_FP)+.23, labels=paste("2-way ANOVA: p=",round(P_VAL2,4),sep="") )
dev.off()


save(FP_ARR, file=paste(PathToSave,"FP_Counts",SI,"_",UF,".Rdata",sep="")) 

#####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####
	#####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####
#####   #####   #####   #####   #####  HEREIAM  #####   #####   #####   #####   #####   #####
	#####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####
#####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####
	#####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####   #####

#############################################################
## PLOT MISSES/FPs vs % EUROPEAN ############################
#############################################################

# Load Table w/ Total % Admixture of each Group
ADM <- read.table("/projects/janssen/Ancestry/ANC_ENV_TEST/Group_Admix_Perc.txt",sep="\t",header=T)
 # ADM goes E0-3, AN0-3, MIX0-3, ASN0-3
ADM_MS <- NORM_MSd[c(13:16,5:8,1:4,9:12),1]
ADM_FP <- NORM_FPd[c(13:16,5:8,1:4,9:12),1]
ADM_DAT <- data.frame(PEUR=ADM$europe/10,MS=ADM_MS,FP=ADM_FP)
rownames(ADM_DAT) <- rownames(ADM)
write.table(ADM_DAT,paste(PathToSave,SI,"_",UF,"_4-Perc_EUR.txt",sep=""), sep="\t",row.names=T,col.names=T,quote=F)

Y_LIM <- c(0.85,1.15)
jpeg( paste(PathToSave,SI,"_",UF,"_4-Perc_EUR.jpeg",sep=""), height=1000, width=2000, pointsize=24)
par(mfrow=c(1,2))
plot( ADM_DAT$MS ~ ADM_DAT$PEUR, main="Normalized Misses vs Percent European", xlab="% European in Group", ylab="Normalized Misses", pch=20, col="tomato",xlim=c(0,1),ylim=Y_LIM)
MOD <- lm(ADM_DAT$MS ~ ADM_DAT$PEUR)
abline(MOD, col="tomato", lwd=3)
text(.2,min(Y_LIM)+.06, labels=paste("R2=",round(summary(MOD)$r.squared,4)))
text(.2,min(Y_LIM)+.03, labels=paste("p=",round(summary(MOD)$coefficients[2,4],4)))
plot( ADM_DAT$FP ~ ADM_DAT$PEUR, main="Normalized False Positives vs Percent European", xlab="% European in Group", ylab="Normalized False Positives", pch=20, col="mediumpurple3",xlim=c(0,1),ylim=Y_LIM)
MOD <- lm(ADM_DAT$FP ~ ADM_DAT$PEUR)
abline(MOD, col="mediumpurple3", lwd=3)
text(.2,min(Y_LIM)+.06, labels=paste("R2=",round(summary(MOD)$r.squared,4)))
text(.2,min(Y_LIM)+.03, labels=paste("p=",round(summary(MOD)$coefficients[2,4],4)))
dev.off()

NORM_ER <- NORM_MS + NORM_FP
NORM_ERp <- t(t(NORM_ER)/LENS)
NORM_ERb <- data.frame( PERC=c(t(NORM_ERp)),GRP=c(rep("M",4),rep("N",4),rep("A",4),rep("E",4)), CHR=rep(c("0","1","2","3"),4) )
MOD <- lm( PERC ~ GRP + as.factor(CHR), data=NORM_ERb )
anova(MOD)
























##############################################################################################
## END OF DOC ################################################################################
##############################################################################################














## Fix Positions Listed Twice #################
V <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	PathToV <- paste(PathToFiles,"ENV_MRGD_NA12878_",IDS[chrom],".",SI,".vcf", sep="")
	print(PathToV)
	V[[chrom]] <- read.table( PathToV, header=F, sep="\t", colClass=c("factor", "numeric", "character", "factor", "factor", "numeric", "factor", "factor", rep("character",5))) 
	colnames(V[[chrom]]) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "H", "N", "M", "A", "E")
	print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
}
names(V) <- IDS

## Fix Positions Listed Twice #################
RPT_V <- RPT_Vi <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	RPT_V[[chrom]] <- V[[chrom]]$POS[which(duplicated(V[[chrom]]$POS))]
	RPT_Vi[[chrom]] <- which( V[[chrom]]$POS %in% RPT_V[[chrom]] ) # Some positions are listed 3 times
	print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
}

## Sets mentioned in the columen "INFO" refer to first REF/ALT value
 # Sets NOT mentioned in "INFO" are AFTER the ":" in REF/ALT
FLAT_V <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	if (length(RPT_Vi[[chrom]])>0) {
		FLAT_V[[chrom]] <- GET_RPTS(V[[chrom]], RPT_V[[chrom]])
		# Take out Rows in RPT_Xi, put in X_RPTs, sort by position
		V[[chrom]] <- rbind(V[[chrom]][-RPT_Vi[[chrom]],],FLAT_V[[chrom]]$Table)
		V[[chrom]] <- V[[chrom]][order(V[[chrom]]$POS),]
		print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
	}
}

###############################################
## Organize into More Convenient Formats ######

# Pull out genotypes at each position (using V tables)
Genos <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	Genos[[chrom]] <- list()	

	Genos[[chrom]]$H <- sapply( strsplit(as.character(V[[chrom]]$H), ":"), "[", 1)
	Genos[[chrom]]$H <- gsub("|", "/", Genos[[chrom]]$H, fixed=T)
	temp <- t(sapply( strsplit(Genos[[chrom]]$H,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$H[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[chrom]]$M <- sapply( strsplit(as.character(V[[chrom]]$M), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[chrom]]$M,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$M[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[chrom]]$N <- sapply( strsplit(as.character(V[[chrom]]$N), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[chrom]]$N,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$N[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[chrom]]$A <- sapply( strsplit(as.character(V[[chrom]]$A), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[chrom]]$A,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$A[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	Genos[[chrom]]$E <- sapply( strsplit(as.character(V[[chrom]]$E), ":"), "[", 1)
	temp <- t(sapply( strsplit(Genos[[chrom]]$E,"/"), "[", c(1,2) ))
	misorder <- which(temp[,1]>temp[,2])
	Genos[[chrom]]$E[misorder] <- paste( temp[misorder,2],temp[misorder,1], sep="/")

	print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
}
names(Genos) <- IDS

# Make array with easy-access genotypes
C_Vr <- SNP <- SNP_1 <- list()
start <- proc.time()
for ( chrom in 1:NUM ) {
	C_Vr[[chrom]] <- data.frame(V[[chrom]][,2:5], H=Genos[[chrom]]$H, M=Genos[[chrom]]$M, N=Genos[[chrom]]$N, A=Genos[[chrom]]$A, E=Genos[[chrom]]$E)	

	# Make table of alt allele count (SNP Array) for each position and group
	temp <- as.matrix(C_Vr[[chrom]][,5:9])
	SNP[[chrom]] <- array( 2, dim(C_Vr[[chrom]][,5:9]) )
	SNP[[chrom]][grep("./.",temp,fixed=T)] <- 0
	SNP[[chrom]][grep("0/",temp,)] <- 1
	SNP[[chrom]][grep("0/0",temp,)] <- 0
	SNP[[chrom]] <- data.frame(SNP[[chrom]])

	# Make table of variant positions
	SNP_1[[chrom]] <- data.frame(ceiling(SNP[[chrom]]/2))
	colnames(SNP[[chrom]]) <- colnames(SNP_1[[chrom]]) <- c("HI", "MIX", "NA", "ASN", "EUR")

	print(paste( "Done with:",IDS[chrom],"-",round((proc.time()-start)[3],2) ))
}
names(SNP) <- names(SNP_1) <- IDS



