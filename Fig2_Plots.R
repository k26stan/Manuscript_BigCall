## Plot ANC_ENV_TEST Results ##
 # Include 
## June 26, 2014 ##
## Kristopher Standish ##

## Plot HC Runtimes vs Group Size
## Plot Summary Statistics for Ancestry Investigation
 # Sensitivity
   # Out of the total TRUE positives, what fraction do you identify?
   # aka, True Positive Rate
 # Specificity
   # Out of the total TRUE negatives, what fraction do you identify?
   # aka, True Negative Rate
 # Accuracy
   # How many calls (either way) did you get right out of all possible calls?
 # PPV
   # Out of the calls you identify as Positive, what fraction actually ARE positive?
   # 1 - False Discovery Rate
 # NPV
   # Out of the calls you identify as Negative, what fraction actually ARE negative?
   # 1 - False Omission Rate



# Plot Sens/Spec vs % EUR all on same 2 plots ##

PLOT_COLS <- c("firebrick2","sienna2","gold2","chartreuse3","cadetblue2","steelblue3","slateblue3","black","grey50","white")
PLOT_COLS.4 <- c("firebrick4","sienna4","gold4","chartreuse4","cadetblue4","steelblue4","slateblue4","black","grey50","white")

#####################################################################
## LOAD SENS/SPEC DATA ##############################################
#####################################################################

## Set Date & Path to Save
DATE <- gsub("-","",Sys.Date())
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/PLOTS/Revisions/",DATE,"_Fig2/",sep="" )
dir.create( PathToSave )

## Set Paths
PathToData <- "/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/DATA/Revisions/Ancestry/"

## Load Summarized Data
SNP <- read.table( paste(PathToData,"Summary_Stats_SNP_FIL.txt",sep=""), sep="\t",header=T,stringsAsFactors=F )
IND <- read.table( paste(PathToData,"Summary_Stats_IND_FIL.txt",sep=""), sep="\t",header=T,stringsAsFactors=F )
SNP$GRP[which(is.na(SNP$GRP))] <- "NAA"
IND$GRP[which(is.na(IND$GRP))] <- "NAA"
 # Calculate Negative Predictive Value
SNP$NPV <- SNP$True_N / SNP$Test_N
IND$NPV <- IND$True_N / IND$Test_N

#####################################################################
## LOAD/ORGANIZE HC DATA ############################################
#####################################################################

## Set Paths
PathToData <- "/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/DATA/Revisions/Updates/HC_vs_GrpSize/ProgressMeter_Files/" ; FILE_TAG <- "Chr_21"
# PathToData <- "/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/DATA/Revisions/Updates/HC_vs_GrpSize/ProgressMeter_Files/Chr_20/" ; FILE_TAG <- "Chr_20"

## Load Files
File_List <- list.files( PathToData )
File_List <- grep( "Progress",File_List, value=T )
File_Names <- gsub( "Progress_HC_","",File_List )
File_Names <- gsub( ".run.oe","",File_Names, fixed=T )
Grp_Size <- as.numeric(gsub("Grp","", lapply(strsplit(File_Names,"_"),function(x)x[1]) ))
Grp_Size.uniq <- sort(unique(Grp_Size))
Iter <- as.numeric(gsub("It","", lapply(strsplit(File_Names,"_"),function(x)x[2]) ))
RAW <- DAT <- list()
for ( f in 1:length(File_List) ) {
	file <- File_List[f]
	name <- File_Names[f]
	# Load File
	RAW[[name]] <- readLines( paste(PathToData,file,sep="") )
	# Get Final Run Stats
	if ( any(grepl("done",RAW[[name]])) ) {
		Incomplete <- F
		Final.site <- strsplit( RAW[[name]][grep("done",RAW[[name]])], " +" )[[1]][6]
		Final.sec <- strsplit( RAW[[name]][grep("done",RAW[[name]])+1], " +" )[[1]][7]
		Final <- as.numeric(c( Final.site, Final.sec)) ; names(Final) <- c("Sites","Secs")
		Split.1 <- strsplit( RAW[[name]][4:(grep("done",RAW[[name]])-1)], " +" )
	}else{
		Incomplete <- T
		Split.1 <- strsplit( RAW[[name]][4:length(RAW[[name]])], " +" )
	}
	# Get Running Progress
	Split.loc <- as.numeric(unlist(lapply( strsplit( unlist(lapply(Split.1,function(x)x[5])), ":"), function(x)x[2] )))
	Split.time.1 <- matrix( unlist(lapply(Split.1,function(x)x[7:8])), ncol=2,byrow=T )
	Split.time <- as.numeric( Split.time.1[,1] )
	Split.time[which(Split.time.1[,2]=="m")] <- 60 * Split.time[which(Split.time.1[,2]=="m")]
	Split.time[which(Split.time.1[,2]=="h")] <- 3600 * Split.time[which(Split.time.1[,2]=="h")]
	Split.time[which(Split.time.1[,2]=="d")] <- 24* 3600 * Split.time[which(Split.time.1[,2]=="d")]
	if ( Incomplete==T ) { Final <- c( Split.loc[length(Split.loc)], Split.time[length(Split.time)] ) ; names(Final) <- c("Sites","Secs") }
	# Combine Stats
	RUN <- data.frame( LOC=Split.loc, SEC=Split.time )
	TAGS <- c( Grp_Size[f], Iter[f] )
	COMPILE <- list( FIN=Final, RUN=RUN, TAG=TAGS, INC=Incomplete )
	DAT[[name]] <- COMPILE
}

#####################################################################
## FCT: PLOT SENS/SPEC STATS ########################################
#####################################################################

PLOT_SS <- function() {

	#########################################################
	## Sensitivity & Specificity: 2 Plots in single window

	## Set universal parameters
	COLS <- PLOT_COLS[c(2,6)] ; names(COLS) <- c("SNP","IND") # PLOT_COLS[c(4,6,4,6)] # PLOT_COLS[c(7,1,7,1)] # PLOT_COLS[c(7,1,4,6)]
	LINE_COL <- PLOT_COLS[9]
	LTYS <- c(1,1) ; names(LTYS) <- c("SNP","IND")
	PCHS <- c(15:18)[factor(SNP$GRP)] # rep( 15:18, rep(4,4) ) # rep( c("*","+","o","x"),c(4,4,4,4) )
	XLIM <- c(0,1)
	## Sensitivity
	 # Calculate Resids vs Chr
	SENS.SNP.chr.res <- resid(lm( SNP$Sensitivity ~ SNP$CHR ))
	SENS.IND.chr.res <- resid(lm( IND$Sensitivity ~ IND$CHR ))
	 # Model Fits & Plot Lines
	MOD.SNP <- lm( Sensitivity ~ europe+CHR, data=SNP )
	P_VAL.SNP <- anova(MOD.SNP)["europe","Pr(>F)"] # summary(MOD.SNP)$coefficients[2,4]
	MOD.IND <- lm( Sensitivity ~ europe+CHR, data=IND )
	P_VAL.IND <- anova(MOD.IND)["europe","Pr(>F)"] # summary(MOD.IND)$coefficients[2,4]
	 # Parameters
	YLIM.SENS <- range( SENS.SNP.chr.res, SENS.IND.chr.res )
	 # Open Plot
	plot(0,0,type="n",xlim=XLIM,ylim=YLIM.SENS,main="Sensitivity vs % European Admixture",xlab="% European",ylab="Sensitivity (Residuals vs Chr)")
	 # Plot Dashed Lines
	ORDER <- floor(log10(diff(YLIM.SENS)))
	LINES <- seq( -3*10^(ORDER), 3*10^(ORDER), 2*10^(ORDER-1) )
	abline( v=seq(0,1,.2), lty=3,col=LINE_COL,lwd=1)
	abline( h=LINES, lty=3,col=LINE_COL,lwd=1) ; abline(h=0,lty=1,col="black",lwd=1)
	 # Plot Points
	points( SNP$europe, SENS.SNP.chr.res, col=COLS["SNP"],pch=PCHS) # points( SNP$europe, SNP$Sensitivity, col=COLS["SNP"],pch=PCHS)
	points( IND$europe, SENS.IND.chr.res, col=COLS["IND"],pch=PCHS)
	 # Plot Lines of Fit
	MOD.SNP.line <- lm( SENS.SNP.chr.res ~ SNP$europe )
	MOD.IND.line <- lm( SENS.IND.chr.res ~ IND$europe )
	abline(MOD.SNP.line, col=COLS["SNP"], lty=LTYS["SNP"], lwd=3)
	abline(MOD.IND.line, col=COLS["IND"], lty=LTYS["IND"], lwd=3)
	 # Add P-Values to Plot
	STAR.SNP <- STAR.IND <- " "
	if ( P_VAL.SNP<0.05 & P_VAL.SNP>.01 ) { STAR.SNP <- "*" } ; if ( P_VAL.SNP<.01 ) { STAR.SNP <- "**" }
	text( quantile(XLIM,.99),quantile(YLIM.SENS,.05),labels=paste(STAR.SNP,"p =",formatC(P_VAL.SNP, digits=2, format="e")),col=COLS["SNP"],pos=2)
	if ( P_VAL.IND<0.05 & P_VAL.IND>.01 ) { STAR.IND <- "*" } ; if ( P_VAL.IND<.01 ) { STAR.IND <- "**" }
	text( quantile(XLIM,.99),quantile(YLIM.SENS,.10),labels=paste(STAR.IND,"p =",formatC(P_VAL.IND, digits=2, format="e")),col=COLS["IND"],pos=2)
	## Specificity
	 # Calculate Resids vs Chr
	SPEC.SNP.chr.res <- resid(lm( SNP$Specificity ~ SNP$CHR ))
	SPEC.IND.chr.res <- resid(lm( IND$Specificity ~ IND$CHR ))
	 # Model Fits & Plot Lines
	MOD.SNP <- lm( Specificity ~ europe+CHR, data=SNP )
	P_VAL.SNP <- anova(MOD.SNP)["europe","Pr(>F)"] # summary(MOD.SNP)$coefficients[2,4]
	MOD.IND <- lm( Specificity ~ europe+CHR, data=IND )
	P_VAL.IND <- anova(MOD.IND)["europe","Pr(>F)"] # summary(MOD.IND)$coefficients[2,4]
	 # Parameters
	YLIM.SPEC <- range( SPEC.SNP.chr.res, SPEC.IND.chr.res )
	 # Open Plot
	plot(0,0,type="n",xlim=XLIM,ylim=YLIM.SPEC,main="Specificity vs % European Admixture",xlab="% European",ylab="Specificity (Residuals vs Chr)")
	 # Plot Dashed Lines
	ORDER <- floor(log10(diff(YLIM.SPEC)))
	LINES <- seq( -3*10^(ORDER), 3*10^(ORDER), 2*10^(ORDER-1) )
	abline( v=seq(0,1,.2), lty=3,col=LINE_COL,lwd=1)
	abline( h=LINES, lty=3,col=LINE_COL,lwd=1) ; abline(h=0,lty=1,col="black",lwd=1)
	 # Plot Points
	points( SNP$europe, SPEC.SNP.chr.res, col=COLS["SNP"],pch=PCHS) # points( SNP$europe, SNP$Specificity, col=COLS["SNP"],pch=PCHS)
	points( IND$europe, SPEC.IND.chr.res, col=COLS["IND"],pch=PCHS)
	 # Plot Lines of Fit
	MOD.SNP.line <- lm( SPEC.SNP.chr.res ~ SNP$europe )
	MOD.IND.line <- lm( SPEC.IND.chr.res ~ IND$europe )
	abline(MOD.SNP.line, col=COLS["SNP"], lty=LTYS["SNP"], lwd=3)
	abline(MOD.IND.line, col=COLS["IND"], lty=LTYS["IND"], lwd=3)
	 # Add P-Values to Plot
	STAR.SNP <- STAR.IND <- " "
	if ( P_VAL.SNP<0.05 & P_VAL.SNP>.01 ) { STAR.SNP <- "*" } ; if ( P_VAL.SNP<.01 ) { STAR.SNP <- "**" }
	text( quantile(XLIM,.99),quantile(YLIM.SPEC,.05),labels=paste(STAR.SNP,"p =",formatC(P_VAL.SNP, digits=2, format="e")),col=COLS["SNP"],pos=2)
	if ( P_VAL.IND<0.05 & P_VAL.IND>.01 ) { STAR.IND <- "*" } ; if ( P_VAL.IND<.01 ) { STAR.IND <- "**" }
	text( quantile(XLIM,.99),quantile(YLIM.SPEC,.10),labels=paste(STAR.IND,"p =",formatC(P_VAL.IND, digits=2, format="e")),col=COLS["IND"],pos=2)
	## Legend
	legend(quantile(XLIM,0),quantile(YLIM.SPEC,1),legend=c("SNPs","Indels"),col=COLS,lty=LTYS,lwd=3,cex=.7)
	legend(quantile(XLIM,0.3),quantile(YLIM.SPEC,1),legend=unique(SNP$GRP),pch=unique(PCHS),lwd=3,cex=.7,ncol=2)

} # Close "PLOT_SS" Function

#####################################################################
## PLOT RESULTS #####################################################
#####################################################################

PLOT_HC <- function(DAT,Grp_Size) {

	##############################################
	## Final Stats vs Group Size #################

	## Compile Data
	TIMES <- unlist(lapply( DAT, function(x) x$FIN["Secs"] ))
	FINAL <- data.frame( SIZE=Grp_Size, SEC=TIMES, SEC_SAMP=TIMES/Grp_Size, SU=round(16*TIMES/3600,0), SU_SAMP=round(16*TIMES/3600,0)/Grp_Size )

	## Plot Raw Numbers (Run Time vs Group Size)
	 # Raw Parameters
	XLIM <- c( 0, max( FINAL$SIZE) ) * c(1,2) # c(1,4) # 
	YLIM <- c( 0, max(FINAL$SU) ) * c(1,3) # c(1,5) # 
	COLS <- PLOT_COLS[c(4,7,8)] # PLOT_COLS[c(2,6,8)]
	 # Make Plot
	plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Group Size",ylab="Total CPU Hours (SU)",main=paste("Computational Cost vs Groupsize -",FILE_TAG),xaxt="n",yaxt="n" )
	abline( h=seq(0,YLIM[2]+500,500), lty=3,col="grey50",lwd=1 )
	axis( 2, at=seq(0,10e3,5e2) )
	# abline( v=c(Grp_Size.uniq,50,100,200,300,400,437), lty=3,col="grey50",lwd=1 )
	# axis( 1, at=c(Grp_Size.uniq,50,100,200,300,400,437) )
	abline( v=seq(0,500,20), lty=3,col="grey50",lwd=1 )
	axis( 1, at=seq(0,500,20) )
	 # Plot Fits
	Grp_Size.seq <- seq( 1, 5*max(Grp_Size.uniq), 1 )
	MOD <- nls( SU ~ a + b^SIZE, data=FINAL, start=list(a=0,b=1.01))
	MOD.q <- lm( SU ~ I(SIZE^2)+SIZE, data=FINAL )
	lines( Grp_Size.seq, predict(MOD.q,list(SIZE=Grp_Size.seq)), col=COLS[1], lwd=3 )
	MOD.l <- lm( SU ~ SIZE, data=FINAL )
	lines( Grp_Size.seq, predict(MOD.l,list(SIZE=Grp_Size.seq)), col=COLS[2], lwd=3,lty=2)
	 # Summarize Fits
	ADJ.q <- round( summary(MOD.q)$adj.r.squared, 4)
	ADJ.l <- round( summary(MOD.l)$adj.r.squared, 4)
	text( quantile(XLIM,.01),quantile(YLIM,.66), label=paste("Quad: Adj.R2=",ADJ.q,sep=""),col=COLS[1], pos=4 )
	text( quantile(XLIM,.01),quantile(YLIM,.60), label=paste("Linr: Adj.R2=",ADJ.l,sep=""),col=COLS[2], pos=4 )
	 # Plot Actual Data
	points( SU ~ SIZE, data=FINAL, col=COLS[3], pch="+",cex=1.5,lwd=2 )
	 # Legend
	legend("topleft",col=COLS[1:2],legend=c("Quad","Lin"),title="Fits",lty=c(1,2),lwd=3 )

	## Plot Per Sample Numbers (Run Time vs Group Size)
	 # Per Sample Parameters
	# XLIM <- range( FINAL$SIZE )
	YLIM <- c( 0, max(FINAL$SU_SAMP) ) # * c(1,3) #
	 # Make Plot
	plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Group Size",ylab="CPU Hours (SU) per Sample",main=paste("Per Sample Computational Cost vs Groupsize -",FILE_TAG),xaxt="n" )
	abline( h=seq(0,YLIM[2]+20,2), lty=3,col="grey50",lwd=1 )
	# abline( v=c(Grp_Size.uniq,50,100,200,300,400,437), lty=3,col="grey50",lwd=1 )
	# axis( 1, at=c(Grp_Size.uniq,50,100,200,300,400,437) )
	abline( v=seq(0,500,20), lty=3,col="grey50",lwd=1 )
	axis( 1, at=seq(0,500,20) )
	 # Plot Fits
	MOD.q.2 <- lm( SU_SAMP ~ I(1/SIZE)+SIZE, data=FINAL )
	lines( Grp_Size.seq, predict(MOD.q.2,list(SIZE=Grp_Size.seq)), col=COLS[1], lwd=3 )
	MOD.l.2 <- lm( SU_SAMP ~ I(1/SIZE), data=FINAL )
	lines( Grp_Size.seq, predict(MOD.l.2,list(SIZE=Grp_Size.seq)), col=COLS[2], lwd=3,lty=2 )
	 # Summarize Fits
	ADJ.q.2 <- round( summary(MOD.q.2)$adj.r.squared, 4)
	ADJ.l.2 <- round( summary(MOD.l.2)$adj.r.squared, 4)
	text( quantile(XLIM,.01),quantile(YLIM,.26), label=paste("Quad: Adj.R2=",ADJ.q.2,sep=""),col=COLS[1], pos=4 )
	text( quantile(XLIM,.01),quantile(YLIM,.20), label=paste("Linr: Adj.R2=",ADJ.l.2,sep=""),col=COLS[2], pos=4 )
	 # Plot Actual Data
	points( SU_SAMP ~ SIZE, data=FINAL, col=COLS[3], pch="+",cex=1.5,lwd=2 )

} # Close "PLOT_HC" Function


#####################################################################
## PLOT ALL TOGETHER ################################################
#####################################################################

## Which HC Tests Finished Running?
DAT.fin <- DAT[ which(unlist(lapply(DAT,function(x)x$INC))==F) ]
Grp_Size.fin <- Grp_Size[ which(unlist(lapply(DAT,function(x)x$INC))==F) ]

## Create File for Entire Figure 2 (At least first parts)
png(paste(PathToSave,"2_A-D.2.png",sep=""),width=2400,height=2000,pointsize=40)
 # Specify Layout
par(mfrow=c(2,2))
 # Create Plots
PLOT_HC(DAT.fin, Grp_Size.fin)
PLOT_SS()

dev.off()





#####################################################################
## SUPPLEMENTAL PLOTS ###############################################
#####################################################################


#########################################################
## Sites vs Time ########################################

## Save Plot of Sites vs Time
png( paste(PathToSave,"Supp-HC_Sites_v_Time.",FILE_TAG,".png",sep=""), height=1000,width=1400,pointsize=30 )
## Set Plotting Parameters
COLS.which <- colorRampPalette(PLOT_COLS[1:7])(length(Grp_Size.uniq))
COLS <- COLS.which[factor(Grp_Size)]
XLIM <- range( lapply(DAT, function(x) range(x$RUN$SEC,na.rm=T) ) )
BY <- 5*3600
XTICK <- seq(0,XLIM[2]+BY,BY)
YLIM <- range( lapply(DAT, function(x) range(x$RUN$LOC,na.rm=T) ) )
## Create Plot
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Walltime (Hrs)",ylab="Sites", main=paste("ProgressMeter -",FILE_TAG), xaxt="n" )
axis( 1, at=XTICK, label=(BY/3600)*0:(length(XTICK)-1) )
abline( v=seq(0,XLIM[2]+1800,1800), lty=3,col="grey50",lwd=1 )
abline( h=seq(0,YLIM[2]+1e7,1e7), lty=3,col="grey50",lwd=1 )
## Add Data
for ( f in 1:length(DAT) ) {
	name <- names(DAT)[f]
	xvals <- DAT[[name]]$RUN$SEC
	yvals <- DAT[[name]]$RUN$LOC
	final <- DAT[[name]]$FIN
	points( xvals, yvals, col=COLS[f], type="l",lwd=2 )
	points( final["Secs"], final["Sites"], col=COLS[f], pch=10, lwd=2,cex=1.5 )
}
 # Legend
legend( "bottomright", fill=COLS.which,legend=levels(factor(Grp_Size)), bg="white",title="Group Size",ncol=8,cex=.8 )
dev.off()


#########################################################
## Accuracy & PPV & NPV: 2 Plots in single window #######

## Set universal parameters
COLS <- PLOT_COLS[c(2,6)] ; names(COLS) <- c("SNP","IND") # PLOT_COLS[c(4,6,4,6)] # PLOT_COLS[c(7,1,7,1)] # PLOT_COLS[c(7,1,4,6)]
LINE_COL <- PLOT_COLS[9]
LTYS <- c(1,1) ; names(LTYS) <- c("SNP","IND")
PCHS <- c(15:18)[factor(SNP$GRP)] # rep( 15:18, rep(4,4) ) # rep( c("*","+","o","x"),c(4,4,4,4) )
XLIM <- c(0,1)
## Open File to Save to
png(paste(PathToSave,"Supp-Perc_EUR_AccuPPV.png",sep=""),height=1000,width=3000,pointsize=36)
par(mfrow=c(1,3))
## Accuracy
 # Calculate Resids vs Chr
ACCU.SNP.chr.res <- resid(lm( SNP$Accuracy ~ SNP$CHR ))
ACCU.IND.chr.res <- resid(lm( IND$Accuracy ~ IND$CHR ))
 # Model Fits & Plot Lines
MOD.SNP <- lm( Accuracy ~ europe+CHR, data=SNP )
P_VAL.SNP <- anova(MOD.SNP)["europe","Pr(>F)"] # summary(MOD.SNP)$coefficients[2,4]
MOD.IND <- lm( Accuracy ~ europe+CHR, data=IND )
P_VAL.IND <- anova(MOD.IND)["europe","Pr(>F)"] # summary(MOD.IND)$coefficients[2,4]
 # Parameters
YLIM.ACCU <- range( ACCU.SNP.chr.res, ACCU.IND.chr.res )
 # Open Plot
plot(0,0,type="n",xlim=XLIM,ylim=YLIM.ACCU,main="Accuracy vs % European Admixture",xlab="% European",ylab="Accuracy (Residuals vs Chr)")
 # Plot Dashed Lines
ORDER <- floor(log10(diff(YLIM.ACCU)))
LINES <- seq( -3*10^(ORDER), 3*10^(ORDER), 2*10^(ORDER-1) )
abline( h=LINES, lty=3,col=LINE_COL,lwd=1) ; abline(h=0,lty=1,col="black",lwd=1)
 # Plot Points
points( SNP$europe, ACCU.SNP.chr.res, col=COLS["SNP"],pch=PCHS) # points( SNP$europe, SNP$Accuracy, col=COLS["SNP"],pch=PCHS)
points( IND$europe, ACCU.IND.chr.res, col=COLS["IND"],pch=PCHS)
 # Plot Lines of Fit
MOD.SNP.line <- lm( ACCU.SNP.chr.res ~ SNP$europe )
MOD.IND.line <- lm( ACCU.IND.chr.res ~ IND$europe )
abline(MOD.SNP.line, col=COLS["SNP"], lty=LTYS["SNP"], lwd=3)
abline(MOD.IND.line, col=COLS["IND"], lty=LTYS["IND"], lwd=3)
 # Add P-Values to Plot
STAR.SNP <- STAR.IND <- " "
if ( P_VAL.SNP<0.05 & P_VAL.SNP>.01 ) { STAR.SNP <- "*" } ; if ( P_VAL.SNP<.01 ) { STAR.SNP <- "**" }
text( quantile(XLIM,.95),quantile(YLIM.ACCU,.05),labels=paste(STAR.SNP,"p =",formatC(P_VAL.SNP, digits=2, format="e")),col=COLS["SNP"],pos=2)
if ( P_VAL.IND<0.05 & P_VAL.IND>.01 ) { STAR.IND <- "*" } ; if ( P_VAL.IND<.01 ) { STAR.IND <- "**" }
text( quantile(XLIM,.95),quantile(YLIM.ACCU,.10),labels=paste(STAR.IND,"p =",formatC(P_VAL.IND, digits=2, format="e")),col=COLS["IND"],pos=2)
## PPV
 # Calculate Resids vs Chr
PPV.SNP.chr.res <- resid(lm( SNP$PPV ~ SNP$CHR ))
PPV.IND.chr.res <- resid(lm( IND$PPV ~ IND$CHR ))
 # Model Fits & Plot Lines
MOD.SNP <- lm( PPV ~ europe+CHR, data=SNP )
P_VAL.SNP <- anova(MOD.SNP)["europe","Pr(>F)"] # summary(MOD.SNP)$coefficients[2,4]
MOD.IND <- lm( PPV ~ europe+CHR, data=IND )
P_VAL.IND <- anova(MOD.IND)["europe","Pr(>F)"] # summary(MOD.IND)$coefficients[2,4]
 # PPV Parameters
YLIM.PPV <- range( PPV.SNP.chr.res, PPV.IND.chr.res )
 # Open Plot
plot(0,0,type="n",xlim=XLIM,ylim=YLIM.PPV,main="PPV vs % European Admixture",xlab="% European",ylab="PPV (Residuals vs Chr)")
 # Plot Dashed Lines
ORDER <- floor(log10(diff(YLIM.PPV)))
LINES <- seq( -3*10^(ORDER), 3*10^(ORDER), 2*10^(ORDER-1) )
abline( h=LINES, lty=3,col=LINE_COL,lwd=1) ; abline(h=0,lty=1,col="black",lwd=1)
 # Plot Points
points( SNP$europe, PPV.SNP.chr.res, col=COLS["SNP"],pch=PCHS) # points( SNP$europe, SNP$PPV, col=COLS["SNP"],pch=PCHS)
points( IND$europe, PPV.IND.chr.res, col=COLS["IND"],pch=PCHS)
 # Plot Lines of Fit
MOD.SNP.line <- lm( PPV.SNP.chr.res ~ SNP$europe )
MOD.IND.line <- lm( PPV.IND.chr.res ~ IND$europe )
abline(MOD.SNP.line, col=COLS["SNP"], lty=LTYS["SNP"], lwd=3)
abline(MOD.IND.line, col=COLS["IND"], lty=LTYS["IND"], lwd=3)
 # Add P-Values to Plot
STAR.SNP <- STAR.IND <- " "
if ( P_VAL.SNP<0.05 & P_VAL.SNP>.01 ) { STAR.SNP <- "*" } ; if ( P_VAL.SNP<.01 ) { STAR.SNP <- "**" }
text( quantile(XLIM,.95),quantile(YLIM.PPV,.05),labels=paste(STAR.SNP,"p =",formatC(P_VAL.SNP, digits=2, format="e")),col=COLS["SNP"],pos=2)
if ( P_VAL.IND<0.05 & P_VAL.IND>.01 ) { STAR.IND <- "*" } ; if ( P_VAL.IND<.01 ) { STAR.IND <- "**" }
text( quantile(XLIM,.95),quantile(YLIM.PPV,.10),labels=paste(STAR.IND,"p =",formatC(P_VAL.IND, digits=2, format="e")),col=COLS["IND"],pos=2)
## Legend
legend(quantile(XLIM,0),quantile(YLIM.PPV,1),legend=c("SNPs","Indels"),col=COLS,lty=LTYS,lwd=3,cex=.7)
legend(quantile(XLIM,0.3),quantile(YLIM.PPV,1),legend=unique(SNP$GRP),pch=unique(PCHS),lwd=3,cex=.7,ncol=2)
## NPV
 # Calculate Resids vs Chr
NPV.SNP.chr.res <- resid(lm( SNP$NPV ~ SNP$CHR ))
NPV.IND.chr.res <- resid(lm( IND$NPV ~ IND$CHR ))
 # Model Fits & Plot Lines
MOD.SNP <- lm( NPV ~ europe+CHR, data=SNP )
P_VAL.SNP <- anova(MOD.SNP)["europe","Pr(>F)"] # summary(MOD.SNP)$coefficients[2,4]
MOD.IND <- lm( NPV ~ europe+CHR, data=IND )
P_VAL.IND <- anova(MOD.IND)["europe","Pr(>F)"] # summary(MOD.IND)$coefficients[2,4]
 # NPV Parameters
YLIM.NPV <- range( NPV.SNP.chr.res, NPV.IND.chr.res )
 # Open Plot
plot(0,0,type="n",xlim=XLIM,ylim=YLIM.NPV,main="NPV vs % European Admixture",xlab="% European",ylab="NPV (Residuals vs Chr)")
 # Plot Dashed Lines
ORDER <- floor(log10(diff(YLIM.NPV)))
LINES <- seq( -3*10^(ORDER), 3*10^(ORDER), 2*10^(ORDER-1) )
abline( h=LINES, lty=3,col=LINE_COL,lwd=1) ; abline(h=0,lty=1,col="black",lwd=1)
 # Plot Points
points( SNP$europe, NPV.SNP.chr.res, col=COLS["SNP"],pch=PCHS) # points( SNP$europe, SNP$NPV, col=COLS["SNP"],pch=PCHS)
points( IND$europe, NPV.IND.chr.res, col=COLS["IND"],pch=PCHS)
 # Plot Lines of Fit
MOD.SNP.line <- lm( NPV.SNP.chr.res ~ SNP$europe )
MOD.IND.line <- lm( NPV.IND.chr.res ~ IND$europe )
abline(MOD.SNP.line, col=COLS["SNP"], lty=LTYS["SNP"], lwd=3)
abline(MOD.IND.line, col=COLS["IND"], lty=LTYS["IND"], lwd=3)
 # Add P-Values to Plot
STAR.SNP <- STAR.IND <- " "
if ( P_VAL.SNP<0.05 & P_VAL.SNP>.01 ) { STAR.SNP <- "*" } ; if ( P_VAL.SNP<.01 ) { STAR.SNP <- "**" }
text( quantile(XLIM,.95),quantile(YLIM.NPV,.05),labels=paste(STAR.SNP,"p =",formatC(P_VAL.SNP, digits=2, format="e")),col=COLS["SNP"],pos=2)
if ( P_VAL.IND<0.05 & P_VAL.IND>.01 ) { STAR.IND <- "*" } ; if ( P_VAL.IND<.01 ) { STAR.IND <- "**" }
text( quantile(XLIM,.95),quantile(YLIM.NPV,.10),labels=paste(STAR.IND,"p =",formatC(P_VAL.IND, digits=2, format="e")),col=COLS["IND"],pos=2)
dev.off()




#####################################################################
## END OF DOC #######################################################
#####################################################################
