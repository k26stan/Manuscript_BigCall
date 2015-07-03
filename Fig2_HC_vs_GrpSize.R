## Check out ProgressMeter Files for HaplotypeCaller ##
## July 1, 2015 ##
## Kristopher Standish ##
PLOT_COLS <- c("firebrick2","sienna2","gold2","chartreuse3","cadetblue2","steelblue3","slateblue3","black","grey50","white")
PLOT_COLS.4 <- c("firebrick4","sienna4","gold4","chartreuse4","cadetblue4","steelblue4","slateblue4","black","grey50","white")

#####################################################################
## LOAD/ORGANIZE DATA ###############################################
#####################################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths
PathToData <- "/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/DATA/Revisions/Updates/HC_vs_GrpSize/ProgressMeter_Files/" ; FILE_TAG <- "Chr_21"
# PathToData <- "/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/DATA/Revisions/Updates/HC_vs_GrpSize/ProgressMeter_Files/Chr_20/" ; FILE_TAG <- "Chr_20"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/PLOTS/Revisions/",DATE,"_HCvsGrpSize/",sep="" )
dir.create( PathToSave )

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
	Final.site <- strsplit( RAW[[name]][grep("done",RAW[[name]])], " +" )[[1]][6]
	Final.sec <- strsplit( RAW[[name]][grep("done",RAW[[name]])+1], " +" )[[1]][7]
	Final <- as.numeric(c( Final.site, Final.sec)) ; names(Final) <- c("Sites","Secs")
	# Get Running Progress
	Split.1 <- strsplit( RAW[[name]][4:(grep("done",RAW[[name]])-1)], " +" )
	Split.loc <- as.numeric(unlist(lapply( strsplit( unlist(lapply(Split.1,function(x)x[5])), ":"), function(x)x[2] )))
	Split.time.1 <- matrix( unlist(lapply(Split.1,function(x)x[7:8])), ncol=2,byrow=T )
	Split.time <- as.numeric( Split.time.1[,1] )
	Split.time[which(Split.time.1[,2]=="m")] <- 60 * Split.time[which(Split.time.1[,2]=="m")]
	Split.time[which(Split.time.1[,2]=="h")] <- 3600 * Split.time[which(Split.time.1[,2]=="h")]
	Split.time[which(Split.time.1[,2]=="d")] <- 24* 3600 * Split.time[which(Split.time.1[,2]=="d")]
	# Combine Stats
	RUN <- data.frame( LOC=Split.loc, SEC=Split.time )
	TAGS <- c( Grp_Size[f], Iter[f] )
	COMPILE <- list( FIN=Final, RUN=RUN, TAG=TAGS )
	DAT[[name]] <- COMPILE
}

#####################################################################
## PLOT RESULTS #####################################################
#####################################################################

##############################################
## Sites vs Time #############################

## Save Plot of Sites vs Time
png( paste(PathToSave,"HC_Sites_v_Time.",FILE_TAG,".png",sep=""), height=1000,width=1000,pointsize=30 )
## Set Plotting Parameters
COLS.which <- colorRampPalette(PLOT_COLS[1:7])(length(Grp_Size.uniq))
COLS <- COLS.which[factor(Grp_Size)]
XLIM <- range( lapply(DAT, function(x) range(x$RUN$SEC,na.rm=T) ) )
XTICK <- seq(0,XLIM[2]+3600,3600)
YLIM <- range( lapply(DAT, function(x) range(x$RUN$LOC,na.rm=T) ) )
## Create Plot
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Walltime (Hrs)",ylab="Sites", main=paste("ProgressMeter -",FILE_TAG), xaxt="n" )
axis( 1, at=XTICK, label=1:length(XTICK) )
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
legend( "bottomright", fill=COLS.which,legend=levels(factor(Grp_Size)), bg="white",title="Group Size" )
dev.off()

##############################################
## Final Stats vs Group Size #################

## Save Plot of Final Stats vs Group Size
png( paste(PathToSave,"HC_Final_v_GrpSize.",FILE_TAG,".png",sep=""), height=1000,width=2000,pointsize=30 )
par(mfrow=c(1,2))

## Compile Data
TIMES <- unlist(lapply( DAT, function(x) x$FIN["Secs"] ))
FINAL <- data.frame( SIZE=Grp_Size, SEC=TIMES, SEC_SAMP=TIMES/Grp_Size, SU=round(16*TIMES/3600,0), SU_SAMP=round(16*TIMES/3600,0)/Grp_Size )

## Plot Raw Numbers (Run Time vs Group Size)
 # Raw Parameters
XLIM <- range( FINAL$SIZE ) * c(1,2)
YLIM <- c( 0, max(FINAL$SU) ) * c(1,3)
COLS <- PLOT_COLS[c(2,6,8)]
 # Make Plot
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Group Size",ylab="Total CPU Hours (SU)",main=paste("Runtime vs Groupsize -",FILE_TAG),xaxt="n",yaxt="n" )
axis( 1, at=c(Grp_Size.uniq,50,100,200,300,400,437) )
axis( 2, at=seq(0,10e3,1e2) )
abline( h=seq(0,YLIM[2]+100,100), lty=3,col="grey50",lwd=1 )
abline( v=c(Grp_Size.uniq,50,100,200,300,400,437), lty=3,col="grey50",lwd=1 )
 # Plot Fits
Grp_Size.seq <- seq( 1, 5*max(Grp_Size.uniq), 1 )
MOD <- nls( SU ~ a + b^SIZE, data=FINAL, start=list(a=0,b=1.01))
MOD.q <- lm( SU ~ I(SIZE^2)+SIZE, data=FINAL )
lines( Grp_Size.seq, predict(MOD.q,list(SIZE=Grp_Size.seq)), col=COLS[1], lwd=2 )
MOD.l <- lm( SU ~ SIZE, data=FINAL )
lines( Grp_Size.seq, predict(MOD.l,list(SIZE=Grp_Size.seq)), col=COLS[2], lwd=2,lty=2)
 # Summarize Fits
ADJ.q <- round( summary(MOD.q)$adj.r.squared, 4)
ADJ.l <- round( summary(MOD.l)$adj.r.squared, 4)
text( quantile(XLIM,.01),quantile(YLIM,.76), label=paste("Quad: Adj.R2=",ADJ.q,sep=""),col=COLS[1], pos=4 )
text( quantile(XLIM,.01),quantile(YLIM,.70), label=paste("Linr: Adj.R2=",ADJ.l,sep=""),col=COLS[2], pos=4 )
 # Plot Actual Data
points( SU ~ SIZE, data=FINAL, col=COLS[3], pch="+",cex=1.5,lwd=2 )
 # Legend
legend("topleft",col=COLS[1:2],legend=c("Quad","Lin"),title="Fits",lty=c(1,2),lwd=2 )

## Plot Per Sample Numbers (Run Time vs Group Size)
 # Per Sample Parameters
# XLIM <- range( FINAL$SIZE )
YLIM <- c( 0, max(FINAL$SU_SAMP) )
 # Make Plot
plot( 0,0,type="n", xlim=XLIM,ylim=YLIM, xlab="Group Size",ylab="CPU Hours (SU) per Sample",main=paste("Per Sample Runtime vs Groupsize -",FILE_TAG),xaxt="n" )
abline( h=seq(0,YLIM[2]+20,2), lty=3,col="grey50",lwd=1 )
abline( v=c(Grp_Size.uniq,50,100,200,300,400,437), lty=3,col="grey50",lwd=1 )
axis( 1, at=c(Grp_Size.uniq,50,100,200,300,400,437) )
 # Plot Fits
MOD.q.2 <- lm( SU_SAMP ~ I(1/SIZE)+SIZE, data=FINAL )
lines( Grp_Size.seq, predict(MOD.q.2,list(SIZE=Grp_Size.seq)), col=COLS[1], lwd=2 )
MOD.l.2 <- lm( SU_SAMP ~ I(1/SIZE), data=FINAL )
lines( Grp_Size.seq, predict(MOD.l.2,list(SIZE=Grp_Size.seq)), col=COLS[2], lwd=2 )
 # Summarize Fits
ADJ.q.2 <- round( summary(MOD.q.2)$adj.r.squared, 4)
ADJ.l.2 <- round( summary(MOD.l.2)$adj.r.squared, 4)
text( quantile(XLIM,.01),quantile(YLIM,.26), label=paste("Quad: Adj.R2=",ADJ.q.2,sep=""),col=COLS[1], pos=4 )
text( quantile(XLIM,.01),quantile(YLIM,.20), label=paste("Linr: Adj.R2=",ADJ.l.2,sep=""),col=COLS[2], pos=4 )
 # Plot Actual Data
points( SU_SAMP ~ SIZE, data=FINAL, col=COLS[3], pch="+",cex=1.5,lwd=2 )

dev.off()
























#####################################################################
## END OF DOC #######################################################
#####################################################################
