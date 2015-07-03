## Plot ANC_ENV_TEST Results ##
 # Include 
## June 26, 2014 ##
## Kristopher Standish ##

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
## LOAD DATA ########################################################
#####################################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths
PathToData <- "/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/DATA/Revisions/Ancestry/"
PathToPlot <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/PLOTS/Revisions/",DATE,"_Anc/",sep="" )
dir.create( PathToPlot )

## Load Summarized Data
SNP <- read.table( paste(PathToData,"Summary_Stats_SNP_FIL.txt",sep=""), sep="\t",header=T,stringsAsFactors=F )
IND <- read.table( paste(PathToData,"Summary_Stats_IND_FIL.txt",sep=""), sep="\t",header=T,stringsAsFactors=F )
SNP$GRP[which(is.na(SNP$GRP))] <- "NAA"
IND$GRP[which(is.na(IND$GRP))] <- "NAA"
 # Calculate Negative Predictive Value
SNP$NPV <- SNP$True_N / SNP$Test_N
IND$NPV <- IND$True_N / IND$Test_N

#####################################################################
## PLOT SUMMARY STATS ###############################################
#####################################################################

#########################################################
## Sensitivity & Specificity: 2 Plots in single window

## Set universal parameters
COLS <- PLOT_COLS[c(2,6)] ; names(COLS) <- c("SNP","IND") # PLOT_COLS[c(4,6,4,6)] # PLOT_COLS[c(7,1,7,1)] # PLOT_COLS[c(7,1,4,6)]
LINE_COL <- PLOT_COLS[9]
LTYS <- c(1,1) ; names(LTYS) <- c("SNP","IND")
PCHS <- c(15:18)[factor(SNP$GRP)] # rep( 15:18, rep(4,4) ) # rep( c("*","+","o","x"),c(4,4,4,4) )
XLIM <- c(0,1)
## Open File to Save to
png(paste(PathToPlot,"4-Perc_EUR.png",sep=""),height=1000,width=2000,pointsize=36)
par(mfrow=c(1,2))
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
text( quantile(XLIM,.95),quantile(YLIM.SENS,.05),labels=paste(STAR.SNP,"p =",formatC(P_VAL.SNP, digits=2, format="e")),col=COLS["SNP"],pos=2)
if ( P_VAL.IND<0.05 & P_VAL.IND>.01 ) { STAR.IND <- "*" } ; if ( P_VAL.IND<.01 ) { STAR.IND <- "**" }
text( quantile(XLIM,.95),quantile(YLIM.SENS,.10),labels=paste(STAR.IND,"p =",formatC(P_VAL.IND, digits=2, format="e")),col=COLS["IND"],pos=2)
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
text( quantile(XLIM,.95),quantile(YLIM.SPEC,.05),labels=paste(STAR.SNP,"p =",formatC(P_VAL.SNP, digits=2, format="e")),col=COLS["SNP"],pos=2)
if ( P_VAL.IND<0.05 & P_VAL.IND>.01 ) { STAR.IND <- "*" } ; if ( P_VAL.IND<.01 ) { STAR.IND <- "**" }
text( quantile(XLIM,.95),quantile(YLIM.SPEC,.10),labels=paste(STAR.IND,"p =",formatC(P_VAL.IND, digits=2, format="e")),col=COLS["IND"],pos=2)
## Legend
legend(quantile(XLIM,0),quantile(YLIM.SPEC,1),legend=c("SNPs","Indels"),col=COLS,lty=LTYS,lwd=3,cex=.7)
legend(quantile(XLIM,0.3),quantile(YLIM.SPEC,1),legend=unique(SNP$GRP),pch=unique(PCHS),lwd=3,cex=.7,ncol=2)
dev.off()

# text(.2,min(YLIM)+.12-.03*i, labels=paste("R2=",round(summary(MOD)$r.squared,4),";","p=",P_VAL),col=COLS[i])
# if (P_VAL<0.05) { text(.05,min(YLIM)+.12-.03*i, labels="*",col=COLS[i]) }


#########################################################
## Accuracy & PPV & NPV: 2 Plots in single window

## Set universal parameters
COLS <- PLOT_COLS[c(2,6)] ; names(COLS) <- c("SNP","IND") # PLOT_COLS[c(4,6,4,6)] # PLOT_COLS[c(7,1,7,1)] # PLOT_COLS[c(7,1,4,6)]
LINE_COL <- PLOT_COLS[9]
LTYS <- c(1,1) ; names(LTYS) <- c("SNP","IND")
PCHS <- c(15:18)[factor(SNP$GRP)] # rep( 15:18, rep(4,4) ) # rep( c("*","+","o","x"),c(4,4,4,4) )
XLIM <- c(0,1)
## Open File to Save to
png(paste(PathToPlot,"4-Perc_EUR_AccuPPV.png",sep=""),height=1000,width=3000,pointsize=36)
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















#######################################
## For UNfiltered Sets
YLIM <- c(0.75,1.25)
WHICH_SETS <- 1:2 # 1:4

png(paste(PathToPlot,"4-Perc_EUR_Unf.png",sep=""),height=1000,width=2000,pointsize=36)
par(mfrow=c(1,2))
# Misses
plot(0,0,type="n",xlim=c(0,1),ylim=YLIM,main="Misses vs Percent European Group",xlab="% European in Group",ylab="Normalized Misses")
abline( h=seq(0,2,.1), lty=2,col=PLOT_COLS[9],lwd=1)
for (i in WHICH_SETS) {
points(DAT[[i]]$PEUR,DAT[[i]]$MS,col=COLS[i],pch=PCH)
MOD <- lm(DAT[[i]]$MS ~ DAT[[i]]$PEUR)
abline(MOD, col=COLS[i], lty=LTYS[i], lwd=3)
P_VAL <- summary(MOD)$coefficients[2,4]
P_VAL_SCI <- formatC(P_VAL, digits=2, format="e") #summary(MOD)$coefficients[2,4]
STAR <- " "
if ( P_VAL<0.05 & P_VAL>.01 ) { STAR <- "*" }
if ( P_VAL<.01 ) { STAR <- "**" }
text(0,min(YLIM)+.12-.03*i, labels=paste(STAR," R2=",round(summary(MOD)$r.squared,3),";","p=",P_VAL_SCI),col=COLS[i],pos=4)
# text(.2,min(YLIM)+.12-.03*i, labels=paste("R2=",round(summary(MOD)$r.squared,4),";","p=",P_VAL),col=COLS[i])
# if (P_VAL<0.05) { text(.05,min(YLIM)+.12-.03*i, labels="*",col=COLS[i]) }
}
legend(0.45,max(YLIM),legend=c("Unfiltered SNPs","Unfiltered Indels","Filtered SNPs","Filtered Indels")[1:2],col=COLS[1:2],lty=LTYS[1:2],lwd=3) #,pch=20)
# FPs
plot(0,0,type="n",xlim=c(0,1),ylim=YLIM,main="False Positives vs Percent European Group",xlab="% European in Group",ylab="Normalized FPs")
abline( h=seq(0,2,.1), lty=2,col=PLOT_COLS[9],lwd=1)
for (i in WHICH_SETS) {
points(DAT[[i]]$PEUR,DAT[[i]]$FP,col=COLS[i],pch=PCH)
MOD <- lm(DAT[[i]]$FP ~ DAT[[i]]$PEUR)
abline(MOD, col=COLS[i], lty=LTYS[i], lwd=3)
P_VAL <- summary(MOD)$coefficients[2,4]
P_VAL_SCI <- formatC(P_VAL, digits=2, format="e")
# P_VAL <- summary(MOD)$coefficients[2,4]
STAR <- " "
if ( P_VAL<0.05 & P_VAL>.01 ) { STAR <- "*" }
if ( P_VAL<.01 ) { STAR <- "**" }
text(0,min(YLIM)+.12-.03*i, labels=paste(STAR," R2=",round(summary(MOD)$r.squared,3),";","p=",P_VAL_SCI),col=COLS[i],pos=4)
}
# legend(0.7,max(YLIM)-.05,legend=c("Unfiltered SNPs","Unfiltered Indels","Filtered SNPs","Filtered Indels"),col=COLS,lty=1,pch=20)
dev.off()

############################################
## CALCULATE CONTINGENCY TABLES ############
############################################

## Gold Standard vs My Calls
 # Calculate TP, FP, TN, FN
 # Sensitivity = TP / ( TP + FN ) = TP / Actual_P
 # Specificity = TN / ( FP + TN ) = TN / Actual_N


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

## Replot Figure w/ Variant Counts for BigCall Paper ##
## July 1, 2014 ##
## Kristopher Standish ##


## Load Data
DAT <- read.table("/Users/kstandis/Dropbox/Schork/JNJ11/Slides/20131029 - JnJ Mtg/Filtered_Variants.txt",header=T,colClasses=c("numeric","character"))

## Set Parameters
COLS <- PLOT_COLS[c(5,7)] # c("cadetblue2","mediumpurple2")
Split1 <- sapply(strsplit(DAT$FILE,"/"),"[",4)
GRP_NMS <- sapply(strsplit(Split1,".",fixed=T),"[",1)

PLOT_ARR <- array(c(DAT[c(seq(1,38,2),seq(2,38,2)),1]),dim=c(19,2))
rownames(PLOT_ARR) <- GRP_NMS[seq(1,38,2)]
colnames(PLOT_ARR) <- c("Indel","SNP")

PLOT_ARR <- PLOT_ARR[c(1:4,15:19,5:14),]

##
png(paste(PathToPlot,"Var_Counts.png",sep=""), width=1600, height=1000, pointsize=30)
xvals <- barplot(t(PLOT_ARR), main="Filtered Variant Calls by Group", ylab="Variant Positions (Millions)", xlab="Group",col=COLS, beside=T, xaxt="n", width=0.8, space=c(0,.25), yaxt="n")
abline( h=seq(0,2e7,1e6), col=PLOT_COLS[9],lty=2,lwd=1)
# axis( 2, at=seq(0,2e7,2e6), labels=formatC(seq(0,2e7,2e6),format="e",digits=1), las=2, cex=.8 )
axis( 2, at=seq(0,2e7,2e6), labels=seq(0,20,2), las=2, cex=.8 )
# text(cex=.9, pos=2, x=colMeans(xvals)+1, y=-.035*max(PLOT_ARR), labels=rownames(PLOT_ARR), xpd=T, srt=35)
text(cex=.9, pos=2, x=colMeans(xvals)+1, y=-.035*max(PLOT_ARR), labels=paste("Group",1:19), xpd=T, srt=35)
xvals <- barplot(t(PLOT_ARR), col=COLS, beside=T, xaxt="n", width=0.8, space=c(0,.25), add=T, yaxt="n")
dev.off()
#axis(1,at=1:19-.5,labels=rownames(PLOT_ARR), srt=45)

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

## Re-plot Ancestral Estimates for Cohort ##
## April 7, 2014 ##
## Kristopher Standish ##

PathToSave <- "/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/PLOTS/20141113/"

##################################################
## LOAD DATA #####################################
##################################################

DR <- read.table("/Users/kstandis/Data/Burn/DRSA436EE.txt",sep="\t",header=T)

GROUP_NAMES <- as.character(read.table("/Users/kstandis/Data/Burn/GROUP_NAMES_ALL.txt")[,1])
LISTS <- list()
for (i in 1:19) {
	PathToFile <- paste("/Users/kstandis/Data/Burn/GRP_NMS/",GROUP_NAMES[i],sep="")
	LISTS[[i]] <- as.character(read.table(PathToFile)[,1])
}
names(LISTS) <- GROUP_NAMES
ORD_LIST <- LISTS[c(3,12:14,7:11,15:19,1:2,4:6)]

##################################################
## PLOT ALL SAMPLES ##############################
##################################################

## Pull out Ondrej's Ancestry Columns
ANC <- DR[,c(14:19)]
rownames(ANC) <- DR$ID

## Sort by eastAsian and AmerInd and Eur
ANC <- ANC[order(ANC$eastAsian),]
ANC <- ANC[order(ANC$amerind),]
ANC <- ANC[order(ANC$europe),]

## Plot
COLS <- c("aquamarine3", "coral2", "darkgoldenrod2", "mediumorchid2", "cadetblue2", "slateblue1")
COLS <- PLOT_COLS[c(4,2,3,7,5,6)]
png(paste(PathToSave,"Admix_Full.png",sep=""), width=2700, height=1000, pointsize=40)
par(mar=c(2,5,4,0))
barplot(t(ANC), col=COLS, border="white", width=1, space=0, main="Admixture Estimates for 437 Patients", xlab="",ylab="% Admixture",legend=F, xaxt="n") #  xlab="Patient", 
legend( 340, .98, legend=c("Europe","Africa","Amer. Ind.","East Asian","Oceanic","Central Asian"), fill=COLS)
text( 220, -.07, labels="Patients", xpd=T )
dev.off()

##################################################
## PLOT GROUPS ###################################
##################################################

png(paste(PathToSave,"Admix_Group.png",sep=""), width=2400, height=1200, pointsize=30)
par(mfrow=c(4,5))
par(mar=c(2,5,2,0))
loop <- 0
for (LIST in ORD_LIST) {
	loop <- loop + 1
	DATA <- ANC[which(rownames(ANC) %in% LIST),]
	barplot(t(DATA), col=COLS, border="white", width=1, space=.05, main=paste("Group",loop),ylab="% Admixture", names.arg=rep("",nrow(DATA)) ) # ,las=2)
	text( 12, -.15, labels="Patients", xpd=T )
}
# Make Key for the last square
plot(0,0,type="n", xaxt="n",yaxt="n", xlab="",ylab="",xlim=c(0,1), ylim=c(0,1))
legend(0.2,1, legend=colnames(ANC), fill=COLS, title="Ancestry", cex=1)
dev.off()




#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################


########## DO THIS ON TSCC ##################

## Make plot of Ancestry Groups for ANC_ENV_TEST (32 grps - 4chr x 4grp x 2rep) ##
## May 01, 2014 ##
## Kristopher Standish ##


############################################################
## LOAD DATA TABLE #########################################
############################################################

TAB <- read.table("/Users/kstandis/Data/Burn/DRSA436EE.txt",sep="\t",header=T)

CHROMS <- c(17,14,22,20)

# PathToPlot <- "/projects/janssen/Ancestry/ANC_ENV_TEST/PLOTS/20140501_Anc_Grps/"

############################################################
## LOAD SAMPLE LISTS #######################################
############################################################

SAMP_LIST_NAMES <- c("EUR0b","EUR1","EUR1b","EUR2","EUR2b","EUR3","EUR3b","AN0b","AN1","AN1b","AN2","AN2b","AN3","AN3b","MIX0b","MIX1","MIX1b","MIX2","MIX2b","MIX3","MIX3b","ASN0b","ASN1","ASN1b","ASN2","ASN2b","ASN3","ASN3b")
SAMP_LISTS <- list()
for (samp in 1:length(SAMP_LIST_NAMES)) {
	PATH <- paste("/Users/kstandis/Data/Burn/ANC_ENV_GRPS/",SAMP_LIST_NAMES[samp],".list",sep="")
	SAMP_LISTS[[samp]] <- as.character(read.table(PATH,sep="")[,1])
}
names(SAMP_LISTS) <- SAMP_LIST_NAMES
SAMP_LISTS$EUR0 <- c("X061041","Q665950","X040662","P179281","X039092","X148910","X478693","U270591","X155365","Q662090")
SAMP_LISTS$ASN0 <- c("X291529","Q633972","Q428373","Q428350","Q428901","Q428900","Q429084","Q429093","V695532","V695540")
SAMP_LISTS$AN0 <- c("G628582","F537279","L687135","M597358","T984023","T984100","S048822","T984080","M613448","L687157")
SAMP_LISTS$MIX0 <- c("Y484678","Q429100","X291529","X061041","F537279","T984100","X040162","S048852","L687168","X842080")

############################################################
## PLOT ####################################################
############################################################

## Set up some parameters
SAMP_LIST_NAMES_PT1 <- SAMP_LIST_NAMES_ALL[seq(1,31,2)]
COLS <- c("aquamarine3", "coral2", "darkgoldenrod2", "mediumorchid2", "cadetblue2", "slateblue1")
COLS <- PLOT_COLS[c(4,2,3,7,5,6)]
LEG <- rep(F,length(SAMP_LISTS)) # c(T,rep(F,length(SAMP_LISTS)-1))
SAMP_LIST_NAMES_ALL <- c("EUR0","EUR0b","EUR1","EUR1b","EUR2","EUR2b","EUR3","EUR3b","AN0","AN0b","AN1","AN1b","AN2","AN2b","AN3","AN3b","MIX0","MIX0b","MIX1","MIX1b","MIX2","MIX2b","MIX3","MIX3b","ASN0","ASN0b","ASN1","ASN1b","ASN2","ASN2b","ASN3","ASN3b")


png(paste(PathToPlot,"Anc_Grps.png",sep=""),width=2400,height=1600,pointsize=36)
par(mfrow=c(4,4))
par(mar=c(2.5,5,2.5,.5)) # rep(2.5,4))
for (i in 1:16) {
	# samp <- SAMP_LIST_NAMES_PT1[i]
	samp <- SAMP_LIST_NAMES_PT1[i]
	DATA <- TAB[which(TAB$ID %in% SAMP_LISTS[[samp]]),14:19]
	barplot(t(DATA), col=COLS, border="white", width=1, space=0, main=paste("Ancestry - Chr",rep(CHROMS,4)[i],"-",SAMP_LIST_NAMES_PT1[i]), ylab="% Admixture", xlab="", legend=LEG[i], xaxt="n")	
	text( 5, -.15, labels="Patients", xpd=T )
}
dev.off()


# ## 2 Separate Plots (half each)
# png(paste(PathToPlot,"Anc_Grps_1.png",sep=""),width=2400,height=2000)
# par(mfrow=c(4,4))
# for (i in 1:16) {
# 	samp <- SAMP_LIST_NAMES_ALL[i]
# 	DATA <- TAB[which(TAB$ID %in% SAMP_LISTS[[samp]]),14:19]
# 	barplot(t(DATA), col=COLS, border="white", width=1, space=0, main=paste("Ancestry -",SAMP_LIST_NAMES[i]), ylab="% Admixture", xlab="Patient", legend=LEG[i], xaxt="n")	
# }
# dev.off()
# png(paste(PathToPlot,"Anc_Grps_2.png",sep=""),width=2400,height=2000)
# par(mfrow=c(4,4))
# for (i in 17:32) {
# 	samp <- SAMP_LIST_NAMES_ALL[i]
# 	DATA <- TAB[which(TAB$ID %in% SAMP_LISTS[[samp]]),14:19]
# 	barplot(t(DATA), col=COLS, border="white", width=1, space=0, main=paste("Ancestry -",SAMP_LIST_NAMES[i]), ylab="% Admixture", xlab="Patient", legend=LEG[i], xaxt="n")	
# }
# dev.off()

# ## Plot all in 1 Plot
# png(paste(PathToPlot,"Anc_Grps.png",sep=""),width=3000,height=1500)
# par(mfrow=c(4,8))
# for (i in 1:32) {
# 	samp <- SAMP_LIST_NAMES_ALL[i]
# 	DATA <- TAB[which(TAB$ID %in% SAMP_LISTS[[samp]]),14:19]
# 	barplot(t(DATA), col=COLS, border="white", width=1, space=0, main=paste("Ancestry -",SAMP_LIST_NAMES[i]), ylab="% Admixture", xlab="Patient", legend=LEG[i], xaxt="n")	
# }
# dev.off()

# ############################################################
# ## PLOT TOTAL ADMIXTURE FOR GROUP ##########################
# ############################################################

# SUMS <- array(,dim=c(32,6))
# rownames(SUMS) <- SAMP_LIST_NAMES_ALL
# colnames(SUMS) <- names(TAB[,14:19])
# for (i in 1:32) {
# 	samp <- SAMP_LIST_NAMES_ALL[i]
# 	DATA <- TAB[which(TAB$ID %in% SAMP_LISTS[[samp]]),14:19]
# 	SUMS[i,] <- colSums(DATA)
# }

# png(paste(PathToPlot,"Total_Anc_Grps.png",sep=""),width=1600,height=800)
# barplot(t(SUMS), col=COLS, border="white", width=1, space=0, main="Total Ancestry By Group", ylab="% Admixture", xlab="Group", legend=T, las=2)	
# dev.off()

# ############################################################
# ## PLOT FIRST 8 GROUPS ONLY ################################
# ############################################################
# SAMP_LIST_NAMES_PT1 <- SAMP_LIST_NAMES_ALL[seq(1,31,2)]
# SUMS_1 <- array(,dim=c(16,6))
# rownames(SUMS_1) <- SAMP_LIST_NAMES_PT1
# colnames(SUMS_1) <- names(TAB[,14:19])
# for (i in 1:16) {
# 	samp <- SAMP_LIST_NAMES_PT1[i]
# 	DATA <- TAB[which(TAB$ID %in% SAMP_LISTS[[samp]]),14:19]
# 	SUMS_1[i,] <- colSums(DATA)
# }

# jpeg(paste(PathToPlot,"Total_Anc_Grps_Pt1.jpeg",sep=""),width=1000,height=800)
# barplot(t(SUMS_1), col=COLS, border="white", width=1, space=0, main="Total Ancestry By Group", ylab="% Admixture", xlab="Group", legend=T, las=2)	
# dev.off()

# write.table(SUMS_1, "/projects/janssen/Ancestry/ANC_ENV_TEST/Group_Admix_Perc.txt",sep="\t",row.names=T,col.names=T,quote=F)








