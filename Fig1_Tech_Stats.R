## Compile Usage Stats for JnJ ##
 # Taken from "Slides/20131001 - JnJ_Report/20131001-JnJ_Report_Rscript.txt"
 # Script based on "Big_Call_Paper/Code/20141113/20141113_1_Tech_Plots.R"
## June 29, 2014 ##
## Kristopher Standish ##

## Compile Table of Statistics for Variant Calling Pipeline
 # 

PLOT_COLS <- c("firebrick2","sienna2","gold2","chartreuse3","cadetblue2","steelblue3","slateblue3","black","grey50","white")
PLOT_COLS.4 <- c("firebrick4","sienna4","gold4","chartreuse4","cadetblue4","steelblue4","slateblue4","black","grey50","white")

#####################################################################
## LOAD DATA ########################################################
#####################################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Path To Data
PathToData <- "/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/DATA/Updates/"
PathToSave <- paste("/Users/kstandis/Dropbox/Schork/JNJ11/Big_Call_Paper/PLOTS/Revisions/",DATE,"_Tech/",sep="" )
dir.create( PathToSave )

## Load Job Info
JB <- read.table(paste(PathToData,"Update_Jobs.txt",sep=""),sep="\t",header=T, stringsAsFactors=F)
# Jobs <- read.table(paste(PathToData,"Update_Jobs.txt",sep=""), sep="\t", header=T, stringsAsFactors=F)

## Load File Info
FL1 <- read.table(paste(PathToData,"Update_Files_ALL_2013-09-12",sep=""),sep="\t",header=T)
FL2 <- read.table(paste(PathToData,"Update_Files_ALL_2013-09-20",sep=""),sep="\t",header=T)
fastq <- read.table(paste(PathToData,"Update_fastq",sep=""), sep="\t", header=T, stringsAsFactors=F)
FILES_914 <- read.table(paste(PathToData,"BackUp/9_14_2_Pre-Delete/Update_Files",sep=""), sep="\t", header=T, stringsAsFactors=F)
FILES_930 <- read.table(paste(PathToData,"BackUp/10_1_1_Pre-Delete/Update_Files",sep=""), sep="\t", header=T, stringsAsFactors=F)
FILES_FIN <- read.table(paste(PathToData,"Update_Files.txt",sep=""), sep="\t", header=T, stringsAsFactors=F)
colnames(FILES_930) <- colnames(FILES_914) <- colnames(FILES_FIN) <- c("FILE_NAME", "f_DATE", "SIZE")

## Load Sample Names
SAMPLE_NAMES <- read.table(paste(PathToData,"SAMPLE.list",sep=""), sep="\t", header=F, stringsAsFactors=F)[,1]

## Specify Step Names
STEP_NAMES <- c("Map", "Bam", "Merge", "Sort", "Mrk_Dps", "Trgt_Crt", "Indl_Rlgn", "Bs_Rcl", "Prnt_Rds")
N.Steps <- length(STEP_NAMES)
STEP_SUFX <- c(".sam",".bam","_m.bam","_s.bam","_md.bam",".intervals","_real.bam","_recal.table","_recal.bam")

#####################################################################
## MERGE/FILTER/ORGANIZE TABLES #####################################
#####################################################################

dim(JB)
# dim(Jobs)
dim(fastq)
dim(FILES_914)
dim(FILES_930)
dim(FL1)
dim(FL2)

#####################################
## Reorganize Files Table ###########

## Throw All Output Files into single data frame (most recent first)
FILES_dups <- rbind( FILES_FIN, FILES_930, FILES_914 )
 # Get rid of duplicated Files
FILES <- FILES_dups[which(!duplicated(FILES_dups$FILE_NAME)),]

## Separate Files by Step
Files <- list()
 # FastQ's
Files$FQ <- fastq[,c("FQ_FILE","FILE_TIME_F","FILE_SIZE_F")] ; colnames(Files$FQ) <- c("FILE_NAME","f_DATE","SIZE")
 # Output Files
for ( s in 1:N.Steps ) {
  step <- STEP_NAMES[s]
  suffix <- STEP_SUFX[s]
  if ( step != "Bam" ) {
    Files[[step]] <- FILES[grep(suffix, FILES$FILE_NAME),]
  }else{
    Files[[step]] <- FILES[intersect(grep("_L", FILES$FILE_NAME),grep(suffix, FILES$FILE_NAME)),]
  }
} ; unlist(lapply(Files,nrow))

## Remove Certain Sample Files????
# Files$Map <- FILES[grep(".sam", FILES$FILE_NAME),]
# Remove_from_Map <- c(grep("Q665673-27_I224_L1.sam", Files$Map$FILE_NAME), grep("T988163-27_I154_L5.sam", Files$Map$FILE_NAME), grep("X040172-27_I649_L3.sam", Files$Map$FILE_NAME))
# Remove_from_Map <- grep("Q665673-27_I224_L1.sam|T988163-27_I154_L5.sam|X040172-27_I649_L3.sam", Files$Map$FILE_NAME)
# Remove_from_Map <- grep("Q665673-27_I224_L1|T988163-27_I154_L5|X040172-27_I649_L3", Files$Map$FILE_NAME)
# Files$Map <- Files$Map[setdiff(1:nrow(Files$Map),Remove_from_Map),]

# Files$Bam <- FILES[intersect(grep("_L", FILES$FILE_NAME),grep(".bam", FILES$FILE_NAME)),]
# Remove_from_Bam <- c(grep("Q665673-27_I224_L1.bam", Files$Bam$FILE_NAME), grep("T988163-27_I154_L5.bam", Files$Bam$FILE_NAME), grep("X040172-27_I649_L3.bam", Files$Bam$FILE_NAME))
# Files$Bam <- Files$Bam[setdiff(1:nrow(Files$Bam),Remove_from_Bam),]

## Names of all Lanes/Bams
# Names_Mapped <- gsub(".sam", "", Files$Map$FILE_NAME)
# Names_Bammed <- gsub(".bam", "", Files$Bam$FILE_NAME)
# setdiff( Names_Mapped, Names_Bammed )
# setdiff( Names_Bammed, Names_Mapped )

#####################################
## Reorganize Jobs Table ############

## Remove Duplicated Files that also have an Error
Error_Files <- which( JB$ERR_FLAG=="Error" )
Dup_Files <- which(duplicated( JB$OUTPUT_FILE ))
Dup_Filenames <- JB$OUTPUT_FILE[ Dup_Files ]
Remove_Files <- which( JB$OUTPUT_FILE %in% Dup_Filenames & JB$ERR_FLAG=="Error" )
JB <- JB[ -Remove_Files, ]

## Remove Jobs that Times Out
Timed_Out <- which( JB$EXIT_STATUS=="Timed_Out" )
JB <- JB[ -Timed_Out, ]

## Check out & Remove Duplicates
Dup_Files <- which(duplicated( JB$OUTPUT_FILE ))
Dup_Filenames <- JB$OUTPUT_FILE[ Dup_Files ]
Remove_Files <- intersect( which( JB$STEP=="Merge" ), Dup_Files )
JB <- JB[ -Remove_Files, ]

# ## Remove Incomplete or Error Files
# MD_Error_Files <- as.character( JB$OUTPUT_FILE[ which( JB$STEP=="MrkDps" & JB$ERR_FLAG=="Error") ] )
# MD_Dup_Files <- as.character( JB$OUTPUT_FILE[ which(duplicated( JB$OUTPUT_FILE)) ] )
# MD_Error_Dup_Files <- intersect( MD_Error_Files, MD_Dup_Files )
# MD_Remove_Files <- which( JB$OUTPUT_FILE %in% MD_Error_Dup_Files & JB$ERR_FLAG=="Error" )
# JB <- JB[ -MD_Remove_Files, ]

## Reformat File Names in Jobs Table
 # Remove "O=" from Filename for JB$STEP=="MrkDups"
JB$OUTPUT_FILE <- gsub( "O=","", JB$OUTPUT_FILE )
 # Add ".bam" to several Sort Output Files
JB$OUTPUT_FILE[ grep("_s$",JB$OUTPUT_FILE) ] <- paste( grep("_s$",JB$OUTPUT_FILE,value=T),".bam",sep="" )

#####################################
## Merge Tables #####################

MG <- merge(x=JB,y=FILES,by.x="OUTPUT_FILE",by.y="FILE_NAME")
MG <- merge(x=JB,y=FILES,by.x="OUTPUT_FILE",by.y="FILE_NAME", all.x=T)

#####################################################################
## COMPILE STATS PER STEP ###########################################
#####################################################################

#####################################
## CALCULATE TIME PER JOB ###########

## Calculate SUs
Split <- strsplit(as.character(MG$WALL_TIME), ":")
Times <- t(sapply(Split,"[",1:3))
Secs <- 3600*as.numeric(Times[,1]) + 60*as.numeric(Times[,2]) + as.numeric(Times[,3])
SUs <- round(Secs*MG$CORES/3600,0)
SUs[which(SUs==0)] <- 1
 # Add to Table
MG <- data.frame( MG, SUs, Wall=Secs )

## Calculate CPU_Time
Split <- strsplit(as.character(MG$CPU_TIME), ":")
Times <- t(sapply(Split,"[",1:3))
Secs <- 3600*as.numeric(Times[,1]) + 60*as.numeric(Times[,2]) + as.numeric(Times[,3])
 # Add to Table
MG <- data.frame( MG, CPU=Secs )

#####################################
## PARSE JOB METRICS ################

## Strip MEM/VMEM of "kb" tag
Strip.MEM <- as.numeric( gsub("kb","",MG$MEM) )
Strip.VMEM <- as.numeric( gsub("kb","",MG$VMEM) )

## Number of Commands per Job
COMMANDS_SPLIT2 <- strsplit(as.character(MG$JOB_NAME), "_")
COMMANDS_SPLIT <- as.numeric(sapply(COMMANDS_SPLIT2,"[",c(2)))
COMMANDS_SPLIT[MG$JOB_NAME=="Bam_Manual_1_20130919"] <- 6
COMMANDS_SPLIT[MG$JOB_NAME=="Bam_Manual_2_20130919"] <- 5
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_5_Test16"] <- 1
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_5b_Test15"] <- 15
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_r3_16a"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_r3_16b"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_6_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_7_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_8_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_9_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_10_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_11_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_12_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_13_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_14_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_15_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_16_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_17_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_18_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_19_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_20_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_21_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_22_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_23_16"] <- 16
COMMANDS_SPLIT[MG$JOB_NAME=="Sort_24_6"] <- 6

 # Add to Table
MG <- data.frame( MG, MEM.kb=Strip.MEM, VMEM.kb=Strip.VMEM, COMMANDS_per_JOB=COMMANDS_SPLIT,
  CORES_per_COMMAND=MG$CORES/COMMANDS_SPLIT, SUs_per_COMMAND=SUs/COMMANDS_SPLIT, CPUtime_per_COMMAND=MG$CPU/COMMANDS_SPLIT,
  VMEM_per_COMMAND=Strip.VMEM/COMMANDS_SPLIT, MEM_per_COMMAND=Strip.MEM/COMMANDS_SPLIT )

#####################################
## COMPILE STATS by SAMPLE ##########
Stats_by_Sample <- c("SUs_per_COMMAND","CPUtime_per_COMMAND","Wall","CPU","SIZE")
frame <- MG[ ,Stats_by_Sample ]
By_Sample <- aggregate( frame, list( STEP=MG$STEP,SAMPLE=MG$SAMPLE), sum )
names(By_Sample)

#####################################
## COMPILE STATS ####################
# Stats_to_Compile <- c("SUs","Wall","CPU","SIZE","CORES_per_COMMAND","COMMANDS")
By_Step <- list()
 # Commands/Sample by Step
Commands_Samp_Tab <- table( MG$SAMPLE, MG$STEP )
By_Step$Command_per_Sample <- data.frame( STEP=colnames(Commands_Samp_Tab), MN=colMeans(Commands_Samp_Tab), SD=apply(Commands_Samp_Tab,2,sd) )
 # Cores/Command by Step
By_Step$Cores_per_Command <- merge( aggregate(CORES_per_COMMAND ~ STEP,data=MG,mean), aggregate(CORES_per_COMMAND ~ STEP,data=MG,sd), by="STEP")
colnames(By_Step$Cores_per_Command) <- c("STEP","MN","SD")
 # Commands/Node by Step
By_Step$Command_per_Node <- merge( aggregate(COMMANDS_per_JOB ~ STEP,data=MG,mean), aggregate(COMMANDS_per_JOB ~ STEP,data=MG,sd), by="STEP")
colnames(By_Step$Command_per_Node) <- c("STEP","MN","SD")
 # Wall_Time/Sample by Step
By_Step$Walltime_per_Sample <- merge( aggregate(Wall ~ STEP,data=By_Sample,mean), aggregate(Wall ~ STEP,data=By_Sample,sd), by="STEP")
By_Step$Walltime_per_Sample[,2:3] <- By_Step$Walltime_per_Sample[,2:3] / 3600
colnames(By_Step$Walltime_per_Sample) <- c("STEP","MN","SD")
 # CPU_Time/Sample by Step
By_Step$CPUtime_per_Sample <- merge( aggregate(CPUtime_per_COMMAND ~ STEP,data=By_Sample,mean), aggregate(CPUtime_per_COMMAND ~ STEP,data=By_Sample,sd), by="STEP")
By_Step$CPUtime_per_Sample[,2:3] <- By_Step$CPUtime_per_Sample[,2:3] / 3600
colnames(By_Step$CPUtime_per_Sample) <- c("STEP","MN","SD")
 # SUs/Sample by Step
By_Step$SUs_per_Sample <- merge( aggregate(SUs_per_COMMAND ~ STEP,data=By_Sample,mean), aggregate(SUs_per_COMMAND ~ STEP,data=By_Sample,sd), by="STEP")
colnames(By_Step$SUs_per_Sample) <- c("STEP","MN","SD")
 # MEM/Command by Step
By_Step$MEM_per_Command <- merge( aggregate(MEM_per_COMMAND ~ STEP,data=MG,mean), aggregate(MEM_per_COMMAND ~ STEP,data=MG,sd), by="STEP")
By_Step$MEM_per_Command[,2:3] <- By_Step$MEM_per_Command[,2:3] / 1e6
colnames(By_Step$MEM_per_Command) <- c("STEP","MN","SD")
 # VMEM/Command by Step
By_Step$VMEM_per_Command <- merge( aggregate(VMEM_per_COMMAND ~ STEP,data=MG,mean), aggregate(VMEM_per_COMMAND ~ STEP,data=MG,sd), by="STEP")
By_Step$VMEM_per_Command[,2:3] <- By_Step$VMEM_per_Command[,2:3] / 1e6
colnames(By_Step$VMEM_per_Command) <- c("STEP","MN","SD")
 # File_Size/Sample by Step
By_Step$FileSize_per_Sample <- merge( aggregate(SIZE ~ STEP,data=By_Sample,mean), aggregate(SIZE ~ STEP,data=By_Sample,sd), by="STEP")
By_Step$FileSize_per_Sample[,2:3] <- By_Step$FileSize_per_Sample[,2:3] / 1e9
colnames(By_Step$FileSize_per_Sample) <- c("STEP","MN","SD")
## Re-order Steps
By_Step <- lapply( By_Step, function(x) x[c(4,1,5,8,6,9,3,2,7),] )

#####################################
## PUT INTO SAMPLE ##################
table_cols <- c("Order","Step","Tool","Commands/Sample","Cores/Command","Commands/16cores","Wall_Time/Sample","CPU_Time/Sample","SUs/Sample","MEM/Command","VMEM/Command","Output_File/Sample")
TABLE <- array( ,c( length(STEP_NAMES)+1,length(table_cols) ) )
colnames(TABLE) <- table_cols ; rownames(TABLE) <- c("FQ",STEP_NAMES)
 # Steps
TABLE[,"Order"] <- 1:nrow(TABLE)-1
TABLE[,"Step"] <- rownames(TABLE)
TABLE[,"Tool"] <- c( "-","BWA","Samtools","Samtools","Samtools","PicardTools","GATK","GATK","GATK","GATK" )
 # Cores/Commands
Print_Jobs <- paste( round(By_Step$Command_per_Sample$MN,2), "+/-", round(By_Step$Command_per_Sample$SD,2), sep="" )
TABLE[,"Commands/Sample"] <- c( NA, Print_Jobs ) ; TABLE[,"Commands/Sample"][4:10] <- 1
# TABLE[,"Cores/Command"] <- c( NA,8,1,1,1,2,2,3,8,8 )
Print_Cores <- paste( round(By_Step$Cores_per_Command$MN,2), "+/-", round(By_Step$Cores_per_Command$SD,2), sep="" )
TABLE[,"Cores/Command"] <- c( NA, Print_Cores )
Print_perNode <- paste( round(By_Step$Command_per_Node$MN,2), "+/-", round(By_Step$Command_per_Node$SD,2), sep="" )
TABLE[,"Commands/16cores"] <- c( NA, Print_perNode )
 # Time/SUs
Print_Wall_Time <- paste( round(By_Step$Walltime_per_Sample$MN,2), "+/-", round(By_Step$Walltime_per_Sample$SD,2), sep="" )
TABLE[,"Wall_Time/Sample"] <- c(NA, Print_Wall_Time)
Print_CPU_Time <- paste( round(By_Step$CPUtime_per_Sample$MN,2), "+/-", round(By_Step$CPUtime_per_Sample$SD,2), sep="" )
TABLE[,"CPU_Time/Sample"] <- c(NA, Print_CPU_Time)
Print_SUs <- paste( round(By_Step$SUs_per_Sample$MN,2), "+/-", round(By_Step$SUs_per_Sample$SD,2), sep="" )
TABLE[,"SUs/Sample"] <- c(NA, Print_SUs)
 # Memory
Print_MEM <- paste( round(By_Step$MEM_per_Command$MN,2), "+/-", round(By_Step$MEM_per_Command$SD,2), sep="" )
TABLE[,"MEM/Command"] <- c( NA, Print_MEM )
Print_VMEM <- paste( round(By_Step$VMEM_per_Command$MN,2), "+/-", round(By_Step$VMEM_per_Command$SD,2), sep="" )
TABLE[,"VMEM/Command"] <- c( NA, Print_VMEM )
 # File Size
Print_Size <- paste( round(By_Step$FileSize_per_Sample$MN,2), "+/-", round(By_Step$FileSize_per_Sample$SD,2), sep="" )
TABLE[,"Output_File/Sample"] <- c(NA, Print_Size )
FQ_SIZE <- aggregate( as.numeric(fastq[,"FILE_SIZE_F"]), list( SAMPLE=fastq$SAMPLE_F), sum )[,2] / 1e9
TABLE["FQ","Output_File/Sample"] <- paste( round(mean(FQ_SIZE),2), "+/-", round(sd(FQ_SIZE),2) )

#############################################
## FIGURE 1B - MAPPING EFFICIENCY ###########
#############################################

MAP <- MG[ MG$STEP=="Map", ]
COLORS <- PLOT_COLS[c(4,6)] # c("slateblue3","tomato2")
COLORS.2 <- PLOT_COLS.4[c(4,6)]
CORE_COLS <- as.numeric(factor(MAP$CORES)) # ; names(CORE_COLS) <- c("8","16")
 # Calculate fit for 8 and 16 cores
C.8 <- which( MAP$CORES==8 )
MOD.8 <- lm( MAP$SUs[C.8] ~ I(MAP$SIZE[C.8]/1e9) ) ; COEF.8 <- round( coef(MOD.8)[2], 3 )
MOD.16 <- lm( MAP$SUs[-C.8] ~ I(MAP$SIZE[-C.8]/1e9) ) ; COEF.16 <- round( coef(MOD.16)[2], 3 )
XLIM <- c(0,max(MAP$SIZE/1e9,na.rm=T))
YLIM <- c(0,max(MAP$SUs,na.rm=T))
png(paste(PathToSave,"1-Mapping_Efficiency.png",sep=""),width=1000,height=1000, pointsize=30)
plot(0,0, type="n", xlim=XLIM, ylim=YLIM, xlab="Sam File Size (GB)", ylab="SUs (Cores-Hours)", main="Read Mapping Efficiency")
abline(h=seq(0,200,20), lty=2, ,lwd=1, col=PLOT_COLS[9] )
abline(v=seq(0,200,20), lty=2, ,lwd=1, col=PLOT_COLS[9] )
points( MAP$SIZE/1e9, MAP$SUs, pch="+", col=COLORS[CORE_COLS] )
abline( MOD.8, col=COLORS.2[1], lty=2, lwd=6 )
abline( MOD.16, col=COLORS.2[2], lty=2, lwd=6 )
legend(20,80,legend=c("8","16"),title="# Cores",pch="+",col=COLORS)
text( 0, c(55,50), labels=paste(c(8,16),"Cores:",c(COEF.8,COEF.16),"SU/GB"), col=COLORS.2, pos=4 )
dev.off()





###############################
## WHICH SAMPLES ARE WHERE? ###
###############################

fastq_Size <- rep(0,length(SAMPLE_NAMES))
Step_Run <- array("N", dim=c(length(SAMPLE_NAMES),9))
Step_Size <- array(, dim=c(length(SAMPLE_NAMES),9))
#Step_Run[,1:3] <- "Y"
colnames(Step_Run) <- colnames(Step_Size) <- c(STEP_NAMES)
#colnames(Step_Run)[1:8] <- c(STEP_NAMES[1:8], "Prnt_Rds_1", "Prnt_Rds_2", "Prnt_Rds_3", "Prnt_Rds_4", "Prnt_Rds_5")
rownames(Step_Run) <- rownames(Step_Size) <- SAMPLE_NAMES

for (i in 1:length(SAMPLE_NAMES)) {

  fastq_Size[i] <- sum(fastq$FILE_SIZE_F[grep(SAMPLE_NAMES[i], fastq$SAMPLE_F)]) 

  Step_Run[i,1] <- length(grep(SAMPLE_NAMES[i], Map_F$FILE_NAME))
  Step_Size[i,1] <- sum(Map_F[grep(SAMPLE_NAMES[i], Map_F$FILE_NAME),3]) 

  Step_Run[i,2] <- length(grep(SAMPLE_NAMES[i], Bam_F$FILE_NAME))
  Step_Size[i,2] <- sum(Bam_F[grep(SAMPLE_NAMES[i], Bam_F$FILE_NAME),3]) 

  if (length(grep(SAMPLE_NAMES[i], Merge_F$FILE_NAME)) > 0) {
    Step_Size[i,3] <- Merge_F[grep(SAMPLE_NAMES[i], Merge_F$FILE_NAME),3] }

  if (length(grep(SAMPLE_NAMES[i], Sort_F$FILE_NAME)) > 0) {
    Step_Run[i,4] <- "Y" 
    Step_Size[i,4] <- Sort_F[grep(SAMPLE_NAMES[i], Sort_F$FILE_NAME),3] }

  if (length(grep(SAMPLE_NAMES[i], MrkDps_F$FILE_NAME)) > 0) {
    Step_Run[i,5] <- "Y"
    Step_Size[i,5] <- MrkDps_F[grep(SAMPLE_NAMES[i], MrkDps_F$FILE_NAME),3] }

  if (length(grep(SAMPLE_NAMES[i], TrgtCrtr_F$FILE_NAME)) > 0) {
    Step_Run[i,6] <- "Y"
    Step_Size[i,6] <- TrgtCrtr_F[grep(SAMPLE_NAMES[i], TrgtCrtr_F$FILE_NAME),3] }

  if (length(grep(SAMPLE_NAMES[i], IndelRlgnr_F$FILE_NAME)) > 0) {
    Step_Run[i,7] <- "Y"
    Step_Size[i,7] <- IndelRlgnr_F[grep(SAMPLE_NAMES[i], IndelRlgnr_F$FILE_NAME),3] }

  if (length(grep(SAMPLE_NAMES[i], BsRecal_F$FILE_NAME)) > 0) {
    Step_Run[i,8] <- "Y"
    Step_Size[i,8] <- BsRecal_F[grep(SAMPLE_NAMES[i], BsRecal_F$FILE_NAME),3] }

  Step_Run[i,9] <- length(grep(SAMPLE_NAMES[i], PrntRds_F$FILE_NAME))
  Step_Size[i,9] <- sum(PrntRds_F[grep(SAMPLE_NAMES[i], PrntRds_F$FILE_NAME),3]) }

Step_Size[which(Step_Size[,9]==0),9] <- NA
Step_Size <- data.frame(FastQ=fastq_Size, Step_Size)

## Show Memory Necessity for each Step

# Show Boxplot of the steps
png(paste(PathToSave,"2-Per_Sample_Storage.png",sep=""),width=1500,height=1000,pointsize=30)
boxplot(Step_Size[,1:10]/1e9, col=PLOT_COLS[2], main="Storage Required Per Sample", ylab="Output File Size (GB)", xlab="Processing Step",pch="+", xaxt="n")
text(x=1:10+.25, par("usr")[3] - 0.05*diff(range(Step_Size/1e9,na.rm=T)),labels=colnames(Step_Size), pos=2, srt=30, cex=1.05,xpd=T)
axis(1, at=1:10, labels=F)
abline( h=seq(0,1000,100), col=PLOT_COLS[9], lty=2,lwd=1)
boxplot(Step_Size[,1:10]/1e9, col=PLOT_COLS[2], pch="+", xaxt="n", add=T)
# for ( i in 1:10 ) { points( rep(i,nrow(Step_Size)), Step_Size[,i]/1e9, pch="+" ) }
dev.off()

# Add up total File Size for each Step and divide by 1000000000 (to convert to GB)
# (To extrapolate for incomplete steps, find mean & multiply by 438)
Step_Means <- c(mean(Step_Size[,1], na.rm=T), mean(Step_Size[,2], na.rm=T), mean(Step_Size[,3], na.rm=T), mean(Step_Size[,4], na.rm=T), mean(Step_Size[,5], na.rm=T), mean(Step_Size[,6], na.rm=T), mean(Step_Size[,7], na.rm=T), mean(Step_Size[,8], na.rm=T), mean(Step_Size[,9], na.rm=T), mean(Step_Size[,10], na.rm=T))
Step_Means_Round <- round(Step_Means/1000000000,3) # GB

Step_Total <- Step_Means*438/1000000000 # GB
Step_Total_Round <- round(Step_Total/1000,3) # Tb

png(paste(PathToSave,"2-Total_Storage.png",sep=""),width=1500,height=1000,pointsize=30)
barplot(Step_Total/1000, ylim=c(0,200), col=PLOT_COLS[2], main="Total (437 Samples) Storage Required by Step", ylab="Output File Size (TB)", xlab="Processing Step", xaxt="n", width=.8,space=.25)
text(x=1:10-.25, par("usr")[3] - 0.05*diff(range(Step_Total/1e3,na.rm=T)),labels=colnames(Step_Size), pos=2, srt=30, cex=1.05,xpd=T)
axis(1, at=1:10-.4, labels=F)
abline( h=seq(0,1000,50), col=PLOT_COLS[9], lty=2,lwd=1)
barplot(Step_Total/1000, col=PLOT_COLS[2], width=.8,space=.25, add=T)
dev.off()

## Demonstrate Correlation of File Size for Troubleshooting
png(paste(PathToSave,"2-Storage_Correlation.png",sep=""),width=2000,height=2000,pointsize=34)
pairs(Step_Size[,1:6]/1e9, col=PLOT_COLS[2], pch="+", main="Correlating File Size (GB) Between Steps")
dev.off()

##############################################################################################################################################################
##############################################################################################################################################################

###########################################################
## COMPILE SOME COMPUTING NUMBERS #########################
###########################################################

STEP_NAMES <- c("Map", "Bam", "Merge", "Sort", "MrkDps", "TrgtCrtr", "IndelRlgnr", "BsRecal", "PrntRds")
SAMPLE_NAMES <- unique(Jobs$SAMPLE)

## Calculate total FQ size by Sample
FQ_SIZE <- array(,dim=c(length(SAMPLE_NAMES),1))
colnames(FQ_SIZE) <- "FQ_SIZE"
rownames(FQ_SIZE) <- SAMPLE_NAMES
#FQ_SIZE <- rep(0,length(SAMPLE_NAMES))
#names(FQ_SIZE) <- SAMPLE_NAMES
for (i in 1:length(SAMPLE_NAMES)) {
  FQ_SIZE[i,1] <- sum(fastq$FILE_SIZE_F[grep(SAMPLE_NAMES[i], fastq$SAMPLE_F)])/1000000000
} # Divide by 1e9 to convert to GB

## Calculate SUs for the Jobs Array
JOBS_WALL_SPLIT2 <- strsplit(as.character(Jobs$WALL_TIME), ":")
JOBS_WALL_SPLIT <- t(sapply(JOBS_WALL_SPLIT2,"[",c(1:3)))
JOBS_WALL_SEC <- 3600*as.numeric(JOBS_WALL_SPLIT[,1]) + 60*as.numeric(JOBS_WALL_SPLIT[,2]) + as.numeric(JOBS_WALL_SPLIT[,3])

JOBS_CPU_SPLIT2 <- strsplit(as.character(Jobs$CPU_TIME), ":")
JOBS_CPU_SPLIT <- t(sapply(JOBS_CPU_SPLIT2,"[",c(1:3)))
JOBS_CPU_SEC <- 3600*as.numeric(JOBS_CPU_SPLIT[,1]) + 60*as.numeric(JOBS_CPU_SPLIT[,2]) + as.numeric(JOBS_CPU_SPLIT[,3])

COMMANDS_SPLIT2 <- strsplit(as.character(Jobs$JOB_NAME), "_")
COMMANDS_SPLIT <- as.numeric(sapply(COMMANDS_SPLIT2,"[",c(2)))
COMMANDS_SPLIT[Jobs$JOB_NAME=="Bam_Manual_1_20130919"] <- 6
COMMANDS_SPLIT[Jobs$JOB_NAME=="Bam_Manual_2_20130919"] <- 5
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_5_Test16"] <- 1
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_5b_Test15"] <- 15
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_r3_16a"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_r3_16b"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_6_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_7_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_8_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_9_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_10_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_11_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_12_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_13_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_14_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_15_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_16_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_17_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_18_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_19_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_20_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_21_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_22_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_23_16"] <- 16
COMMANDS_SPLIT[Jobs$JOB_NAME=="Sort_24_6"] <- 16

J <- data.frame(Jobs[c(1,2,3,4,7,8,9,14)], WALL_SEC=JOBS_WALL_SEC, WALLxCORES=JOBS_WALL_SEC*Jobs$CORES, CPU_SEC=JOBS_CPU_SEC, SU=round(JOBS_WALL_SEC*Jobs$CORES/3600,0), COMMperJOB=COMMANDS_SPLIT, SUperCOMM=round(JOBS_WALL_SEC*Jobs$CORES/3600,0)/COMMANDS_SPLIT)

#write.table(J, file="/home/kris/Desktop/Link to Dropbox/Schork/JNJ11/Usage/Jobs.xls", sep="\t", col.names=T, row.names=F, quote=F)

#########################################
## WHAT DO I WANT TO SHOW?? #############
#########################################

#1) Total Jobs & SUs broken up by Step (so far)
J1 <- J[which(J$STEP==STEP_NAMES[1]),]
J2 <- J[which(J$STEP==STEP_NAMES[2]),]
J3 <- J[which(J$STEP==STEP_NAMES[3]),]
J4 <- J[which(J$STEP==STEP_NAMES[4]),]
J5 <- J[which(J$STEP==STEP_NAMES[5]),]
J6 <- J[which(J$STEP==STEP_NAMES[6]),]
J7 <- J[which(J$STEP==STEP_NAMES[7]),]
J8 <- J[which(J$STEP==STEP_NAMES[8]),]
J9 <- J[which(J$STEP==STEP_NAMES[9]),]

# How Many Samples Done per Step So Far
SampByStep <- c(length(unique(J1$SAMPLE)), length(unique(J2$SAMPLE)), length(unique(J3$SAMPLE)), length(unique(J4$SAMPLE)), length(unique(J5$SAMPLE)), length(unique(J6$SAMPLE)), length(unique(J7$SAMPLE)), length(unique(J8$SAMPLE)), length(unique(J9$SAMPLE)))
# barplot(SampByStep, main="Samples Complete By Step", xlab="Processing Step", ylab="# Samples", names=STEP_NAMES, col="orange2", ylim=c(0,450))

# How Many SU's Per Step So Far
StepSUM <- c(sum(J1$SUperCOMM), sum(J2$SUperCOMM,na.rm=T), sum(J3$SUperCOMM), sum(J4$SUperCOMM,na.rm=T), sum(J5$SUperCOMM), sum(J6$SUperCOMM), sum(J7$SUperCOMM), sum(J8$SUperCOMM), sum(J9$SUperCOMM))
png(paste(PathToSave,"3-Total_Computing.png",sep=""),width=1500,height=1000,pointsize=30)
barplot(StepSUM, main="Total Computational Cost", xlab="Processing Step", ylab="SUs (Core-Hours)", xaxt="n", ylim=c(0,100000), col=PLOT_COLS[2], width=.8,space=.25 )
text(x=1:9-.25, par("usr")[3] - 0.05*diff(range(StepSUM,na.rm=T)),labels=colnames(Step_Size)[2:10], pos=2, srt=30, cex=1.05,xpd=T)
# axis(2, at=seq(0,1e5,1e4), labels=as.character(formatC( seq(0,1e5,1e4), format="e", digits=0)), las=2)
axis(1, at=1:9-.4, labels=F)
abline(h=seq(0,1e5,1e4), lty=2, ,lwd=1, col=PLOT_COLS[9] )
barplot(StepSUM, ylim=c(0,1e5), col=PLOT_COLS[2], add=T, xaxt="n", width=.8,space=.25 )
dev.off()
sum(StepSUM)

#2) How many SUs per sample per Step
By_Sample1 <- array(NA,dim=c(length(SAMPLE_NAMES),9))
colnames(By_Sample1) <- c("MAP_SU", "BAM_SU", "MERGE_SU", "SORT_SU", "DUPS_SU", "TRGT_SU", "INDEL_SU", "RECAL_SU", "PRINT_SU")
rownames(By_Sample1) <- SAMPLE_NAMES
for (i in 1:length(SAMPLE_NAMES)) {
  sample_info <- J[which(J$SAMPLE==SAMPLE_NAMES[i]),]
  By_Sample1[i,1] <- sum(sample_info$SUperCOMM[which(sample_info$STEP==STEP_NAMES[1])],na.rm=T)
  By_Sample1[i,2] <- sum(sample_info$SUperCOMM[which(sample_info$STEP==STEP_NAMES[2])],na.rm=T)
  By_Sample1[i,3] <- sum(sample_info$SUperCOMM[which(sample_info$STEP==STEP_NAMES[3])],na.rm=T)
  By_Sample1[i,4] <- sum(sample_info$SUperCOMM[which(sample_info$STEP==STEP_NAMES[4])],na.rm=T)
  By_Sample1[i,5] <- sum(sample_info$SUperCOMM[which(sample_info$STEP==STEP_NAMES[5])],na.rm=T)
  By_Sample1[i,6] <- sum(sample_info$SUperCOMM[which(sample_info$STEP==STEP_NAMES[6])],na.rm=T)
  By_Sample1[i,7] <- sum(sample_info$SUperCOMM[which(sample_info$STEP==STEP_NAMES[7])],na.rm=T)
  By_Sample1[i,8] <- sum(sample_info$SUperCOMM[which(sample_info$STEP==STEP_NAMES[8])],na.rm=T)
  By_Sample1[i,9] <- sum(sample_info$SUperCOMM[which(sample_info$STEP==STEP_NAMES[9])],na.rm=T) }
for (i in 1:ncol(By_Sample1)) {
  By_Sample1[which(By_Sample1[,i]==0),i] <- NA }
By_Sample <- data.frame(By_Sample1)
By_Sample_w_FQ <- data.frame(FQ_SIZE=FQ_SIZE, By_Sample1)

# Plot SUs per Sample by Step
png(paste(PathToSave,"3-Per_Sample_Computing.png",sep=""),width=1500,height=1000,pointsize=30)
boxplot(By_Sample, main="Computational Cost Per Sample", xlab="Processing Step", ylab="SUs (Core-Hours)", xaxt="n", col=PLOT_COLS[2],pch="+")
text(x=1:9+.1, par("usr")[3] - 0.05*diff(range(By_Sample,na.rm=T)),labels=colnames(Step_Size)[2:10], pos=2, srt=30, cex=1.05,xpd=T)
axis(1, at=1:9, labels=F)
abline(h=seq(0,1e3,50), lty=2, ,lwd=1, col=PLOT_COLS[9] )
boxplot(By_Sample, col=PLOT_COLS[2],pch="+", add=T, xaxt="n")
for (i in 1:9) {
  points(rep(i,438),By_Sample[,i], pch="+") }
dev.off()

# Boxplot of PrintReads step by # Cores
# boxplot(J9$SU ~ J9$CORES)

# Plot file size vs SUs for Print Reads Step
PRNT_RDS_MG <- merge(Jobs, PrntRds_F, by="FILE_NAME")
  # Calculate SUs
Split <- strsplit(as.character(PRNT_RDS_MG$WALL_TIME), ":")
Times <- t(sapply(Split,"[",1:3))
Secs <- 3600*as.numeric(Times[,1]) + 60*as.numeric(Times[,2]) + as.numeric(Times[,3])
SUs <- round(Secs*PRNT_RDS_MG$CORES/3600,0)
SUs[which(SUs==0)] <- 1
COLORS <- PLOT_COLS[c(4,6)] # c("slateblue3","tomato2")
COLORS.2 <- PLOT_COLS.4[c(4,6)]
CORE_COLS <- as.numeric(factor(MG$CORES)) # ; names(CORE_COLS) <- c("8","16")
 # Calculate fit for 8 and 16 cores
C.8 <- which( PRNT_RDS_MG$CORES==8 )
MOD.8 <- lm( SUs[C.8] ~ I(PRNT_RDS_MG$SIZE[C.8]/1e9) ) ; COEF.8 <- round( coef(MOD.8)[2], 3 )
MOD.16 <- lm( SUs[-C.8] ~ I(PRNT_RDS_MG$SIZE[-C.8]/1e9) ) ; COEF.16 <- round( coef(MOD.16)[2], 3 )
CORE_COLS <- as.numeric(factor(PRNT_RDS_MG$CORES)) # ; names(CORE_COLS) <- c("8","16")
png(paste(PathToSave,"4-PrintReads_Efficiency.png",sep=""),width=1000,height=1000,pointsize=30)
plot(PRNT_RDS_MG$SIZE/1e9, SUs, col=COLORS[CORE_COLS], pch="+", xlab="Recalibrated Bam File Size (GB)", ylab="SUs (Cores-Hours)", main="Print Reads Efficiency", xaxt="n")
axis(1, at=seq(0,70,10))
abline(h=seq(0,200,10), lty=2, ,lwd=1, col=PLOT_COLS[9] )
abline(v=seq(0,200,10), lty=2, ,lwd=1, col=PLOT_COLS[9] )
points( PRNT_RDS_MG$SIZE/1e9, SUs, pch="+", col=COLORS[CORE_COLS] )
abline( MOD.8, col=COLORS.2[1], lty=2, lwd=6 )
abline( MOD.16, col=COLORS.2[2], lty=2, lwd=6 )
legend(10,100,legend=c("8","16"),title="# Cores",pch="+",col=COLORS)
text( 0, c(55,50), labels=paste(c(8,16),"Cores:",c(COEF.8,COEF.16),"SU/GB"), col=COLORS.2, pos=4 )
dev.off()

#3) Extrapolate to Estimate How Many SUs will be required for each step
MEANS <- apply(By_Sample, 2, mean, na.rm=T)
EXTRAP <- MEANS*438
# x11()
# barplot(EXTRAP, main="Extrapolation of Mean Values for Total CPU Usage", xlab="Processing Step", ylab="SUs", names=STEP_NAMES, col="green4")
# sum(EXTRAP)












###############################
## TROUBLE SHOOTING!!! ########
###############################

plot(Step_Size[,1:2])
PERC_10 <- Step_Size[,2]/Step_Size[,1]
Temp <- Step_Size[which(PERC_10<3.5),]
CHECK_MAP <- c()
for (i in rownames(Temp)) {
CHECK_MAP <- rbind(CHECK_MAP, FILES_930[grep(i, FILES_930$FILE_NAME),], FILES_914[grep(i, FILES_914$FILE_NAME),]) }

plot(Step_Size[,2:3])
PERC_21 <- Step_Size[,3]/Step_Size[,2]
Temp <- Step_Size[which(PERC_21<.3),]
CHECK_BAM <- c()
for (i in rownames(Temp)) {
CHECK_BAM <- rbind(CHECK_BAM, FILES_930[grep(i, FILES_930$FILE_NAME),], FILES_914[grep(i, FILES_914$FILE_NAME),]) }

plot(Step_Size[,3:4])
PERC_32 <- Step_Size[,4]/Step_Size[,3]
Temp <- Step_Size[which(PERC_32<.3),]
CHECK_MERGE <- c()
for (i in rownames(Temp)) {
CHECK_MERGE <- rbind(CHECK_MERGE, FILES_930[grep(i, FILES_930$FILE_NAME),], FILES_914[grep(i, FILES_914$FILE_NAME),]) }

plot(Step_Size[,c(1,4)])
PERC_30 <- Step_Size[,4]/Step_Size[,1]
Temp <- Step_Size[which(PERC_30<.3),]
CHECK_MERGE_FQ <- c()
for (i in rownames(Temp)) {
CHECK_MERGE_FQ <- rbind(CHECK_MERGE_FQ, FILES_930[grep(i, FILES_930$FILE_NAME),], FILES_914[grep(i, FILES_914$FILE_NAME),]) }

