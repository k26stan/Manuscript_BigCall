## Call Variants on JnJ ART3001 Cohort w/ Various Group Sizes ##
 # Only Chr 21
 # Repeat a few times w/ the smaller groups
sort -R input | head -n 100 >output

#####################################################
## WRITE SHELL FILE #################################
#####################################################

#!/bin/bash
#PBS -N HC_PERMNAME
#PBS -q shared
#PBS -v QOS=10
#PBS -l walltime=72:00:00
#PBS -m abe
#PBS -l nodes=1:ppn=16:native:flash
#PBS -o HC_PERMNAME.run.oe
#PBS -j oe
#PBS -M k26stan@yahoo.com
#PBS -V
#PBS -A jan100
cd /projects/janssen/scripts/20150701_HC_Time_Test/
echo "<startTime>"`date`"</startTime>"
echo "<output>"

## Specify Paths to Tools
GATK_JAR=/projects/janssen/Tools/gatk3.3-0/GenomeAnalysisTK.jar
REF_FA=/projects/janssen/ref/ref.fa
OUTPUT_DIR=/projects/janssen/scripts/20150701_HC_Time_Test/VCFs
SCRATCH_DIR=/scratch/$USER/$PBS_JOBID
CHR_LIST=/projects/janssen/scripts/20150701_HC_Time_Test/Chr_21.list
DB_SNP=/projects/janssen/gatk_bundle/dbsnp_137.b37.vcf


## Write HaplotypeCaller Command
java -Xmx52g \
-Djava.io.tmpdir=${SCRATCH_DIR} \
-jar ${GATK_JAR} \
-T HaplotypeCaller \
-R ${REF_FA} \
INSERT_BAMS_HERE \
-nct 16 \
--dbsnp ${DB_SNP} \
-L ${CHR_LIST} \
-dcov 200 \
-stand_call_conf 40 \
-stand_emit_conf 10 \
-o ${OUTPUT_DIR}/${PERMNAME}.raw.vcf

echo "</output>"
echo "<exitStatus>"$?"</exitStatus>"
echo "<stopTime>"`date`"</stopTime>"
qstat -f $PBS_JOBID | grep Job
qstat -f $PBS_JOBID | grep resources

#####################################################
## WRITE/RUN w/ VARYING GROUP SIZES #################
#####################################################

## Set up Paths
HC_SHELL=/projects/janssen/scripts/20150701_HC_Time_Test/HC_Shell.run
HC_RUNFILES=/projects/janssen/scripts/20150701_HC_Time_Test/RUN_FILES/
BAM_FILE_LOCS=/projects/janssen/scripts/FINAL_BAM_FILENAMES_PT4

## Loop Through Various Group Sizes
 # Smaller groups...5 times each
for g in 1 {5..25..5} ; do
	for i in `seq 5` ; do
		tag=Grp${g}_It${i}
		echo $tag
		## Rename File/Permname
		sed "s/PERMNAME/$tag/g" ${HC_SHELL} > ${HC_RUNFILES}/HC_${tag}.run
		## Sample g Patients Randomly
		PATIENTS=`sort -R $BAM_FILE_LOCS | head -${g}| tr '\n' ' '`
		sed -i "s|INSERT_BAMS_HERE|${PATIENTS}|g" ${HC_RUNFILES}/HC_${tag}.run
	done
done
 # Slightly Larger...4 times each
for g in 1 {30..50..10} ; do
	for i in `seq 5` ; do
		tag=Grp${g}_It${i}
		echo $tag
		## Rename File/Permname
		sed "s/PERMNAME/$tag/g" ${HC_SHELL} > ${HC_RUNFILES}/HC_${tag}.run
		## Sample g Patients Randomly
		PATIENTS=`sort -R $BAM_FILE_LOCS | head -${g}| tr '\n' ' '`
		sed -i "s|INSERT_BAMS_HERE|${PATIENTS}|g" ${HC_RUNFILES}/HC_${tag}.run
	done
done
 # Large Groups...2 times each
for g in {100..200..100} ; do
	for i in `seq 2` ; do
		tag=Grp${g}_It${i}
		echo $tag
		## Rename File/Permname
		sed "s/PERMNAME/$tag/g" ${HC_SHELL} > ${HC_RUNFILES}/HC_${tag}.run
		## Sample g Patients Randomly
		PATIENTS=`sort -R $BAM_FILE_LOCS | head -${g}| tr '\n' ' '`
		sed -i "s|INSERT_BAMS_HERE|${PATIENTS}|g" ${HC_RUNFILES}/HC_${tag}.run
	done
done
 # Largest Group...1 time
for g in 300 400 437 ; do
	for i in `seq 1` ; do
		tag=Grp${g}_It${i}
		echo $tag
		## Rename File/Permname
		sed "s/PERMNAME/$tag/g" ${HC_SHELL} > ${HC_RUNFILES}/HC_${tag}.run
		## Sample g Patients Randomly
		PATIENTS=`sort -R $BAM_FILE_LOCS | head -${g}| tr '\n' ' '`
		sed -i "s|INSERT_BAMS_HERE|${PATIENTS}|g" ${HC_RUNFILES}/HC_${tag}.run
	done
done

#####################################################
## COMPILE RUN STATS FROM OUTPUT FILE ###############
#####################################################
HC_RUNFILES=/projects/janssen/scripts/20150701_HC_Time_Test/RUN_FILES/

cd $HC_RUNFILES
for file in `ls HC*oe`; do
	grep ProgressMeter $file > Progress_${file}
done




#####################################################
## END OF DOC #######################################
#####################################################