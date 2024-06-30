#!/bin/bash
#WES Inhouse pipeline for pre-processing
#Edited in Server #3

bwa="/DATA07/home/gangtl22/Tools/bwa-0.7.17"
gatk="python3.8 /DATA07/home/gangtl22/Tools/GATK/gatk-4.3.0.0"
trim="java -jar /DATA07/home/gangtl22/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
picard="java -jar /DATA07/home/gangtl22/Tools/picard-2.21.1/picard.jar"
ref=/DATA07/home/gangtl22/PIPE/script/ref/INHOUSE

echo "1. Sample Prefix: ${1}"
echo "2. Normal Prefix: ${2}"
echo "3. /Path/To/Normal/Input: ${3}"
echo "4. Tumor Prefix: ${4}"
echo "5. /Path/To/Tumor/Input: ${5}"
echo "6. Reference Genome: ${6}"
echo "7. EGFT_MT: ${7}"
echo "8. Output directory: ${8}"
if [ ${9} -eq 0 ];then
	echo "9. Analysis whole region"
else
	echo "9. Anlaysis Chromosome ${9} region"
fi
echo "10. Number of threads: ${10}"
echo "11. Library Kit Name: ${11}"
if [ ${12} == 'Y' ];then
	echo "12. Force to process this step: Yes"
else
	echo "12. Force to process this step: No"
fi
sleep 3s
mkdir -p ${8}/${1}/report/QC/
mkdir -p ${8}/log

#--------------------------------------------------------------------------------------------
# STEP00.Indexing Reference
#--------------------------------------------------------------------------------------------
IndexRef(){
	echo "## STEP00.Indexing Reference - `date` ## " >> ${8}/log/${1}_Tlog.txt
	${bwa}/bwa index ${ref}/${6}/${6}.fa

	##Index reference .fai file before running haplotype caller
	samtools faidx ${ref}/${6}/${6}.fa

	## Ref dict before running haplotype caller 
	${picard} CreateSequenceDictionary R=${ref}/${6}/${6}.fa O=${ref}/${6}/${6}.dict
}


#--------------------------------------------------------------------------------------------
# STEP01.Trimming and Quality Check
#--------------------------------------------------------------------------------------------

Trimming(){
	echo "## STEP01.Trimming and Quality Check ##" >> ${8}/log/${1}_Tlog.txt
	echo "## Trimming & QC START - `date` ##" >> ${8}/log/${1}_Tlog.txt


	if [ ! -e ${8}/${1}/${4}_2.trim.fastq.gz ] && [ ! -e ${8}/${1}/${4}.recal.bam ];then
		if [ ${4} = 'NA' ];then
			echo "## No Tumor Sample ##" >> ${8}/log/${1}_Tlog.txt
			echo "## No Tumor Sample ##"
		else
			${trim} PE \
				-phred33 \
	       			-threads ${10} \
	       			${5}/${4}_1.fastq.gz ${5}/${4}_2.fastq.gz \
				${8}/${1}/${4}_1.trim.fastq.gz ${8}/${1}/report/QC/${4}_1.untrim.fastq.gz \
	       			${8}/${1}/${4}_2.trim.fastq.gz ${8}/${1}/report/QC/${4}_2.untrim.fastq.gz \
	      			ILLUMINACLIP:/DATA07/home/gangtl22/Tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
			fastqc ${8}/${1}/${4}_1.trim.fastq.gz -o ${8}/${1}/report/QC/
			fastqc ${8}/${1}/${4}_2.trim.fastq.gz -o ${8}/${1}/report/QC/
		fi
	else
		echo "## Tumor trimming pass ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Tumor trimming pass ##"
	fi


	if [ ! -e ${8}/${1}/${2}_2.trim.fastq.gz ] && [ ! -e ${8}/${1}/${2}.recal.bam ];then
		if [ ${2} = 'NA' ];then
			echo "## No Normal Sample ##" >> ${8}/log/${1}_Tlog.txt
			echo "## No Normal Sample ##"

		else
			${trim} PE \
				-phred33 \
				-threads ${10} \
				${3}/${2}_1.fastq.gz ${3}/${2}_2.fastq.gz \
				${8}/${1}/${2}_1.trim.fastq.gz ${8}/${1}/report/QC/${2}_1.untrim.fastq.gz \
				${8}/${1}/${2}_2.trim.fastq.gz ${8}/${1}/report/QC/${2}_2.untrim.fastq.gz \
				ILLUMINACLIP:/DATA07/home/gangtl22/Tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

			fastqc ${8}/${1}/${2}_1.trim.fastq.gz -o ${8}/${1}/report/QC/
			fastqc ${8}/${1}/${2}_2.trim.fastq.gz -o ${8}/${1}/report/QC/
		fi
	else
		echo "## Normal trimming pass ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Normal trimming pass ##"
	fi
	chmod 755 ${8}/${1}/*.trim.fastq.gz
	echo "## Trimming & QC END - `date` ##" >> ${8}/log/${1}_Tlog.txt
	echo "## Trimming & QC END - `date` ##"
}

#--------------------------------------------------------------------------------------------
# STEP02.Map to reference using BWA-MEM
#--------------------------------------------------------------------------------------------

Mapping(){
	if [ ${9} = 0 ];then
		echo "## STEP02. Map to reference using BWA-MEM ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Mapping START - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## STEP02. Map to reference using BWA-MEM ##"
		echo "## Mapping START - `date` ##"
		#Tumor
		if [ ! -e ${8}/${1}/${4}.aln.bam ] && [ ! -e ${8}/${1}/${4}.recal.bam ];then
			if [ ${4} = 'NA' ];then
				echo "## Tumor Pass Cuz No Sample ##" >> ${8}/log/${1}_Tlog.txt 
				echo "## Tumor Pass Cuz No Sample ##"
			else
				bwa mem \
					-t ${10} \
					-R "@RG\tID:${1}\tPL:NA\tSM:${4}" \
					${ref}/${6}/${6}.fa \
					${8}/${1}/${4}_1.trim.fastq.gz \
					${8}/${1}/${4}_2.trim.fastq.gz |\
					samtools view -@ ${10} -Shb -o ${8}/${1}/${4}.aln.bam

				# Add read group information to the header
#				${picard} AddOrReplaceReadGroups \
#					I=${8}/${1}/${4}.nrg.aln.bam \
#					O=${8}/${1}/${4}.aln.bam \
#					RGID=${1} \
#					RGLB=NA \
#					RGPL=NA \
#					RGPU=NA \
#					RGSM=${2} \
#					CREATE_INDEX=true
				samtools flagstat ${8}/${1}/${4}.aln.bam > ${8}/${1}/report/${4}.flags.results
				samtools sort ${8}/${1}/${4}.aln.bam -o ${8}/${1}/${4}.sorted.bam
			fi

		else
			echo "## Tumor Mapping Is Passed Cuz Already Existed ##"  >> ${8}/log/${1}_Tlog.txt
			echo "## Tumor Mapping Is Passed Cuz Already Existed ##"
		fi

		#Normal
		if [ ! -e ${8}/${1}/${2}.aln.bam ] && [ ! -e ${8}/${1}/${2}.recal.bam ];then
			if [ ${2} = 'NA' ];then
				echo "## Normal Pass Cuz No Sample ##" >> ${8}/log/${1}_Tlog.txt
				echo "## Normal Pass Cuz No Sample ##"
			else
				bwa mem \
					-t ${10} \
					-R "@RG\tID:${1}\tPL:NA\tSM:${2}" \
					${ref}/${6}/${6}.fa \
					${8}/${1}/${2}_1.trim.fastq.gz \
					${8}/${1}/${2}_2.trim.fastq.gz |\
					samtools view -@ ${10} -Shb -o ${8}/${1}/${2}.aln.bam

				# Add read group information to the header
#				${picard} AddOrReplaceReadGroups \
#			    		I=${8}/${1}/${2}.nrg.aln.bam \
#			    		O=${8}/${1}/${2}.aln.bam \
#			    		RGID=${1} \
#					RGLB=NA \
#					RGPL=NA \
#					RGPU=NA \
#					RGSM=${2} \
#			    		CREATE_INDEX=true

				samtools flagstat ${8}/${1}/${2}.aln.bam > ${8}/${1}/report/${2}.flags.results
				samtools sort ${8}/${1}/${2}.aln.bam -o ${8}/${1}/${2}.sorted.bam
			fi
		else
			echo "## Normal Mapping Is Passed Cuz Already Existed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Normal Mapping Is Passed Cuz Already Existed ##"
		fi
		echo "## Mapping END - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Mapping END - `date` ##"
	else
		if [ ! -e ${8}/${1}/CHR${9} ];then
			mkdir -p ${8}/${1}/CHR${9}
		else
			echo "Pass"
		fi
		echo "## STEP02. Map Only CHR${9} region to reference using BWA-MEM ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Mapping START - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## STEP02. Map Only CHR${9} region to reference using BWA-MEM ##"
		echo "## Mapping START - `date` ##"
		#Tumor
		if [ ${4} != 'NA' ] && [ ! -e ${8}/${1}/CHR${9}/${4}.aln.bam ] && [ ! -e ${8}/${1}/CHR${9}/${4}.recal.bam ];then
			if [ ${4} = 'NA' ];then
				echo "## Tumor Pass Cuz No Sample ##" >> ${8}/log/${1}_Tlog.txt
				echo "## Tumor Pass Cuz No Sample ##"
			else
				bwa mem \
					-t ${10} \
					-R "@RG\tID:${1}\tPL:NA\tSM:${4}" \
					${ref}/${6}/fasta_chr${9}.fa \
					${8}/${1}/${4}_1.trim.fastq.gz \
					${8}/${1}/${4}_2.trim.fastq.gz |\
					samtools view -@ ${10} -Shb chr${9} -o ${8}/${1}/CHR${9}/${4}.nrg.aln.bam

				# Add read group information to the header
				${picard} AddOrReplaceReadGroups \
					I=${8}/${1}/CHR${9}/${4}.nrg.aln.bam \
					O=${8}/${1}/CHR${9}/${4}.aln.bam \
					RGID=${1} \
					RGLB=NA \
					RGPL=NA \
					RGPU=NA \
					RGSM=${4} \
					CREATE_INDEX=true
				samtools flagstat ${8}/${1}/CHR${9}/${4}.aln.bam > ${8}/${1}/report/${4}.chr${9}.flags.results
				samtools sort ${8}/${1}/CHR${9}/${4}.aln.bam -o ${8}/${1}/CHR${9}/${4}.sorted.bam
			fi
		else
			echo "## Tumor Mapping ON CHR${9} Region Is Passed Cuz Already Existed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Tumor Mapping ON CHR${9} Region Is Passed Cuz Already Existed ##"
		fi
		#Normal
		if [ ! -e ${8}/${1}/CHR${9}/${2}.aln.bam ] && [ ! -e ${8}/${1}/CHR${9}/${2}.recal.bam ];then
			if [ ${2} = 'NA' ];then
				echo "## Normal Pass Cuz No Sample ##" >> ${8}/log/${1}_Tlog.txt
				echo "## Normal Pass Cuz No Sample ##"
			else
				bwa mem \
					-t ${10} \
					-R "@RG\tID:${1}\tPL:NA\tSM:${2}" \
					${ref}/${6}/fasta_chr${9}.fa \
					${8}/${1}/${2}_1.trim.fastq.gz \
					${8}/${1}/${2}_2.trim.fastq.gz |\
					samtools view -@ ${10} -Shb chr${9} -o ${8}/${1}/CHR${9}/${2}.nrg.aln.bam
	
				# Add read group information to the header
				${picard} AddOrReplaceReadGroups \
					I=${8}/${1}/CHR${9}/${2}.nrg.aln.bam \
					O=${8}/${1}/CHR${9}/${2}.aln.bam \
					RGID=${1} \
					RGLB=NA \
					RGPL=NA \
					RGPU=NA \
					RGSM=${2} \
					CREATE_INDEX=true
				samtools flagstat ${8}/${1}/CHR${9}/${2}.aln.bam > ${8}/${1}/report/${2}.chr${9}.flags.results
				samtools sort ${8}/${1}/CHR${9}/${2}.aln.bam -o ${8}/${1}/CHR${9}/${2}.sorted.bam
			fi
		else
			echo "## Normal Mapping ON CHR${9} Region Is Passed Cuz Already Existed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Tumor Mapping ON CHR${9} Region Is Passed Cuz Already Existed ##"
		fi
		echo "## Mapping ON CHR${9} Region END - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Mapping ON CHR${9} Region END - `date` ##"
	fi
}

#--------------------------------------------------------------------------------------------
# STEP03.Mark Duplicates and Sort - GATK4
#--------------------------------------------------------------------------------------------

MarkDup(){
	if [ ${9} = 0 ];then
		echo "## STEP03.Mark Duplicates and Sort - GATK4 ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Duplicate Marking START - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## STEP03.Mark Duplicates and Sort - GATK4 ##"
		echo "## Duplicate Marking START - `date` ##"		
		#Tumor
		if [ ! -e ${8}/${1}/${4}.mark_dup.bam ] && [ ${4} != 'NA' ] && [ ! -e ${8}/${1}/${4}.recal.bam ];then
			${gatk}/gatk MarkDuplicatesSpark \
				-I ${8}/${1}/${4}.sorted.bam \
				-O ${8}/${1}/${4}.mark_dup.bam \
				-M ${8}/${1}/report/${4}.mark_dup_metrics.txt \
				--conf 'spark.executor.cores=3'
				#--remove-sequencing-duplicates
		else
			echo "## Marking Duplicates From Tumor Sample Is Passed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Marking Duplicates From Tumor Sample Is Passed ##"
		fi

		#Normal
		if [ ! -e ${8}/${1}/${2}.mark_dup.bam ] && [ ${2} != 'NA' ] && [ ! -e ${8}/${1}/${2}.recal.bam ];then
			${gatk}/gatk MarkDuplicatesSpark \
				-I ${8}/${1}/${2}.sorted.bam \
				-O ${8}/${1}/${2}.mark_dup.bam \
				-M ${8}/${1}/report/${2}.mark_dup_metrics.txt \
				--conf 'spark.executor.cores=3'
				#--remove-sequencing-duplicates
		else
			echo "## Marking Duplicates From Normal Sample Is Passed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Marking Duplicates From Normal Sample Is Passed ##"

		fi
		echo "## Duplicate Marking END - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Duplicate Marking END - `date` ##"
	else
		echo "## STEP03.Mark Duplicates and Sort on CHR${9} region - GATK4 ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Duplicate Marking START - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## STEP03.Mark Duplicates and Sort on CHR${9} region - GATK4 ##"
		echo "## Duplicate Marking START - `date` ##"

		#Tumor
		if [ ! -e ${8}/${1}/CHR${9}/${4}.mark_dup.bam ] && [ ${4} != 'NA' ] && [ ! -e ${8}/${1}/CHR${9}/${4}.recal.bam ];then
			${gatk}/gatk MarkDuplicatesSpark \
				-I ${8}/${1}/CHR${9}/${4}.sorted.bam \
				-O ${8}/${1}/CHR${9}/${4}.mark_dup.bam \
				-M ${8}/${1}/report/${4}.chr${9}.mark_dup_metrics.txt \
				--conf 'spark.executor.cores=3'
				#--remove-sequencing-duplicates
		else
			echo "## Marking Duplicates From Tumor CHR${9} Region Is Passed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Marking Duplicates From Tumor CHR${9} Region Is Passed ##"
		fi

		#Normal
		if [ ! -e ${8}/${1}/CHR${9}/${2}.mark_dup.bam ] && [ ${2} != 'NA' ] && [ ! -e ${8}/${1}/CHR${9}/${2}.recal.bam ];then
			${gatk}/gatk MarkDuplicatesSpark \
				-I ${8}/${1}/CHR${9}/${2}.sorted.bam \
				-O ${8}/${1}/CHR${9}/${2}.mark_dup.bam \
				-M ${8}/${1}/report/${2}.chr${9}.mark_dup_metrics.txt \
				--conf 'spark.executor.cores=3'
				#--remove-sequencing-duplicates
		else
			echo "## Marking Duplicates From Normal CHR${9} Region Is Passed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Marking Duplicates From Normal CHR${9} Region Is Passed ##"
		fi
		echo "## Duplicate Marking ON CHR${9} Region END - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Duplicate Marking ON CHR${9} Region END - `date` ##"
	fi
}
#--------------------------------------------------------------------------------------------
# STEP04.BQSR -GATK4
#
# BASE Quality Score Recalibration
# A BQSR step is then performed using BaseRecalibrator. This step adjusts base quality scores 
# based on detectable and systematic errors. This step also increases the accuracy of downstream 
# variant calling algorithms. 
# Note that the original quality scores are kept in the QQ field of co-cleaned BAM files. 
# These scores should be used if conversion of BAM files to FASTQ format is desired.
#
# GATK4 does not support Indel Realignment home/gangtl22/Tools anymore (RealignerTargetCreator and IndelRealigner)
#--------------------------------------------------------------------------------------------
Recal(){
	if [ ${6} = 'hg19' ];then
		knownSites=${ref}/hg19_vcf/dbsnp_138.hg19.vcf
		knownSites2=${ref}/hg19_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
	else
		knownSites=${ref}/hg38_vcf/Homo_sapiens_assembly38.dbsnp138.vcf
		knownSites2=${ref}/hg38_vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
	fi
	if [ ${9} = 0 ];then
		echo "## BASE Quality Score Recalibration START - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## BASE Quality Score Recalibration START - `date` ##"
		
		
		if [ ! -e ${8}/${1}/${4}.recal.bam ] && [ ${4} != 'NA' ];then
			#1. Build the model - Tumor
			${gatk}/gatk BaseRecalibrator \
				-R ${ref}/${6}/${6}.fa \
				-I ${8}/${1}/${4}.mark_dup.bam \
				--known-sites ${knownSites}\
				--known-sites ${knownSites2}\
				-O ${8}/${1}/${4}_recal.table
			#2. Apply the model to adjust the base quality scores -Tumor
			${gatk}/gatk ApplyBQSR\
				-R ${ref}/${6}/${6}.fa \
				-I ${8}/${1}/${4}.mark_dup.bam \
				--bqsr-recal-file ${8}/${1}/${4}_recal.table \
				-O ${8}/${1}/${4}.recal.bam
		else
			echo "## Recalibration Step -Tumor Is Passed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Recalibration Step -Tumor Is Passed ##"

		fi
		
		if [ ! -e ${8}/${1}/${2}.recal.bam ] && [ ${2} != 'NA' ];then
			#1. Build the model - Normal
			${gatk}/gatk BaseRecalibrator \
				-R ${ref}/${6}/${6}.fa \
				-I ${8}/${1}/${2}.mark_dup.bam \
				--known-sites ${knownSites}\
				--known-sites ${knownSites2}\
				-O ${8}/${1}/${2}_recal.table
			#2. Apply the model to adjust the base quality scores -Normal
			${gatk}/gatk ApplyBQSR\
				-R ${ref}/${6}/${6}.fa \
				-I ${8}/${1}/${2}.mark_dup.bam \
				--bqsr-recal-file ${8}/${1}/${4}_recal.table \
				-O ${8}/${1}/${2}.recal.bam
		else
			echo "## Recalibration Step -Normal Is Passed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Recalibration Step -Normal Is Passed ##"
		fi

		echo "## BASE Quality Score Recalibration END - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## BASE Quality Score Recalibration END - `date` ##"
	else
		echo "## BASE Quality Score Recalibration on CHR${9} region START - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## BASE Quality Score Recalibration on CHR${9} region START - `date` ##"
		if [ ! -e ${8}/${1}/CHR${9}/${4}.recal.bam ] && [ ${4} != 'NA' ];then
			#1. Build the model - Tumor
			${gatk}/gatk BaseRecalibrator \
				-R ${ref}/${6}/fasta_chr${9}.fa \
				-I ${8}/${1}/CHR${9}/${4}.recal.bam \
				--known-sites ${knownSites}\
				--known-sites ${knownSites2}\
				-O ${8}/${1}/CHR${9}/${4}_recal.table
			#2. Apply the model to adjust the base quality scores -Tumor
			${gatk}/gatk ApplyBQSR\
				-R ${ref}/${6}/fasta_chr${9}.fa \
				-I ${8}/${1}/CHR${9}/${4}.mark_dup.bam \
				--bqsr-recal-file ${8}/${1}/CHR${9}/${4}_recal.table \
				-O ${8}/${1}/CHR${9}/${4}.recal.bam
		else
			echo "## Recalibration Step -Tumor CHR${9} Region Is Passed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Recalibration Step -Tumor CHR${9} Region Is Passed ##"
		fi
			
		if [ ! -e ${8}/${1}/CHR${9}/${2}.recal.bam ] && [ ${2} != 'NA' ];then
			#1. Build the model - Normal
			${gatk}/gatk BaseRecalibrator \
					-R ${ref}/${6}/fasta_chr${9}.fa \
					-I ${8}/${1}/CHR${9}/${2}.mark_dup.bam \
					--known-sites ${knownSites}\
					--known-sites ${knownSites2}\
					-O ${8}/${1}/CHR${9}/${2}_recal.table
				
			#2. Apply the model to adjust the base quality scores -Normal
			${gatk}/gatk ApplyBQSR\
				-R ${ref}/${6}/fasta_chr${9}.fa \
				-I ${8}/${1}/CHR${9}/${2}.mark_dup.bam \
				--bqsr-recal-file ${8}/${1}/CHR${9}/${2}_recal.table \
				-O ${8}/${1}/CHR${9}/${2}.recal.bam
			
		else
			echo "## Recalibration Step -Normal CHR${9} Region Is Passed ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Recalibration Step -Noraml CHR${9} Region Is Passed ##"


		fi

		echo "## BASE Quality Score Recalibration END - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## BASE Quality Score Recalibration END - `date` ##"
	fi
	if [ ${9} = 0 ];then
		if [ ${2} != 'NA' ] && [ ${4} != 'NA' ];then
			echo -e "${1}\t${2}\t${8}/${1}/${2}.recal.bam\t${4}\t${8}/${1}/${4}.recal.bam\t${6}\t${7}\t${11}"
			echo -e "${1}\t${2}\t${8}/${1}/${2}.recal.bam\t${4}\t${8}/${1}/${4}.recal.bam\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.output.txt
		
		elif [ ${2} = 'NA' ] && [ ${4} != 'NA' ];then
			echo -e "${1}\t${2}\tNA\t${4}\t${8}/${1}/${4}.recal.bam\t${6}\t${7}\t${11}"
			echo -e "${1}\t${2}\tNA\t${4}\t${8}/${1}/${4}.recal.bam\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.output.txt
	
		elif [ ${2} != 'NA' ] && [ ${4} = 'NA' ];then
			echo -e "${1}\t${2}\t${8}/${1}/${2}.recal.bam\t${4}\tNA\t${6}\t${7}\t${11}"
			echo -e "${1}\t${2}\t${8}/${1}/${2}.recal.bam\t${4}\tNA\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.output.txt
		else
			echo -e "${1}\t${2}\tNA\t${4}\tNA\t${6}\t${7}\t${11}"
			echo -e "${1}\t${2}\tNA\t${4}\tNA\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.output.txt
		fi
	else
		if [ ${2} != 'NA' ] && [ ${4} != 'NA' ];then
			echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/${2}.recal.bam\t${4}\t${8}/${1}/CHR${9}/${4}.recal.bam\t${6}\t${7}\t${11}"
			echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/${2}.recal.bam\t${4}\t${8}/${1}/CHR${9}/${4}.recal.bam\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.CHR${9}.output.txt
		
                elif [ ${2} = 'NA' ] && [ ${4} != 'NA' ];then
			echo -e "${1}\t${2}\tNA\t${4}\t${8}/${1}/CHR${9}/${4}.recal.bam\t${6}\t${7}\t${11}"
			echo -e "${1}\t${2}\tNA\t${4}\t${8}/${1}/CHR${9}/${4}.recal.bam\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.CHR${9}.output.txt

		elif [ ${2} != 'NA' ] && [ ${4} = 'NA' ];then
			echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/${2}.recal.bam\t${4}\tNA\t${6}\t${7}\t${11}"
			echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/${2}.recal.bam\t${4}\tNA\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.CHR${9}.output.txt
 		else
			echo -e "${1}\t${2}\tNA\t${4}\tNA\t${6}\t${7}\t${11}"
			echo -e "${1}\t${2}\tNA\t${4}\tNA\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.CHR${9}.output.txt
		fi
	fi

}

ExtractChr(){
	if [ ${9} == 'egfr' ];then
		echo "## Extrating EGFR Location Start - `date` ##" >> ${8}/log/${1}_Tlog.txt
		if [ ! -e ${8}/EGFR/ ];then
			mkdir -p ${8}/EGFR
		else
			echo "pass"
		fi
		#EGFR (hg19) : chr7:55086710-55279321
		#E19 deletion loucs : chr7:55,242,118-55,242,826 (hg19)
		#L858R locus : chr7:55,259,473-55,259,561 (hg19)
		#EGFR (hg38) : chr7:55019017-55211628

		if [ -e ${5}/${4}*recal.bam ] && [ ${4} != 'NA' ];then
			if [ ${6} = 'hg19' ];then
				samtools view -bh ${5}/${4}*recal.bam \
					chr7:55086710-55279321 > ${8}/EGFR/${4}.egfr.bam
				samtools index ${8}/EGFR/${4}.egfr.bam
				if [ ! -e ${3}/${2}*recal.bam ];then
					echo "Normal is not existed"
					echo -e "${1}\t${2}\tNA\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}"
					echo -e "${1}\t${2}\tNA\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}" >> ${8}/EGFR/prep.inhouse.output.txt

				else
					samtools view -bh ${3}/${2}*recal.bam \
						chr7:55086710-55279321 > ${8}/EGFR/${2}.egfr.bam
					samtools index ${8}/EGFR/${2}.egfr.bam
					echo -e "${1}\t${2}\t${8}/EGFR/${2}.egfr.bam\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}"
 					echo -e "${1}\t${2}\t${8}/EGFR/${2}.egfr.bam\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}" >> ${8}/EGFR/prep.inhouse.output.txt

				fi
			else
				samtools view -bh ${5}/${4}*recal.bam \
					chr7:55019017-55211628 > ${8}/EGFR/${4}.egfr.bam
				samtools index ${8}/EGFR/${2}.egfr.bam
				if [ ! -e ${3}/${2}*recal.bam ];then
					echo "Normal is not existed"
					echo -e "${1}\t${2}\tNA\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}"
					echo -e "${1}\t${2}\tNA\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}" >> ${8}/EGFR/prep.inhouse.output.txt
				else
					samtools view -bh ${3}/${2}*recal.bam \
						chr7:55019017-55211628 > ${8}/EGFR/${2}.egfr.bam
					samtools index ${8}/EGFR/${2}.egfr.bam
					echo -e "${1}\t${2}\t${8}/EGFR/${2}.egfr.bam\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}"
					echo -e "${1}\t${2}\t${8}/EGFR/${2}.egfr.bam\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}" >> ${8}/EGFR/prep.inhouse.output.txt
				fi
			fi
		else
			if [ ${6} = 'hg19' ];then
				samtools view -bh ${8}/${1}/${4}*recal.bam \
					chr7:55086710-55279321 > ${8}/EGFR/${4}.egfr.bam
				samtools index ${8}/EGFR/${4}.egfr.bam
				if [ ! -e ${8}/${2}*recal.bam ];then
					echo "Normal is not existed"
					echo -e "${1}\t${2}\tNA\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}"
					echo -e "${1}\t${2}\tNA\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}" >> ${8}/EGFR/prep.inhouse.output.txt
				else
					samtools view -bh ${8}/${1}/${2}*recal.bam \
						chr7:55086710-55279321 > ${8}/EGFR/${2}.egfr.bam
					samtools index ${8}/EGFR/${2}.egfr.bam
					echo -e "${1}\t${2}\t${8}/EGFR/${2}.egfr.bam\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}"
 					echo -e "${1}\t${2}\t${8}/EGFR/${2}.egfr.bam\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}" >> ${8}/EGFR/prep.inhouse.output.txt
				fi

			else
				samtools view -bh ${8}/${1}/${4}*recal.bam \
					chr7:55019017-55211628 > ${8}/EGFR/${4}.egfr.bam
				samtools index ${8}/EGFR/${2}.egfr.bam
				if [ ! -e ${8}/${2}*recal.bam ];then
					echo "Normal is not existed"
 					echo -e "${1}\t${2}\tNA\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}"
					echo -e "${1}\t${2}\tNA\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}" >> ${8}/EGFR/prep.inhouse.output.txt
				else
					samtools view -bh ${8}/${1}/${2}*recal.bam \
						chr7:55019017-55211628 > ${8}/EGFR/${2}.egfr.bam
					samtools index ${8}/EGFR/${2}.egfr.bam
					echo -e "${1}\t${2}\t${8}/EGFR/${2}.egfr.bam\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}"
 					echo -e "${1}\t${2}\t${8}/EGFR/${2}.egfr.bam\t${4}\t${8}/EGFR/${4}.egfr.bam\t${6}\t${7}\t${11}" >> ${8}/EGFR/prep.inhouse.output.txt
				fi
			fi
		fi
		echo "## Extracting EGFR Location End - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Extracting EGFR Location End - `date` ##"
	else
		echo "## Extract EGFR IS Skipped ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Extract EGFR IS Skipped ##"
	fi
}

#--------------------------------------------------------------------------------------------
# STEP05.Collect Alignment and Insert Size Metrics
#
# CollectAlignmentSummaryMetrics: Produces a summary of alignment metics from a SAM or BAM file.
# This tool produces metics detiling the qulity of the rad alignments as well as the proportion 
# of the reads that passed machine signal-to-noise threshold quality filters. Note that these 
# quality filters are specific to Illumina data.
#--------------------------------------------------------------------------------------------
CollectAlgn(){
	if [ ${9} = 0 ];then
		echo "## STEP05.Collect Alignment and Insert Size Metrics ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Alignment collect & Insert Size Check START - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## STEP05.Collect Alignment and Insert Size Metrics ##"
		echo "## Alignment collect & Insert Size Check START - `date` ##"
		if [ ${4} != 'NA' ] && [ -e ${8}/${1}/${4}.recal.bam ];then
			#Tumor
			${gatk}/gatk CollectAlignmentSummaryMetrics \
				R=${ref}/${6}/${6}.fa \
				I=${8}/${1}/${4}.recal.bam \
				O=${8}/${1}/report/${4}.aln_metrics.txt
			
			${gatk}/gatk CollectInsertSizeMetrics \
				INPUT=${8}/${1}/${4}.recal.bam \
				OUPUT=${8}/${1}/${4}.insert_size.txt \
				HISTOGRAM_FILE=${8}/${1}/report/${4}.insert_size_hist.pdf
		else
			echo "## Collect Alignment From Tumor Sample Is Skipped ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Collect Alignment From Tumor Sample Is Skipped ##"
		fi
	
		if [ ${2} != 'NA' ] && [ -e ${8}/${1}/${2}.recal.bam ];then
			#Normal
			${gatk}/gatk CollectAlignmentSummaryMetrics \
				R=${ref}/${6}/${6}.fa \
				I=${8}/${1}/${2}.recal.bam \
				O=${8}/${1}/report/${2}.aln_metrics.txt
			
			${gatk}/gatk CollectInsertSizeMetrics \
				INPUT=${8}/${1}/${2}.recal.bam \
				OUPUT=${8}/${1}/${2}.insert_size.txt \
				HISTOGRAM_FILE=${8}/${1}/report/${2}.insert_size_hist.pdf
		else
			echo "## Collect Alignment From Normal Sample Is Skipped ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Collect Alignment From Normal Sample Is Skipped ##"
		fi

		echo "## Alignment collect & Insert Size Check END - `date` ##" >> ${8}/log/${1}_Tlog.txt

	else
		echo "## STEP05.Collect Alignment and Insert Size Metrics on CHR${9} region ##" >> ${8}/log/${1}_Tlog.txt
		echo "## Alignment collect & Insert Size Check START on CHR${9} region - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo "## STEP05.Collect Alignment and Insert Size Metrics on CHR${9} region ##"
		echo "## Alignment collect & Insert Size Check START on CHR${9} region - `date` ##"


		if [ ${4} != 'NA' ] && [ -e ${8}/${1}/CHR${9}/${4}.recal.bam ];then
			#Tumor
			${gatk}/gatk CollectAlignmentSummaryMetrics \
				R=${ref}/${6}/fasta_chr${9}.fa \
				I=${8}/${1}/CHR${9}/${4}.recal.bam \
				O=${8}/${1}/report/${4}.chr${9}.aln_metrics.txt
			${gatk}/gatk CollectInsertSizeMetrics \
				INPUT=${8}/${1}/CHR${9}/${4}.recal.bam \
				OUPUT=${8}/${1}/CHR${9}/${4}.insert_size.txt \
				HISTOGRAM_FILE=${8}/${1}/report/${4}.chr${9}.insert_size_hist.pdf
		else
			echo "## Collect Alignment From Tumor Sample CHR${9} Region Is Skipped ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Collect Alignment From Tumor Sample CHR${9} Region Is Skipped ##"

		fi

		if [ ${2} != 'NA' ] && [ -e ${8}/${1}/CHR${9}/${2}.recal.bam ];then
			#Normal
			${gatk}/gatk CollectAlignmentSummaryMetrics \
				R=${ref}/${6}/fasta_chr${9}.fa \
				I=${8}/${1}/CHR${9}/${2}.recal.bam \
				O=${8}/${1}/report/${2}.chr${9}.aln_metrics.txt
			
			${gatk}/gatk CollectInsertSizeMetrics \
				INPUT=${8}/${1}/CHR${9}/${2}.recal.bam \
				OUPUT=${8}/${1}/CHR${9}/${2}.insert_size.txt \
				HISTOGRAM_FILE=${8}/${1}/report/${2}.chr${9}.insert_size_hist.pdf
		else
			echo "## Collect Alignment From Normal Sample CHR${9} Region Is Skipped ##" >> ${8}/log/${1}_Tlog.txt
			echo "## Collect Alignment From Normal Sample CHR${9} Region Is Skipped ##"


		fi
		
		echo "## Alignment collect & Insert Size Check on CHR${9} region END - `date` ##" >> ${8}/log/${1}_Tlog.txt
	fi
}

#--------------------------------------------------------------------------------------------
# Running
#--------------------------------------------------------------------------------------------


echo "## Whole CHR will be mapped ##" >> ${8}/log/${1}_Tlog.txt
if [ -e ${ref}/${6}/${6}.dict ];then
	echo "## Indexed Reference exists ##" >> ${8}/log/${1}_Tlog.txt
	echo "Reference Indexing is skipped"
else
	IndexRef ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11}
fi

if [ ${4} != 'NA' ] && [ ! -e ${5}/${4}_1.fastq.gz ] && [ -e ${5}/${4}.tar ];then
	echo "## Uncompressing Sample ${4}.tar is required ##" >> ${8}/log/${1}_Tlog.txt
	echo "## Uncompressing Sample ${4}.tar is required ##"

elif [ ${2} != 'NA' ] && [ ! -e ${3}/${2}_1.fastq.gz ] && [ -e ${3}/${2}.tar ];then
	echo "## Uncompressing Sample ${2}.tar is required ##" >> ${8}/log/${1}_Tlog.txt
	echo "## Uncompressing Sample ${2}.tar is required ##"

else
	if [ ${12} == 'Y' ];then
		Trimming ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
		Mapping ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
		MarkDup ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
		Recal ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
		ExtractChr ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
		#CollectAlgn ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
	else
#		if [ -e ${5}/${4}.recal.bam ] && [ -e ${3}/${2}.recal.bam ];then
#			echo "## WES Preprocess Step is not required - `date` ##" >> ${8}/log/${1}_Tlog.txt
#			echo -e "${1}\t${2}\t${3}/${2}.recal.bam\t${4}\t${5}/${4}.recal.bam\t${6}\t${7}\t${11}"
#			echo -e "${1}\t${2}\t${3}/${2}.recal.bam\t${4}\t${5}/${4}.recal.bam\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.output.txt
		if [ -e ${8}/${1}/${4}.recal.bam ] && [ -e ${8}/${1}/${2}.recal.bam ];then
			echo "## WES Preprocessing Step is not required - `date` ##" >> ${8}/log/${1}_Tlog.txt
			echo "Preprocessing is skipped"
			echo -e "${1}\t${2}\t${8}/${1}/${2}.recal.bam\t${4}\t${8}/${1}/${4}.recal.bam\t${6}\t${7}\t${11}"
			echo -e "${1}\t${2}\t${8}/${1}/${2}.recal.bam\t${4}\t${8}/${1}/${4}.recal.bam\t${6}\t${7}\t${11}" >> ${8}/prep.inhouse.output.txt
		else
			Trimming ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
			Mapping ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
			MarkDup ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
			Recal ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
			ExtractChr ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
			#CollectAlgn ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}
		fi
	fi
	if [ -e ${8}/EGFR/${4}.egfr.bam ] && [ -e ${8}/EGFR/${2}.recal.bam ];then
		echo "## Extraction already is done ##" >> ${8}/log/${1}_Tlog.txt
	else
		ExtractChr ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10}
	fi
	
#	rm ${8}/${1}/*.trim.fastq.gz
#	rm ${8}/${1}/*.aln.bam
#	rm ${8}/${1}/*.sorted.bam
#	rm ${8}/${1}/*.mark_dup.bam

	chmod 755 ${8}/${1}/*
	
	if [ -d ${8}/${1}/CHR${9} ];then
#		rm ${8}/${1}/*.trim.fastq.gz
#		rm ${8}/${1}/CHR${9}/*.nrg.aln.bam
#		rm ${8}/${1}/CHR${9}/*.aln.bam
#		rm ${8}/${1}/CHR${9}/*.sorted.bam
#		rm ${8}/${1}/CHR${9}/*.mark_dup.bam
		chmod 755 ${8}/${1}/CHR${9}/*
	else
		chmod 755 ${8}/${1}/*
	fi
fi

chmod 755 ${8}/${1}
chmod 755 ${8}/${1}/*
