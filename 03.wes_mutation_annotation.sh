#!/bin/bash
bwa="/DATA03/Tools/bwa-0.7.17"
gatk="python3.8 /DATA03/Tools/GATK/gatk-4.3.0.0"
gatk_i="java -jar /DATA03/Tools/GATK/gatk-4.1.4.0/gatk-package-4.1.4.0-local.jar"
trim="java -jar /DATA03/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar"
picard="java --java-options "-Xmx4G" -jar /DATA03/Tools/picard-2.21.1/picard.jar"
ref=/DATA07/home/gangtl22/PIPE/script/ref/INHOUSE

echo "This is WES Somatic Variant Annotating Process by Funcotator"
echo "1. Sample Prefix: ${1}"
echo "2. Normal Prefix: ${2}"
echo "3. /Path/To/Filtered/VCF: ${3}"
echo "4. Tumor Prefix: ${4}"
echo "5. /Path/To/Filtered/VCF: ${5}"
echo "6. Reference Genome: ${6}"
echo "7. EGFT_MT: ${7}"
echo "8. Output directory: ${8}"
if [ ${9} -eq 0 ];then
	echo "9. Analysis whole region"
else
	echo "9. Anlaysis Chromosome ${9} region"
fi
echo "10. Number of threads: ${10}"
echo "11. Librayr Kit Name: ${11}"
if [ ${12} == 'Y' ];then
	echo "12. Force to process this step: Yes"
else
	echo "12. Force to process this step: No"
fi
sleep 3s
if [ ! -e ${8}/${1} ];then
	mkdir -p ${8}/${1}/
	mkdir -p ${8}/log
fi

if [ ! -e ${8}/funcotator.maf.output.txt ];then
	echo -e "Patient_Num(Sample_Prefix)\tNormal_ID\tNormal_Path\tTumor_ID\tTumor_Path\tRef_Genome\tEGFR_MT\tLib_Kit" >> ${8}/funcotator.maf.output.txt
        chmod 755 ${8}/funcotator.maf.output.txt
fi

if [ ! -e ${8}/funcotator.vcf.output.txt ];then
	echo -e "Patient_Num(Sample_Prefix)\tNormal_ID\tNormal_Path\tTumor_ID\tTumor_Path\tRef_Genome\tEGFR_MT\tLib_Kit" >> ${8}/funcotator.vcf.output.txt
	chmod 755 ${8}/funcotator.vcf.output.txt
fi

#For somatic data sources:
#${gatk}/gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download

#For germline data sources:
#${gatk}/gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download

#-------------------------------------------------------------------------------------------
# Annoate Variants - GATK4 Funcotator
#--------------------------------------------------------------------------------------------
Funcotator(){
	if [ ! -e ${8}/${1}/VCF ];then
		mkdir -p ${8}/${1}/VCF
	fi
	if [ ! -e ${8}/${1}/MAF ];then
		mkdir -p ${8}/${1}/MAF
	fi

	vf="${3%.gz}"
	
	echo "${vf}"
	
	sleep 3s
	
	if [ ! -e ${vf} ];then
		gunzip -dk ${3}
	fi
	if [ ! -e ${vf}.idx ];then
		${gatk_i} IndexFeatureFile -F ${vf}
	fi
 
	if [ ${9} = 0 ];then
		#Annotate using Funcotator
		echo "Variants Annotation - gatk4 Funcotator"
		echo "## Varant Annotation by Funcotator Starts - `date` ##" >> ${8}/log/${1}_Tlog.txt
		if [ -e ${8}/${1}/VCF/${1}.funcotator.vcf ] && [ ${12} == 'N' ];then
			
			echo "## Funcotator Annotation Getting VCF is done ##"

		elif [ ${12} == 'Y' ];then
			
			${gatk}/gatk --java-options '-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' Funcotator \
				--variant ${vf} \
				--reference ${ref}/${6}/${6}.fa \
				--ref-version ${6} \
				--data-sources-path /DATA03/Tools/funcotator/funcotator_dataSources.v1.7.20200521s \
				--output ${8}/${1}/VCF/${1}.funcotator.vcf \
				--output-file-format VCF
		
		else
			${gatk}/gatk --java-options '-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' Funcotator \
				--variant ${vf} \
				--reference ${ref}/${6}/${6}.fa \
				--ref-version ${6} \
				--data-sources-path /DATA03/Tools/funcotator/funcotator_dataSources.v1.7.20200521s \
				--output ${8}/${1}/VCF/${1}.funcotator.vcf \
				--output-file-format VCF
		fi

		if [ -e ${8}/${1}/MAF/${1}.funcotator.maf ] && [ ${12} == 'N' ];then
			
			echo "## Funcotator Annotation Getting MAF is done ##"

		elif [ ${12} == 'Y' ];then
	
			${gatk}/gatk --java-options '-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' Funcotator \
				--variant ${vf} \
				--reference ${ref}/${6}/${6}.fa \
				--ref-version ${6} \
				--data-sources-path /DATA03/Tools/funcotator/funcotator_dataSources.v1.7.20200521s \
				--output ${8}/${1}/MAF/${1}.funcotator.maf \
				--output-file-format MAF
		else
			${gatk}/gatk --java-options '-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' Funcotator \
				--variant ${vf} \
				--reference ${ref}/${6}/${6}.fa \
				--ref-version ${6} \
				--data-sources-path /DATA03/Tools/funcotator/funcotator_dataSources.v1.7.20200521s \
				--output ${8}/${1}/MAF/${1}.funcotator.maf \
				--output-file-format MAF
		fi
		
		echo "## Varant Annotation (VCF) by  Funcotator Ends - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo -e "${1}\t${2}\t${8}/${1}/VCF/${1}.funcotator.vcf\t${4}\t${8}/${1}/VCF/${1}.funcotator.vcf\t${6}\t${7}\t${11}"
		echo -e "${1}\t${2}\t${8}/${1}/VCF/${1}.funcotator.vcf\t${4}\t${8}/${1}/VCF/${1}.funcotator.vcf\t${6}\t${7}\t${11}" >> ${8}/funcotator.vcf.output.txt
		echo "## Varant Annotation (MAF) by  Funcotator Ends - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo -e "${1}\t${2}\t${8}/${1}/MAF/${1}.funcotator.maf\t${4}\t${8}/${1}/MAF/${1}.funcotator.maf\t${6}\t${7}\t${11}"
		echo -e "${1}\t${2}\t${8}/${1}/MAF/${1}.funcotator.maf\t${4}\t${8}/${1}/MAF/${1}.funcotator.maf\t${6}\t${7}\t${11}" >> ${8}/funcotator.maf.output.txt
	else
		#Annotate using Funcotator
		echo "Variants Annotation - gatk4 Funcotator"
		echo "## Varant on CHR${9} region Annotation by Funcotator Starts - `date` ##" >> ${8}/log/${1}_Tlog.txt
	
		if [ -e ${8}/${1}/CHR${9}/VCF/${1}.funcotator.vcf ] && [ ${12} == 'N' ];then
	
			echo "## VEP Annotation Getting MAF is done ##"

		elif [ ${12} == 'Y' ];then
			
			${gatk}/gatk --java-options '-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' Funcotator \
				--variant ${vf} \
				--reference ${ref}/${6}/${6}.fa \
				--ref-version ${6} \
				--data-sources-path /DATA03/Tools/funcotator/funcotator_dataSources.v1.7.20200521s \
				--output ${8}/${1}/CHR${9}/VCF/${1}.funcotator.vcf \
				--output-file-format VCF
		
		else

			${gatk}/gatk --java-options '-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' Funcotator \
				--variant ${vf} \
				--reference ${ref}/${6}/${6}.fa \
				--ref-version ${6} \
				--data-sources-path /DATA03/Tools/funcotator/funcotator_dataSources.v1.7.20200521s \
				--output ${8}/${1}/CHR${9}/VCF/${1}.funcotator.vcf \
				--output-file-format VCF
		
		fi
		
		if [ -e ${8}/${1}/CHR${9}/MAF/${1}.funcotator.vcf ] && [ ${12} == 'N' ];then
	
			echo "## VEP Annotation Getting MAF is done ##"

		elif [ ${12} == 'Y' ];then
			
			${gatk}/gatk --java-options '-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' Funcotator \
 				--variant ${vf} \
				--reference ${ref}/${6}/${6}.fa \
				--ref-version ${6} \
				--data-sources-path /DATA03/Tools/funcotator/funcotator_dataSources.v1.7.20200521s \
				--output ${8}/${1}/CHR${9}/MAF/${1}.funcotator.maf \
				--output-file-format MAF
		else

			${gatk}/gatk --java-options '-Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true' Funcotator \
				--variant ${vf} \
				--reference ${ref}/${6}/${6}.fa \
				--ref-version ${6} \
				--data-sources-path /DATA03/Tools/funcotator/funcotator_dataSources.v1.7.20200521s \
				--output ${8}/${1}/CHR${9}/MAF/${1}.funcotator.maf \
				--output-file-format MAF
		fi

		echo "## Varant on CHR${9} region Annotation (VCF) by Funcotator Ends - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/VCF/${1}.funcotator.vcf\t${4}\t${8}/${1}/CHR${9}/VCF/${1}.funcotator.vcf\t${6}\t${7}\t${11}"
		echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/VCF/${1}.funcotator.vcf\t${4}\t${8}/${1}/CHR${9}/VCF/${1}.funcotator.vcf\t${6}\t${7}\t${11}" >> ${8}/funcotator.vcf.output.txt
		echo "## Varant on CHR${9} region Annotation (VCF) by Funcotator Ends - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/MAF/${1}.funcotator.maf\t${4}\t${8}/${1}/CHR${9}/MAF/${1}.funcotator.maf\t${6}\t${7}\t${11}"
		echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/MAF/${1}.funcotator.maf\t${4}\t${8}/${1}/CHR${9}/MAF/${1}.funcotator.maf\t${6}\t${7}\t${11}" >> ${8}/funcotator.maf.output.txt
	fi

}

Funcotator ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11}

chmod 755 ${8}/${1}/VCF/
chmod 755 ${8}/${1}/VCF/*
chmod 755 ${8}/${1}/MAF/
chmod 755 ${8}/${1}/MAF/*
if [ -e ${8}/${1}/CHR${9} ];then
	chmod 755 ${8}/${1}/CHR${9}/VCF/
	chmod 755 ${8}/${1}/CHR${9}/VCF/*
	chmod 755 ${8}/${1}/CHR${9}/MAF/
	chmod 755 ${8}/${1}/CHR${9}/MAF/*
fi
