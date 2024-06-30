#!/bin/bash
#script by GT.Lee Jan.2.2023 
#In clinvar vcf file, chromoseom regions are represented as integers without chromosome letters
#FACETS is required to fit region representation between bam files and vcf files
#Bam and vcf files should be sorted before running
#Before start, conda activate facets should be rquired
ref=/DATA07/home/gangtl22/PIPE/script/ref
konda=/usr/anaconda3/bin

echo "This is WES Somatic Variant Calling &Filtration Process by FACETS"
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
echo "11. Librayr Kit Name: ${11}"
if [ ${12} == 'Y' ];then
	echo "12. Force to process this step: Yes"
else
	echo "12. Force to process this step: No"
fi
echo "SCNV Calling - FACETS"
sleep 3s


vcf=${ref}/INHOUSE/${6}_vcf/rev.clinvar.vcf.gz

echo "conda activate facets"
source ${konda}/activate cnv_facets

if [ ! -e ${8}/${1}/CNV/ ];then
	mkdir -p ${8}/${1}/CNV/
fi

cnv_dir=${8}/${1}/CNV

if [ ! -e ${8}/log/${1}_Tlog.txt ];then
	echo "## CNV Detection by FACETS Start -`date` ##"
else
	echo "## CNV Detection by FACETS Start -`date` ##" >> ${8}/log/${1}_Tlog.txt
fi

if [ -e ${8}/${1}/CNV/${1}.vcf.gz ] && [ ${12} == 'N' ];then
	echo "## CNV Calling File is already done  ##"
	echo "## CNV Calling File is already done  ##" >> ${8}/log/${1}_Tlog.txt
else
	if [ -e ${5} ] && [ -e ${3} ];then
		cnv_facets.R -t ${5} \
			-n ${3} \
			-vcf ${vcf}\
			-o ${cnv_dir}/${1} \
			--gbuild ${6}

	else
		echo "## Either Tumor.recal.bam or Normal.recal.bam at raw data directory and output directory ##"
		echo "## Missing Files exists. Check Input Files ##"
	fi
fi
if [ -e ${8}/${1}/CNV/${1}.vcf.gz ];then
	echo -e "${1}\t${2}\t${cnv_dir}/${1}.vcf.gz\t${4}\t${cnv_dir}/${1}.vcf.gz\t${6}\t${7}\t${11}"
	echo -e "${1}\t${2}\t${cnv_dir}/${1}.vcf.gz\t${4}\t${cnv_dir}/${1}.vcf.gz\t${6}\t${7}\t${11}" >> ${8}/facets.output.txt
	echo "## CNV Detection by FACETS End -`date` ##"
	echo "## CNV Detection by FACETS End -`date` ##" >> ${8}/log/${1}_Tlog.txt
else
	echo "## Error occurred while CNV Calling Step##"
	echo "## Error occurred while CNV Calling Step##" >> ${8}/log/${1}_Tlog.txt
fi

chmod 755 ${8}/${1}/CNV/
chmod 755 ${8}/${1}/CNV/*
source ${konda}/deactivate
