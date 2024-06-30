#!/bin/bash
gatk="python3.8 /DATA03/Tools/GATK/gatk-4.3.0.0"
ref=/DATA07/home/gangtl22/PIPE/script/ref/INHOUSE

echo "This is WES Somatic Variant Calling Process by GATK4"
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
sleep 3s

#--------------------------------------------------------------------------------------------
# STEP06. Somatic Variant Calling
# Variant calls are reported by each pipeline in a VCF formatted file. 
#--------------------------------------------------------------------------------------------
#MuTect2
echo "Somatic Variant Calling - MuTect2"
if [ ! -e ${8}/${1} ];then
	mkdir -p ${8}/${1}/
	chmod 755 ${8}/${1}/
	mkdir -p ${8}/log
fi

if [ ! -e ${8}/filter.output.txt ];then
	echo -e "Patient_Num(Sample_Prefix)\tNormal_ID\tNormal_Path\tTumor_ID\tTumor_Path\tRef_Genome\tEGFR_MT\tLib_Kit" >> ${8}/filter.output.txt
	chmod 755 ${8}/filter.output.txt
fi


Mut2(){
	if [ ${6} = hg19 ];then
		gnome=${ref}/hg19_vcf/conv_af-only-gnomad.raw.sites.vcf.gz
		pon=${ref}/hg19_vcf/conv_Mutect2-exome-panel.vcf.gz
		common=${ref}/hg19_vcf/rev_small_exac_common_3.vcf
	else
		gnome=${ref}/hg38_vcf/af-only-gnomad.hg38.vcf.gz
		pon=${ref}/hg38_vcf/1000g_pon.hg38.vcf.gz
		common=${ref}/hg38_vcf/small_exac_common_3.hg38.vcf.gz
	fi

	if [ ${9} = 0 ];then
		echo "STEP06.Somatic Mutation Calling"
		echo "## Somatic Mutation Calling START - `date` ##" >> ${8}/log/${1}_Tlog.txt
		if [ -e ${5} ];then
			if [ ${2} = 'NA' ];then 
				echo "## Tumor-only mode ##"
				sleep 2s

				if [ -e ${8}/${1}/${1}.raw.vcf.gz ] && [ ${12} = N ];then
					echo "## Raw VCF file already exists ##"
				else	
					${gatk}/gatk Mutect2 \
						-R ${ref}/${6}/${6}.fa \
						-I ${5} \
						--germline-resource ${gnome} \
						--panel-of-normals ${pon} \
						-O ${8}/${1}/${1}.raw.vcf.gz
				fi
				
				##Summarizes counts of rads that support reference (Tumor)
				if [ -e ${8}/${1}/tumor.pileups.table ] && [ ${12} = N ];then
					echo "## Getting Pileup Summary from Tumor step is passed ##"
				else
					${gatk}/gatk GetPileupSummaries \
						-I ${5} \
						-V ${common} \
						-L ${common} \
						-O ${8}/${1}/tumor.pileups.table
				fi
				
				##Calculate Contamination
				if [ -e ${8}/${1}/${1}.contamination.table ] && [ ${12} = N ];then
					echo "## Calculating Contamination step is passed ##"
				else
					${gatk}/gatk CalculateContamination \
						-I ${8}/${1}/tumor.pileups.table \
						-O ${8}/${1}/${1}.contamination.table
				fi

				##Filtering VCF file
				if [ -e ${8}/${1}/${8}/${1}/${1}.filtered.vcf.gz ] && [ ${12} = N ];then
					echo "## Filtering step is passed ##"
				else
					${gatk}/gatk FilterMutectCalls \
						-R ${ref}/${6}/${6}.fa \
						-V ${8}/${1}/${1}.raw.vcf.gz \
						--contamination-table ${8}/${1}/${1}.contamination.table \
						-O ${8}/${1}/${1}.filtered.vcf.gz
				fi

			else
				echo "## Tumor with matched normal ##"
				sleep 2s

				if [ -e ${8}/${1}/${1}.raw.vcf.gz ] && [ ${12} = N ];then
					echo "## Raw VCF file already exists ##"
				else
					${gatk}/gatk Mutect2 \
						-R ${ref}/${6}/${6}.fa \
						-I ${5} \
						-I ${3} \
						--tumor-sample ${4} \
						--normal-sample ${2} \
						--germline-resource ${gnome} \
						--panel-of-normals ${pon} \
						-O ${8}/${1}/${1}.raw.vcf.gz
				fi

				##Summarizes counts of rads that support reference (Tumor)
				if [ -e ${8}/${1}/tumor.pileups.table ] && [ ${12} = N ];then
					echo "## Getting Pileup Summary from Tumor step is passed ##"
				else
					${gatk}/gatk GetPileupSummaries \
						-I ${5} \
						-V ${common} \
						-L ${common} \
						-O ${8}/${1}/tumor.pileups.table
				fi

				##Summarizes counts of rads that support reference (Normal)
				if [ -e ${8}/${1}/normal.pileups.table ] && [ ${12} = N ];then
					echo "## Getting Pileup Summary from Normal step is passed ##"
				else
					${gatk}/gatk GetPileupSummaries \
					     	-I ${3} \
					     	-V ${common} \
					     	-L ${common} \
					     	-O ${8}/${1}/normal.pileups.table
				fi

				##Calculate Contamination
				if [ -e ${8}/${1}/${1}.contamination.table ] && [ ${12} = N ];then
					echo "## Calculating Contamination step is passed ##"
				else
					${gatk}/gatk CalculateContamination \
						-I ${8}/${1}/tumor.pileups.table \
						-matched ${8}/${1}/normal.pileups.table \
						-O ${8}/${1}/${1}.contamination.table
				fi
				
				##Filtering VCF file
				if [ -e ${8}/${1}/${8}/${1}/${1}.filtered.vcf.gz ] && [ ${12} = N ];then
					echo "## Filtering step is passed ##"
				else
					${gatk}/gatk FilterMutectCalls \
						-R ${ref}/${6}/${6}.fa \
						-V ${8}/${1}/${1}.raw.vcf.gz \
						--contamination-table ${8}/${1}/${1}.contamination.table \
						-O ${8}/${1}/${1}.filtered.vcf.gz
				fi
			fi
		fi
		echo "## Somatic Mutation Calling END - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo -e "${1}\t${2}\t${8}/${1}/${1}.raw.vcf.gz\t${4}\t${8}/${1}/${1}.raw.vcf.gz\t${6}\t${7}\t${11}"
		echo -e "${1}\t${2}\t${8}/${1}/${1}.raw.vcf.gz\t${4}\t${8}/${1}/${1}.raw.vcf.gz\t${6}\t${7}\t${11}" >> ${8}/mutcall.output.txt
		echo -e "${1}\t${2}\t${8}/${1}/${1}.filtered.vcf.gz\t${4}\t${8}/${1}/${1}.filtered.vcf.gz\t${6}\t${7}\t${11}"
		echo -e "${1}\t${2}\t${8}/${1}/${1}.filtered.vcf.gz\t${4}\t${8}/${1}/${1}.filtered.vcf.gz\t${6}\t${7}\t${11}" >> ${8}/filter.output.txt

	else
		echo "STEP06.Somatic Mutation Calling"
		echo "## Somatic Mutation calling START on CHR${9} region - `date` ##" >> ${8}/log/${1}_Tlog.txt
		mkdir -p ${8}/${1}/CHR${9}
		if [ ${2} = 'NA' ];then
			echo "## Tumor-only mode ##"
			sleep 2s
			${gatk}/gatk Mutect2 \
				-R ${ref}/${6}/fasta_chr${9}.fa \
				-I ${5} \
				--germline-resource ${ref}/${6}_vcf/conv_af-only-gnomad.raw.sites.vcf.gz \
				--panel-of-normals ${ref}/${6}_vcf/conv_Mutect2-exome-panel.vcf.gz \
				-O ${8}/${1}/CHR${9}/${1}.raw.vcf.gz
			
			##Summarizes counts of rads that support reference (Tumor)
			${gatk}/gatk GetPileupSummaries \
				-I ${5} \
				-V ${common} \
				-L ${common} \
				-O ${8}/${1}/CHR${9}/tumor.pileups.table
			
			##Calculate Contamination
			${gatk}/gatk CalculateContamination \
				-I ${8}/${1}/CHR${9}/tumor.pileups.table \
				-O ${8}/${1}/CHR${9}/${1}.contamination.table
			
			##Filtering VCF file
			${gatk}/gatk FilterMutectCalls \
				-R ${ref}/${6}/fasta_chr${9}.fa \
				-V ${8}/${1}/CHR${9}/${1}.raw.vcf.gz \
				--contamination-table ${8}/${1}/CHR${9}/${1}.contamination.table \
				-O ${8}/${1}/CHR${9}/${1}.filtered.vcf.gz
		else
			echo "## Tumor with matched normal ##"
			sleep 2s

			${gatk}/gatk Mutect2 \
				-R ${ref}/${6}/fasta_chr${9}.fa \
				-I ${8}/${1}/CHR${9}/${5} \
				-I ${8}/${1}/CHR${9}/${3} \
				--tumor-sample ${4} \
				--normal-sample ${2} \
				--germline-resource ${ref}/${6}_vcf/conv_af-only-gnomad.raw.sites.vcf.gz \
				--panel-of-normals ${ref}/${6}_vcf/conv_Mutect2-exome-panel.vcf.gz \
				-O ${8}/${1}/CHR${9}/${1}.raw.vcf.gz

			##Summarizes counts of rads that support reference (Tumor)
			${gatk}/gatk GetPileupSummaries \
				-I ${5} \
				-V ${common} \
				-L ${common} \
				-O ${8}/${1}/CHR${9}/tumor.pileups.table
	
			##Summarizes counts of rads that support reference (Normal)
			${gatk}/gatk GetPileupSummaries \
				-I ${3} \
				-V ${common} \
				-L ${common} \
				-O ${8}/${1}/CHR${9}/normal.pileups.table
			
			##Calculate Contamination
			${gatk}/gatk CalculateContamination \
				-I ${8}/${1}/CHR${9}/tumor.pileups.table \
				-matched ${8}/${1}/CHR${9}/normal.pileups.table \
				-O ${8}/${1}/CHR${9}/${1}.contamination.table

			##Filtering VCF file
			${gatk}/gatk FilterMutectCalls \
				-R ${ref}/${6}/fasta_chr${9}.fa \
				-V ${8}/${1}/CHR${9}/${1}.raw.vcf.gz \
				--contamination-table ${8}/${1}/CHR${9}/${1}.contamination.table \
				-O ${8}/${1}/CHR${9}/${1}.filtered.vcf.gz	

		fi

		echo "## Somatic Mutation Calling END - `date` ##" >> ${8}/log/${1}_Tlog.txt
		echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/${1}.raw.vcf.gz\t${4}\t${8}/${1}/CHR${9}/${1}.raw.vcf.gz\t${6}\t${7}\t${11}"
		echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/${1}.raw.vcf.gz\t${4}\t${8}/${1}/CHR${9}/${1}.raw.vcf.gz\t${6}\t${7}\t${11}" >> ${8}/mutcall.output.txt
		echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/${1}.filtered.vcf.gz\t${4}\t${8}/${1}/CHR${9}/${1}.filtered.vcf.gz\t${6}\t${7}\t${11}"
		echo -e "${1}\t${2}\t${8}/${1}/CHR${9}/${1}.filtered.vcf.gz\t${4}\t${8}/${1}/CHR${9}/${1}.filtered.vcf.gz\t${6}\t${7}\t${11}" >> ${8}/filter.output.txt

	fi
}



#--------------------------------------------------------------------------------------------
# Running
#--------------------------------------------------------------------------------------------

Mut2 ${1} ${2} ${3} ${4} ${5} ${6} ${7} ${8} ${9} ${10} ${11} ${12}

chmod 755 ${8}/${1}
chmod 755 ${8}/${1}/*
