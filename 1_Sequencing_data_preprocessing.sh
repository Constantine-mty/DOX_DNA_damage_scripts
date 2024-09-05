##------------------------------------------------RNA-seq-----------------------------------------------------------------

file_path="/home/user"
file_name_1="treatment_1"
file_name_2="treatment_2"
file_name_3="treatment_3"

#hisat2
ENSEMBL_RELEASE=84
ENSEMBL_GRCh38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

HISAT2_BUILD_EXE=./hisat2-build
if [ ! -x "$HISAT2_BUILD_EXE" ] ; then
	if ! which hisat2-build ; then
		echo "Could not find hisat2-build in current directory or in PATH"
		exit 1
	else
		HISAT2_BUILD_EXE=`which hisat2-build`
	fi
fi

rm -f genome.fa
F=Homo_sapiens.GRCh38.dna.primary_assembly.fa
if [ ! -f $F ] ; then
	get ${ENSEMBL_GRCh38_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	mv $F genome.fa
fi

CMD="${HISAT2_BUILD_EXE} -p 4 genome.fa genome"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi

#mapping
for i in file_name_1 file_name_2 file_name_3
do
	#hisat2做序列比对
	hisat2 -p 10 -t --dta --summary-file ${i}.txt -x /home/ze/data/hg38/grch38/genome \
	-1 ${file_path}/${i}_1.fq.gz \
	-2 ${file_path}/${i}_2.fq.gz \
	-S ${o}.sam
	#samtools转换格式和质量过滤
	samtools view -@ 10 -S ${i}.sam -b > ${i}.bam
	samtools sort -@ 10 ${i}.bam -o ${i}_sorted.bam
	samtools index -@ 10 ${file_path}/${i}_sorted.bam
	#stringtie为转录本拼接和定量的软件
	stringtie  -p 10 ${i}_sorted.bam\
	-A ${i}.tab\
	-C ${i}_cov_refs.gtf\
	-G /home/ze/data/gtf/Homo_sapiens.GRCh38.105.gtf\
	-o ${i}.gtf
done

##bamCoverage
for i in file_name_1 file_name_2 file_name_3
do
	bamCoverage -p 10 \
	--binSize 50 \
	--normalizeUsing CPM \
	--smoothLength 60 \
	-b ${file_path}/${i}_sorted.bam \
	-o ${i}_RNA.bw
done

#featurecounts
featurecounts -T 10 -p --countReadPairs -t exon -g gene_id -O \
 -a /home/ze/data/gtf/Homo_sapiens.GRCh38.105.gtf \
 -o all.id.txt \
 *_sorted.bam

##------------------------------------------------ATAC-seq-----------------------------------------------------------------

file_path="/home/user"
file_name_1="treatment_1"
file_name_2="treatment_2"
file_name_3="treatment_3"

for i in file_name_1 file_name_2 file_name_3
do
	#bowtie2
	bowtie2 -q -p 10 --no-unal -x /home/ze/data/hg38/bowtie_index/hg38 \
	-1 ${file_path}/${i}_1.fq.gz \
	-2 ${file_path}/${i}_2.fq.gz \
	-S ${file_path}/${i}.sam
	#samtools
	samtools view -@ 10 -S ${file_path}/${i}.sam -b > ${file_path}/${i}.bam
	samtools sort -@ 10 ${file_path}/${i}.bam -o ${file_path}/${i}_sorted.bam
	samtools view -@ 10 -b -f 2 -q 30 ${file_path}/${i}_sorted.bam -o ${file_path}/${i}_filtered.bam
	picard MarkDuplicates I=${file_path}/${i}_filtered.bam O=${file_path}/${i}_dedup.bam M=${file_path}/${i}_dedup_metrics.txt REMOVE_DUPLICATES=true
	samtools index ${file_path}/${i}_dedup.bam
	#bam to tagAlign
	samtools sort -@ 10 -n ${file_path}/${i}_dedup.bam -o ${file_path}/${i}_dedup_namesorted.bam
	bedtools bamtobed -bedpe -mate1 -i ${file_path}/${i}_dedup_namesorted.bam | gzip -nc > ${file_path}/${i}_dedup.bedpe.gz
	zcat ${file_path}/${i}_dedup.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${file_path}/${i}_dedup.tagAlign.gz
done

#macs2 callpeaks
for i in file_name_1 file_name_2 file_name_3
do
	macs2 callpeak -t ${file_path}/${i}_dedup.tagAlign.gz \
	-f BED -q 0.05 -B -g hs --nomodel --shift -100 --extsize 200 --keep-dup all --SPMR -n ${i}-ATAC
	cut -f 1-6 ${file_path}/${i}-ATAC_peaks.narrowPeak > ${file_path}/${i}-ATAC_peaks.bed
done
##bamCoverage
for i in file_name_1 file_name_2 file_name_3
do
	bamCoverage -p 10 \
	--normalizeUsing CPM \
	--binSize 50 \
	--smoothLength 60 \
	-b ${file_path}/${i}_dedup.bam \
	-o ${i}_ATAC.bw
done

#MAnorm
manorm --pe --wa \
 --p1 ${file_path}/${file_name_2}-ATAC_peaks.bed \
 --p2 ${file_path}/${file_name_1}-ATAC_peaks.bed \
 --pf bed \
 --r1 ${file_path}/${file_name_2}_dedup.bam \
 --r2 ${file_path}/${file_name_1}_dedup.bam \
 --rf bam \
 -o ${file_path}/${file_name_1}_compare_${file_name_2}

manorm --pe --wa \
 --p1 ${file_path}/${file_name_3}-ATAC_peaks.bed \
 --p2 ${file_path}/${file_name_1}-ATAC_peaks.bed \
 --pf bed \
 --r1 ${file_path}/${file_name_3}_dedup.bam \
 --r2 ${file_path}/${file_name_1}_dedup.bam \
 --rf bam \
 -o ${file_path}/${file_name_1}_compare_${file_name_3}

manorm --pe --wa \
 --p1 ${file_path}/${file_name_3}-ATAC_peaks.bed \
 --p2 ${file_path}/${file_name_2}-ATAC_peaks.bed \
 --pf bed \
 --r1 ${file_path}/${file_name_3}_dedup.bam \
 --r2 ${file_path}/${file_name_2}_dedup.bam \
 --rf bam \
 -o ${file_path}/${file_name_2}_compare_${file_name_3}

##------------------------------------------------ChIP-seq-----------------------------------------------------------------

file_path="/home/user"
file_name_1="treatment_1"
file_name_2="treatment_2"
file_name_3="treatment_3"

for i in file_name_1 file_name_2 file_name_3
do
	#bowtie2
	bowtie2 -q -p 10 --no-unal -x /home/ze/data/hg38/bowtie_index/hg38 \
	-1 ${file_path}/${i}_Input_1.fq.gz \
	-2 ${file_path}/${i}_Input_2.fq.gz \
	-S ${file_path}/${i}_Input.sam
	bowtie2 -q -p 10 --no-unal -x /home/ze/data/hg38/bowtie_index/hg38 \
	-1 ${file_path}/${i}_IP_1.fq.gz \
	-2 ${file_path}/${i}_IP_2.fq.gz \
	-S ${file_path}/${i}_IP.sam
	#samtools
	samtools view -@ 10 -S ${file_path}/${i}_Input.sam -b > ${file_path}/${i}_Input.bam
	samtools sort -@ 10 ${file_path}/${i}_Input.bam -o ${file_path}/${i}_Input_sorted.bam
	samtools view -@ 10 -b -f 2 -q 30 ${file_path}/${i}_Input_sorted.bam -o ${file_path}/${i}_Input_filtered.bam
	picard MarkDuplicates I=${file_path}/${i}_Input_filtered.bam O=${file_path}/${i}_Input_dedup.bam M=${file_path}/${i}_Input_dedup_metrics.txt REMOVE_DUPLICATES=true
	samtools index ${file_path}/${i}_Input_dedup.bam
	samtools view -@ 10 -S ${file_path}/${i}_IP.sam -b > ${file_path}/${i}_IP.bam
	samtools sort -@ 10 ${file_path}/${i}_IP.bam -o ${file_path}/${i}_IP_sorted.bam
	samtools view -@ 10 -b -f 2 -q 30 ${file_path}/${i}_IP_sorted.bam -o ${file_path}/${i}_IP_filtered.bam
	picard MarkDuplicates I=${file_path}/${i}_IP_filtered.bam O=${file_path}/${i}_IP_dedup.bam M=${file_path}/${i}_IP_dedup_metrics.txt REMOVE_DUPLICATES=true
	samtools index ${file_path}/${i}_IP_dedup.bam
done
##samtools merge
for i in file_name_1 file_name_2 file_name_3
do
	for q in Input IP
	do
		#bam to tagAlign
		samtools sort -@ 10 -n ${file_path}/${i}_${q}_merged.bam -o ${file_path}/${i}_${q}_dedup_namesorted.bam
		bedtools bamtobed -bedpe -mate1 -i ${file_path}/${i}_${q}_dedup_namesorted.bam | gzip -nc > ${file_path}/${i}_${q}_dedup.bedpe.gz
		zcat ${file_path}/${i}_${q}_dedup.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${file_path}/${i}_${q}_dedup.tagAlign.gz
	done
done
rm *_namesorted.bam
#macs2 callpeaks
for i in file_name_1 file_name_2 file_name_3
do      
	macs2 callpeak -t ${file_path}/${i}_IP_dedup.tagAlign.gz \
	-c ${file_path}/${i}_Input_dedup.tagAlign.gz \
	-f BED -B -g hs -q 0.05 --keep-dup all --SPMR -n ${i}-h3k27ac
	cut -f 1-6 ${file_path}/${i}-h3k27ac_peaks.narrowPeak > ${file_path}/${i}-H3K27ac_peaks.bed
done

##bowtie2
for i in file_name_1 file_name_2 file_name_3
do
	bamCompare -p 10 -b1 ${file_path}/${i}_IP_merged.bam \
	-b2 ${file_path}/${i}_Input_merged.bam \
	--binSize 50 \
	--scaleFactorsMethod None \
	--normalizeUsing CPM \
	--operation subtract \
	--smoothLength 60 \
	-o ${i}-H3K27ac.bw
done

#MAnorm
manorm --pe --wa \
 --p1 ${file_path}/${file_name_2}-H3K27ac_peaks.bed \
 --p2 ${file_path}/${file_name_1}-H3K27ac_peaks.bed \
 --pf bed \
 --r1 ${file_path}/${file_name_2}_dedup.bam \
 --r2 ${file_path}/${file_name_1}_dedup.bam \
 --rf bam \
 -o ${file_path}/${file_name_1}_compare_${file_name_2}

manorm --pe --wa \
 --p1 ${file_path}/${file_name_3}-H3K27ac_peaks.bed \
 --p2 ${file_path}/${file_name_1}-H3K27ac_peaks.bed \
 --pf bed \
 --r1 ${file_path}/${file_name_3}_dedup.bam \
 --r2 ${file_path}/${file_name_1}_dedup.bam \
 --rf bam \
 -o ${file_path}/${file_name_1}_compare_${file_name_3}

manorm --pe --wa \
 --p1 ${file_path}/${file_name_3}-H3K27ac_peaks.bed \
 --p2 ${file_path}/${file_name_2}-H3K27ac_peaks.bed \
 --pf bed \
 --r1 ${file_path}/${file_name_3}_dedup.bam \
 --r2 ${file_path}/${file_name_2}_dedup.bam \
 --rf bam \
 -o ${file_path}/${file_name_2}_compare_${file_name_3}

##------------------------------------------------Enhancer and gene filtering-----------------------------------------------------------------

file_path="/home/user"
file_name_1="treatment_1"
file_name_2="treatment_2"
file_name_3="treatment_3"

#dELS（distal Enhancer-like States）；pELS（proximal Enhancer-like States）
zcat ENCFF287MQD.bed.zip | awk '$10=="dELS"||$10=="dELS,CTCF-bound"||$10=="pELS"||$10=="pELS,CTCF-bound" {print}' ENCFF287MQD.bed > encode_ELS_hg38.bed
zcat ENCFF287MQD.bed.zip | awk '$10=="PLS"||$10=="PLS,CTCF-bound" {print}' > encode_PLS_hg38.bed

#gene TSS/TES  
sed 's/"/\t/g' /newdisk/jinyuan/p53/rna_seq/Homo_sapiens.GRCh38.109.gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="gene") {if($7=="+") {start=$4-1000; end=$5+1000;} else {if($7=="-") start=$4-1000; end=$5+1000; } if(start<0) start=0; print "chr"$1,start,end,$14,$10,$7;}}' >geneExtension.bed

sort -k1,1 -k2,2n geneExtension.bed > sorted_geneExtension.bed

bedtools intersect -a encode_ELS_hg38.bed -b sorted_geneExtension.bed -v > encode_enhancer_rmgeneoverlap.bed

#TSS
awk 'BEGIN{OFS="\t"} $6 == "+" {start=$2; end=$2+1} $6 == "-" {end=$3; start=$3-1} {print $1,start,end,$4,$5,$6}' ensembl_gene_region.bed > TSS.bed
bedtools intersect -a encode_PLS_hg38.bed -b TSS.bed -wa -wb > a.bed
awk 'BEGIN{OFS="\t"}  {strand=$17} {print $1,$2,$3,$4,$5,strand,$13,$16}' a.bed > look_promoter.bed
awk 'BEGIN{OFS="\t"}  {strand=$17} {print $1,$2,$3,$4,$5,strand}' a.bed > promoter.bed

bedtools intersect -a encode_ELS_hg38.bed -b ./GRCh38_unified_blacklist.bed -v > encode_eRNA_rmgeneoverlap_rmblacklist.bed

bedtools random -g hg38.genome.txt -n 10000 -l 4000 > random_region.bed
bedtools sort -i random_region.bed > sorted_random_region.bed

#multiBigwigSummary
for name in intergenic_enhancer promoter distal_enhancer proximal_enhancer random_region
do
	multiBigwigSummary BED-file -p 15 \
	--bwfiles ${file_path}/${file_name_1}_ATAC.bw ${file_path}/${file_name_2}_ATAC.bw ${file_path}/${file_name_3}_ATAC.bw\
	--BED ${name}.bed \
	--labels ${file_name_1} ${file_name_2} ${file_name_3} \
	--outFileName ATAC_${name}.npz --outRawCounts ATAC_${name}.tab
done
for name in intergenic_enhancer promoter distal_enhancer proximal_enhancer random_region
do
	multiBigwigSummary BED-file -p 15 \
	--bwfiles  ${file_path}/${file_name_1}-H3K27ac.bw  ${file_path}/${file_name_2}-H3K27ac.bw  ${file_path}/${file_name_3}-H3K27ac.bw  \
	--BED ${name}.bed \
	--labels ${file_name_1} ${file_name_2} ${file_name_3} \
	--outFileName H3K27ac_${name}.npz --outRawCounts H3K27ac_${name}.tab
done
for name in intergenic_enhancer promoter distal_enhancer proximal_enhancer random_region
do
	multiBigwigSummary BED-file -p 15 \
	--bwfiles  ${file_path}/${file_name_1}_RNA.bw ${file_path}/${file_name_2}_RNA.bw ${file_path}/${file_name_3}_RNA.bw \
	--BED ${name}.bed \
	--labels ${file_name_1} ${file_name_2} ${file_name_3} \
	--outFileName RNA_${name}.npz --outRawCounts RNA_${name}.tab
done
rm -rf *.npz


