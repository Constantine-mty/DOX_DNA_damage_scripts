#bedfile
bedfile1=gene.bed
regionsLabel1="gene"
color="#003049"

##RNAseq
computeMatrix scale-regions -p 10 \
 -S ${file_path}/${file_name_1}_RNA.bw ${file_path}/${file_name_2}_RNA.bw ${file_path}/${file_name_3}_RNA.bw \
 -R ${bedfile1} \
 --binSize 50 \
 --regionBodyLength 4000 \
 -a 1000 -b 1000 \
 -out RNA_expression.tab.gz
plotProfile \
 -m RNA_expression.tab.gz \
 -out RNA_expression.png \
 --colors "#8CABFF" "#4477CE" "#003049" \
 --regionsLabel ${regionsLabel1} \
 --plotHeight 8 --plotWidth 10 \
 --perGroup \
 --samplesLabel "1-RNA" "2-RNA" "3-RNA"\
 --plotTitle 'RNA expression'

##ATACseq
computeMatrix scale-regions -p 10 \
 -S ${file_path}/${file_name_1}_ATAC.bw ${file_path}/${file_name_2}_ATAC.bw ${file_path}/${file_name_3}_ATAC.bw \
 -R ${bedfile1} \
 --binSize 50 \
 --regionBodyLength 4000 \
 -a 1000 -b 1000 \
 -out ATAC.tab.gz
plotProfile \
 -m ATAC.tab.gz \
 -out ATAC.png \
 --colors "#8CABFF" "#4477CE" "#003049" \
 --regionsLabel ${regionsLabel1} \
 --plotHeight 8 --plotWidth 10 \
 --perGroup \
 --samplesLabel "1-ATAC" "2-ATAC" "3-ATAC"\
 --plotTitle 'ATAC signal'

##Chipseq
computeMatrix reference-point -p 10 \
 -S ${file_path}/${file_name_1}-H3K27ac.bw  ${file_path}/${file_name_2}-H3K27ac.bw  ${file_path}/${file_name_3}-H3K27ac.bw \
 -R ${bedfile1} \
 --binSize 50 \
 --regionBodyLength 4000 \
 -a 1000 -b 1000 \
 -out H3K27ac.tab.gz
plotProfile \
 -m H3K27ac.tab.gz \
 -out H3K27ac.png \
 --colors "#8CABFF" "#4477CE" "#003049" \
 --regionsLabel ${regionsLabel1} \
 --plotHeight 8 --plotWidth 10 \
 --perGroup \
 --samplesLabel "1-H3K27ac" "1-H3K27ac" "3-H3K27ac" \
 --plotTitle 'H3K27ac signal'


#plotProfile
for i in promoter intergenic_enhancer distal_enhancer proximal_enhancer random_region
do
	computeMatrix reference-point -p 15 \
	-S ${file_path}/${file_name_1}_RNA.bw ${file_path}/${file_name_2}_RNA.bw ${file_path}/${file_name_3}_RNA.bw \
	-R ../${i}.bed \
	--binSize 50 \
	--referencePoint center \
	-a 500 -b 500 \
	-out ${i}.tab.gz
	cd ./${i}
	bash plotprofile_pipeline.sh
	cd ../
done
for i in promoter intergenic_enhancer distal_enhancer proximal_enhancer random_region
do
	computeMatrix reference-point -p 15 \
	-S ${file_path}/${file_name_1}_ATAC.bw ${file_path}/${file_name_2}_ATAC.bw ${file_path}/${file_name_3}_ATAC.bw \
	-R ./${i}.bed \
	--binSize 50 \
	--referencePoint center \
	-a 1000 -b 1000 \
	-out ${i}.tab.gz
	cd ./${i}
	bash plotprofile_pipeline.sh
	cd ../
done
for i in promoter intergenic_enhancer distal_enhancer proximal_enhancer random_region
do
	computeMatrix reference-point -p 15 \
	-S ${file_path}/${file_name_1}-H3K27ac.bw  ${file_path}/${file_name_2}-H3K27ac.bw  ${file_path}/${file_name_3}-H3K27ac.bw  \
	-R ../${i}.bed \
	--binSize 50 \
	--referencePoint center \
	-a 1000 -b 1000 \
	-out ${i}.tab.gz
	cd ./${i}
	bash plotprofile_pipeline.sh
	cd ../
done


name=distal_enhancer
##RNAseq
plotProfile \
 -m ../${name}.tab.gz \
 -out ${name}.png \
 --colors "#FFCCD5" "#FF4D6D" "#800F2F" \
 --regionsLabel ${name} \
 --plotHeight 8 --plotWidth 10 \
 --perGroup \
 --yMax 0.065 \
 --samplesLabel "1-RNA" "2-RNA" "3-RNA"\
 --plotTitle 'RNA signal'

name=promoter
##RNAseq
plotProfile \
 -m ../${name}.tab.gz \
 -out ${name}.png \
 --colors "#B7EFC5" "#2DC653" "#155D27" \
 --regionsLabel ${name} \
 --plotHeight 8 --plotWidth 10 \
 --perGroup \
 --samplesLabel "1-RNA" "2-RNA" "3-RNA"\
 --plotTitle 'RNA signal'

name=proximal_enhancer
##RNAseq
plotProfile \
 -m ../${name}.tab.gz \
 -out ${name}.png \
 --colors "#FFCCD5" "#FF4D6D" "#800F2F" \
 --regionsLabel ${name} \
 --plotHeight 8 --plotWidth 10 \
 --perGroup \
 --yMax 0.065 \
 --samplesLabel "1-RNA" "2-RNA" "3-RNA"\
 --plotTitle 'RNA signal'

name=random_region
##RNAseq
plotProfile \
 -m ../${name}.tab.gz \
 -out ${name}.png \
 --colors "#eae2b7" "#FFD700" "#B8860B" \
 --regionsLabel ${name} \
 --plotHeight 8 --plotWidth 10 \
 --perGroup \
 --samplesLabel "1-RNA" "2-RNA" "3-RNA"\
 --plotTitle 'RNA signal'

#plotHeatmap
for i in promoter intergenic_enhancer distal_enhancer proximal_enhancer random_region
do
	computeMatrix reference-point -p 15 \
	-S ${file_path}/${file_name_1}_ATAC.bw ${file_path}/${file_name_2}_ATAC.bw ${file_path}/${file_name_3}_ATAC.bw \
	-R ../${i}.bed \
	--binSize 50 \
	--referencePoint center \
	-a 1000 -b 1000 \
	-out ${i}.tab.gz
	mkdir -p ${i}
done
for i in promoter intergenic_enhancer distal_enhancer proximal_enhancer random_region
do
	computeMatrix reference-point -p 15 \
	-S ${file_path}/${file_name_1}-H3K27ac.bw  ${file_path}/${file_name_2}-H3K27ac.bw  ${file_path}/${file_name_3}-H3K27ac.bw \
	-R ../${i}.bed \
	--binSize 50 \
	--referencePoint center \
	-a 1000 -b 1000 \
	-out ${i}.tab.gz
	cd ./${i}
	bash plotHeatmap.sh
	cd ../
done
#ATACseq
name=distal_enhancer
plotHeatmap \
 -m ../${name}.tab.gz \
 -out ${name}.png \
 --sortRegions ascend \
 --colorList "black,#FFACAC,white" \
 --missingDataColor 1 \
 --heatmapHeight 7 \
 --heatmapWidth 2.5 \
 --whatToShow "heatmap and colorbar" \
 --boxAroundHeatmaps no \
 --refPointLabel "ATAC-seq peak center" \
 --samplesLabel "1-ATAC" "2-ATAC" "3-ATAC" \
 --regionsLabel ${name}
name=promoter
plotHeatmap \
 -m ../${name}.tab.gz \
 -out ${name}.png \
 --sortRegions ascend \
 --colorList "black,#CBFFA9,white" \
 --missingDataColor 1 \
 --heatmapHeight 7 \
 --heatmapWidth 2.5 \
 --whatToShow "heatmap and colorbar" \
 --boxAroundHeatmaps no \
 --refPointLabel "ATAC-seq peak center" \
 --samplesLabel "1-ATAC" "2-ATAC" "3-ATAC" \
 --regionsLabel ${name}
name=proximal_enhancer
plotHeatmap \
 -m ../${name}.tab.gz \
 -out ${name}.png \
 --sortRegions ascend \
 --colorList "black,#FFACAC,white" \
 --missingDataColor 1 \
 --heatmapHeight 7 \
 --heatmapWidth 2.5 \
 --whatToShow "heatmap and colorbar" \
 --boxAroundHeatmaps no \
 --refPointLabel "ATAC-seq peak center" \
 --samplesLabel "1-ATAC" "2-ATAC" "3-ATAC" \
 --regionsLabel ${name}
name=random_region
plotHeatmap \
 -m ../${name}.tab.gz \
 -out ${name}.png \
 --sortRegions ascend \
 --colorList "black,#FFFEC4,white" \
 --missingDataColor 1 \
 --heatmapHeight 7 \
 --heatmapWidth 2.5 \
 --whatToShow "heatmap and colorbar" \
 --boxAroundHeatmaps no \
 --refPointLabel "ATAC-seq peak center" \
 --samplesLabel "1-ATAC" "2-ATAC" "3-ATAC" \
 --regionsLabel ${name}