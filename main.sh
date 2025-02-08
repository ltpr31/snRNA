#!/bin/bash
workdir=/snRNA/
scripts=/snRNA/scripts
datadir=/snRNA/matrix

cd $workdir
mkdir 01.data_QC/
cd 01.data_QC/

########################################## quality_control #########################################

Rscript $scripts/quality_control1.r \
-d $datadir/ck0h1/ \
-p ck0h1 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value ck0h1 mock 0h ck0

Rscript $scripts/quality_control1.r \
-d $datadir/ck0h2/ \
-p ck0h2 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value ck0h2 mock 0h ck0

Rscript $scripts/quality_control1.r \
-d $datadir/ck12h1/ \
-p ck12h1 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value ck12h1 ck 12h ck12

Rscript $scripts/quality_control1.r \
-d $datadir/ck12h2/ \
-p ck12h2 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value ck12h2 ck 12h ck12


Rscript $scripts/quality_control1.r \
-d $datadir/ck24h1/ \
-p ck24h1 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value ck24h1 ck 24h ck24

Rscript $scripts/quality_control1.r \
-d $datadir/ck24h2/ \
-p ck24h2 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value ck24h2 ck 24h ck24


Rscript $scripts/quality_control1.r \
-d $datadir/ck48h1/ \
-p ck48h1 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value ck48h1 ck 48h ck48

Rscript $scripts/quality_control1.r \
-d $datadir/ck48h2/ \
-p ck48h2 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value ck48h2 ck 48h ck48



Rscript $scripts/quality_control1.r \
-d $datadir/tr0h1/ \
-p tr0h1 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value tr0h1 mock 0h tr0

Rscript $scripts/quality_control1.r \
-d $datadir/tr0h2/ \
-p tr0h2 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value tr0h2 mock 0h tr0

Rscript $scripts/quality_control1.r \
-d $datadir/tr12h1/ \
-p tr12h1 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value tr12h1 tr 12h tr12

Rscript $scripts/quality_control1.r \
-d $datadir/tr12h2/ \
-p tr12h2 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value tr12h2 tr 12h tr12


Rscript $scripts/quality_control1.r \
-d $datadir/tr24h1/ \
-p tr24h1 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value tr24h1 tr 24h tr24

Rscript $scripts/quality_control1.r \
-d $datadir/tr24h2/ \
-p tr24h2 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value tr24h2 tr 24h tr24


Rscript $scripts/quality_control1.r \
-d $datadir/tr48h1/ \
-p tr48h1 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value tr48h1 tr 48h tr48

Rscript $scripts/quality_control1.r \
-d $datadir/tr48h2/ \
-p tr48h2 \
--nUMI.min 500 \
--nGene.min 250 \
--log10GenesPerUMI 0.8 \
--metadata.col.name sname stim time stti \
--metadata.value tr48h2 tr 48h tr48



######### doublet removing ########################


for i in  ck0h1 ck0h2 tr0h1 tr0h2 ck12h1 ck12h2 tr12h1 tr12h2 ck24h1 ck24h2 tr24h1 tr24h2 ck48h1 ck48h2 tr48h1 tr48h2; do
Rscript $scripts/single_sample_clustering.r \
--rds $workdir/01.data_QC/${i}.afterQC.rds \
-p ${i} \
--resolution 0.5 \
--high.variable.genes 2000

Rscript $scripts/DoubletFinder.r \
 -i ${i}.rds \
 -p ${i} 
 -o $workdir/01.data_QC

Rscript $scripts/quality_control2.r \
-i $workdir/01.data_QC/${i}.doubletFinder.rds \
-p ${i} \
--nUMI.min 500 \
--nGene.min 700
done


###################################### samples merge #####################################

cd $workdir
mkdir 02.merge
cd 02.merge
Rscript $scripts/merge_rds.r \
 -i $workdir/01.data_QC/ck0h1.afterQC.rds \
    $workdir/01.data_QC/ck0h2.afterQC.rds \
    $workdir/01.data_QC/tr0h1.afterQC.rds \
    $workdir/01.data_QC/tr0h2.afterQC.rds \
    $workdir/01.data_QC/ck12h1.afterQC.rds \
    $workdir/01.data_QC/ck12h2.afterQC.rds \
    $workdir/01.data_QC/tr12h1.afterQC.rds \
    $workdir/01.data_QC/tr12h2.afterQC.rds \
    $workdir/01.data_QC/ck24h1.afterQC.rds \
    $workdir/01.data_QC/ck24h2.afterQC.rds \
    $workdir/01.data_QC/tr24h1.afterQC.rds \
    $workdir/01.data_QC/tr24h2.afterQC.rds \
    $workdir/01.data_QC/ck48h1.afterQC.rds \
    $workdir/01.data_QC/ck48h2.afterQC.rds \
    $workdir/01.data_QC/tr48h1.afterQC.rds \
    $workdir/01.data_QC/tr48h2.afterQC.rds \
    -p all.sample.merged


#################################### integration ######################################

Rscript $scripts/integrate.r \
--rds $workdir/02.merge/all.sample.merged.rds  \
--batch.id sname \
--integration.reduction rpca \
-p merge \
-o ./

############################# dimensionality reduction #################################

Rscript $scripts/dimensionality_reduction.r \
--rds $workdir/02.merge/merge.integrated.rds  \
--resolution 0.5 -d 33  --high.variable.genes 2000 \
-p reduced \
-o ./

######################################## data scale ####################################

Rscript $scripts/scale_data.r \
-r $workdir/02.merge/reduced.rds \
-o ./   \
-p scaled


#################################  trajectory construction  #############################

cd $workdir
mkdir 03.trajectory
cd 03.trajectory

#### data subset #################################

Rscript $scripts/subset_rds.r \
--rds $workdir/02.merge/all.sample.merged.rds \
--subset 'time %in% c("0h","24h","48h")' \
-p no_12h \
-o $workdir/03.trajectory

Rscript $scripts/subset_rds.r \
--rds $workdir/03.trajectory/no_12h.rds \
--subset 'stim %in% c("mock","ck")' \
-p ck \
-o $workdir/03.trajectory

Rscript $scripts/subset_rds.r \
--rds $workdir/03.trajectory/no_12h.rds \
--subset 'stim %in% c("mock","tr")' \
-p tr \
-o $workdir/03.trajectory


Rscript $scripts/subset_rds.r \
--rds $workdir/03.trajectory/ck.rds \
--subset 'Celltype %in% c("Epidermis")' \
-p Epidermis.ck \
-o $workdir/03.trajectory

Rscript $scripts/subset_rds.r \
--rds $workdir/03.trajectory/tr.rds \
--subset 'Celltype %in% c("Epidermis")' \
-p Epidermis.tr \
-o $workdir/03.trajectory


Rscript $scripts/subset_rds.r \
--rds $workdir/03.trajectory/ck.rds \
--subset 'Celltype %in% c("Procambium")' \
-p Procambium.ck \
-o $workdir/03.trajectory

Rscript $scripts/subset_rds.r \
--rds $workdir/03.trajectory/tr.rds \
--subset 'Celltype %in% c("Procambium")' \
-p Procambium.tr \
-o $workdir/03.trajectory


for i in  Epidermis Procambium ; do
#### trajectory constructed by all genes ###########
Rscript $scripts/trajectory_construction.r \
  -i $workdir/03.trajectory/${i}.ck.rds \
  -p ${i}.ck \
  --gene.select.method monocle \
  -o $workdir/03.trajectory

Rscript $scripts/trajectory_construction.r \
  -i $workdir/03.trajectory/${i}.tr.rds \
  -p ${i}.tr \
  --gene.select.method monocle \
  -o $workdir/03.trajectory


#### gene trim #####################################


Rscript $scripts/gene_pick.r  \
-i $workdir/03.trajectory/${i}.tr.Time_diff_all.tsv  \
-m $workdir/03.trajectory/${i}.ck.Time_diff_all.tsv  \
-o $workdir/03.trajectory/${i}_trimmed.txt  \
--ntr 2000 --nck 500


#### infection-related trajectory construction######

Rscript $scripts/trajectory_construction.r \
  -i $workdir/03.trajectory/${i}.ck.rds \
  -p ${i}.ck.trimmed \
  --gene.select.method Customization \
  -o $workdir/03.trajectory  \
  --gene.list $workdir/03.trajectory/${i}_trimmed.txt

Rscript $scripts/trajectory_construction.r \
  -i $workdir/03.trajectory/${i}.tr.rds \
  -p ${i}.tr.trimmed \
  --gene.select.method Customization \
  -o $workdir/03.trajectory  \
  --gene.list $workdir/03.trajectory/${i}_trimmed.txt
done













