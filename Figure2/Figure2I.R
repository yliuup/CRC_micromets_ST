1.	Bam & Raw depth (50Kb mean)
file/
  */*_rmdup.bam 
*/depth_binned.bed

2.	Copy Number
Inferred CN: 
  *_allelicCN/CN/as_cna_profile_* _rmdup*.txt
<total_copy_number> for CN plot
<nA> <nB> for tree construction

Inferred ploidy and purity:
  *_allelicCN/CN/summary_ascat_sc_*.txt

Code for purity0.8:
  library(ASCAT.sc)
run_sc_sequencing(tumour_bams=bams,
                  allchr=paste0("chr",c(1:22)),
                  sex=sex_list,
                  chrstring_bam="chr",
                  purs = rep(list(seq(0.8, 1, 0.01)), length(bams)),
                  ploidies = rep(list(seq(1.7,5, 0.01)), length(bams)),
                  maxtumourpsi=5,
                  binsize=500000,
                  build="hg38",
                  MC.CORES=8,
                  outdir="12.ascat.sc.multipcf.purity/",
                  projectname=paste0("ascat_sc_", donor),
                  path_to_phases=path_to_phases,
                  list_ac_counts_paths=list_ac_counts_paths,
                  multipcf=TRUE, predict_refit = FALSE)

Code for refit:
  library(ASCAT.sc)
run_sc_sequencing(tumour_bams=bams,
                  allchr=paste0("chr",c(1:22)),
                  sex=sex_list,
                  chrstring_bam="chr",
                  purs = rep(list(seq(0.1, 1, 0.01)), length(bams)),
                  ploidies = rep(list(seq(1.7,5, 0.01)), length(bams)),
                  maxtumourpsi=5,
                  binsize=500000,
                  build="hg38",
                  MC.CORES=8,
                  outdir="12.ascat.sc.multipcf/",
                  projectname=paste0("ascat_sc_", donor),
                  path_to_phases=path_to_phases,
                  list_ac_counts_paths=list_ac_counts_paths,
                  multipcf=TRUE)

3.	Tree
Inferred tree:
  *_allelicCN/Tree/*_as_cna_profile_*.output_allelicCN/*allelicCN_final_tree.new

Input:
  *_allelicCN/Tree/* _as_cna_profile_*.allelicCN.txt

Code: 
  medicc2 \
${mode}_as_cna_profile_${donor}.allelicCN.txt \
${mode}_as_cna_profile_${donor}.output_allelicCN
