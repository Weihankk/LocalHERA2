#!/usr/bin/Rscript
library(parallel)

# How many daligner run at the same time ?
cl <- makeCluster(30)

Use_Working_Script <- function(x){
  use.script <- paste0(Working_Script,"/",x)
  return(use.script)
}

# args[1]: LocalHERA2_Parameters.R
#args <- commandArgs(T)
args <- c("/store/whzhang/LocalHERA2/LocalHERA2_Parameters.R")
source(args[1])

print("--------------------------------------------")
print("              Start LocalHERA2              ")
print("--------------------------------------------")


# Make the working dirs
#----------------------------------
system("mkdir 01-Pacbio_And_NonScaffold")
setwd("./01-Pacbio_And_NonScaffold")
system(Use_Working_Script("Check"))
setwd("../")

system("mkdir 02-Pacbio-Alignment")
setwd("./02-Pacbio-Alignment")
system(Use_Working_Script("Check"))
setwd("../")

system("mkdir 03-Pacbio-SelfAlignment")
setwd("./03-Pacbio-SelfAlignment")
system(Use_Working_Script("Check"))
setwd("../")

system("mkdir 04-Graphing")
setwd("./04-Graphing")
system(Use_Working_Script("Check"))
setwd("../")

system("mkdir 05-PathContig")
setwd("./05-PathContig")
system(Use_Working_Script("Check"))
setwd("../")

system("mkdir 06-Daligner")
setwd("./06-Daligner")
system(Use_Working_Script("Check"))
setwd("../")

system("mkdir 07-FilledGap")
setwd("./07-FilledGap")
system(Use_Working_Script("Check"))
setwd("../")

system("mkdir 08-PathContig_Consensus")
setwd("./08-PathContig_Consensus")
system(Use_Working_Script("Check"))
setwd("../")

system("mkdir 09-ReAssembly")
system(Use_Working_Script("Check"))
#----------------------------------

#convert the fasta to lines
#----------------------------------
print("Convert the fasta to lines ...")
system(paste(Use_Working_Script("readstoline"), genome_seq, paste0(genome_name,"-Genome.fasta"),"C"))
#----------------------------------

#split the sequences into two files with large contigs and small contigs
#----------------------------------
print("Split the sequences into two files with large contigs and small contigs ...")
system(paste(Use_Working_Script("01-Filter_Raw_Contig_By_Length"), paste0(genome_name,"-Genome.fasta"), "Large_Contig.fasta", "Small_Contig.fasta", "50000", "15000"))
#----------------------------------

#covert the fasta formate to lines
#----------------------------------
print("Covert the fasta formate to lines ...")
system(paste(Use_Working_Script("readstoline"), Corrected_Pacbio, paste0(genome_name,"-CorrectedPacbio.fasta"), "P"))
Corrected_Pacbio <- paste0(genome_name, "-CorrectedPacbio.fasta")
#----------------------------------

#Merge the non-scaffolded contig with corrected pacbio and they are all used to construct overlaping graph
#----------------------------------
print("Merge the non-scaffolded contig with corrected pacbio and they are all used to construct overlaping graph ..")
system(paste("cat", Bionano_NonScaffolded_Contig, Corrected_Pacbio ,"> Query_Merged.fasta"))
setwd("./01-Pacbio_And_NonScaffold/")
#----------------------------------

#Split the corrected pacbios and non-scaffolded contigs into parts
#----------------------------------
print("Split the corrected pacbios and non-scaffolded contigs into parts ...")
system(paste(Use_Working_Script("03-fasta-splitter"),"--n-parts 100 ../Query_Merged.fasta"))
#----------------------------------

#Make the list of split sequence
#----------------------------------
print("Make the list of split sequence ...")
setwd("../")
system("ls ./01-Pacbio_And_NonScaffold/*.fasta >list_Split.txt")
#----------------------------------

#Make the index of Contig
#----------------------------------
print("Make the index of Contig ...")
system(paste("/store/whzhang/tools/bwa-0.7.10/bwa index",Bionano_Scaffolded_Contig))
#----------------------------------

#Align the corrected pacbios and non-scaffolded contigs to scaffolded contigs
#----------------------------------
print("Align the corrected pacbios and non-scaffolded contigs to scaffolded contigs ...")
list_Split <- read.table("list_Split.txt", header = F)
list_Split$V1 <- as.character(list_Split$V1)

RunBWA_1 <- function(i.n = "./01-Pacbio_And_NonScaffold/Query_Merged.part-001.fasta"){
  i.n.count <- gsub(pattern = "Query_Merged.part-",replacement = "", x = basename(i.n))
  i.n.count <- as.character(as.integer(gsub(pattern = ".fasta", replacement = "", x = i.n.count)))
  system(paste("/store/whzhang/tools/bwa-0.7.10/bwa mem -a -t 10 Large_Contig.fasta", i.n, paste0("> ./02-Pacbio-Alignment/Part_Alignment_",i.n.count,".sam")))
}

clusterExport(cl, "RunBWA_1")
a <- parLapply(cl = cl, X = list_Split$V1, fun = RunBWA_1)

RunSAM2BLASR_1 <- function(i.n = "./01-Pacbio_And_NonScaffold/Query_Merged.part-001.fasta", sam2blasr = Use_Working_Script("sam2blasr.pl")){
  i.n.count <- gsub(pattern = "Query_Merged.part-",replacement = "", x = basename(i.n))
  i.n.count <- as.character(as.integer(gsub(pattern = ".fasta", replacement = "", x = i.n.count)))
  system(paste("perl",sam2blasr,paste0("./02-Pacbio-Alignment/Part_Alignment_",i.n.count,".sam"),paste0("./02-Pacbio-Alignment/Part_Alignment_",i.n.count,".txt")))
}

clusterExport(cl, "RunSAM2BLASR_1")
a <- parLapply(cl = cl, X = list_Split$V1, fun = RunSAM2BLASR_1, sam2blasr = Use_Working_Script("sam2blasr.pl"))
#----------------------------------

#Merge all alignment into an single file
#----------------------------------
print("Merge all alignment into an single file ...")
system("cat ./02-Pacbio-Alignment/Part_Alignment_*.txt > ./02-Pacbio-Alignment/Total_Alignment.txt")
#----------------------------------

#Remove the pacbios and non-scaffolded contigs aligned to the internal scaffolded contigs
#----------------------------------
print("Remove the pacbios and non-scaffolded contigs aligned to the internal scaffolded contigs ...")
system(paste(Use_Working_Script("05-Filtered_InterIncluded_Pacbio"),"./02-Pacbio-Alignment/Total_Alignment.txt ./02-Pacbio-Alignment/InterIncluded_Pacbio.txt", InterIncluded_Identity, InterIncluded_Coverage, InterIncluded_Side))
#----------------------------------

#Record the pacbio alignment of contig's head and end
#----------------------------------
print("Record the pacbio alignment of contig's head and end ...")
system(paste(Use_Working_Script("06-Extract_Contig_Head_Tail_Pacbio_Alignment"), "-Align=./02-Pacbio-Alignment/Total_Alignment.txt", paste0("-MinIden=",MinIdentity), paste0("-MinCov=",MinCoverage), paste0("-HTLen=",InterIncluded_Side), paste0("-MinLen=",MinLength)))
#----------------------------------

#Change the aligned positions into positive chain
#----------------------------------
print("Change the aligned positions into positive chain ...")
system(paste(Use_Working_Script("10-Switch_Locus_To_Positive"),"Contig_Head_Tail_Pacbio.txt ./04-Graphing/Contig_Head_Tail_Pacbio_Pos.txt"))
#----------------------------------

#Extract the sequence of corrected pacbio and non-scaffoled contigs which are nonaligned or aligned to the start or end of the contigs
#----------------------------------
print("Extract the sequence of corrected pacbio and non-scaffoled contigs which are nonaligned or aligned to the start or end of the contigs ...")
system(paste(Use_Working_Script("07-extract_fasta_seq_by_name"),"./02-Pacbio-Alignment/InterIncluded_Pacbio.txt ./Query_Merged.fasta ./02-Pacbio-Alignment/Both_Side_Pacbio.fasta"))
#----------------------------------

#Split the remained pacbio or contigs into parts
#----------------------------------
print("Split the remained pacbio or contigs into parts")
setwd("./03-Pacbio-SelfAlignment/")
system(paste(Use_Working_Script("03-fasta-splitter"),"--n-parts 30 ../02-Pacbio-Alignment/Both_Side_Pacbio.fasta"))
#----------------------------------

#Make index for every part of the pacbios and non-scaffolded contigs
#----------------------------------
setwd("../")
system("ls ./03-Pacbio-SelfAlignment/*.fasta >list_outer_pacbio.txt")
list_outer <- read.table("list_outer_pacbio.txt", header = F)
list_outer$V1 <- as.character(list_outer$V1)

RunBWA_INDEX <- function(i.n = "./01-Pacbio_And_NonScaffold/Query_Merged.part-001.fasta"){
  system(paste("/store/whzhang/tools/bwa-0.7.10/bwa index",i.n))
}

clusterExport(cl, "RunBWA_INDEX")
a <- parLapply(cl = cl, X = list_outer$V1, fun = RunBWA_INDEX)
#----------------------------------

#Align the corrected pacbios and non-scaffolded contigs to each other for finding overlaps
#----------------------------------
print("Align the corrected pacbios and non-scaffolded contigs to each other for finding overlaps ...")
i.v <- c()
j.v <- c()
for (i in seq(length(list_outer$V1))){
  part.i <- list_outer$V1[i]
  j <- i
  for (j in seq(j,length(list_outer$V1))){
    part.j <- list_outer$V1[j]
    #print(c(i,j))
    i.v <- c(i.v, i)
    j.v <- c(j.v, j)
  }
}

bwa.pair.list <- paste(i.v, j.v, sep = ",")

RunBWA_2 <- function(x = "1,1", fasta.list = list_outer$V1){
  x1 <- strsplit(x,",")[[1]][1]
  x2 <- strsplit(x,",")[[1]][2]
  if (x1 == x2){
    system(paste("/store/whzhang/tools/bwa-0.7.10/bwa mem -a -e -t 8",fasta.list[as.integer(x1)], fasta.list[as.integer(x2)], paste0("> ./03-Pacbio-SelfAlignment/Part_SelfAlignment_",as.character(x1),"-",as.character(x2),".sam"),paste0("2> ./03-Pacbio-SelfAlignment/Part_SelfAlignment_",as.character(x1),"-",as.character(x2),".bwalog")))
  }else{
    system(paste("/store/whzhang/tools/bwa-0.7.10/bwa mem -a -t 8",fasta.list[as.integer(x1)], fasta.list[as.integer(x2)], paste0("> ./03-Pacbio-SelfAlignment/Part_SelfAlignment_",as.character(x1),"-",as.character(x2),".sam"),paste0("2> ./03-Pacbio-SelfAlignment/Part_SelfAlignment_",as.character(x1),"-",as.character(x2),".bwalog")))
  }
}

clusterExport(cl, "RunBWA_2")
a <- parLapply(cl = cl, X = bwa.pair.list, fun = RunBWA_2, fasta.list = list_outer$V1)

RunSAM2BLASR_2 <- function(x = "1,1", fasta.list = list_outer$V1, Working_Script = Working_Script){
  Use_Working_Script <- function(x){
    use.script <- paste0(Working_Script,"/",x)
    return(use.script)
  }
  
  x1 <- strsplit(x,",")[[1]][1]
  x2 <- strsplit(x,",")[[1]][2]
  system(paste(Use_Working_Script("sam2blasr.pl"), paste0("./03-Pacbio-SelfAlignment/Part_SelfAlignment_",as.character(x1),"-",as.character(x2),".sam"), paste0("./03-Pacbio-SelfAlignment/Part_SelfAlignment_",as.character(x1),"-",as.character(x2),".txt")))
}

clusterExport(cl, "RunSAM2BLASR_2")
a <- parLapply(cl = cl, X = bwa.pair.list, fun = RunSAM2BLASR_2, fasta.list = list_outer$V1, Working_Script = Working_Script)

FilterAlign <- function(x = "1,1", fasta.list = list_outer$V1, Working_Script = Working_Script, MaxOverhang_Overlap = MaxOverhang_Overlap, MinIdentity_Overlap = MinIdentity_Overlap, MinOverlap_Overlap = MinOverlap_Overlap, MinExtend_Overlap = MinExtend_Overlap ){
  Use_Working_Script <- function(x){
    use.script <- paste0(Working_Script,"/",x)
    return(use.script)
  }
  
  x1 <- strsplit(x,",")[[1]][1]
  x2 <- strsplit(x,",")[[1]][2]
  system(paste(Use_Working_Script("11-PacbioAlignmentFilter"), paste0("./03-Pacbio-SelfAlignment/Part_SelfAlignment_",as.character(x1),"-",as.character(x2),".txt"), MaxOverhang_Overlap, MinIdentity_Overlap, MinOverlap_Overlap, MinExtend_Overlap, paste0("> ./04-Graphing/PacbioAlignmentFiltered_",as.character(x1),"-",as.character(x2),".txt")))
  
}

clusterExport(cl, "FilterAlign")
a <- parLapply(cl = cl, X = bwa.pair.list, fun = FilterAlign, fasta.list = list_outer$V1, Working_Script = Working_Script, MaxOverhang_Overlap = MaxOverhang_Overlap, MinIdentity_Overlap = MinIdentity_Overlap, MinOverlap_Overlap = MinOverlap_Overlap, MinExtend_Overlap = MinExtend_Overlap)
#----------------------------------

#Merge all filtered alignment into an single file
#----------------------------------
setwd("./04-Graphing/")
print("Merge all filtered alignment into an single file ...")
system(paste("cat PacbioAlignmentFiltered_* > PacbioAlignmentFiltered.txt"))
setwd("../")
#----------------------------------

#Find the proper overlap for constructing the graph
#----------------------------------
print("Find the proper overlap for constructing the graph ...")
system(paste(Use_Working_Script("12-PacbioAlignmentLinker"), "./04-Graphing/PacbioAlignmentFiltered.txt", MaxOverhang_Overlap, MinExtend_Overlap, "> ./04-Graphing/PacbioAlignmentLinked.txt"))
#----------------------------------

#Constrct graph by the alignment of pacbios, and the nodes are pacbios and the edges are overlaps.
#Then Finding Contigs Pathway with the Correct Orientatios
#----------------------------------
print("Constrct graph by the alignment of pacbios, and the nodes are pacbios and the edges are overlaps.")
print("Then Finding Contigs Pathway with the Correct Orientatios")
setwd("./04-Graphing/")
system(paste(Use_Working_Script("Selected_Best_Pairs"), "PacbioAlignmentLinked.txt PacbioAlignmentLinked_BestMatch.txt"))
system(paste(Use_Working_Script("13-Graph_By_Finding_Best_MaxExtending_Random_Path"), "PacbioAlignmentLinked_BestMatch.txt >check"))
#----------------------------------

#Output the uniq path
#----------------------------------
print("Output the uniq path ...")
system("cat ctg_clusters.txt |sort |uniq > ../05-PathContig/ctg_clusters_uniq.txt")
system("cat cluster_ori.txt |sort |uniq > ../05-PathContig/cluster_ori_uniq.txt")
setwd("../")
setwd("./05-PathContig/")
#----------------------------------

#Make the corrected pacbios and non-scaffolded contigs into a line
#----------------------------------
print("Make the corrected pacbios and non-scaffolded contigs into a line ...")
system(paste(Use_Working_Script("14-make_ctg_line"), "cluster_ori_uniq.txt cluster_ori_same_chain.txt"))
system(paste(Use_Working_Script("18-compute_fasta_file_len"), "../Query_Merged.fasta Query_Len.txt"))
#----------------------------------

#Change the path into the same chain of bionano scaffolds
#----------------------------------
print("Change the path into the same chain of bionano scaffolds")
system(paste(Use_Working_Script("15-make_junction_by_pos"), "../04-Graphing/ctg_pairs.txt Query_Len.txt cluster_ori_same_chain.txt cluster_ori_same_chain_pos.txt"))
#----------------------------------

#Extract the aligned information of pacbios for final pathcontigs
#----------------------------------
print("Extract the aligned information of pacbios for final pathcontigs")
system(paste(Use_Working_Script("16-extract_ctg_infor_for_seq"), "cluster_ori_same_chain_pos.txt cluster_ori_same_chain_pos_for_seq.txt"))
system("echo '>NA' >NA.fasta")
system("echo 'ATCG' >>NA.fasta")
#----------------------------------

#Output the final contigs of path used to fill the gap of bionano
#----------------------------------
print("Output the final contigs of path used to fill the gap of bionano")
system(paste(Use_Working_Script("17-extract_seq_by_pos"), "cluster_ori_same_chain_pos_for_seq.txt ../Query_Merged.fasta NA.fasta PathContig.fasta"))
#----------------------------------

#Compute the length of pathcontigs
#----------------------------------
print("Compute the length of pathcontigs")
system(paste(Use_Working_Script("18-compute_fasta_file_len"), "PathContig.fasta ../06-Daligner/PathContig_Len.txt"))
setwd("../")
#----------------------------------

#Make the working dirs
#----------------------------------
print("Make the working dirs : 10-Contig_Pairs ...")
system("mkdir 10-Contig_Pairs")
setwd("./10-Contig_Pairs/")
system(Use_Working_Script("Check"))
system("touch overlap.txt")
#----------------------------------

#Formating the contig pairs based on the paths
#----------------------------------
print("Formating the contig pairs based on the paths ...")
system(paste(Use_Working_Script("03-Formate_Contig_Pairs_By_Paths"), "overlap.txt ../05-PathContig/ctg_clusters_uniq.txt Contig_Pairs.txt"))
system(paste0("cat Contig_Pairs.txt |awk '{if(($5+$6/3+$7/6)>='",MinPathNum,"'){$8=$5+$6/3+$7/6;print $0;}}' >Contig_Pairs_Filtered.txt"))
#----------------------------------

#Selecting the final contig pairs with clustering based on scores
#----------------------------------
print("Selecting the final contig pairs with clustering based on scores ...")
system(paste(Use_Working_Script("05-Merge_With_HighestScore_To_Sequence_By_Path"), "Contig_Pairs_Filtered.txt ../Large_Contig.fasta SuperContig.fasta >Selected_Pairs.txt"))
setwd("../")
#----------------------------------

#Extract the paths which connects the final selected contigs
#----------------------------------
setwd("./06-Daligner/")
print("Extract the paths which connects the final selected contigs ...")
system(paste(Use_Working_Script("19-Path2Scaffold_NoBioNano"), "../10-Contig_Pairs/Selected_Pairs.txt ../05-PathContig/ctg_clusters_uniq.txt PathContig_Len.txt Path_Scaffold.txt"))
#----------------------------------

#Rename the path contigs
#----------------------------------
print("Rename the path contigs ...")
system(paste(Use_Working_Script("20-PathContig-Rename_NoBioNano"), "Path_Scaffold.txt ../05-PathContig/PathContig.fasta PathContig_Rename.fasta >log"))
system(paste(Use_Working_Script("Rename1"), "../10-Contig_Pairs/SuperContig.fasta  SuperContig_Rename.fasta >Rename_Pairs.txt"))
system(paste(Use_Working_Script("Rename2"), "Rename_Pairs.txt PathContig_Rename.fasta PathContig_Rename2.fasta"))
system("mv -f PathContig_Rename2.fasta PathContig_Rename.fasta")
#----------------------------------

#Formating the connected scaffold
#----------------------------------
print("Formating the connected scaffold ...")
system(paste(Use_Working_Script("01-Gap_Count"), "SuperContig_Rename.fasta", Enzyme, "Gap.txt"))
system(paste(Use_Working_Script("01-Finding_Contigs_Gap"), "Gap.txt Scaffold2Ctg_Gap.txt"))
system(paste(Use_Working_Script("02-Split_Scaffold_To_Contigs"), "SuperContig_Rename.fasta Prosudo_ScaffoldNonEnzyme2Contig.fasta", Enzyme))
#----------------------------------

#Aligning the path-contigs to scaffold
#----------------------------------
print("Aligning the path-contigs to scaffold ...")
system(paste("perl", Use_Working_Script("21-Daligner_New.pl"), "Scaffold2Ctg_Gap.txt Prosudo_ScaffoldNonEnzyme2Contig.fasta PathContig_Rename.fasta", queue, "qsub", Working_Script, genome_name, DAZZ_DB, DALIGNER))
ss <- dir(pattern = "^Super.*pbs")
ss <- gsub(pattern = ".pbs", replacement = "", x = ss)

RunDalign <- function(Gap_Info = "Super-Scaffold_1-1-2",Working_Script = Working_Script, threads = "8", wd = getwd()){
  Use_Working_Script <- function(x){
    use.script <- paste0(Working_Script,"/",x)
    return(use.script)
  }
  setwd(paste(wd,Gap_Info,sep = "/"))
  print(Gap_Info)
  print(getwd())
  system(paste("perl", Use_Working_Script("lines_to_split.pl"), paste0(Gap_Info,".fasta"), paste0(Gap_Info,"-formated.fasta")))
  system(paste("/store/whzhang/anaconda2/envs/denovo_asm/bin/fasta2DAM", Gap_Info, paste0(Gap_Info,"-formated.fasta")))
  system(paste("/store/whzhang/anaconda2/envs/denovo_asm/bin/DBdust", paste0(Gap_Info,".dam")))
  system(paste("/store/whzhang/anaconda2/envs/denovo_asm/bin/DBsplit", "-x1000 -s50", paste0(Gap_Info,".dam")))
  system(paste("/store/whzhang/anaconda2/envs/denovo_asm/bin/HPC.daligner", paste0(Gap_Info,".dam"),paste0("-T",threads),  ">", paste0(Gap_Info,".sh")))
  system(paste("time sh",paste0(Gap_Info,".sh")))
  
  system(paste("rm -f",paste0(Gap_Info,".*.",Gap_Info,".*.?*.las")))
  system(paste("/store/whzhang/anaconda2/envs/denovo_asm/bin/DBdump -rh", paste0(Gap_Info,".dam"),"| perl",Use_Working_Script("ParseDAZZDB.pl"),">ParseDAZZDB.txt"))
  system(paste("cat", paste0(Gap_Info,"*.las > All.las")))
  system(paste("/store/whzhang/anaconda2/envs/denovo_asm/bin/LAdump -cd", paste0(Gap_Info,".dam All.las | perl"), Use_Working_Script("ParseLA.pl"),">",paste0(Gap_Info,"-Final.txt")))
  system(paste("perl", Use_Working_Script("Daligner_Reformate.pl"), paste0(Gap_Info,"-Final.txt"), paste0(Gap_Info,"-Final_Reformated.txt")))
}

cl <- makeCluster(10)
clusterExport(cl, "RunDalign")
a <- parLapply(cl = cl, X = ss, fun = RunDalign, Working_Script = Working_Script, threads = "8", wd = getwd())
#----------------------------------

#Filling the gaps with the path-contigs and the final sequence of Supercontig is SuperContig.fasta
#----------------------------------
print("Filling the gaps with the path-contigs and the final sequence of Supercontig is SuperContig.fasta ...")
system(paste(Use_Working_Script("22-Filling-Gap"), "Scaffold2Ctg_Gap.txt Prosudo_ScaffoldNonEnzyme2Contig.fasta PathContig_Rename.fasta SuperContig.fasta"))
#----------------------------------

#formating the final genome
#warnings: In general, users need to filter the used small contigs in "Small_Contig.fasta" rather than simply merge files. Users can map the small contigs to SuperContigs by the tools of bwa and filter the small contigs by coverage
#----------------------------------
print("Formating the final genome ...")
system(paste(Use_Working_Script("Final_Formating.sh"), Bionano_NonScaffolded_Contig, genome_name))
#----------------------------------


