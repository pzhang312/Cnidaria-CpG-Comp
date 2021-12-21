library(dplyr)
library(stringr)
library(ggplot2)
setwd("/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/02_eta")
AuDE = read.table("/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/RNA-seq/DE_genes_renamed.txt")
data1 = read.delim("/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Orthogroups-1.tsv",na.strings="")
colnames(data1) = c("Orthogroup", "Acropora_digitifera", "Alatina_alata", "Aurelia_coerulia", "Calvadosia_cruxmelitensis", "Clytia_hemisphaerica", "Hydra_vulgaris", "Aiptasia_pallida","Kudoa_iwatai", "Morbakka_virulenta", "Nematostella_vectensis", "Physalia_physalis")


#11way orthologs#
present_all = na.omit(data1)
colnames(present_all) = c("Orthogroup", "Acropora_digitifera", "Alatina_alata", "Aurelia_coerulia", "Calvadosia_cruxmelitensis", "Clytia_hemisphaerica", "Hydra_vulgaris", "Aiptasia_pallida","Kudoa_iwatai", "Morbakka_virulenta", "Nematostella_vectensis", "Physalia_physalis")
#single copy orthogroups#
SingleCopyOrtholog = c("OG0009698", "OG0009770", "OG0009783", "OG0009815", "OG0009844", "OG0009852", "OG0009859", "OG0009862") 
write.table(present_all, file = "11_way_orthologs_protein_seqs.txt", quote = F, sep = "\t")
#changed some seq IDs, reload present_all#
present_all = read.table("/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11_way_orthologs_protein_seqs_renamed.txt", sep = "\t", header = T)



# Aurelia # 
aurelia = unlist(strsplit(present_all$Aurelia.Genome_v1.2_Protein_Models_12.28.18, ", "))
AureliaCpG = read.table("Aurelia_eta_whole_seqs2.txt", sep = " ", header = T)
AureliaCpG = filter(AureliaCpG, Length > 300 & !is.na(eta))
AureliaCpG$Ortho = "F"
AureliaCpG[AureliaCpG$ID %in% aurelia,]$Ortho = "T"
test = t.test(CGo_e ~ Ortho, data = AureliaCpG)

alatina = unlist(strsplit(present_all$Alatina_alata, ", "))
AlatinaCpG = read.table("Alatina_alata_eta_whole_seqs2.txt", sep = " ", header = T)
AlatinaCpG = filter(AlatinaCpG, Length > 300 & !is.na(eta))
AlatinaCpG$Ortho = "F"
AlatinaCpG[AlatinaCpG$ID %in% alatina,]$Ortho = "T"
AlatinaCpG[48021,]
which(alatina1 == "c19626_g1_i1")
test = t.test(CGo_e ~ Ortho, data = AlatinaCpG)

### function for t test ### 

OrthoCpG1 = function(input, species, output){
  OrthoGenes = unlist(strsplit(present_all[,which(colnames(present_all) == species)], ", "))
  CpGinput = read.table(input, sep = " ", header = T)
  CpGinput = filter(CpGinput, Length > 300 & !is.na(eta))
  CpGinput$Ortho = "F"
  CpGinput[CpGinput$ID %in% OrthoGenes,]$Ortho = "T"
  test = t.test(CGo_e ~ Ortho, data = CpGinput)
  cat(paste(species, "\t", test$estimate[1],"\t", test$estimate[2],"\t", test$statistic,"\t",test$p.value, "\n"), file=output, append = TRUE)
}


cat(paste("species", "\t", "mean_CpG_of_non_Ortho", "\t", "mean_CpG_of_Ortho","\t", "t_stat", "\t", "p-value", "\n"), file='/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt', append = TRUE)

OrthoCpG1(input = "Aurelia_eta_whole_seqs2.txt", species = "Aurelia_coerulia", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")
OrthoCpG1(input = "MOR05_r06_mRNA_whole_seqs_eta.txt", species = "Morbakka_virulenta", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")
OrthoCpG1(input = "Alatina_alata_eta_whole_seqs2.txt", species = "Alatina_alata", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")
OrthoCpG1(input = "Calvadosia_cruxmelitensis_eta_whole_seqs2.txt", species = "Calvadosia_cruxmelitensis", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")
OrthoCpG1(input = "Clytia_hemisphaerica_HAMU01_eta_whole_seqs2.txt", species = "Clytia_hemisphaerica", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")
OrthoCpG1(input = "Kudoa_iwatai_eta_whole_seqs2.txt", species = "Kudoa_iwatai", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")
OrthoCpG1(input = "Physalia_physalis_eta_whole_seqs2.txt", species = "Physalia_physalis", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")


### function for t test for the NCBI protein sources ###
AcroporaOrtho = unlist(strsplit(present_all[,"Acropora_digitifera"], ", "))
write.table(AcroporaOrtho, file = "Acropora_digitifera_11_way_orthologs_protein_seqs.txt", quote = F, sep = "\n", row.names = F, col.names = F)
HydraOrtho = unlist(strsplit(present_all[,"Hydra_vulgaris"], ", "))
write.table(HydraOrtho, file = "Hydra_vulgaris_11_way_orthologs_protein_seqs.txt", quote = F, sep = "\n", row.names = F, col.names = F)
AiptasiaOrtho = unlist(strsplit(present_all[,"Aiptasia_pallida"], ", "))
write.table(AiptasiaOrtho, file = "Aiptasia_pallida_11_way_orthologs_protein_seqs.txt", quote = F, sep = "\n", row.names = F, col.names = F)
NematostellaOrtho = unlist(strsplit(present_all[,"Nematostella_vectensis"], ", "))
write.table(NematostellaOrtho, file = "Nematostella_vectensis_11_way_orthologs_protein_seqs.txt", quote = F, sep = "\n", row.names = F, col.names = F)

OrthoCpG2 = function(Ortholist, CpGdata, species, output){
  Ortho = read.table(Ortholist, sep = "\t", header = T)
  CpGinput = read.table(CpGdata, sep = " ", header = T)
  CpGinput = filter(CpGinput, Length > 300 & !is.na(eta))
  CpGinput$Ortho = "F"
  CpGinput$ID = str_replace(CpGinput$ID,"\\..*","")
  CpGinput[CpGinput$ID %in% Ortho$RefSeq.mRNA.Accession,]$Ortho = "T"
  test = t.test(CGo_e ~ Ortho, data = CpGinput)
  cat(paste(species, "\t", test$estimate[1],"\t", test$estimate[2],"\t", test$statistic,"\t",test$p.value, "\n"), file=output, append = TRUE)
}
OrthoCpG2(Ortholist = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Acropora_giditifera_protein2mRNA_accession_mapping_11_way_ortholog.txt", CpGdata = "Acropora_digitefera_GCF_000222465_rna_eta_whole_seqs.txt", species = "Acropora_digitifera", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")
OrthoCpG2(Ortholist = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Aiptasia_pallida_protein2mRNA_accession_mapping_11way_ortholog.txt", CpGdata = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/03_plots by group from 02_eta/extras/Aiptasia_pallida_GCF_001417965_rna_eta_whole_seqs.txt", species = "Aiptasia_pallida", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")
OrthoCpG2(Ortholist = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Hvu_protein2mRNA_accession_mapping_11way_ortholog.txt", CpGdata = "Hydra_vulgaris_GCF_000004095_rna_eta_whole_seqs.txt", species = "Hydra_vulgaris", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")
OrthoCpG2(Ortholist = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Nve_protein2mRNA_accession_mapping_11way_ortholog.txt", CpGdata = "Nematostella_vectensis_GCF_000209225.1_ASM20922v1_rna_eta_whole_seqs2.txt", species = "Nematostella_vectensis", output = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/11-way-orthologs-stats.txt")

### Pair-wise single ortholog CpG Comparison - Aurelia, Nematostella, Hydra, Kudoa###
Aco_single = data1[str_count(data1$Aurelia_coerulia, pattern = "\\.") == 1, c(1,4)]
Aco_single = na.omit(Aco_single)

Nve_single = data1[str_count(data1$Nematostella_vectensis, pattern = "\\.") == 1, c(1,11)]
Nve_single = na.omit(Nve_single)
write.table(Nve_single$Nematostella_vectensis, file = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Nematostella_vectensis_single_orthologs_protein_seqs.txt", quote = F, sep = "\n", row.names = F, col.names = F)

Hvu_single = data1[str_count(data1$Hydra_vulgaris, pattern = "\\.") == 1, c(1,7)]
Hvu_single = na.omit(Hvu_single)
write.table(Hvu_single$Hydra_vulgaris, file = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Hydra_vulgaris_single_orthologs_protein_seqs.txt", quote = F, sep = "\n", row.names = F, col.names = F)
Kiwa_single = data1[str_count(data1$Kudoa_iwatai, pattern = "\\.") == 1, c(1,9)]
Kiwa_single = na.omit(Kiwa_single)
Kiwa_single$Kudoa_iwatai = str_replace(Kiwa_single$Kudoa_iwatai,"\\..*","")
Adi_single = data1[str_count(data1$Acropora_digitifera, pattern = "\\.") == 1, c(1,2)]
Adi_single = na.omit(Adi_single)
Adi_single$Acropora_digitifera = str_replace(Adi_single$Acropora_digitifera,"\\..*","")
write.table(Adi_single$Acropora_digitifera, file = "/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Acropora_digitifera_single_orthologs_protein_seqs.txt", quote = F, sep = "\n", row.names = F, col.names = F)
Clytia_single = data1[str_count(data1$Clytia_hemisphaerica, pattern = "\\.") == 2, c(1,6)]
Clytia_single = na.omit(Clytia_single)
Clytia_single$Clytia_hemisphaerica = str_replace(Clytia_single$Clytia_hemisphaerica,".p1","")



Aco_CpG = read.table("Aurelia_eta_whole_seqs2.txt", sep = " ", header = T)
Aco_CpG = Aco_CpG[,c(1,2,7)]
Aco_CpG = filter(Aco_CpG, Length > 300 & !is.na(CGo_e))
Kudoa_CpG = read.table("Kudoa_iwatai_eta_whole_seqs2.txt", sep = " ", header = T)
Kudoa_CpG = Kudoa_CpG[,c(1,2,7)]
Kudoa_CpG = filter(Kudoa_CpG, Length > 300 & !is.na(CGo_e))
Hvu_CpG = read.table("Hydra_vulgaris_GCF_000004095_rna_eta_whole_seqs.txt", sep = " ", header = T)
Hvu_CpG = Hvu_CpG[,c(1,2,7)]
Hvu_CpG = filter(Hvu_CpG, Length > 300 & !is.na(CGo_e))
Nve_CpG = read.table("Nematostella_vectensis_GCF_000209225.1_ASM20922v1_rna_eta_whole_seqs2.txt", sep = " ", header = T)
Nve_CpG = Nve_CpG[,c(1,2,7)]
Nve_CpG = filter(Nve_CpG, Length > 300 & !is.na(CGo_e))


#Aurelia vs. Kudoa#
Aco_Kudoa = merge(Aco_single, Kiwa_single, by = "Orthogroup", all = F )
Aco_Kudoa = merge(Aco_Kudoa, Aco_CpG, by.x = "Aurelia_coerulia", by.y = "ID", all = F)
Aco_Kudoa = merge(Aco_Kudoa, Kudoa_CpG, by.x = "Kudoa_iwatai", by.y = "ID", all = F)
colnames(Aco_Kudoa)[5] = "Aurelia_CpG"
colnames(Aco_Kudoa)[7] = "Kudoa_CpG"
ggplot(Aco_Kudoa,aes(x = Aurelia_CpG, y = Kudoa_CpG))+ geom_point()
Aco_Kudoa %>% summarise(low_low = sum(Aurelia_CpG>MeanCpG[1,2] & Kudoa_CpG>MeanCpG[4,2]), low_high = sum(Aurelia_CpG>MeanCpG[1,2] & Kudoa_CpG<MeanCpG[4,2]), high_low = sum(Aurelia_CpG<MeanCpG[1,2] & Kudoa_CpG>MeanCpG[4,2]), high_high = sum(Aurelia_CpG<MeanCpG[1,2] & Kudoa_CpG<MeanCpG[4,2])) 
chisq.test(Aco_Kudoa %>% summarise(low_low = sum(Aurelia_CpG>MeanCpG[1,2] & Kudoa_CpG>MeanCpG[4,2]), low_high = sum(Aurelia_CpG>MeanCpG[1,2] & Kudoa_CpG<MeanCpG[4,2]), high_low = sum(Aurelia_CpG<MeanCpG[1,2] & Kudoa_CpG>MeanCpG[4,2]), high_high = sum(Aurelia_CpG<MeanCpG[1,2] & Kudoa_CpG<MeanCpG[4,2])))


#Aurelia vs. Nematostella #
Aco_Nve = merge(Aco_single, Nve_single, by = "Orthogroup", all = F)
Nve_single_mRNA = read.table("/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Nve_single_orthologs_protein2mRNA_accession_mapping.txt", sep = "\t", header = T)
Nve_CpG$ID = str_replace(Nve_CpG$ID,"\\..*","")
Aco_Nve$Nematostella_vectensis = str_replace(Aco_Nve$Nematostella_vectensis,"\\..*","")

Nve_single_mRNA = merge(Nve_single_mRNA[,c(1,2)], Nve_CpG[,c(1,3)], by.x = 'RefSeq.mRNA.Accession', by.y = 'ID', all = F)
Aco_Nve = merge(Aco_Nve, Nve_single_mRNA[,c(2,3)], by.x = 'Nematostella_vectensis', by.y = 'RefSeq.Protein.Accession', all = F)
Aco_Nve = merge(Aco_Nve, Aco_CpG[,c(1,3)], by.x = 'Aurelia_coerulia', by.y = 'ID', all = F)
colnames(Aco_Nve)[4] = 'Nve_CpG'
colnames(Aco_Nve)[5] = 'Aco_CpG'
ggplot(Aco_Nve,aes(x = Aco_CpG, y = Nve_CpG))+ geom_point()
Aco_Nve %>% summarise(low_low = sum(Aco_CpG>MeanCpG[1,2] & Nve_CpG>MeanCpG[3,2]), low_high = sum(Aco_CpG>MeanCpG[1,2] & Nve_CpG<MeanCpG[3,2]), high_low = sum(Aco_CpG<MeanCpG[1,2] & Nve_CpG>MeanCpG[3,2]), high_high = sum(Aco_CpG<MeanCpG[1,2] & Nve_CpG<MeanCpG[3,2])) 


#Aurelia vs. Hydra # 
Aco_Hvu = merge(Aco_single, Hvu_single, by = "Orthogroup", all = F)
Hvu_single_mRNA = read.table("/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/04_Ortholog_analyses/Hvu_single_orthologs_protein2mRNA_accession_mapping.txt", sep = "\t", header = T)
Hvu_CpG$ID = str_replace(Hvu_CpG$ID,"\\..*","")
Aco_Hvu$Hydra_vulgaris = str_replace(Aco_Hvu$Hydra_vulgaris,"\\..*","")
Hvu_single_mRNA = merge(Hvu_single_mRNA[,c(1,2)], Hvu_CpG[,c(1,3)], by.x = 'RefSeq.mRNA.Accession', by.y = 'ID', all = F)
Aco_Hvu = merge(Aco_Hvu, Hvu_single_mRNA[,c(2,3)], by.x = 'Hydra_vulgaris', by.y = 'RefSeq.Protein.Accession', all = F)
Aco_Hvu = merge(Aco_Hvu, Aco_CpG[,c(1,3)], by.x = 'Aurelia_coerulia', by.y = 'ID', all = F)
colnames(Aco_Hvu)[4] = 'Hvu_CpG'
colnames(Aco_Hvu)[5] = 'Aco_CpG'
ggplot(Aco_Hvu,aes(x = Aco_CpG, y = Hvu_CpG))+ geom_point()
Aco_Hvu %>% summarise(low_low = sum(Aco_CpG>MeanCpG[1,2] & Hvu_CpG>MeanCpG[2,2]), low_high = sum(Aco_CpG>MeanCpG[1,2] & Hvu_CpG<MeanCpG[2,2]), high_low = sum(Aco_CpG<MeanCpG[1,2] & Hvu_CpG>MeanCpG[2,2]), high_high = sum(Aco_CpG<MeanCpG[1,2] & Hvu_CpG<MeanCpG[2,2])) 

#Nematostella vs. Hydra#
Nve_Hvu = merge(Nve_single, Hvu_single, by = "Orthogroup", all = F)
Nve_Hvu$Hydra_vulgaris = str_replace(Nve_Hvu$Hydra_vulgaris,"\\..*","")
Nve_Hvu$Nematostella_vectensis = str_replace(Nve_Hvu$Nematostella_vectensis,"\\..*","")
Nve_Hvu = merge(Nve_Hvu, Hvu_single_mRNA[,c(2,3)], by.x = 'Hydra_vulgaris', by.y = 'RefSeq.Protein.Accession', all = F)
Nve_Hvu = merge(Nve_Hvu, Nve_single_mRNA[,c(2,3)], by.x = 'Nematostella_vectensis', by.y = 'RefSeq.Protein.Accession', all = F)
colnames(Nve_Hvu)[4] = 'Hvu_CpG'
colnames(Nve_Hvu)[5] = 'Nve_CpG'
ggplot(Nve_Hvu,aes(x = Nve_CpG, y = Hvu_CpG))+ geom_point()
Nve_Hvu %>% summarise(low_low = sum(Nve_CpG>MeanCpG[3,2] & Hvu_CpG>MeanCpG[2,2]), low_high = sum(Nve_CpG>MeanCpG[3,2] & Hvu_CpG<MeanCpG[2,2]), high_low = sum(Nve_CpG<MeanCpG[3,2] & Hvu_CpG>MeanCpG[2,2]), high_high = sum(Nve_CpG<MeanCpG[3,2] & Hvu_CpG<MeanCpG[2,2])) 

#Nematostella vs. Kudoa#
Nve_Kiwa = merge(Nve_single, Kiwa_single, by = "Orthogroup", all = F)
Nve_Kiwa$Kudoa_iwatai = str_replace(Nve_Kiwa$Kudoa_iwatai,"\\..*","")
Nve_Kiwa$Nematostella_vectensis = str_replace(Nve_Kiwa$Nematostella_vectensis,"\\..*","")
Nve_Kiwa = merge(Nve_Kiwa, Kudoa_CpG[,c(1,3)], by.x = 'Kudoa_iwatai', by.y = 'ID', all = F)
Nve_Kiwa = merge(Nve_Kiwa, Nve_single_mRNA[,c(2,3)], by.x = 'Nematostella_vectensis', by.y = 'RefSeq.Protein.Accession', all = F)
colnames(Nve_Kiwa)[4] = 'Kiwa_CpG'
colnames(Nve_Kiwa)[5] = 'Nve_CpG'
ggplot(Nve_Kiwa,aes(x = Nve_CpG, y = Kiwa_CpG))+ geom_point()
Nve_Kiwa %>% summarise(low_low = sum(Nve_CpG>MeanCpG[3,2] & Kiwa_CpG>MeanCpG[4,2]), low_high = sum(Nve_CpG>MeanCpG[3,2] & Kiwa_CpG<MeanCpG[4,2]), high_low = sum(Nve_CpG<MeanCpG[3,2] & Kiwa_CpG>MeanCpG[4,2]), high_high = sum(Nve_CpG<MeanCpG[3,2] & Kiwa_CpG<MeanCpG[4,2])) 

#Hydra vs. Kudoa#
Hvu_Kiwa = merge(Hvu_single, Kiwa_single, by = "Orthogroup", all = F)
Hvu_Kiwa$Kudoa_iwatai = str_replace(Hvu_Kiwa$Kudoa_iwatai,"\\..*","")
Hvu_Kiwa$Hydra_vulgaris = str_replace(Hvu_Kiwa$Hydra_vulgaris,"\\..*","")
Hvu_Kiwa = merge(Hvu_Kiwa, Kudoa_CpG[,c(1,3)], by.x = 'Kudoa_iwatai', by.y = 'ID', all = F)
Hvu_Kiwa = merge(Hvu_Kiwa, Hvu_single_mRNA[,c(2,3)], by.x = 'Hydra_vulgaris', by.y = 'RefSeq.Protein.Accession', all = F)
colnames(Hvu_Kiwa)[4] = 'Kiwa_CpG'
colnames(Hvu_Kiwa)[5] = 'Hvu_CpG'
ggplot(Hvu_Kiwa,aes(x = Hvu_CpG, y = Kiwa_CpG))+ geom_point()
Hvu_Kiwa %>% summarise(low_low = sum(Hvu_CpG>MeanCpG[2,2] & Kiwa_CpG>MeanCpG[4,2]), low_high = sum(Hvu_CpG>MeanCpG[2,2] & Kiwa_CpG<MeanCpG[4,2]), high_low = sum(Hvu_CpG<MeanCpG[2,2] & Kiwa_CpG>MeanCpG[4,2]), high_high = sum(Hvu_CpG<MeanCpG[2,2] & Kiwa_CpG<MeanCpG[4,2])) 








### Processing for GO enrichment analyses ### 
Hvu_Kiwa_low = Hvu_Kiwa[Hvu_Kiwa$Kiwa_CpG > MeanCpG[4,2] & Hvu_Kiwa$Hvu_CpG>MeanCpG[2,2],]
 

