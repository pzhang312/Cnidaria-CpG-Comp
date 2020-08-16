library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggpubr)
data1 = read.table("/Users/zhangxinhui/Box Sync/Aurelia_Genetics_Epigenetics/Methylation/CGcontent/Inverts_CpG/CpG_data/outputs/cnidaria/01_stats/cnidarian_sum1_with_new_eta_script.txt", header = T, sep = "\t", quote = "")
data1$Class = factor(data1$Class, levels = c("Anthozoa","Polypodiozoa", "Myxozoa", "Staurozoa", "Cubozoa", "Scyphozoa", "Hydrozoa"))
data1$Order = factor(data1$Order, levels = c("Alcyonacea","Helioporacea", "Actiniaria","Corallimorpharia","Scleractinia","Polypodiidea", "Multivalvulida", "Bivalvulida", "Stauromedusae", "Carybdeida", "Chirodropoda", "Coronatae", "Semaeostomeae", "Rhizostomeae", "Limnomedusae", "Narcomedusae", "Siphonophorae", "Leptothecata", "Anthoathecata"), ordered = TRUE)

bp1 = ggplot(data1, aes(x=Order, y=Mean_GB_CpG, color=Class)) + geom_boxplot(outlier.shape = NA)
bp1 + theme(axis.text.x = element_text(angle = 90), axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 12,face="bold")) + xlab("Order") + ylab("Average Gene Body CpG o/e") + geom_dotplot(binaxis='y', stackdir='center', dotsize=45, binwidth = 0.0002, aes(fill=factor(genes))) + labs(color = "Class", fill = "genes") + scale_fill_manual(values = c("red","gray", "green", "blue"))

#fig 1#
bp2 = ggplot(data1, aes(x=Order, y=Mean_GB_CpG))
fig1 = bp2 + theme(axis.text.x = element_text(angle = 90), axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 14,face="bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 12), legend.position = c(0.21,0.8)) + xlab("Order") + ylab("Average Gene Body CpG o/e") + geom_dotplot(binaxis='y', stackdir='center', dotsize=100, binwidth = 0.0001, aes(fill=factor(genes))) + labs(fill = "Relevant genes", element_text(size = 10)) + scale_fill_manual(values = c("#D55E00","#000000", "#56B4E9", "#F0E442"))
tiff("fig1 v1.tiff",width = 10, height = 8, units = 'in', res = 300)
fig1
dev.off()

#no boxes, with species labels: 
bp2 = ggplot(data1, aes(x=Order, y=Mean_GB_CpG))
bp2 + theme(axis.text.x = element_text(angle = 90), axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 12,face="bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 11.5)) + xlab("Order") + ylab("Average Gene Body CpG o/e") + geom_dotplot(binaxis='y', stackdir='center', dotsize=60, binwidth = 0.0001, aes(fill=factor(genes))) + labs(fill = "Relevant enes", element_text(size = 10)) + scale_fill_manual(values = c("red","dark gray", "blue", "green")) + geom_text(aes(label = Species), size = 1.5, vjust = 0, hjust = 0, color = "black") + scale_x_discrete(expand = c(0,2))


### plot TpG~CpG slope
bp3 = ggplot(data1, aes(x=Order, y=lm.TpG.CpG._Slope, color=Class)) +
  geom_boxplot()
bp3 + xlab("Order") + ylab("TpG~CpG slope") + labs(fill = "Class") + theme(axis.text.x = element_text(angle = 90), axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 12,face="bold"))  + geom_dotplot(binaxis='y', stackdir='center', dotsize=25, binwidth = 0.0005, aes(fill=factor(genes))) + labs(color = "Class", fill = "genes") + scale_fill_manual(values = c("red","gray", "green", "blue"))

### plots for repeats ### 
cbPalette <- c("#E69F00", "#D55E00","#56B4E9", "#009E73", "#CC79A7", "#0072B2","#F0E442")


dp1 = ggplot(data1, aes(x=genome_size, y=mean_repeats_CpG, color = Class)) + geom_point(size = 3)
dp1 + scale_colour_manual(values=cbPalette) + scale_x_continuous(trans='log10', limits = c(30,3000)) + theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 14,face="bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 12)) + xlab("Genome Size (Mbp)") + ylab("Average Repeat CpG o/e") + geom_text(aes(label = abbriviation), size =3.5, vjust = 0, hjust = 0, color = "black")



dp2 = ggplot(data1, aes(x=genome_size, y=Mean_GB_CpG, color = Class)) + geom_point(size = 3)
dp2 + scale_colour_manual(values=cbPalette) + scale_x_continuous(trans='log10') + theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 14,face="bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 12), legend.position = c(0.8,0.75)) + xlab("Genome Size") + ylab("Average Gene Body CpG o/e")  + geom_text(aes(label = Species), size =3.5, vjust = 0, hjust = 0, color = "black")


dp3 = ggplot(data1, aes(x=repeat_content, y=mean_repeats_CpG, color = Class)) + geom_point(size = 3)
dp3 + scale_colour_manual(values=cbPalette) + scale_x_continuous(trans='log10') + theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 14,face="bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 12)) + xlab("Repeat Content (%)") + ylab("Average Repeat CpG o/e") + geom_text(aes(label = Species), size =2.5, vjust = 0, hjust = 0, color = "black")
, legend.position = c(0.8,0.75)

dp4 = ggplot(data1, aes(x=repeat_content_me, y=mean_repeats_CpG, color = Class)) + geom_point(size = 3)
dp4 + scale_colour_manual(values=cbPalette) + scale_x_continuous(trans='log10') + theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 14,face="bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 12)) + xlab("Repeat Content (%) \n from this study") + ylab("Average Repeat CpG o/e")+ geom_text(aes(label = Species), size =2.5, vjust = 0, hjust = 0, color = "black")




# fig2a: genome size and repeat and CpG in one graph # 
data1$Class = factor(data1$Class, levels = c("Anthozoa","Cubozoa", "Hydrozoa", "Myxozoa", "Scyphozoa"))

dp6 = ggplot(data1, aes(x=repeat_content, y=mean_repeats_CpG))+ geom_point(aes(size = genome_size, color = Class))
fig2a = dp6 + scale_colour_manual(values=cbPalette, na.translate = F) + scale_size_continuous(limits = c(50,1500))+ theme(axis.title = element_text(size=11,face="bold"), axis.text = element_text(size = 11,face="bold"), legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + xlab("Repeat Content (%)") + ylab("Average Repeat CpG o/e") + geom_text_repel(aes(label=abbriviation), size=5, box.padding = unit(0.15, "lines")) + labs(size = "Genome Size", color = "Class") + expand_limits(x = c(10, 60))

dp7 = ggplot(data1, aes(x=repeat_content, y=Mean_GB_CpG))+ geom_point(aes(size = genome_size, color = Class))
fig2b = dp7 + scale_colour_manual(values=cbPalette, na.translate = F) + scale_size_continuous(limits = c(50,1500))+ theme(axis.title = element_text(size=11,face="bold"), axis.text = element_text(size = 11,face="bold"), legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + xlab("Repeat Content (%)") + ylab("Average Gene Body CpG o/e") + geom_text_repel(aes(label=abbriviation), size=5, box.padding = unit(0.15, "lines")) + labs(size = "Genome Size", color = "Class") + expand_limits(x = c(10, 60))
#fig2c: GB CpG changes with repeat CpG#
dp5 = ggplot(data1, aes(x=mean_repeats_CpG, y=Mean_GB_CpG)) + geom_point(aes(size = genome_size, color = Class))
fig2c = dp5 + scale_colour_manual(values=cbPalette, na.translate = F) + scale_size_continuous(limits = c(50,1500)) + theme(axis.title = element_text(size=11,face="bold"), axis.text = element_text(size = 11,face="bold"), legend.title = element_text(size = 10), legend.text = element_text(size = 10)) + xlab("Average Repeat CpG o/e") + ylab("Average Gene Body CpG o/e") + geom_text_repel(aes(label=abbriviation), size=5, box.padding = unit(0.15, "lines"))+ labs(size = "Genome Size", color = "Class")

tiff("Fig2v7.tiff",width = 8, height = 10, units = 'in', res = 300)
ggarrange(fig2a, fig2b, fig2c, labels = c("A", "B", "C"),common.legend = TRUE, nrow = 3, legend = 'right')
dev.off()
jpeg("fig2.jpeg",width = 8, height = 10, units = 'in', res = 300)
ggarrange(fig2a, fig2b, fig2c, labels = c("A", "B", "C"),common.legend = TRUE, nrow = 3, legend = 'right')
dev.off()
# genome size with repeat and GB CpG #

sp1 = ggplot(data1, aes(x=genome_size, y=value, color = variable))+ geom_point(aes(y =mean_repeats_CpG, col = abbriviation, shape=23)) + geom_point(aes(y =Mean_GB_CpG, col = abbriviation, shape=16))
sp1 + scale_x_continuous(trans='log10') + theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 14,face="bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 12)) + xlab("Genome Size") + ylab("CpG o/e")+ labs(shape = "Methylation context", col = "species", element_text(size = 10))
sp1

#repeat content and GB CpG#
sp2 = ggplot(data1, aes(x=repeat_content, y=Mean_GB_CpG))+ geom_point(aes(size = genome_size, color = Class))
sp2 + scale_colour_manual(values=cbPalette) + theme(axis.title = element_text(size=14,face="bold"), axis.text = element_text(size = 14,face="bold"), legend.title = element_text(size = 12), legend.text = element_text(size = 12)) + xlab("Repeat Content (%)") + ylab("Average Gene Body CpG o/e")+ geom_text_repel(aes(label=abbriviation), size=4, box.padding = unit(0.35, "lines"))+ scale_x_continuous(limits = c(0,60)) + scale_size_continuous(limits = c(50,1500))  + labs(fill = "Class", size = "Genome Size")

### stat tests ###
cor.test(data1$repeat_content, data1$mean_repeats_CpG, method = c("pearson", "kendall", "spearman"), use = "complete.obs")
cor.test(data1$repeat_content, data1$mean_repeats_CpG, method = c("kendall", "spearman"), use = "complete.obs")




