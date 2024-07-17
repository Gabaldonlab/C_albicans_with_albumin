library(DESeq2)
library("RColorBrewer")
library(reshape2)
library(ggplot2)
library(ggrepel)
library(tibble)   ### for inserting coulmn in dataframe
library(VennDiagram)
library(openxlsx)
library(patchwork)
library(GO.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)


yeast<- read.table("~/Dropbox (Personal)/calb_albumin/yeast_data.tsv", sep = "\t", check.names = F)

### PCA plot for all data
metadata<-NULL
metadata<-cbind(metadata,colnames(yeast))
metadata<-as.data.frame(metadata)
metadata$Time <- ifelse(grepl("0h", metadata$V1),"0.5h",
                              ifelse(grepl("3h", metadata$V1),"3h","24h"))
metadata$Condition <- ifelse(grepl("HSA", metadata$V1),"+ albumin","control (- albumin)")
metadata$Media <- ifelse(grepl("plastic", metadata$V1),"WT only","WT on host cells")

metadata_final <- metadata[,-c(1)]


colData_pca <- DataFrame(metadata_final)

pre_dds_pca <- DESeqDataSetFromMatrix(yeast, colData_pca, design = ~Time)#just a fake desing, does not affect the data normalization and transforamtion
dds_pca <- DESeq(pre_dds_pca)

### PCA
for_pca_vst<-vst(dds_pca)
for_pca_vst_counts<-assay(for_pca_vst)

pca <- prcomp(t(for_pca_vst_counts))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
df_plotting<-data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],colData_pca) 

df_plotting$Time <- factor(df_plotting$Time, levels = c("0.5h","3h","24h"))


df_plotting <- df_plotting[df_plotting$Media ==  "WT on host cells",]

ggplot(df_plotting, aes(PC1,PC2,label=rownames(df_plotting)))+
  stat_ellipse(geom = "polygon",aes(fill = Time),alpha = 0.1,level = 0.90)+
  stat_ellipse(aes(fill = Time), level = 0.90)+
  scale_fill_manual(values = c("red", "darkblue", "yellow"))+
  #geom_point(aes(shape=Media), size=4.5)+
  #geom_point(aes(color=Condition,shape=Media), size=4)+
  geom_point(aes(color=Condition), size=4)+
  
  #geom_text_repel(label=rownames(df_plotting), size = 3, segment.color = "black",arrow = arrow(length = unit(0, 'npc')))+
  xlab(paste0("PC1: ",round(percentVar[1]*100),"% variance")) +
  ylab(paste0("PC2: ",round(percentVar[2]*100),"% variance")) + 
  theme_bw()+
  scale_shape_manual(values=c(17, 15))+
  scale_color_manual(breaks = c("+ albumin", "control (- albumin)"),values=c("#addead", "#a0a0a4"))+
guides(fill=guide_legend(title="Time point"))

### DE analysis
################################################################################
######################     separate comparisons       ##########################
################################################################################

colData_separate<-metadata_final


colData_separate <- as.data.frame(apply(colData_separate, 2, function(x) gsub(" ", "_", x)))



colData_separate$Albumin <- ifelse(grepl("control", colData_separate$Condition),"noAlbumin","withAlbumin")

colData_separate$Time <- factor(colData_separate$Time, levels = c("0.5h","3h","24h"))

colData_separate$Albumin <- factor(colData_separate$Albumin, levels = c("noAlbumin","withAlbumin"))

#colData_separate$Media <-  ifelse(grepl("only", colData_separate$Media),"Plastic","Host")
#colData_separate$Media <- factor(colData_separate$Media, levels = c("Plastic","Host"))
colData_separate$Group <- as.factor(paste0(colData_separate$Time,"_", colData_separate$Albumin))

colData_separate <- colData_separate[,-c(2,3)]

pre_dds_alb <- DESeqDataSetFromMatrix(yeast, colData_separate, design = ~Group)
dds_alb <- DESeq(pre_dds_alb)
resultsNames(dds_alb)

res_0h_host_alb_vs_nonalb_separate<-results(dds_alb, contrast = c("Group","0.5h_withAlbumin", "0.5h_noAlbumin"))

res_3h_host_alb_vs_nonalb_separate<-results(dds_alb, contrast = c("Group","3h_withAlbumin", "3h_noAlbumin"))

res_24h_host_alb_vs_nonalb_separate<-results(dds_alb, contrast = c("Group","24h_withAlbumin", "24h_noAlbumin"))





###

library(ggbreak)

plot_volcano <- function(res, title) {
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$significant <- "Not Significant"
  res_df$significant[res_df$padj < 0.01 & res_df$log2FoldChange > 1] <- "Upregulated"
  res_df$significant[res_df$padj < 0.01 & res_df$log2FoldChange < -1] <- "Downregulated"
  
  
  # Select top 5 upregulated and downregulated genes for labeling
  top5_up <- res_df[res_df$significant == "Upregulated", ] %>%
    arrange(padj) %>%
    head(5)
  top5_down <- res_df[res_df$significant == "Downregulated", ] %>%
    arrange(padj) %>%
    head(5)
  top_genes <- rbind(top5_up, top5_down)
  
  
  num_upregulated <- sum(res_df$significant == "Upregulated")
  num_downregulated <- sum(res_df$significant == "Downregulated")
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point() +
    scale_color_manual(values = c("blue", "black", "red")) +
    theme_minimal() +
    labs(title = title, x = "log2(Fold Change)", y = "-log10(Adjusted p-value)") +
    geom_label_repel(data = top_genes,
                     aes(label = gene), size = 3, force = 5, max.overlaps = Inf,
                     box.padding = 0.5, fill = "white") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray") +  # Horizontal cutoff line
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +     # Vertical cutoff lines
    scale_y_break(c(40, 70), scales = c(1, 50))+ # Adjust these values based on your data
    theme(legend.position = "none",  # Remove the legend
          axis.line = element_line(size = 1, color = "black"),  # Thicker axis lines
          axis.line.x.bottom = element_line(size = 1, color = "black"),
          axis.line.y.left = element_line(size = 1, color = "black"),
          axis.text = element_text(size = 12),  # Larger axis numbers
          axis.ticks = element_line(size = 1),
          axis.ticks.length = unit(0.2, "cm")) +  # Thicker ticks
    
    annotate("text", x = -8, y = 0, size=8,label = paste(num_downregulated), color = "blue", hjust = 0) +
    annotate("text", x = 8, y = 0, size=8,label = paste(num_upregulated), color = "red", hjust = 1)
}

volcano_0.5h <- plot_volcano(res_0h_host_alb_vs_nonalb_separate, "0.5 h")
volcano_0.5h
volcano_3h <- plot_volcano(res_3h_host_alb_vs_nonalb_separate, "3 h")
volcano_3h
volcano_24h <- plot_volcano(res_24h_host_alb_vs_nonalb_separate, "24 h")


# Arrange plots side by side
library(gridExtra)
grid_plot <- arrangeGrob(volcano_0.5h, volcano_3h,volcano_24h, ncol = 3)
ggsave("~/Dropbox/calb_albumin/fixed_results_050524/volcano_plots.png",grid_plot, device = "png",units="in", width=20, height=6, dpi=600) 
ggsave("~/Dropbox/calb_albumin/fixed_results_050524/volcano_plots.pdf",grid_plot, device = "pdf",units="in", width=20, height=6, dpi=600) 








############## Excel

wb_DE <- createWorkbook("calb_albumin")
header<-c("gene", "baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
gene_descr<-read.csv("~/Dropbox (Personal)/calb_albumin/C_albicans_SC5314_A22_current_chromosomal_feature.tab", sep = "\t",skip = 8, header = F)
gene_descr<-gene_descr[grepl("_A$", gene_descr$V1),]



gene_descr_short <- gene_descr[,c(1,3,11)]
b<-data.frame(do.call('rbind', strsplit(as.character(gene_descr_short$V3),'|',fixed=TRUE)))
gene_descr_short$V12 <- b$X1
gene_descr_short <- gene_descr_short[,c(1,3,4)]

colnames(gene_descr_short)<-c("genes","gene_description","ORF_id")

for (file in file_names_sep){
  assign("fung",as.data.frame(eval(parse(text=paste0(file)))))
  fung <- cbind("genes"=rownames(fung), data.frame(fung, row.names=NULL))
  
  fung$regulation <- ifelse(
    fung$log2FoldChange >1 & 
      fung$padj < 0.01 &
      !is.na(fung$padj),"Up_regulated",
    ifelse(
      fung$log2FoldChange < -1 &
        fung$padj < 0.01 & 
        !is.na(fung$padj),"Down_regulated", "Non-significant"))
  
  fung_with_gene_descr <- merge(fung,gene_descr_short,by="genes")
  sheet_name<-gsub("_separate","",file)
  addWorksheet(wb_DE, sheet_name)
  writeData(wb_DE, sheet_name, fung_with_gene_descr,keepNA = T)
}

#saveWorkbook(wb_DE, "~/Dropbox/calb_albumin/fixed_results_050524/DE_data_with_ORF.xlsx", overwrite = TRUE)


########## GO TERMS

### FUNGAL GO TERMS
### Finding MF,BP and CC terms
all_go_tems<-as.data.frame(GOTERM)

all_go_tems<-all_go_tems[!(all_go_tems$go_id %in% c("GO:0003674","GO:0008150","GO:0005575")),]


MF_terms<-all_go_tems$go_id[all_go_tems$Ontology=="MF"]
BP_terms<-all_go_tems$go_id[all_go_tems$Ontology=="BP"]
CC_terms<-all_go_tems$go_id[all_go_tems$Ontology=="CC"]


### Association of GOID with description
goterms <- Term(GOTERM)
a<-as.data.frame(goterms)
go_names<-cbind(row.names(a),a)


### CALB GOIDS
calb_go<-read.table("~/Dropbox (Personal)/calb_albumin/calb_go.txt")


MF_universe<-calb_go[calb_go$V1 %in% MF_terms,]
BP_universe<-calb_go[calb_go$V1 %in% BP_terms,]
CC_universe<-calb_go[calb_go$V1 %in% CC_terms,]









#for (res in file_names_sep){

for (res in c("res_3h_host_alb_vs_nonalb_separate","res_24h_host_alb_vs_nonalb_separate")){
  df <- get(res)
  
  up <- row.names(df[!is.na(df$padj) & df$padj<0.01 & df$log2FoldChange > 1,])
  down <- row.names(df[!is.na(df$padj) & df$padj<0.01 & df$log2FoldChange < -1,])
  
  for (reg in c("up","down")){
    genes <- get(reg)
    if (!is.null(genes)){
      ego<-enricher(genes, pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH", universe = as.character(BP_universe$V2), minGSSize = 2, 
                    maxGSSize = 10000, TERM2GENE = BP_universe,TERM2NAME = go_names)
      
      
      if (!is.null(ego)){
        p<-dotplot(ego,showCategory = 15)
        ggplot2::ggsave(sprintf("~/Dropbox/calb_albumin/fixed_results_050524/%s_%s_BP.png",res,reg),
                        units="in", width=10, heigh=15, dpi=600)
        write.table(ego,sprintf("~/Dropbox/calb_albumin/fixed_results_050524/%s_%s_BP.txt",res,reg), sep="\t")
        ego_df <- as.data.frame(ego)
        ego_df <- cbind.data.frame("Cluster"=res,ego_df)
        assign(paste0("GO_df_",res,"_",reg),as.data.frame(ego_df), envir = .GlobalEnv)
        
      }
    }
  }
}




#### composite GO

### Up-regualted composite

all_up<-NULL
for (df in grep("GO_df_res_.*host_alb_vs_nonalb_separate_up",names(.GlobalEnv),value=TRUE)){
  df <- get(df)
  all_up <- rbind(all_up,df)
}



all_up$Media <- ifelse(grepl("host",all_up$Cluster), "WT on host","WT only")


all_up$Cluster <- gsub("_alb_vs_nonalb_separate","",all_up$Cluster)
all_up$Cluster <- gsub("res_","",all_up$Cluster)
all_up$Cluster <- factor(all_up$Cluster, levels = c("0h_host","3h_host","24h_host","0h_plastic","3h_plastic","24h_plastic"))

mydf <- data.frame(Entrez=c('1', '103', '1000', '100101467',
                            '100127206', '100128071'),
                   logFC = c(1.1, -0.5, 5, 2.5, -3, 3),
                   group = c('A', 'A', 'A', 'B', 'B', 'B'),
                   othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))

to_hack <- compareCluster(Entrez~group, data=mydf,
                          fun='enrichGO',OrgDb='org.Hs.eg.db')


to_hack_up <- to_hack
to_hack_up@compareClusterResult <- all_up

simpl_up<-  simplify(to_hack_up, cutoff=0.7, by="p.adjust", select_fun=min)



# Assuming `simpl_up` is your compareClusterResult object
# Extract the compareClusterResult data frame
simpl_up_data <- simpl_up@compareClusterResult

# Remove "_host" from the Cluster column
simpl_up_data$Cluster <- gsub("h_host", " h", simpl_up_data$Cluster)

# Set the order of factors in the Cluster column
simpl_up_data$Cluster <- factor(simpl_up_data$Cluster, levels = c("3 h", "24 h"))

# Create a new compareClusterResult object with modified labels and ordered factors
simpl_up_modified <- simpl_up
simpl_up_modified@compareClusterResult <- simpl_up_data


{
  p <- dotplot(simpl_up_modified, showCategory = 8) + 
    ggplot2::facet_grid(~Media, scales = "free")+
        xlab("") + 
    scale_y_discrete(guide = guide_axis(n.dodge = 1)) +
    theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      strip.text.x = element_blank()

    )
  
  ggplot2::ggsave("~/Dropbox/calb_albumin/fixed_results_050524/FIG_3C_composite_GO_UP.png", 
                  plot = p, units = "in", width = 15, height = 8, dpi = 600)
  ggplot2::ggsave("~/Dropbox/calb_albumin/fixed_results_050524/FIG_3C_composite_GO_UP.pdf", 
                  plot = p, units = "in", width = 15, height = 8, dpi = 600)
}
    
 
#######################################
####### Down-regualted composite
  
  all_down<-NULL
  for (df in grep("GO_df_res_.*separate_down",names(.GlobalEnv),value=TRUE)){
    df <- get(df)
    all_down <- rbind(all_down,df)
  }

  all_down$Media <- ifelse(grepl("host",all_down$Cluster), "WT on host","WT only")
   all_down$Cluster <- gsub("_alb_vs_nonalb_separate","",all_down$Cluster)
  all_down$Cluster <- gsub("res_","",all_down$Cluster)
  all_down$Cluster <- factor(all_down$Cluster, levels = c("0h_host","3h_host","24h_host","0h_plastic","3h_plastic","24h_plastic"))

  mydf <- data.frame(Entrez=c('1', '100', '1000', '100101467',
                              '100127206', '100128071'),
                     logFC = c(1.1, -0.5, 5, 2.5, -3, 3),
                     group = c('A', 'A', 'A', 'B', 'B', 'B'),
                     othergroup = c('good', 'good', 'bad', 'bad', 'good', 'bad'))
  
  to_hack <- compareCluster(Entrez~group+othergroup, data=mydf,
                            fun='enrichGO', OrgDb='org.Hs.eg.db')
  
  to_hack_down <- to_hack
  to_hack_down@compareClusterResult <- all_down
  
  
  
  simpl_down<-  simplify(to_hack_down, cutoff=0.7, by="p.adjust", select_fun=min)
  
  
  # Assuming `simpl_up` is your compareClusterResult object
  # Extract the compareClusterResult data frame
  simpl_down_data <- simpl_down@compareClusterResult
  
  # Remove "_host" from the Cluster column
  simpl_down_data$Cluster <- gsub("h_host", " h", simpl_down_data$Cluster)
  
  # Set the order of factors in the Cluster column
  simpl_down_data$Cluster <- factor(simpl_down_data$Cluster, levels = c("3 h", "24 h"))
  
  # Create a new compareClusterResult object with modified labels and ordered factors
  simpl_down_modified <- simpl_down
  simpl_down_modified@compareClusterResult <- simpl_down_data
  
  
  
  
  {
    p <- dotplot(simpl_down_modified, showCategory = 8) + 
      ggplot2::facet_grid(~Media, scales = "free")+
      xlab("") + 
      scale_y_discrete(guide = guide_axis(n.dodge = 1)) +
      theme(
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        strip.text.x = element_blank()
        
      )
    
    ggplot2::ggsave("~/Dropbox/calb_albumin/fixed_results_050524/FIG_3D_composite_GO_DOWN.png", 
                    plot = p, units = "in", width = 10, height = 8, dpi = 600)
    ggplot2::ggsave("~/Dropbox/calb_albumin/fixed_results_050524/FIG_3D_composite_GO_DOWN.pdf", 
                    plot = p, units = "in", width = 10, height = 8, dpi = 600)
  }
  
  
  
  