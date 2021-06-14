#######################################
## Differential Expression Functions ##
#######################################

## This function is to write out the files of differentially expressed genes
## It takes in a DESeqResults object, and two strings -- the numerator and denominator used in the analysis -- and writes out csv files

write_files <- function(results, numerator, denominator){
  # these are all the genes that are differentially expressed between the two conditions, not just the significant ones
  write.csv(results, paste0("analysis/output/DE/",numerator,"_",denominator,"_all.csv"), row.names = TRUE, col.names = TRUE, sep = " ")
}

## This function plots the volcano plot
## It takes in a data frame and two strings which are used for the title of the plot

generate_volcano <- function(data_frame, numerator, denominator){
  lfc.threshold = log2(1.5)
  tmp = as.data.frame(data_frame)
  tmp$signif = ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.01, "U1", 
                      ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.05, "U2",
                             ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.1, "U3",
                                    ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.01, "D1", 
                                           ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.05, "D2",
                                                  ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.1, "D3",                                                  "N"))))))
  tmp$signif = factor(tmp$signif, c("N", "U1", "U2", "U3", "D3", "D2", "D1"))
  
  x = ggplot(data=tmp, aes(x=log2FoldChange, y=-log10(padj), colour= signif)) + geom_point(alpha=1.0, size=2.00) +
    ggtitle(paste("Volcano Plot:", numerator, "vs.", denominator)) + scale_x_continuous("log2(fold change)", limits=c(-15, 15)) +    
    scale_y_continuous("-log10(FDR)") + geom_vline(xintercept = lfc.threshold, linetype="dotdash") +
    geom_vline(xintercept = -1*(lfc.threshold), linetype="dotdash") +
    geom_hline(yintercept = -log10(0.1), colour="gray40", linetype="dotdash") +   
    geom_hline(yintercept = -log10(0.05), colour="gray40", linetype="dotdash") + 
    geom_hline(yintercept = -log10(0.01), colour="gray40", linetype="dotdash") + 
    scale_colour_manual("", values=c("#666666", "#d73027", "#f46d43", "#fdae61", "#abd9e9", "#74add1", "#4575b4" ), labels = c("N", "U1", "U2", "U3", "D3", "D2", "D1")) + theme_classic() + theme(legend.position = "none", plot.title = element_text(size = 20), axis.title=element_text(size=16,face="bold"))
  print(x)
  print(table(tmp$signif))
}

## This function generates the MA plots with significant changes above the threshold coloured in red and significant changes below the threshold coloured in blue
## It takes in a DESeqResults object, uses the plotMA function from DESeq2 to obtain the necessary data frame to plot

generate_ma <- function(results){
  df <- DESeq2::plotMA(results, ylim = c(-10,10), colSig = "red", returnData = TRUE)
  plot <- df %>%
    mutate(signif = ifelse(lfc > lfc.threshold & isDE == TRUE, "U", 
                           ifelse(lfc < -lfc.threshold & isDE == TRUE, "D", "N"))) %>%
    ggplot(aes(x=mean, y=lfc, colour = signif)) + 
    geom_point(size = 1.5, alpha = 0.8) + 
    theme_classic() + 
    geom_hline(yintercept=0, colour="grey40", lwd = 1) + 
    #stat_smooth(se = FALSE, method = "loess", color = "red3") + 
    theme_classic() + 
    scale_colour_manual(values=c("#4575b4","#a3a3a3","#d73027"), labels = c("D", "N", "U")) +
    ylim(c(-10,10)) +
    theme(legend.position = "none") +
    ylab("Log fold change") +
    xlab("Mean of normalized counts") +
    scale_x_log10()
  return(plot)
}

############################################
## Visualisation and Enrichment Functions ##
############################################

## This function makes volcano plots with labels of certain genes
## It takes in a DESeqResults object, Flybase Gene IDs (FBgn IDs) of genes we're interested in labelling, and the title of the plot

make_volcano_with_labels <- function(data, labels, title){
  tmp = as.data.frame(data)
  tmp$genelabels = FALSE
  tmp[labels,]$genelabels = TRUE
  tmp$signif = ifelse(tmp$log2FoldChange > lfc.threshold & tmp$padj< 0.1, "U", 
                      ifelse(tmp$log2FoldChange < -1*lfc.threshold & tmp$padj< 0.1, "D","N"))
  tmp$signif = factor(tmp$signif, c("N", "U", "D"))
  myCols = c("#666666", "#d73027","#4575b4")
  names(myCols) = levels(tmp$signif)
  
  ggplot(data=tmp, aes(x=log2FoldChange, y=-log10(padj), colour = signif)) + 
    geom_point(alpha=1.0) + geom_point(data = tmp[tmp$genelabels == TRUE,], color = "black") +
    ggtitle(title) + scale_x_continuous("log2(fold change)", limits=c(-15, 15)) +    
    scale_y_continuous("-log10(FDR)", limits = c(-30, max(-log10(tmp[["padj"]])))) + 
    geom_vline(xintercept = lfc.threshold, linetype="dotdash") +
    geom_vline(xintercept = -1*(lfc.threshold), linetype="dotdash") +
    geom_hline(yintercept = -log10(0.1), colour="gray40", linetype="dotdash") +   
    scale_colour_manual(name = "signif", values = myCols) + theme_classic() + theme(legend.position = "none", plot.title = element_text(size = 20), axis.title=element_text(size=16,face="bold")) +  geom_label_repel(aes(log2FoldChange, -log10(padj)),label = ifelse(tmp$genelabels, as.character(tmp$external_gene_name),""))
}

## This function gets the zscores table for our dorsal targets of interest
## It takes in the full zscores matrix and a data frame of gene_ids that we're interested in. The data frame needs to contain the external gene names of the genes we're interested in because we use that to relabel the rownames of the zscore matrix

get_zscores = function(zscores_mat, gene_id_df){
  zscores = zscores_mat[rownames(zscores_mat) %in% rownames(gene_id_df),] %>%
    data.frame()
  zscores$gene_name = gene_id_df$external_gene_name[match(rownames(zscores), rownames(gene_id_df))]
  rownames(zscores) = zscores$gene_name
  zscores = zscores[,1:6]
}


## This function returns the counts of the genes we're interested in
## It takes in the FlyBase gene IDs (FBgn IDs) and returns a data frame of counts obtained from normalised dds counts

get_df <- function(gene_id){
  foo = as.data.frame(nc[rownames(nc) %in% gene_id,])
  foo$gene_names = ensembl_genes$external_gene_name[match(rownames(foo), ensembl_genes$gene_id)]
  rownames(foo) = foo$gene_names
  foo = foo[,1:6]
  return(foo)
}

## This function plots the barplots using results from EnrichR 
## It takes in an enrichR result, name of the plot, and the number of terms to show (default being 20)

plot_enrichr <- function(data.frame, name, showCategory = 20){
  plot = data.frame %>%
    mutate(Term = gsub("\\([^()]*\\)", "", Term),
           Term = factor(Term, levels = rev(Term))) %>%
    mutate(Annotated = as.numeric(str_extract(as.character(Overlap), "\\d+$")),
           Significant = as.numeric(str_extract(as.character(Overlap), "^\\d+")),
           Ratio = Significant/Annotated) %>%
    arrange(Adjusted.P.value) %>%
    head(showCategory) %>%
    ggplot(aes(x = reorder(Term,desc(Adjusted.P.value)), y = Ratio, fill = Adjusted.P.value)) +
    geom_bar(stat = "identity") +
    ggpubr::rotate() +
    xlab(NULL) + 
    ylab("Gene Ratio") +
    scale_fill_continuous(low = "red", high = "blue") + 
    scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
    labs(title = name,
         fill = "Adjusted p-value") +
    guides(fill = guide_colorbar(title = "Adjusted p-value", reverse = TRUE)) +
    theme(axis.text.x = element_text(size=8))
  return(plot)
}

## This function is purely for extracting the legend of a plot to be used in plot_grid
extracting_legend <- function(df, name){
  fig <- df[name,] %>%
    melt() %>%
    `colnames<-`(c("Sample", "Value")) %>%
    mutate(Condition = ifelse(grepl("WT", Sample), "WT", "LOF_Mutant")) %>%
    ggplot(aes(x=Condition, y=Value, fill = Condition)) + 
    geom_boxplot() + 
    theme_classic() + 
    theme(axis.text.x = element_blank(), legend.position="top") +
    xlab(NULL) + 
    ylab(NULL) +
    scale_x_discrete(limits = rev) +
    expand_limits(y=0)
  
  legend <- get_legend(fig)
  return(legend)
}

## This function plots the boxplots for gene expression profiles
## It takes in a data frame from which to get the normalised counts, the name of the gene to plot, the lower y-axis limit and upper y-axis limit
gene_boxplot <- function(df, name, lower_y_limit = 0, upper_y_limit = NA){
  fig <- df[name,] %>%
    melt() %>%
    `colnames<-`(c("Sample", "Value")) %>%
    mutate(Condition = ifelse(grepl("WT", Sample), "WT", "LOF_Mutant")) %>%
    ggplot(aes(x=Condition, y=Value, fill = Condition)) + 
    geom_boxplot() + 
    theme_classic() + 
    theme(axis.text.x = element_blank(), legend.position="none") +
    labs(title = name) + 
    xlab(NULL) + 
    ylab(NULL) +
    scale_x_discrete(limits = rev) +
    expand_limits(y=0) +
    scale_y_continuous(limits = c(lower_y_limit, upper_y_limit), breaks = seq(lower_y_limit,upper_y_limit,by = (upper_y_limit-lower_y_limit)/5))
  return(fig)
}