

save_plot <- function(plot, filepath, save_format = pdf) {
    if (!is.null(filepath) && !is.null(save_format)) {
        save_format(filepath)
        show(plot)
        dev.off()
    }
}

plot_surrogates <- function(svseq, dds, filepath=NULL, save_format=pdf) {
    
    par(mfrow=c(2,1),mar=c(3,5,3,1))
    stripchart(svseq$sv[,1] ~ dds$patient.id,vertical=TRUE,main="SV1")
    abline(h=0)
    stripchart(svseq$sv[,2] ~ dds$patient.id,vertical=TRUE,main="SV2")
    abline(h=0)
    res <-  recordPlot()
    if (!is.null(filepath) && !is.null(save_format)) {
        save_format(filepath)
    }
    show(res)
    if (!is.null(filepath) && !is.null(save_format)) {
        dev.off()
    }
    return(res)
}

plot_meanSdPlot <- function(dds, vsd=NULL, filepath=NULL, save_format=pdf) {
    if (is.null(vsd)) {vsd <- vst(dds, blind=FALSE)}
    
    if (!is.null(filepath) && !is.null(save_format)) {
        save_format(filepath)
    }
    res <- meanSdPlot(assay(vsd))
    if (!is.null(filepath) && !is.null(save_format)) {
        dev.off()
    }
    return(res)
}

plot_pheatmap <- function(dds, dds_res, vsd=NULL, dds_shrink_PL_H=NULL, filepath=NULL, save_format=pdf) {
    if (is.null(vsd)) {vsd <- vst(dds, blind=FALSE)}
    if (is.null(dds_shrink_PL_H)) {
        index = which(resultsNames(dds_res) == "condition_PSO_lesional_vs_CTRL_healthy")
        dds_shrink_PL_H <- lfcShrink(dds_res,
                                 coef = index,
                                 type = "apeglm")
    }
    dds <- estimateSizeFactors(dds)
    select <- row.names(head(dds_shrink_PL_H[order(dds_shrink_PL_H$padj), ], 50))
    
    df <- as.data.frame(colData(dds)[, c("disease", "sample.type")])
    rownames(df) <- row.names(df)
    
    custom_color <- list(
        sample.type = c('non.lesional' = 'lightgrey', 'lesional' = '#F44336'),
        disease = c('healthy' = 'lightgrey', 'psoriasis' = 'steelblue', 'AD' = '#FFCC80', 'chronic_AD' = '#D32F2F'))
    
    if (!is.null(filepath) && !is.null(save_format)) {
        save_format(filepath)
    }
    res <- pheatmap(assay(vsd)[select,],
        cluster_rows=TRUE, 
        cluster_cols=TRUE, 
        annotation_col=df,
        fontsize_col = 2,
        fontsize_row = 7,
        show_colnames = FALSE,
        show_rownames = TRUE,
        treeheight_row = 0,
        main = "clustering - top 50 DEGs of PSO_L_vs_H",
        annotation_colors = custom_color,
        color = colorRampPalette(c("#4A148C", "#FFF9C4", "#FF7043"))(100)
    )
    if (!is.null(filepath) && !is.null(save_format)) {
        dev.off()
    }
    
    return(res)
}

plot_pca <- function(dds, vsd=NULL, filepath=NULL, save_format=pdf) {
    if (is.null(vsd)) {vsd <- vst(dds, blind=FALSE)}  

    pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData = TRUE, ntop = 5000)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    if (!is.null(filepath) && !is.null(save_format)) {
        save_format(filepath)
    }
    res <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggtitle("PCA - top 5000 genes") +
      coord_fixed() +
      #stat_ellipse(aes(group = condition), type = "norm", level = 0.95, linewidth = 1) +
      theme(
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        axis.line.x.bottom=element_line(color="black"),
        axis.line.y.left=element_line(color="black"),
        plot.title = element_text(hjust = 0.5)
      ) +
      scale_color_manual(values = c("grey", "#FFB300", "#FFECB3", "#B71C1C", "#EF9A9A", "#283593", "#9FA8DA"))
    
    if (!is.null(filepath) && !is.null(save_format)) {
        dev.off()
    }

    return(res)
}













