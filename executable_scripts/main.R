


###############################################################################
###---------------------------import packages-------------------------------###
###############################################################################

library("pacman")
library("BiocParallel")
library("apeglm")
library("DESeq2")
library("sva")
library("vsn")
library("dplyr")
library("stats")
library("dfidx")
library("ggplot2")
library("ggpubr")
library("ggfortify")
library("pheatmap")
library("rmarkdown")
library("hexbin")


# by Aaron McDaid, from https://stackoverflow.com/questions/8835426/get-filename-and-path-of-sourced-file
get.full.path.to.this.sourced.script = function() {    
    for(i in sys.nframe():1) {  # Go through all the call frames,
                                # in *reverse* order.
        x = sys.frame(i)$ofile
        if(!is.null(x))               # if $ofile exists,
            return(normalizePath(x))  #  then return the full absolute path
    }
}

pth <- get.full.path.to.this.sourced.script()
source(file.path(pth,"..","..", "definitions", "process_data.R"))
source(file.path(pth,"..","..", "definitions", "visualize_data.R"))

###---------------cleaning input data-------###
print("Step 1/5: Cleaning input data...")
meta_data <- read.csv(file.path(pth,"..","..","input","meta_data_subset_object.csv"), header = TRUE) |>
                rename_rows() |> 
                columns_to_factor_levels()
count_matrix<- read.delim(file.path(pth,"..","..","input","GSE121212_readcount.txt"), header = TRUE) |>
                deduplicate_count_matrix()

cleaned_data = clean_data(meta_data, count_matrix)

meta_data = cleaned_data$meta_data
count_matrix = cleaned_data$count_matrix



count_matrix <- count_matrix[, order(colnames(count_matrix))]
meta_data <- meta_data[order(meta_data[,"X"]), ]

dds <- create_SeqDataSet(count_matrix, meta_data)

###############################################################################
###----------------------------DE gene analysis-----------------------------###
###############################################################################

print("Step 2/5: Starting differential gene anaylsis...")
#dds_res <- commit_diff_expr_analysis(dds)
dds_res <- DESeq(dds, parallel = TRUE, BPPARAM = SnowParam(4))
dds_shrink_PL_H <- lfcShrink(dds,
                             coef = 43,
                             type = "apeglm")


# for removing batch effects with large sample sizes > 50
print("Step 3/5: Starting differential gene anaylsis with surrogate variables...")
svseq <- create_surrogate_variables(dds, dds_res, meta_data)

dds_sva_res <- commit_diff_expr_analysis_with_surrogates(dds, svseq)

#dds_sva_res <- DESeq(ddssva)

print("Step 4/5: Writing results to directory output...")

outpth_resOrdered <- file.path(pth,'..','..', 'output', 'resOrdered.csv')
print(glue::glue("Writing results to {outpth_resOrdered}..."))
results <- extract_results(dds_res, outpth_resOrdered)

outpth_resOrderedCorrected <- file.path(pth,"..","..", "output", "resOrdered_corrected.csv")
print(glue::glue('Writing rectified results to {outpth_resOrderedCorrected}...'))
results_sva <- extract_results(dds_sva_res, outpth_resOrderedCorrected)

outpth_diff_expr <- file.path(pth,"..","..", "output", "normalized_counts_dds.csv")
print(glue::glue('Writing differential expression dataset to {outpth_diff_expr}...'))
write.csv(assay(dds, normalized = TRUE), file = outpth_diff_expr)

outpth_meta_data <- file.path(pth,"..","..", "output", "meta_data_subset_object.csv")
print(glue::glue('Writing metadata to {outpth_meta_data}...'))
write.csv(meta_data, file = outpth_meta_data)

################################################################################
###---------------------------Making Visualizations--------------------------###
################################################################################

print("Step 5/5: Saving visualizations to directory output...")

vsd <- vst(dds, blind=FALSE)

outdir <- file.path(pth,"..","..", "output")
plot_meanSdPlot(dds, vsd, filepath=file.path(outdir, "meanSdPlot.pdf"))
plot_pheatmap(dds, vsd, filepath=file.path(outdir, "pheatmap.pdf"))
plot_pca(dds, vsd, filepath=file.path(outdir, "pcs_plot.pdf"))