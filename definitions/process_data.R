


library(glue)


create_meta_data <- function(srarun) {
    # must be done by hand
    return(read.csv("input/meta_data_subset_object.csv", header = TRUE))
}

new_deduped_name <- function(count_matrix, k, separator="_duplicate") {
    nm <- count_matrix$X[k]
    parts <- strsplit("nm", separator)[[1]]
    if (length(parts) == 1) {
        n <- 1
    }
    else if (length(parts) == 2) {
        n <- as.integer(parts[2])
        nm <- parts[1]
    }
    else {
        warning("Error while parsing gene name in function deduplicate_count_matrix!")
        n <- as.integer(parts[2])
        nm <- parts[1]
    }
    
    if (is.na(n)) {n <- 1}
    
    return(glue::glue("{nm}{separator}{n}"))
}

deduplicate_count_matrix <- function(count_matrix) {
    
    k <- anyDuplicated(count_matrix$X)
    while (k != 0) {
        count_matrix$X[k] <- new_deduped_name(count_matrix, k)
        k <- anyDuplicated(count_matrix$X)
    }
    return(count_matrix)
}


rename_rows <- function(meta_data) {
    meta_data$X <- gsub("non_lesional", "non.lesional", meta_data$X)
    meta_data$sample.type <- gsub("non_lesional", "non.lesional", meta_data$sample.type)
    row.names(meta_data) <- meta_data[,"X"]
    return(meta_data)
}

columns_to_factor_levels <- function(meta_data) {
    meta_data$condition <- factor(meta_data$condition)
    meta_data$patient.ID <- factor(meta_data$patient.ID)
    meta_data$patient.id <- factor(meta_data$patient.id)
    meta_data$disease <- factor(meta_data$disease)
    meta_data$sample.type <- factor(meta_data$sample.type)
    
    return(meta_data)
}

clean_data <- function(meta_data, count_matrix) {
    subset_lesional <- (meta_data |> filter(sample.type == "lesional"))$patient.id
    subset_non_lesional <- (meta_data |> filter(sample.type == "non.lesional"))$patient.id
    subset_healthy <- (meta_data |> filter(disease=="healthy"))$patient.id
    ids = intersect(subset_lesional, subset_non_lesional)
    ids = union(ids, subset_healthy)
    
    meta_data_subset <- meta_data |> filter(patient.id %in% ids)
    droplevels(meta_data_subset)
    
    #TODO: give meta_data a column sample.name
    count_matrix_subset <- count_matrix[, colnames(count_matrix) %in% meta_data_subset[,"X"]]
    
    return(list(meta_data=meta_data_subset, count_matrix = count_matrix_subset))
}

create_SeqDataSet <- function(count_matrix, meta_data) {
    count_matrix_subset <- count_matrix[, order(colnames(count_matrix))]
    meta_data_subset <- arrange(meta_data,"X")
    dds <- DESeqDataSetFromMatrix(countData = count_matrix_subset,
                                  colData = meta_data_subset,
                                  design = ~ patient.ID + condition)
    keep <- rowSums(counts(dds)) > 1 #keeps row with > 0 counts
    dds <- dds[keep,]
    dds$condition <- relevel(dds$condition, ref = "CTRL_healthy")
    
    return(dds)
}


commit_diff_expr_analysis <- function(dds) {
    dds_res <- DESeq(dds, parallel = TRUE, BPPARAM = SnowParam(4)) #runs the DESeq2 analysis
    return(dds_res)
}

create_surrogate_variables <- function(dds, dds_res, meta_data_subset) {
    dat <- counts(dds_res, normalized=TRUE)
    idx <- rowMeans(dat) > 1
    dat <- dat[idx,]
    mod <- model.matrix(~ patient.ID + condition, meta_data_subset)
    mod0 <- model.matrix(~ 1, colData(dds_res))
    svseq <- svaseq(dat, mod, mod0, n.sv=2)
    return(svseq)
}


commit_diff_expr_analysis_with_surrogates <- function(dds_res, svseq) {
    ddssva <- dds_res
    ddssva$SV1 <- svseq$sv[,1]
    ddssva$SV2 <- svseq$sv[,2]
    dds_sva_res <- DESeq(ddssva)
    
    return(dds_sva_res)
}

extract_results <- function(dds_res, outputpath=NULL) {
    res <- results(dds_res)
    res <- res[order(res$padj),]
    
    if (!is.null(outputpath)) {
        write.csv(as.data.frame(res),
            row.names = TRUE,
            outputpath)
    }
    
    return(res)
}


extract_psoriasis_genes <- function(dds_res,outdir = NULL) {
    result_PSO <- results(dds_res, contrast = c('condition', 'PSO_lesional', 'PSO_non_lesional'), independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH", parallel = TRUE)
    #filter out NAs from the padj columns as well as all padj > 0.1
    filtered_genes <- result_PSO[!is.na(result_PSO$padj) & result_PSO$padj <= 0.1, ]
    #Order and save the filtered results to a CSV file
    filtered_genes_ordered <- filtered_genes[order(filtered_genes$padj), ]
    if (!is.null(outdir)) {
        #write as a csv file
        write.csv(as.data.frame(filtered_genes_ordered), row.names = TRUE, file.path(outdir, "DEGs_PSO_NL_vs_PSO_L.csv"))
        #save the object
        saveRDS(result_PSO, "DEG_Analysis_result_PSO.rds")
    }
}

extract_AD_genes <- function(dds_res,outdir = NULL) {
    result_AD <- results(dds_res, contrast = c('condition', 'AD_lesional', 'AD_non_lesional'), independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH", parallel = TRUE)
    filtered_genes <- result_AD[!is.na(result_AD$padj) & result_AD$padj < 0.1, ]
    filtered_genes_ordered <- filtered_genes[order(filtered_genes$padj), ]
    if (!is.null(outdir)) {
        write.csv(as.data.frame(filtered_genes_ordered), row.names = TRUE, file.path(outdir, "DEGs_AD_NL_vs_AD_L.csv"))
        saveRDS(result_AD, "DEG_Analysis_result_AD.rds")
    }
}

extract_CAD_genes <- function(dds_res,outdir = NULL) {
    result_CAD <- results(dds_res, contrast = c('condition', 'CAD_chronic_lesion', 'CAD_non_lesional'), independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH", parallel = TRUE)
    filtered_genes <- result_CAD[!is.na(result_CAD$padj) & result_CAD$padj < 0.1, ]
    filtered_genes_ordered <- filtered_genes[order(filtered_genes$padj), ]
    if (!is.null(outdir)) {
        write.csv(as.data.frame(filtered_genes_ordered), row.names = TRUE, file.path(outdir, "DEGs_CAD_NL_vs_CAD_L.csv"))
        saveRDS(result_CAD, "DEG_Analysis_result_CAD.rds")
    }
}

extract_gene_differences_psoriasis_healthy <- function(dds_res,outdir = NULL) {
    result_PSO_H <- results(dds_res, contrast = c('condition', 'PSO_lesional', 'CTRL_healthy'), independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH", parallel = TRUE)
    filtered_genes <- result_PSO_H[!is.na(result_PSO_H$padj) & result_PSO_H$padj < 0.1, ]
    filtered_genes_ordered <- filtered_genes[order(filtered_genes$padj), ]
    if (!is.null(outdir)) {
        write.csv(as.data.frame(filtered_genes_ordered), row.names = TRUE, file.path(outdir, "DEGs_CTRL_H_vs_PSO_L.csv"))
        saveRDS(result_PSO_H, "DEG_Analysis_result_PSO_H.rds")
    }
}

extract_gene_differences_AD_healthy <- function(dds_res,outdir = NULL) {
    result_AD_H <- results(dds_res, contrast = c('condition', 'AD_lesional', 'CTRL_healthy'), independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH", parallel = TRUE)
    filtered_genes <- result_AD_H[!is.na(result_AD_H$padj) & result_AD_H$padj < 0.1, ]
    filtered_genes_ordered <- filtered_genes[order(filtered_genes$padj), ]
    if (!is.null(outdir)) {
        write.csv(as.data.frame(filtered_genes_ordered), row.names = TRUE, file.path(outdir, "DEGs_CTRL_H_vs_AD_L.csv"))
        saveRDS(result_AD_H, "DEG_Analysis_result_AD_H.rds")
    }
}

extract_gene_differences_CAD_healthy <- function(dds_res,outdir = NULL) {
    result_CAD_H <- results(dds_res, contrast = c('condition', 'CAD_chronic_lesion', 'CTRL_healthy'), independentFiltering = TRUE, alpha = 0.05, pAdjustMethod = "BH", parallel = TRUE)
    filtered_genes <- result_CAD_H[!is.na(result_CAD_H$padj) & result_CAD_H$padj < 0.1, ]
    filtered_genes_ordered <- filtered_genes[order(filtered_genes$padj), ]
    if (!is.null(outdir)) {
        write.csv(as.data.frame(filtered_genes_ordered), row.names = TRUE, file.path(outdir, "DEGs_CTRL_H_vs_CAD_L.csv"))
        saveRDS(result_CAD_H, "DEG_Analysis_result_CAD_H.rds")
    }
}

