

library(glue)

# by Aaron McDaid, from https://stackoverflow.com/questions/8835426/get-filename-and-path-of-sourced-file
get.full.path.to.this.sourced.script = function() {    
    for(i in sys.nframe():1) {  # Go through all the call frames,
                                # in *reverse* order.
        x = sys.frame(i)$ofile
        if(!is.null(x))               # if $ofile exists,
            return(normalizePath(x))  #  then return the full absolute path
    }
}

packageList <- c("pacman", "BiocParallel", "DESeq2", "apeglm", "sva", "vsn", "dplyr", "stats",  "dfidx", "ggplot2", "ggpubr", "ggfortify", "pheatmap", "rmarkdown")
citationList <- list()

pth <- get.full.path.to.this.sourced.script()
out_pth <- file.path(pth, "..", "..", "output", "citations.txt")
# Sinking will fail for too large text chunks
sink(out_pth)
for (package in packageList){
    package |> cat()
    package |> citation() |> cat()
}
sink()

for (package in packageList){
  citationList[[package]] <- citation(package)
}

cat(citationList)
cat("")
cat(glue::glue("Citations have been written to {out_pth}"))