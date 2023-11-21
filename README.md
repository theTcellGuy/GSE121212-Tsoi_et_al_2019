# GSE121212-Tsoi_et_al_2019

This is an R script to conduct differential expression analysis of the public data set GSE121212 on the GEO data base repository. The data were created by Tsoi _et al._ in their study from 2019.
It starts from the meta data and processed reads.

doi: https://doi.org/10.1016/j.jid.2018.12.0183
GEO data set link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121212
GEO accession number: GSE121212


The processed read count matrix can be downloaded from the linked repository directly.
The raw counts as well as the run meta data can be downloaded from the SRA run selector https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA496323

The meta data table built by trimming the original file to hold only necessary information. The sample names from the column 'Sample Names' were matched with the sample names listed in the GEO page (link 'Geo data set link') to obtain condition information and all other information was added by hand to obtain the final meta data file provided as an object.

Key steps are saved as .rds objects and .csv files. The normalized counts and differential expression analysis was visulaized in Python using seaborn, and matplotlib.
