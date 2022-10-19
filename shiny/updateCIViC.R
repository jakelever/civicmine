library(data.table)

# URL and filenames for CIViC data
civicURL <- 'https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv'
tempFilename <- 'nightly-ClinicalEvidenceSummaries.tsv.tmp'
outFilename <- 'nightly-ClinicalEvidenceSummaries.tsv'

# Download the CIViC file and load it
download.file(civicURL, tempFilename, method='curl')
civicdb <- fread(tempFilename,sep='\t',header=T)

# Move PubMed IDs to their own column (from citation_id column) to match old CIViC format
civicdb$pubmed_id <- 0
civicdb[civicdb$source_type=='PubMed','pubmed_id'] <- civicdb[civicdb$source_type=='PubMed','citation_id']
civicdb$pubmed_id <- as.integer(civicdb$pubmed_id)

# Check that all the expected columns are there
expectedColumns <- c('pubmed_id','entrez_id','doid','drugs','evidence_type','variant','phenotypes')
foundColumns <- expectedColumns %in% colnames(civicdb)
stopifnot(all(foundColumns)==T)
stopifnot(nrow(civicdb) > 0)

# Write out the smaller version of the CIViC data with only the columns we need
civicdb <- civicdb[, expectedColumns, with=F]
write.table(civicdb, outFilename,row.names=F,sep='\t',quote=F)

# Clean up
file.remove(tempFilename)

print("Success. Downloaded loadable version of CIViCdb")
