library(data.table)

civicURL <- 'https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv'
tempFilename <- 'nightly-ClinicalEvidenceSummaries.tsv.tmp'
outFilename <- 'nightly-ClinicalEvidenceSummaries.tsv'

download.file(civicURL, tempFilename, method='curl')

expectedColumns <- c('pubmed_id','entrez_id','doid','drugs','evidence_type','variant','phenotypes')
civicdb <- fread(tempFilename,sep='\t',header=T)

civicdb$pubmed_id <- 0
civicdb[civicdb$source_type=='PubMed','pubmed_id'] <- civicdb[civicdb$source_type=='PubMed','citation_id']
civicdb$pubmed_id <- as.integer(civicdb$pubmed_id)

foundColumns <- expectedColumns %in% colnames(civicdb)
stopifnot(all(foundColumns)==T)
stopifnot(nrow(civicdb) > 0)

civicdb <- civicdb[, expectedColumns, with=F]
#file.rename(tempFilename,outFilename)
write.table(civicdb, outFilename,row.names=F,sep='\t',quote=F)
file.remove(tempFilename)

print("Success. Downloaded loadable version of CIViCdb")
