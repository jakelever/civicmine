library(data.table)

download <- function(url) {
  tempFilename <- 'table.tmp'
  
  if(.Platform$OS.type == "unix") {
    download.file(url, tempFilename, method='wget', extra='--no-check-certificate')
  } else {
    download.file(url, tempFilename, method='curl')
  }
  
  table <- fread(tempFilename,sep='\t',header=T)
  
  # Clean up
  file.remove(tempFilename)
  
  return(table)
}

outFilename <- 'civicdb.tsv'

# URL and filenames for CIViC data
clinicalEvidenceUrl <- 'https://civicdb.org/downloads/nightly/nightly-ClinicalEvidenceSummaries.tsv'
molecularProfileUrl <- 'https://civicdb.org/downloads/nightly/nightly-MolecularProfileSummaries.tsv'
variantUrl <- 'https://civicdb.org/downloads/nightly/nightly-VariantSummaries.tsv'

clinicalEvidenceTable <- download(clinicalEvidenceUrl)
molecularProfileTable <- download(molecularProfileUrl)
variantTable <- download(variantUrl)

setnames(molecularProfileTable, "variant_ids", "variant_id")
#variantTable['variant_id'] <- as.character(variantTable['variant_id'])
molecularProfileTable[, variant_id := suppressWarnings(as.numeric(variant_id))]

civicdb <- merge(molecularProfileTable, variantTable, by="variant_id")
civicdb <- merge(clinicalEvidenceTable, civicdb, by="molecular_profile_id")


# Move PubMed IDs to their own column (from citation_id column) to match old CIViC format
civicdb$pubmed_id <- 0
civicdb[civicdb$source_type=='PubMed','pubmed_id'] <- civicdb[civicdb$source_type=='PubMed','citation_id']
civicdb$pubmed_id <- as.integer(civicdb$pubmed_id)

# Check that all the expected columns are there
expectedColumns <- c('pubmed_id','entrez_id','doid','therapies','evidence_type','variant','phenotypes')
foundColumns <- expectedColumns %in% colnames(civicdb)
stopifnot(all(foundColumns)==T)
stopifnot(nrow(civicdb) > 0)

# Write out the smaller version of the CIViC data with only the columns we need
civicdb <- civicdb[, expectedColumns, with=F]
write.table(civicdb, outFilename,row.names=F,sep='\t',quote=F)

print("Success. Downloaded loadable version of CIViCdb")
