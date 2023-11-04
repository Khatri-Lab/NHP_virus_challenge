## checking human MVS -> monkey analogs
## then seeing if any of the modules are the problem
path_0 = "/Users/kalani/Desktop/Stanford_Classes/Projects/KR04_MVS-NHP/"
datPath_0 = paste0(path_0, "data/")
mvs = fread(paste0(datPath_0, "gene_mods/mvs/mvs.csv"))
mvs$es = 0
View(mvs)

library(biomaRt)
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")

datasets <- listDatasets(useMart('ensembl'))
rhesus <- useMart('ensembl', dataset = 'mmulatta_gene_ensembl')
human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

table <- getLDS(
  mart = rhesus,
  attributes = c('ensembl_gene_id', 'vgnc', 'external_gene_name', 'chromosome_name'),
  martL = human,
  attributesL = c('ensembl_gene_id','hgnc_symbol','gene_biotype', 'chromosome_name'))

x = mvs$gene
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",  host="https://uswest.ensembl.org")
  rhesus = useMart("ensembl", dataset = "mmulatta_gene_ensembl",  host="https://uswest.ensembl.org")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("vgnc"), martL = rhesus, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
genes <- convertHumanGeneList(mvs$gene)

