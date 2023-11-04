library(babelgene)
require(biomaRt)

# https://support.bioconductor.org/p/9144001/
# listDatasets(useMart('ensembl'))
require(biomaRt)
datasets <- listDatasets(useMart('ensembl'))
# datasets[grep('mfasciculari', datasets[,1]),]
# datasets[grep('mnemestrina', datasets[,1]),]

cynomolgus <- useMart('ensembl', dataset = 'mfascicularis_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
pigtailed <- useMart('ensembl', dataset = 'mnemestrina_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")

table_cyno <- getLDS(
  mart = cynomolgus,
  attributes = c('ensembl_gene_id', 'vgnc', 'external_gene_name', 'chromosome_name'),
  martL = human,
  attributesL = c('ensembl_gene_id','hgnc_symbol','gene_biotype', 'chromosome_name'))

table_pig <- getLDS(
  mart = pigtailed,
  attributes = c('ensembl_gene_id', 'vgnc', 'external_gene_name', 'chromosome_name'),
  martL = human,
  attributesL = c('ensembl_gene_id','hgnc_symbol','gene_biotype', 'chromosome_name'))

inputDat = dat
dataset_name = "Baize"
### FUNCTIONS
convert_monkey_to_human_orthologs = function(inputDat, dataset_name) {
  species = unique(inputDat$originalData[[dataset_name]]$pheno$species_alignment)
  message(unique(inputDat$originalData[[dataset_name]]$pheno$dataset2))
  if (grepl("human|Human|sapien", species)) {
    return(inputDat)
  }
  if (grepl("rhesus|mulatta", species)) {
    message("rhesus data")
    return(convert_rhesus_to_human_orthologs(inputDat, dataset_name))
  }
  else {
    if(grepl("pig|nemestrina", species)) {
      message("pigtailed data")
      keys_convert = table_pig
    }
    if(grepl("crab|cynomolgus|fascicularis", species)) {
      message("cynomolgus data")
      keys_convert = table_cyno
    }
    object_one = inputDat$originalData[[dataset_name]]
    keys_one = as.data.frame(object_one$keys)
    colnames(keys_one) = "symbol"
    # rownames(keys_one) = keys_one$symbol ## remove
    keys_one$probe = rownames(keys_one)
    keys_convert = keys_convert[!(keys_convert$Gene.name == ""),]
    keys_convert = keys_convert[!is.na(keys_convert$HGNC.symbol),]
    keys_convert = keys_convert[,c("Gene.name","HGNC.symbol")] %>% unique()
    keys_new = merge(keys_one, keys_convert[,c("Gene.name","HGNC.symbol")], by.x = "symbol", by.y = "Gene.name",  all.x = TRUE, all.y = FALSE)
    keys_new$keep = ifelse(!is.na(keys_new$HGNC.symbol), keys_new$HGNC.symbol, keys_new$symbol)
    keys_new <- keys_new[match(names(object_one$keys), keys_new$probe),]
    # keys_new <- keys_new[match(object_one$keys, keys_new$probe),] ## remove
    keys_new_list = keys_new$keep
    setequal(rownames(object_one$expr), keys_new$probe); setequal(names(object_one$keys), keys_new$probe)
    names(keys_new_list) = keys_new$probe
    # if(unique(object_one$pheno$platform) == "Illumina") {
    #   identical(rownames(object_one$expr), keys_new$probe)
    #   rownames(object_one$expr) = keys_new$keep
    #   names(keys_new_list) = keys_new$keep
    # }
    object_one$keys = keys_new_list
    if (checkDataObject(object_one, "Dataset") ){
      print(checkDataObject(object_one, "Dataset") )
      inputDat$originalData[[dataset_name]] = object_one
      return(inputDat)
    }
    else {
      message("something is off")
    }
  }
}


convert_rhesus_to_human_orthologs = function(inputDat, dataset_name) {
  object_one = inputDat$originalData[[dataset_name]]
  keys_one = as.data.frame(object_one$keys)
  colnames(keys_one) = "symbol"
  keys_one$probe = rownames(keys_one)
  genes_dat = keys_one$symbol
  genes_dat = genes_dat[!(genes_dat %in% c("", NA))]
  keys_convert = orthologs(genes = genes_dat, species = "rhesus macaque", human = FALSE)
  keys_new = merge(keys_one, keys_convert[,c("symbol","human_symbol")], by = "symbol", all.x = TRUE, all.y = FALSE)
  keys_new$keep = ifelse(!is.na(keys_new$human_symbol), keys_new$human_symbol, keys_new$symbol)
  keys_new <- keys_new[match(names(object_one$keys), keys_new$probe),]
  keys_new_list = keys_new$keep
  names(keys_new_list) = keys_new$probe
  if(unique(object_one$pheno$platform) == "Illumina") {
    identical(rownames(object_one$expr), keys_new$probe)
    rownames(object_one$expr) = keys_new$keep
    names(keys_new_list) = keys_new$keep
  }
  object_one$keys = keys_new_list
  if (checkDataObject(object_one, "Dataset") ){
    print(checkDataObject(object_one, "Dataset") )
    inputDat$originalData[[dataset_name]] = object_one
    return(inputDat)
  }
  else {
    message("something is off")
  }
}
