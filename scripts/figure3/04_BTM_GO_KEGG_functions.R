btm_path = paste0(path,"scripts/BTMS/")


## take deg subset, path_category ("KEGG","Process")
runGOPathAnalysis = function(genes, path_category = "Process", all_genes_measured, species_val = 9606){
  string_db <- STRINGdb$new(version="11", species=as.numeric(species_val), score_threshold=200, input_directory="")
  
  allGenesMeasured = all_genes_measured
  allGenesMeasured = data.frame(allGenesMeasured)
  
  allGenes_mapped <- string_db$map(allGenesMeasured, "allGenesMeasured", removeUnmappedRows = TRUE )
  length(allGenesMeasured$allGenesMeasured)
  backgroundV <- allGenes_mapped$STRING_id
  string_db$set_background(backgroundV)
  genes = as.data.frame(genes)
  colnames(genes) = "gene"
  res_mapped <- string_db$map(genes, "gene", removeUnmappedRows = TRUE )
  hits_resg = res_mapped$STRING_id[1:length(res_mapped$gene)]
  
  
  stringdb_enrichment = string_db$get_enrichment(hits_resg, category = path_category)
  stringdb_enrichment$proportion_genes_in_path = stringdb_enrichment$number_of_genes/length(genes)
  return(stringdb_enrichment)
}

btm_col_table = readRDS(paste0(path,"scripts/BTMS/btm_name_to_SUBGROUP_colors.rds"))
source(paste0(btm_path, '/setEnrichment.R'))
load(paste0(btm_path, 'bloodtranscriptionalmodules.RData'))


get_modgenes_table = function(in_dat = modGenes) {
  all_genes_in_path = lapply(in_dat, function(x) as.data.frame(paste(x, collapse = ",")))
  mod_df <- rbindlist(all_genes_in_path)
  mod_df$path_name = names(in_dat)
  colnames(mod_df)[1] <- "all_genes_in_path"
  return(mod_df)
}

modgenes_df = get_modgenes_table(modGenes)

GetBTMs_withgenes = function(gene_list, all_Genes,  input_Genes = modGenes){
  input_Genes = input_Genes[!(grepl("TBA",names(input_Genes)))]
  pathways_sig = setEnrichment(toupper(gene_list), input_Genes, toupper(all_Genes), names(input_Genes))
  return(pathways_sig)
}

generate_heaty_input = function(pos_frame = NULL, neg_frame = NULL, virus = "", path.name = "", path.p.name = "") {
  out_df = NULL
  if (!(is.null(pos_frame) | nrow(pos_frame) == 0)) {
    out_df = pos_frame %>% as.data.frame()
  }
  if (!(is.null(neg_frame) | nrow(neg_frame) == 0)) {
    neg_frame = as.data.frame(neg_frame)
    neg_frame[path.p.name] = -1*neg_frame[path.p.name]
    if (is.null(out_df)){
      out_df = neg_frame
    }
    else {
      out_df = rbind(out_df, neg_frame)
    }
  }
  out_df = as.data.frame(out_df)
  out_df$virus = virus
  out_df <- out_df[order(out_df[path.name], abs(out_df[path.p.name]) ), ]
  out_df = out_df[ !duplicated(out_df[path.name]), ] 
  rownames(out_df) = NULL
  return(out_df)
}

## cleaning pooled results for deg count
threholding_deg_data_2 = function(dat_obj, padj.threshold = padj_thresh, es.threshold = es_thresh, studies.threshold = min_studies) {
  pooled_results_df = dat_obj$metaAnalysis$pooledResults
  pooled_results_df_og = pooled_results_df
  pooled_results_df$direction = ifelse(pooled_results_df$effectSize > 0, "pos", "neg")
  pooled_results_df$gene = rownames(pooled_results_df)
  pooled_results_df = pooled_results_df[pooled_results_df$effectSizeFDR < padj.threshold & abs(pooled_results_df$effectSize) >= es.threshold & pooled_results_df$numStudies >= studies.threshold,]
  print(table(pooled_results_df$direction))
  out_list = list(pooled_results_df, pooled_results_df_og)
  names(out_list) = c("subset","all")
  return(out_list)
}

get_input = function(dat, p.thresh = padj_thresh, es.thresh = es_thresh, dataset.thresh = min_studies) {
  out_dfs = threholding_deg_data_2(dat, p.thresh, es.thresh, dataset.thresh) 
  all_genes = rownames(out_dfs[[2]])
  pos_genes = rownames(out_dfs[[1]][out_dfs[[1]]$direction == "pos",])
  neg_genes = rownames(out_dfs[[1]][out_dfs[[1]]$direction != "pos",])
  out_list = list(pos_genes, neg_genes, all_genes)
  names(out_list) = c("pos_genes", "neg_genes", "all_genes") 
  return(out_list)
}

get_paths_as_list = function(pos_genes = "", neg_genes = "", all_genes = "", list = FALSE) {
  if (list == TRUE) {
    all_genes = pos_genes$all_genes
    neg_genes = pos_genes$neg_genes
    pos_genes = pos_genes$pos_genes
  }
  btm_pos = GetBTMs_withgenes(pos_genes, all_genes)@enrichmentTable
  btm_neg = GetBTMs_withgenes(neg_genes, all_genes)@enrichmentTable
  go_pos = runGOPathAnalysis(pos_genes, path_category = "Process", all_genes)
  go_neg = runGOPathAnalysis(neg_genes, path_category = "Process", all_genes)
  kegg_pos = runGOPathAnalysis(pos_genes, path_category = "KEGG", all_genes)
  kegg_neg = runGOPathAnalysis(neg_genes, path_category = "KEGG", all_genes)
  out_list = list(btm_pos, btm_neg, go_pos, go_neg, kegg_pos, kegg_neg)
  names(out_list) = c("btm_pos", "btm_neg", "go_pos", "go_neg", 'kegg_pos', 'kegg_neg')
  return(out_list)
}


### adding gene score
add_score_to_btm = function(btm_path_df, dat_obj, modgenes_table = modgenes_df, fdr.thresh = .1, es.thresh = 0){
  fdr.thresh = as.numeric(fdr.thresh)
  es.thresh = as.numeric(es.thresh)
  genemat = dat_obj$metaAnalysis$pooledResults
  score_all = c()
  genemat$gene = rownames(genemat)
  btm_path_df$score = 0
  for(pathway_one in btm_path_df$set.name){
    pathgenes <- as.character(modgenes_table$all_genes_in_path[which(modgenes_table$path_name == pathway_one)])
    pathgenes <- unlist(strsplit(pathgenes, ","))
    score = getGeneScores_version2(df = genemat, genes_in_path = pathgenes, fdr_thresh = fdr.thresh, es_thresh = es.thresh)
    score_all = c(score_all, score)
  }
  btm_path_df$score = score_all   
  return(btm_path_df)
}

getGeneScores_version2 = function(df, genes_in_path, fdr_thresh, es_thresh) {
  df = df[df$gene %in% genes_in_path,]
  dim(df)
  df = df[abs(df$effectSize) >= es_thresh & df$effectSizeFDR < fdr_thresh,]
  pos_frame = df[df$effectSize > 0,]
  neg_frame = df[df$effectSize < 0,]
  
  if(nrow(pos_frame) == 0 && nrow(neg_frame) == 0) {
    return(NA)
  } 
  if (nrow(pos_frame) > 0) {
    if (nrow(neg_frame) > 0) {
      return(geomMean(pos_frame$effectSize) - geomMean(-1*(neg_frame$effectSize)))
    }
    else {
      return(geomMean(pos_frame$effectSize))
    }
  }
  return(-geomMean(-1*(neg_frame$effectSize)))
}


geomMean <- function(x, na.rm = T){
  if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(as.numeric(NA))
  }
  if (na.rm)
    x <- x[!is.na(x)]
  if (any(x < 0))
    stop("'x' contains negative value(s)")
  ### direct method causes overflow errors, use log method instead
  ### return(prod(x)^(1/length(x)))
  return(exp(sum(log(x[x>0]))/length(x)))
}

# df = genemat
# genes_in_path = genes_included
getGeneScores_version2 = function(df, genes_in_path, fdr_thresh, es_thresh) {
  df = df[df$gene %in% genes_in_path,]
  df = df[abs(df$effectSize) >= es_thresh & df$effectSizeFDR < fdr_thresh,]
  pos_frame = df[df$effectSize > 0,]
  neg_frame = df[df$effectSize < 0,]
  
  if(nrow(pos_frame) == 0 && nrow(neg_frame) == 0) {
    return(NA)
  } 
  if (nrow(pos_frame) > 0) {
    if (nrow(neg_frame) > 0) {
      return(geomMean(pos_frame$effectSize) - geomMean(-1*(neg_frame$effectSize)))
    }
    else {
      return(geomMean(pos_frame$effectSize))
    }
  }
  return(-geomMean(-1*(neg_frame$effectSize)))
}