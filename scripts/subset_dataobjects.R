## has the datasets ID'd by group to subset the actual datasets

## by virus
arena = c("GSE41752",
          "GSE49838",
          "GSE49838",
          "GSE5790",
          "Baize")
corona = c("GSE44542",
           "GSE156701",
           "GSE155363",
           "Sedar_1",
           "Sedar_2",
           'GSE184949',
           "GSE184949",
           "GSE184949_GPL29319",
           "GSE184949_GPL29319")
filo = c("GSE58287",
         "GSE103825",
         "GSE24943",
         "PRJNA718880",
         "PRJNA398558",
         "GSE99463")
flavi = c("GSE185797",
          "GSE185797",
          "GSE90868",
          "GSE72430",
          "GSE185797_KFDV",
          "GSE185797_ALKV")
ortho = c("GSE152406",
          "GSE60009")
pox = c("GSE4013_GPL3093",
        "GSE4013_GPL3346",
        "GSE4013_GPL3347")

## by virus family


##subject object by dataset IDs
reformat = function(object, dataset_list, reMetanalaysis = FALSE) {
  for (name in (names(object$originalData))) {
    if (!(name %in% dataset_list)) {
      object$originalData[[name]] = NULL
    }
  }
  if (reMetanalaysis){
    object = fix_classes(object)
    object = runMetaAnalysis(object)
    object = filterGenes(object, isLeaveOneOut = TRUE, FDRThresh = 0.001)
  }
  return(object)
}


## remove empty datasets from subset
remove_empties = function(object, time) {
  for (name in (names(object$originalData))) {
    if (sum(object$originalData[[name]]$pheno$time_cat == time) <2) {
      object$originalData[[name]] = NULL
    }
  }
  return(object)
}

## pick formatted name metadata (check if column has one unique value)
fix_formattedname = function(object, metadata_col = "obj_name") {
  if (length(unique(object$originalData[[1]]$pheno[,metadata_col])) > 1) {
    return(object)
  }
  for (name in (names(object$originalData))) {
    object$originalData[[name]]$formattedName = unique(object$originalData[[name]]$pheno[,metadata_col])
  }
  return(object)
}



get_virus_name = function(object) {
  for (name in (names(object$originalData))) {
    object$originalData[[name]]$formattedName = object$originalData[[name]]$name
  }
  return(object)
}


