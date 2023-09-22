library(MetaIntegrator)
library(COCONUT)

bind_pheno_tables = function(object, columns_of_interest) {
  for (i in 1:length(names(object$originalData))) {
    if (i == 1) {
      # print(names(object$originalData)[i])
      # print(columns_of_interest[!(columns_of_interest %in% colnames(object$originalData[[i]]$pheno))])
      object$originalData[[i]]$pheno = as.data.frame(object$originalData[[i]]$pheno)
      out_df = object$originalData[[i]]$pheno[, columns_of_interest]
    }
    else {
      # print(names(object$originalData)[i])
      # print(columns_of_interest[!(columns_of_interest %in% object$originalData[[i]]$pheno)])
      object$originalData[[i]]$pheno = as.data.frame(object$originalData[[i]]$pheno)
      out_df = rbind(out_df, object$originalData[[i]]$pheno[, columns_of_interest])
    }
  }
  out_df$day = as.numeric(as.character(out_df$day))
  out_df$disease = factor(out_df$disease, levels = c("Healthy", "Disease"))
  return(out_df)
}

fix_classes = function(object) {
  for (dat_names in names(object$originalData)) {
    obj = object$originalData[[dat_names]]
    obj$class = ifelse(obj$pheno$disease == "Healthy",0,1)
    names(obj$class) = rownames(obj$pheno)
    object$originalData[[dat_names]] = obj
  }
  return(object)
}

add_score_to_pheno = function(object, col_name = "score", genes_up = "", genes_down = "") {
  if (!is.null(genes_up)) {
    object$filterResults$FDR0.001_es0_nStudies1_looaTRUE_hetero0$posGeneNames = genes_up
  }
  if (!is.null(genes_down)) {
    object$filterResults$FDR0.001_es0_nStudies1_looaTRUE_hetero0$negGeneNames = genes_down
  }
  scoring = object$filterResults$FDR0.001_es0_nStudies1_looaTRUE_hetero0
  for (item in names(object$originalData)) {
    object$originalData[[item]]$pheno[col_name] <- calculateScore(scoring, object$originalData[[item]])
  }
  return(object)
}


fix_time_cat = function(object, breaks = c(1,3,6,9,14), col_name = NA) {
  for (name in (names(object$originalData))) {
    object$originalData[[name]]$pheno = add_time_cat_hardcode(object$originalData[[name]]$pheno, breaks)
  }
  return(object)
}

add_time_cat_hardcode = function(table, breaks = c(1,3,6,9,14)) {
  table = as.data.frame(table)
  table$day = as.numeric(as.character(table$day))
  for (i in 1:length(breaks)) {
    if (i == 1) {
      table$time_cat = ifelse(table$day >= as.numeric(breaks[i]), paste0("t",i), "t0")
    }
    else {
      table$time_cat = ifelse(table$day >= as.numeric(breaks[i]), paste0("t",i), table$time_cat)
    }
  }
  return(table)
}

add_time_cat = function(table, breaks = c(1,3,6,9,14), col_name) {
  table = as.data.frame(table)
  table$day = as.numeric(as.character(table$day))
  for (i in 1:length(breaks)) {
    if (i == 1) {
      table[col_name] = ifelse(table$day >= breaks[i], paste0("t",i), "t0")
    }
    else {
      table[col_name] = ifelse(table$day >= breaks[i], paste0("t",i), table[,col_name])
    }
  }
  return(table)
}

subset_data_with_dfkey = function(object, df, colname = "obj_name") {
  data_out = list()
  object = fix_time_cat(object)
  for (i in 1:length(names(object$originalData))) {
    dataset = object$originalData[[i]]
    pheno = dataset$pheno
    name = unname(unlist(unique(pheno[colname])))
    time_cat_max = (df[df[colname] == name,]$time_category)
    ## subset pheno by time_cat
    pheno_sub = pheno[pheno$time_cat %in% c("t0", time_cat_max),]
    print(paste0(name,"_",nrow(pheno),"_",nrow(pheno_sub)))
    ## subset by dataset by pheno
    dataset <- subsetOriginalData(dataset, 
                                  keepMe= c(rownames(pheno_sub)))
    # checkDataObject(dataset, "Dataset") 
    object$originalData[[i]] = dataset
    }
  return(object)
}



subset_by_timecat = function(object, time_val) {
  data_out = list()
  for (i in 1:length(names(object$originalData))) {
    dataset = object$originalData[[i]]
    pheno = dataset$pheno
    name = unname(unlist(unique(pheno$dataset)))
    ## subset pheno by time_cat
    pheno_sub = pheno[pheno$time_cat %in% c("t0", time_val),]
    print(paste0(name,"_",nrow(pheno),"_",nrow(pheno_sub)))
    ## subset by dataset by pheno
    dataset <- subsetOriginalData(dataset, 
                                  keepMe= c(rownames(pheno_sub)))
    # checkDataObject(dataset, "Dataset") 
    object$originalData[[i]] = dataset
  }
  return(object)
}

baseline_sub_jobs_MVS_score = function(object) {
  for (dat_names in names(object$originalData)) {
    print(dat_names)
    pheno_df = object$originalData[[dat_names]]$pheno
    pheno_df$animal_id = paste0(pheno_df$animal,"_",pheno_df$obj_name)
    pheno_df$rowname = rownames(pheno_df)
    pheno_df = setDT(pheno_df)[, scores_normed := MVS_score - mean(MVS_score[day == 0]), by = animal_id]
    rownames(pheno_df) = pheno_df$rowname
    object$originalData[[dat_names]]$pheno = pheno_df
  }
  return(object)
}

baseline_sub_jobs_mvs_score = function(object) {
  for (dat_names in names(object$originalData)) {
    pheno_df = object$originalData[[dat_names]]$pheno
    pheno_df$animal_id = paste0(pheno_df$animal,"_",pheno_df$obj_name)
    pheno_df$rowname = rownames(pheno_df)
    pheno_df = setDT(pheno_df)[, mvs_scores_normed := mvs_score - mean(mvs_score[day == 0]), by = animal_id]
    pheno_df = as.data.frame(pheno_df)
    rownames(pheno_df) = pheno_df$rowname
    object$originalData[[dat_names]]$pheno = pheno_df
  }
  return(object)
}



#5
get_datasets_subset = function(object, colname, item) {
  data_out = list()
  i = 1
  for (dat_name in names(object$originalData)) {
    dataset = object$originalData[[dat_name]]
    if (unique(dataset$pheno[,colname]) %in% item) {
      print(dat_name)
      data_out$originalData[[i]] = dataset
      names(data_out$originalData)[[i]] = dat_name
      i = i + 1
    }
  }
  return(data_out)
}



coconut_subset = function(object, column, item, acute = TRUE, formatted_name = NA, version = "a") {
  object_out = get_datasets_subset(object, column, item)
  object_out = coconutMetaIntegrator(object_out)
  object_out = combineCOCOoutput(object_out)
  object_out$class = object_out$class.cntl0.dis1
  if (acute) {
    if (version == "a") {
      to_keep = readRDS(paste0(datPath, "times_to_keep_a.rds"))[,c("obj_name","time_category")]
    }
    if (version == "b") {
      to_keep = readRDS(paste0(datPath, "times_to_keep_b.rds"))[,c("obj_name","time_category")]
    }
    to_keep$time_category = ifelse(grepl("EBOV_2", to_keep$obj_name), "t3", to_keep$time_category) ## why
    colnames(to_keep)[2] = "peak_time"
    object_out$pheno = merge(object_out$pheno, to_keep, by = "obj_name")
    object_out$pheno$acute_time = ifelse(object_out$pheno$time_category == object_out$pheno$peak_time, "yes", "no")
    object_out$pheno$acute_time = ifelse(object_out$pheno$time_category == "t0", "yes", object_out$pheno$acute_time)
    # object_out = get_datasets_subset(object_out, "acute_time", "yes")
    object_out = removeSamples(object_out,which(!object_out$pheno$acute_time %in% c("yes")),expr=F)
  }
  if (is.na(formatted_name)) {
    object_out$formattedName = item
  }
  else {
    object_out$formattedName = formatted_name
  }
  return(object_out)
}

coconut_subset_day = function(object, column, item, acute = TRUE, formatted_name = NA, version = "a") {
  object_out = get_datasets_subset(object, column, item)
  object_out = coconutMetaIntegrator(object_out)
  object_out = combineCOCOoutput(object_out)
  object_out$class = object_out$class.cntl0.dis1
  if (acute) {
    if (version == "a") {
      to_keep = readRDS(paste0(datPath, "times_to_keep_a.rds"))[,c("obj_name","time_category")]
    }
    if (version == "b") {
      to_keep = readRDS(paste0(datPath, "times_to_keep_b.rds"))[,c("obj_name","time_category")]
    }
    to_keep$time_category = ifelse(grepl("EBOV_2", to_keep$obj_name), "t3", to_keep$time_category) ## why
    colnames(to_keep)[2] = "peak_time"
    object_out$pheno = merge(object_out$pheno, to_keep, by = "obj_name")
    object_out$pheno$acute_time = ifelse(object_out$pheno$day == object_out$pheno$peak_time, "yes", "no")
    object_out$pheno$acute_time = ifelse(object_out$pheno$day == 0, "yes", object_out$pheno$acute_time)
    # object_out = get_datasets_subset(object_out, "acute_time", "yes")
    object_out = removeSamples(object_out,which(!object_out$pheno$acute_time %in% c("yes")),expr=F)
  }
  if (is.na(formatted_name)) {
    object_out$formattedName = item
  }
  else {
    object_out$formattedName = formatted_name
  }
  return(object_out)
}
