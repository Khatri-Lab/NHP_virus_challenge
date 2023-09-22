## scripts for metanalysis/forrestplots
get_meta_dat = function(meta_df, timepoint = "", virus = "") {
  tp = timepoint
  vrs = virus
  dat_num = meta_df$k
  case_num = meta_df$n.e.pooled
  control_num = meta_df$n.c.pooled
  Q = meta_df$Q %>% as.numeric()%>% round(1)
  Q_p = meta_df$pval.Q %>% as.numeric()
  I = meta_df$I2 %>% as.numeric()
  I_low =meta_df$lower.I2 %>% as.numeric()
  I_up =meta_df$upper.I2 %>% as.numeric()
  smd_common = meta_df$TE.common  %>% as.numeric()# SMD
  smd_common_se = meta_df$seTE.common  %>% as.numeric()
  smd_common_p = meta_df$pval.common %>% as.numeric()
  lower_common = meta_df$lower.common %>% as.numeric()
  upper_common = meta_df$upper.common %>% as.numeric()
  smd_random = meta_df$TE.random  %>% as.numeric()# SMD
  smd_random_se = meta_df$seTE.random  %>% as.numeric()
  smd_random_p = meta_df$pval.random %>% as.numeric()
  lower_random = meta_df$lower.random %>% as.numeric()
  upper_random = meta_df$upper.random %>% as.numeric()
  out_list = c(tp, vrs, dat_num, case_num, control_num, Q, Q_p, I, I_low, I_up, smd_common, smd_common_se, smd_common_p, lower_common, upper_common, smd_random, smd_random_se, smd_random_p, lower_random, upper_random)
  names(out_list) = c("timepoint","virus","n_datasets","n_cases","n_controls","Q","Q_p","I","I_low","I_up",
                      "smd_common", "smd_common_se", "smd_common_p","common_lower","common_upper", "smd_random", "smd_random_se","smd_random_p", "random_lower","random_upper")
  return(out_list)
}

make_epi_out = function(dataframe, virus, timepoint) {
  dataframe = dataframe[dataframe$dataset5 == virus,]
  dataframe$animal_day =paste0(dataframe$animal,"_",dataframe$day)
  
  df.comp = dataframe[dataframe$time_cat == timepoint,]
  df.comp = df.comp[,c("dataset","MVS_score","animal_day")]
  df.comp = df.comp %>%
    group_by(dataset) %>% 
    summarise(
      n = n(),
      mean = mean(MVS_score,na.rm = T),
      min=min(MVS_score, na.rm=TRUE),
      max=max(MVS_score,na.rm=TRUE), 
      sd= sd(MVS_score, na.rm=TRUE))
  
  df.cont = dataframe[dataframe$time_cat == "t0",]
  df.cont = df.cont[,c("dataset","MVS_score","animal_day")]
  df.cont = df.cont %>%
    group_by(dataset) %>% 
    summarise(
      n.cont = n(),
      mean.cont = mean(MVS_score,na.rm = T),
      min.cont=min(MVS_score, na.rm=TRUE),
      max.cont=max(MVS_score,na.rm=TRUE), 
      sd.cont= sd(MVS_score, na.rm=TRUE))
  
  dataframe = merge(df.comp, df.cont, by = "dataset", all.x = TRUE, all.y = FALSE)
  return(dataframe)
}

make_epi_out2 = function(dataframe, timepoint) {
  dataframe$animal_day =paste0(dataframe$animal,"_",dataframe$day)
  
  df.comp = dataframe[dataframe$time_cat == timepoint,]
  df.comp = df.comp[,c("dataset","MVS_score","animal_day")]
  df.comp = df.comp %>%
    group_by(dataset) %>% 
    summarise(
      n = n(),
      mean = mean(MVS_score,na.rm = T),
      min=min(MVS_score, na.rm=TRUE),
      max=max(MVS_score,na.rm=TRUE), 
      sd= sd(MVS_score, na.rm=TRUE))
  
  df.cont = dataframe[dataframe$time_cat == "t0",]
  df.cont = df.cont[,c("dataset","MVS_score","animal_day")]
  df.cont = df.cont %>%
    group_by(dataset) %>% 
    summarise(
      n.cont = n(),
      mean.cont = mean(MVS_score,na.rm = T),
      min.cont=min(MVS_score, na.rm=TRUE),
      max.cont=max(MVS_score,na.rm=TRUE), 
      sd.cont= sd(MVS_score, na.rm=TRUE))
  
  dataframe = merge(df.comp, df.cont, by = "dataset", all.x = TRUE, all.y = FALSE)
  return(dataframe)
}

out_meta = function(dataframe) {
  res.out =  metacont(n, mean, sd, 
                      n.cont, mean.cont, sd.cont,
                      comb.fixed = T, comb.random = T, studlab = dataset,
                      data = dataframe, sm = "SMD") 
  return(res.out)
}

meta_analysis_out = function(input_df, virus, timepoint) {
  df = make_epi_out(input_df, virus, timepoint)
  df = out_meta(df)
  output = get_meta_dat(df, timepoint, virus)
  return(output)
}

get_meta_dat = function(meta_df, timepoint = "", virus = "") {
  tp = timepoint
  vrs = virus
  dat_num = meta_df$k
  case_num = meta_df$n.e.pooled
  control_num = meta_df$n.c.pooled
  Q = meta_df$Q %>% as.numeric()%>% round(1)
  Q_p = meta_df$pval.Q %>% as.numeric()
  I = meta_df$I2 %>% as.numeric()
  I_low =meta_df$lower.I2 %>% as.numeric()
  I_up =meta_df$upper.I2 %>% as.numeric()
  smd_common = meta_df$TE.common  %>% as.numeric()# SMD
  smd_common_se = meta_df$seTE.common  %>% as.numeric()
  smd_common_p = meta_df$pval.common %>% as.numeric()
  lower_common = meta_df$lower.common %>% as.numeric()
  upper_common = meta_df$upper.common %>% as.numeric()
  smd_random = meta_df$TE.random  %>% as.numeric()# SMD
  smd_random_se = meta_df$seTE.random  %>% as.numeric()
  smd_random_p = meta_df$pval.random %>% as.numeric()
  lower_random = meta_df$lower.random %>% as.numeric()
  upper_random = meta_df$upper.random %>% as.numeric()
  out_list = c(tp, vrs, dat_num, case_num, control_num, Q, Q_p, I, I_low, I_up, smd_common, smd_common_se, smd_common_p, lower_common, upper_common, smd_random, smd_random_se, smd_random_p, lower_random, upper_random)
  names(out_list) = c("timepoint","virus","n_datasets","n_cases","n_controls","Q","Q_p","I","I_low","I_up",
                      "smd_common", "smd_common_se", "smd_common_p","common_lower","common_upper", "smd_random", "smd_random_se","smd_random_p", "random_lower","random_upper")
  return(out_list)
}

prep_for_table = function(data_frame, time_of_interest="t1", days="") {
  out_df = data_frame[data_frame$timepoint == time_of_interest,]
  out_df$` ` = "                                       "
  colnames(out_df)[2:5] = c("Virus","# Datasets","# Cases","# Controls")
  out_df$SMD = round(as.numeric(out_df$smd_common), 2)
  out_df$low = round(as.numeric(out_df$common_lower), 2)
  out_df$high = round(as.numeric(out_df$common_upper), 2)
  out_df$Q_p = round(as.numeric(out_df$Q_p), 2)
  out_df$`SMD (95% CI)` = paste0(out_df$SMD," [",out_df$low,", ",out_df$high,"]")
  tm <- forest_theme(core=list(fg_params=list(hjust = 1, x = 0.9)),
                     colhead=list(fg_params=list(hjust=0.5, x=0.5)))
  out_df$se = as.numeric(out_df$smd_common_se)
  out_df$smd_common_p = as.numeric(out_df$smd_common_p)
  out_df$p = ifelse(out_df$smd_common_p < 0.05, "*","") 
  out_df$p = ifelse(out_df$smd_common_p < 0.01, "**",out_df$p) 
  out_df$p = ifelse(out_df$smd_common_p < 0.001, "***",out_df$p)
  out_df$title = paste0(toupper(time_of_interest),days,": Summary Stats")
  return(out_df)
}