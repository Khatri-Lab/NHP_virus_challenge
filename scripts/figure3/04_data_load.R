## load all genes data
arena_all_peak = rownames(readRDS(file = paste0(datPath,"arena","_timesplit_dat.rds"))$t3$metaAnalysis$pooledResults)
corona_all_peak = rownames(readRDS(file = paste0(datPath,"corona","_timesplit_dat.rds"))$t1$metaAnalysis$pooledResults)
filo_all_peak = rownames(readRDS(file = paste0(datPath,"filo","_timesplit_dat.rds"))$t3 $metaAnalysis$pooledResults)
flavi_all_peak = rownames(readRDS(file = paste0(datPath,"flavi","_timesplit_dat.rds"))$t3$metaAnalysis$pooledResults) 
ortho_all_peak = rownames(readRDS(file = paste0(datPath,"ortho","_timesplit_dat.rds"))$t2$metaAnalysis$pooledResults) 
all = readRDS(file = paste0(datPath, "16a_timecat_peaksubset_unbiased.rds"))
all_peak = rownames(all$metaAnalysis$pooledResults)

## load timepoint data for degs
virus_dat_list_t1 = readRDS(paste0(datPath, "15a_viruslist_t1.rds"))
virus_dat_list_t2 = readRDS(paste0(datPath, "15a_viruslist_t2.rds"))
virus_dat_list_t3 = readRDS(paste0(datPath, "15a_viruslist_t3.rds"))
virus_dat_list_t4 = readRDS(paste0(datPath, "15a_viruslist_t4.rds"))
virus_dat_list_t5 = readRDS(paste0(datPath, "15a_viruslist_t5.rds"))

arena_t1_dat = virus_dat_list_t1$arena
arena_t2_dat = virus_dat_list_t2$arena
arena_t3_dat = virus_dat_list_t3$arena
arena_t4_dat = virus_dat_list_t4$arena
arena_t5_dat = virus_dat_list_t5$arena

corona_t1_dat = virus_dat_list_t1$corona
corona_t2_dat = virus_dat_list_t2$corona
corona_t3_dat = virus_dat_list_t3$corona
corona_t4_dat = virus_dat_list_t4$corona
corona_t5_dat = virus_dat_list_t5$corona

filo_t1_dat = virus_dat_list_t1$filo
filo_t2_dat = virus_dat_list_t2$filo
filo_t3_dat = virus_dat_list_t3$filo
filo_t4_dat = virus_dat_list_t4$filo
filo_t5_dat = virus_dat_list_t5$filo

flavi_t1_dat = virus_dat_list_t1$flavi
flavi_t2_dat = virus_dat_list_t2$flavi
flavi_t3_dat = virus_dat_list_t3$flavi
flavi_t4_dat = virus_dat_list_t4$flavi
flavi_t5_dat = virus_dat_list_t5$flavi

ortho_t1_dat = virus_dat_list_t1$ortho
ortho_t2_dat = virus_dat_list_t2$ortho
ortho_t3_dat = virus_dat_list_t3$ortho
ortho_t4_dat = virus_dat_list_t4$ortho
ortho_t5_dat = virus_dat_list_t5$ortho
