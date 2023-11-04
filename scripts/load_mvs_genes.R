library(data.table)
library(NatParksPalettes)

viral_family_cols = c(natparks.pals("DeathValley")[1],natparks.pals("DeathValley")[2],natparks.pals("DeathValley")[3],natparks.pals("DeathValley")[5],natparks.pals("DeathValley")[7])
viral_families = c("Arenaviridae", "Coronaviridae", "Filoviridae", "Flaviviridae", "Orthomyxoviridae")
viral_family_colors = viral_family_cols
names(viral_family_colors) = viral_families

disease_cols = c(natparks.pals("Acadia")[2], natparks.pals("Acadia")[9])
disease_labels = c("Healthy", "Disease")
disease_names_cols = disease_cols
names(disease_names_cols) = disease_labels


infection_status_cols = c(natparks.pals("Acadia")[2], natparks.pals("Acadia")[9])
infection_status_labels =  c("Healthy", "Infected")
infection_status_label_cols = infection_status_cols
names(infection_status_label_cols) = infection_status_labels


mvs = fread(paste0(paste0(here::here(),"/") ,"data/mvs.csv"))
mvsDirections = mvs[, updown := ifelse(es>0, 'up','down')]
mvs.up = mvs[es>0,gene]
mvs.down = mvs[es<0,gene]
mvs = mvsDirections$gene

lim_mvs = c("ACSL1", "ADM", "ANXA2", "ANXA3", "AQP9", "ARHGAP45", "ATG3", 
            "ATP8B4", "AZU1", "BANF1", "BCAT1", "BCL2A1", "BCL2L11", "BCL6", 
            "BTBD7", "BUB3", "CAMP", "CASP7", "CCL2", "CCR7", "CDT1", "CEACAM8", 
            "CEP55", "CHMP7", "CREG1", "CTSG", "DDB1", "DEFA4", "DOK2", "DOK3", 
            "ELL2", "EPHX2", "EXOC2", "EZH1", "FAM8A1", "FBLN5", "FURIN", 
            "GRN", "H1-0", "HLA-DPB1", "HMMR", "IFITM1", "IFITM2", "IFITM3", 
            "IGFBP2", "IL7R", "ITGB7", "KIF15", "KIF23", "KLHL2", "KLRB1", 
            "KLRD1", "KLRG1", "LAPTM4A", "LCN2", "LRBA", "LTBP3", "MAFB", 
            "MAP3K4", "NAPA", "NQO2", "NUCB1", "OASL", "OLR1", "ORM1", "PFKFB4", 
            "PIK3R1", "PITPNC1", "POMP", "PRC1", "PRF1", "PRSS23", "RAD23B", 
            "RBM15B", "SCAND1", "SIDT1", "SLPI", "SMYD2", "SOCS6", "SRGN", 
            "SSR2", "TCEAL9", "TLN1", "TLR2", "TMEM123", "TRAF3IP3", "TRAF5", 
            "TRIB2", "TRIM28", "TRIP13", "TYK2", "UBE2L6", "USP11", 
            "VAMP5", "VRK2")

Mod1 = c("TXN", "ORM1","PFKFB4", "SLPI","BCL6", "NQO2", "AQP9", "ANXA3","DOK2","KLHL2", "TYK2", "TLN1", "ACSL1","SRGN", "GRN","ADM",  "NUCB1", "TLR2",  "BCL2A1")
Mod2 = c("CAMP","LCN2", "DEFA4", "CTSG","CEACAM8", "BCAT1","BTBD7", "KIF15","AZU1","CEP55", "PRC1","HMMR",  "BCL2L11", "CDT1", "TCEAL9", "OLR1","TRIP13", "ELL2", "SOCS6", "IGFBP2","ATP8B4", "KIF23")
Mod3 = c( "UBE2L6","CASP7","OASL",  "TMEM123", "VRK2","NAPA",  "CCL2", "MAFB","VAMP5","ATG3", "FAM8A1","LAPTM4A", "ANXA2","SSR2","IFITM3", "POMP", "IFITM1", "CREG1","SCAND1","RAD23B",  "H1-0", "IFITM2","FURIN")
Mod4 = c("HLA-DPB1","SMYD2", "SIDT1","TRIB2", "DOK2","KLRB1", "EXOC2", "BUB3", "KLRG1","PRSS23", "KLRD1","PRF1","USP11",  "BANF1","CHMP7","RBM15B",  "MAP3K4",  "ITGB7",  "EPHX2", "IL7R",  "DDB1",  "LRBA",  "TRAF5", "ARHGAP45", "CCR7", "PITPNC1",  "TRAF3IP3","LTBP3","FBLN5",  "TRIM28", "PIK3R1",   "EZH1")

MVSup = c(Mod1,Mod2)
MVSdown = c(Mod3, Mod4)