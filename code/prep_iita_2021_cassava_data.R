### From IITA 2021 GS

library(tidyverse); library(magrittr);
blups<-readRDS(file=here::here("data/IITA_2021GS",
                               "IITA_blupsForModelTraining_twostage_asreml_2021Aug09.rds"))

dosages<-readRDS(file=here::here("data/IITA_2021GS","dosages_IITA_2021Aug09.rds"))

reduced_snp_set<-read.table(here::here("data/IITA_2021GS/",
                                       "samples2keep_IITA_MAFpt01_prune1Mb_50kb_pt6.prune.in"),
                            header = F, stringsAsFactors = F)$V1

# full list
# all_traits<-c("logDYLD","logFYLD","logRTNO","logTOPYLD","MCMDS","DM","BCHROMO",
#           "PLTHT","BRLVLS","BRNHT1","HI")
traits<-c("MCMDS","DM","BCHROMO","PLTHT")

blups %>%
  select(Trait,blups) %>%
  unnest(blups) %>%
  distinct(GID) %$% GID -> gidWithBLUPs

genotypedWithBLUPs<-gidWithBLUPs[gidWithBLUPs %in% rownames(dosages)]
length(genotypedWithBLUPs)
# [1] 8669

blups %<>%
  filter(Trait %in% traits) %>%
  select(Trait,blups,varcomp) %>%
  mutate(blups=map(blups,~filter(.,GID %in% genotypedWithBLUPs)))
blups %<>%
  select(Trait,blups) %>%
  unnest(blups) %>%
  mutate(Cohort=NA,
         Cohort=ifelse(grepl("TMS20",GID,ignore.case = T),"TMS20",
                       ifelse(grepl("TMS19",GID,ignore.case = T),"TMS19",
                              ifelse(grepl("TMS18",GID,ignore.case = T),"TMS18",
                                     ifelse(grepl("TMS17",GID,ignore.case = T),"TMS17",
                                            ifelse(grepl("TMS16",GID,ignore.case = T),"TMS16",
                                                   ifelse(grepl("TMS15",GID,ignore.case = T),"TMS15",
                                                          ifelse(grepl("TMS14",GID,ignore.case = T),"TMS14",
                                                                 ifelse(grepl("TMS13|2013_",GID,
                                                                              ignore.case = T),"TMS13","GGetc")))))))))

## IITA Germplasm Ages
ggcycletime<-readxl::read_xls(here::here("data/IITA_2021GS","PedigreeGeneticGainCycleTime_aafolabi_01122020.xls")) %>%
  mutate(Year_Accession=as.numeric(Year_Accession))


rawdata<-readRDS(here::here("data/IITA_2021GS","IITA_ExptDesignsDetected_2021Aug08.rds"))
blups %>%
  distinct(GID) %>%
  left_join(rawdata %>%
              distinct(germplasmName,GID)) %>%
  group_by(GID) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(ggcycletime %>%
              rename(germplasmName=Accession) %>%
              mutate(Year_Accession=as.numeric(Year_Accession))) %>%
  mutate(Year_Accession=case_when(grepl("2013_|TMS13",germplasmName)~2013,
                                  grepl("TMS14",germplasmName)~2014,
                                  grepl("TMS15",germplasmName)~2015,
                                  grepl("TMS16",germplasmName)~2016,
                                  grepl("TMS17",germplasmName)~2017,
                                  grepl("TMS18",germplasmName)~2018,
                                  grepl("TMS19",germplasmName)~2019,
                                  grepl("TMS20",germplasmName)~2020,
                                  !grepl("2013_|TMS13|TMS14|TMS15|TMS16|TMS17|TMS18|TMS19|TMS20",germplasmName)~Year_Accession)) %>%
  select(GID,germplasmName,Year_Accession) %>%
  left_join(blups,.) -> blups
rm(rawdata)
blups %<>%
  relocate(c(Cohort,Year_Accession),.before=everything()) %>%
  select(-germplasmName)

saveRDS(blups,file=here::here("data/IITA_2021GS","blups_forGWAS_FullPop.rds"))

# Subset and simplify for teaching

# Just the "Genetic Gain" population, pre-2013
gg_drgblups_wide<-blups %>%
  filter(!is.na(Year_Accession),
         Year_Accession<2013) %>%
  select(Cohort,Year_Accession,Trait,GID,drgBLUP) %>%
  spread(Trait,drgBLUP)

saveRDS(gg_drgblups_wide,file=here::here("data/IITA_2021GS","drgblups_forGWAS_geneticgainpop.rds"))

snps<-dosages[rownames(dosages) %in% blups$GID,]
saveRDS(snps,file=here::here("data/IITA_2021GS","snps_forGWAS_FullPop.rds"))

snpsets<-readRDS(file=here::here("data/IITA_2021GS","snpsets.rds"))

reduced_snp_set<-snpsets %>%
  filter(Set=="reduced_set") %>%
  unnest(snps2keep) %$% FULL_SNP_ID

snps_small<-snps[,colnames(snps) %in% reduced_snp_set]

snps_gg<-snps[rownames(snps) %in% gg_drgblups_wide$GID,]
snps_small_gg<-snps_gg[,colnames(snps_gg) %in% reduced_snp_set]
saveRDS(snps_gg,file=here::here("data/IITA_2021GS","snps_forGWAS_geneticgainpop.rds"))
saveRDS(snps_small_gg,file=here::here("data/IITA_2021GS","snps_forGWAS_geneticgainpop_reducedset.rds"))





