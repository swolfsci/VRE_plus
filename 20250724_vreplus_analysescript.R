library(tidyverse)
library(survival)
library(survminer)


vre_plus_cohort <- readxl::read_excel("VRE_plus/VREplus_Basisdoku.xlsx", sheet = "Basisdoku_VRE")
vre_nonplus_cohort <- readxl::read_excel("VRE_plus/VREplus_Basisdoku.xlsx", sheet = "Vergleichskohorte")
datum_erstbefunde <- readxl::read_excel("Kopie von ECFCR.xlsx", sheet = "ECFCR")

vre_plus_cohort <- vre_plus_cohort %>% 
  janitor::clean_names()

vre_nonplus_cohort <- vre_nonplus_cohort %>% 
  janitor::clean_names()

vre_plus_cohort <- vre_plus_cohort %>% 
  mutate(cohort = "vre_plus") %>% 
  dplyr::select(geb_datum, patientennummer, vre_type, cohort, datum_ed, letztes_fu_datum, tod, al, aml)

full_cohort_al <- vre_nonplus_cohort %>% 
  mutate(cohort = "vre_nonplus",
         pat_patientennr = as.numeric(pat_patientennr)) %>% 
  dplyr::rename("patientennummer" = "pat_patientennr", 
                "tod" = "tod_ja_1_nein_2") %>% 
  dplyr::select(geb_datum, patientennummer, cohort, datum_ed, letztes_fu_datum, tod, al, aml) %>% 
  bind_rows(vre_plus_cohort) %>% 
  filter(al == T) 

# calculate age at first diagnosis
full_cohort_al <- full_cohort_al %>% 
  mutate(alter_bei_ed = interval(geb_datum, datum_ed) %/% years(1))

# append first VRE positivity dates
full_cohort_al <- full_cohort_al %>% 
  left_join(datum_erstbefunde, by=c("patientennummer" = "PAT_PATIENTENNR"))

# calculate survtime at interval from first diagnosis and infection_time as interval from first VRE positivity
full_cohort_al <- full_cohort_al %>% 
  mutate(survtime = interval(datum_ed, letztes_fu_datum) %/% weeks(1),
         infection_time = interval(ORD_DATUM, letztes_fu_datum) %/% weeks(1),
         tod = ifelse(tod == 1, 1, 0)) 

# next we calculate both intervals censored at 5 years (260 weeks) 
full_cohort_al <- full_cohort_al %>% 
  mutate(survtime_censored_five_years = case_when(survtime <= 260 ~ survtime,
                                                  TRUE ~ 260),
         infection_time_censored_five_years = case_when(infection_time <= 260 ~ infection_time,
                                                        TRUE ~ 260),
         tod_survtime_censored_five_years = case_when(survtime <= 260 ~ tod,
                                                      TRUE ~ 0),
         tod_infection_time_censored_five_years = case_when(infection_time <= 260 ~ tod,
                                                            TRUE ~ 0)) 

# and mutate VRE status as either having VRE with Dapto, Line or Teico resistance (VRE_plus_mit_teico) or pan-sensible VRE (VRE_nonplus)
full_cohort_al <- full_cohort_al %>% 
  mutate(vre_type = ifelse(is.na(vre_type) & patientennummer %in% vre_mit_teico, "vre_teico", vre_type),
         cohort = ifelse(vre_type == "vre_teico" & !is.na(vre_type), "vre_plus", cohort))

full_cohort_al %>% 
  mutate(aml = ifelse(is.na(aml), F, T)) %>% 
  mutate(vre_type = case_when(AB_Teicoplani == "R" ~ "VRE_teico",
                              TRUE ~ vre_type)) %>% 
  filter(aml == T & alter_bei_ed <70 & !vre_type %in% c("VRE_tige", "VRE_line", "VRE_dapto")) %>% 
  ggsurvplot(survfit(Surv(infection_time_censored_five_years, tod_infection_time_censored_five_years) ~ cohort,.),., pval =  T, risk.table = T, rho = 1)

full_cohort_al %>% 
  mutate(aml = ifelse(is.na(aml), F, T)) %>% 
  mutate(vre_type = case_when(AB_Teicoplani == "R" ~ "VRE_teico",
                              TRUE ~ vre_type)) %>% 
  filter(aml == T & alter_bei_ed <70 & !vre_type %in% c("VRE_tige", "VRE_line", "VRE_dapto")) %>% 
  survdiff(Surv(survtime_censored_five_years, tod_survtime_censored_five_years) ~ cohort,., rho = 1)


vre_mit_teico <- c(8645349L, 8714580L, 9023694L, 9089612L, 9867328L, 9525091L, 8374123L, 9882419L, 9517393L, 8696111L, 8637537L, 9840940L, 9860344L, 9835260L, 8065979L, 8985597L, 9546485L, 9326984L, 9810514L, 9385849L, 8461522L, 9911622L, 9959929L, 8843941L, 9802823L, 8193859L, 8949394L, 8326142L, 9823628L, 9902893L, 8518478L, 9158153L, 9858923L, 9874864L, 8525029L, 9501983L, 9084476L, 9949916L, 9799394L, 9876010L, 9043366L, 9794937L, 9349730L, 9862150L, 9328442L, 9843499L, 8801980L, 9073788L, 9910724L, 9512562L, 8756298L, 9883400L, 8021478L, 9863938L, 9895981L, 8857250L, 9903405L, 9614892L, 8579193L, 9077246L, 9073339L, 9793106L, 9846084L, 8747138L, 8775224L, 8709373L, 9625527L, 9463213L, 9934127L, 8512984L, 9894659L, 9901216L, 9276852L, 9361491L, 8214406L, 9096539L, 9870877L, 9749882L, 9869077L, 8195979L, 9954704L, 8476946L, 9912323L, 9833588L, 9678486L, 9840773L, 9765727L, 9953395L, 9094865L, 8709241L, 9075045L, 9883187L, 9061082L, 8784306L, 8291474L, 8543055L, 9860618L, 8394746L, 9455174L, 9150731L, 8822816L, 9843721L, 9891194L, 9827945L, 8864457L, 9868049L, 9340968L, 9050556L, 9864594L, 9649460L, 9940042L, 9298309L, 8727286L, 8342098L, 8932340L, 8261955L, 9396892L, 9494891L, 9831832L, 9915543L, 9864473L, 8746365L, 9052406L, 9829327L, 9551463L, 9959134L, 9028979L, 9295303L, 8752704L, 9867137L, 9074886L, 9219798L, 8680276L, 9836599L, 8921007L, 9852474L, 8260924L, 9609213L, 8171671L, 9932197L, 9743476L, 9010145L, 9885001L, 9817418L, 9525938L, 8662740L, 9930940L, 9896199L, 8555302L, 9908651L, 9818093L, 9879884L, 9560994L, 9881053L, 3822680L, 9058747L, 9907715L, 9656842L, 8280020L, 8574758L, 9883442L, 8437958L, 8301048L, 8789503L, 8237866L, 8789262L, 9534802L, 9879900L, 8974624L, 9360462L, 9760811L)

