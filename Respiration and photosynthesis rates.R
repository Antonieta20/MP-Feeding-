
## Coral Photosynthesis and Respiration ---------------------------------
# Entire Data analysis related to the photosymbiont response
      

## Upload data tables ---------------------------------------------------

# Read table with infos to IDs, tanks, colonies and treatments

tab_info <- read.csv2("Data/info.csv",sep=";", dec=".",header=T) %>% 


tab_info$treatment <- factor(tab_info$treatment, 
                             levels = c("Control", "MP + HF", "MP"))

# Read scanning data

tab_size <- read.csv2("Data/3D_Scanning_Data_Table.csv",sep=";", dec=".",header=T) %>% 
  na.omit() %>% 
  select(ID, Volume_t0, Surface_t0, Volume_t1, Surface_t1)

# Read oxygen incubation data

tab_incubations <- read.csv("Data/Incubations_Data_Table.csv",sep=";", dec=".",header=T) %>% 
   na.omit() %>% 
 

tab_incubations <-merge(tab_incubations, tab_size, by="ID", all = TRUE)


# Read PAM data

tab_PAM <- read.csv("Data/PAM_Data table.csv",sep=";", dec=".",header=T) %>% 
  rename(ID = Sample, Species = Spi) %>% 
  select(c(-t, -Time, -No, -ML, -ETR, -Tank_info, -X, -X.1)) %>% 
  na.omit()  

names(tab_PAM)

## Data processing for Photosynthesis rates -------------
#Here we calculate net photosynthesis, respiration and gross 
#photosynthesis in µg O2 per cm² per h derrived from oxygen 
#incubations. The net photosynthesis is calculated from the 
#measurements in the light, the respiration is calculated from 
#the measurements in the dark. The gross photosynthesis is the sum 
#of the net photosynthesis plus the amount that is respired. For this, 
#first we calculate the mean oxygen production in the controls. This 
#will be substracted from the incubations with the corals as this is 
#the background production/respiration (depending on if it's positive 
#or negative) of the microorganisms in the water. In the next steps the
#differences in oxygen before and after each incubation will be 
#standardized to the Surface_t1 area of the fragment, the water 
#Volume_t1 in the glass (glass Volume_t1-Volume_t1 of the fragment) 
#and the exact time of incubation.

#First we calculate and extract necessary information 
#from the table to make the data analysis easier

tab_incubations <- tab_incubations %>% 
  mutate(day = as.numeric(as.factor(date)),  #adds column for day as number, converted from the date
         light_diff = (O2_stop_light - O2_start_light)/time_light*60,  #difference in O2 per h
         dark_diff = (O2_stop_dark - O2_start_dark)/time_dark*60,
         waterVolume_t1 = 1005.011 - (Volume_t1 / 1000), #calculate water Volume_t1 in jars
         waterVolume_t1 = replace_na(waterVolume_t1,1005.011), #replace NAs in controls with Volume_t1 of jar
         light_O2_perL = light_diff/waterVolume_t1*1000,
         dark_O2_perL = dark_diff/waterVolume_t1*1000)


#Calculate means from control incubations separated by 
#light and dark incubation    

control_rates_light <-
  tab_incubations %>%
  filter(ID=="Control") %>% 
  group_by(day) %>%
  get_summary_stats(light_O2_perL, type = "mean_sd") %>% 
  rename(con_light_O2_perL = mean)  %>% 
  select(day, con_light_O2_perL)  

control_rates_dark <-
  tab_incubations %>%
  filter(ID=="Control") %>% 
  group_by(day) %>%
  get_summary_stats(dark_O2_perL, type = "mean_sd") %>% 
  rename(con_dark_O2_perL = mean) %>% 
  select(day, con_dark_O2_perL)

# Merge light and dark controls together

control <- merge(control_rates_light, control_rates_dark, by="day", all = TRUE)


# Clean datatable

tab_psrate <- tab_incubations %>% 
  filter(ID!="Control") %>% 
  select(ID, sampling_point, day, position_glass, Surface_t1, light_O2_perL, dark_O2_perL)

# Add controls as row behind data from fragments

tab_psrate <- merge(tab_psrate, control, by="day", all = TRUE)


# Subtract controls from incubations and standardize to 
#Surface_t1 area of the fragment

tab_psrate <- tab_psrate %>% 
  mutate(light_O2_perL_corrected = light_O2_perL - con_light_O2_perL, #subtract controls
         dark_O2_perL_corrected = dark_O2_perL - con_dark_O2_perL, 
         light_O2_perL_percm2 = light_O2_perL_corrected/(Surface_t1/100), #correct for Surface_t1 area
         dark_O2_perL_percm2 = dark_O2_perL_corrected/(Surface_t1/100), 
         net_ps = light_O2_perL_percm2*1000, #convert to µmO2
         respiration = (dark_O2_perL_percm2*1000)*-1,# convert to consumed oxygen to positive values
         gross_ps = net_ps + respiration) %>% 
  select(ID, sampling_point, day, position_glass, net_ps, respiration, gross_ps)

tab_psrate <- merge(tab_info, tab_psrate, by = "ID", all = TRUE)

write.table(tab_psrate, sep = ";", "Output/Tab_photosynthesisrates_processed.csv")


## Data processing for PAM data -----------------------------------------
#Here we calculate the means from all light measurements

tab_light_mean <-
  tab_PAM %>%
  filter(Meas_type=="Light") %>% 
  group_by(ID, Date, Species, Treatment, Tank, Colony) %>%
  get_summary_stats(YII, type = "mean_sd") %>% 
  rename(YII = mean)%>% 
  select(c(-variable,-n,-sd))

str(tab_light_mean)

tab_dark <-
  tab_PAM %>%
  filter(Meas_type=="Dark") %>%
  rename(FvFm = YII)

str(tab_dark)


write.table(tab_light_mean, sep = ";", "Output/Tab_PAMlight_processed.csv")
datatable(tab_light_mean, caption = "Table 3: Effective quantum yield (Y(II)).")  

write.table(tab_dark, sep = ";", "Output/Tab_PAMdark_processed.csv")
datatable(tab_dark, caption = "Table 4: Maximum quantum yield (Fv/Fm).")  

## Data analysis photosynthesis rates  ------------------------------------------------
#Here we analyze photosynthesis rates of the corals. 
#These were recorded performing oxygen incubations at two 
#timepoints of the experiment. We analyze the three parameters 
#net photosynthesis, respiration and gross photosynthesis.

# Descriptive statistics

### Check Outliers


outliers_net_ps <- tab_psrate %>%
  drop_na() %>%
  group_by(treatment, Species) %>%
  identify_outliers(net_ps) %>% 
  filter(is.extreme==TRUE) 

outliers_net_ps
# 1  extreme outlier, removed

outliers_net_ps$ID

outliers_resp <- tab_psrate %>%
  drop_na() %>%
  group_by(treatment, Species) %>%
  identify_outliers(respiration) %>% 
  filter(is.extreme==TRUE) 

outliers_resp
# 2  extreme outliers, removed
outliers_resp$ID

outliers_gross_ps <- tab_psrate %>%
  drop_na() %>%
  group_by(treatment, Species) %>%
  identify_outliers(gross_ps) %>% 
  filter(is.extreme==TRUE) 

outliers_gross_ps

outliers_gross_ps$ID


tab_psrate <- tab_psrate %>%
  filter(ID != c("B-T1-2", "A-T9-2", "B-T1", "B-T5"))
### All removed from the biginin (during the upload tables)

## Overview mean and SD values -----------------------------------------

skim_1 <- tab_psrate %>% 
  skim()

summy_1 <- tab_psrate %>%
  drop_na() %>%
  filter(sampling_point=="T1") %>% 
  group_by(Species, treatment) %>%
  get_summary_stats(net_ps, respiration, gross_ps) 


#skim_1
datatable(summy_1, caption = "Table 5: 
          Descriptive summary statistics of photosynthesis 
          rates (µg O2/cm²/h). ")

##Values of the means of the tab_psrate ---------------------------------

#Means values of the tab_psrate

net_ps_means <- tab_psrate %>% 
  
  drop_na() %>%
  
  filter(sampling_point=="T1") %>% 
  
  group_by(treatment, Species) %>% summarize(net_ps_mean = 
    mean(net_ps), respiration_mean = mean(respiration),
    gross_ps_mean = mean(gross_ps))

## Test null hypotheses -------------------------------------------------
#Now we want to test if our treatments have a signifiant influence. 
#We use linear mixed effects models that include the coral colony 
#as random factor. We construct the model and then test the 
#distribution of the model residuals for normality and 
#heteroscedasticity. We had to log transform the data to meet the 
#model criteria. This might change if the outliers are corrected 
#or removed, so go back here and check again if this is correct. 
#We find significant differences (p.adjusted <0.05) for both 
#species.


# Null hypothesis -T1 -----------------------------------------------------
### All analysis of T1 (at the end of 6 weeks)

##P. verrucosa ------------------------------------------------------------

#Results for P. verrucosa, net photosynthesis 

model_net_ps_Pve <- tab_psrate %>% 
  filter(Species=="Pve", sampling_point=="T1") %>% 
  lmer(scale(net_ps)~treatment + (1|colony), data=.)

#qqPlot(residuals(model_net_ps_Pve))  #not so good model fit, but best one found
check_normality(model_net_ps_Pve)
check_heteroscedasticity(model_net_ps_Pve)

results <- tidy(glht(model_net_ps_Pve, linfct = 
  mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

#Results for P. verrucosa, respiration

model_respiration_Pve <- tab_psrate %>% 
  filter(Species=="Pve", sampling_point=="T1") %>% 
  lmer(log(respiration)~treatment + (1|colony), data=.)

#qqPlot(residuals(model_respiration_Pve))  #good model fit
check_normality(model_respiration_Pve)
check_heteroscedasticity(model_respiration_Pve)

tidy(glht(model_respiration_Pve, linfct = 
  mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

#Results for P. verrucosa, gross photosynthesis 

model_gross_ps_Pve <- tab_psrate %>% 
  filter(Species=="Pve", sampling_point=="T1") %>% 
  lmer(log(gross_ps)~treatment + (1|colony), data=.)

#qqPlot #not so good model fit, but best one found
check_normality(model_gross_ps_Pve)
check_heteroscedasticity(model_gross_ps_Pve)

model_gross_ps_Pve <- tidy(glht(model_gross_ps_Pve, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

model_gross_ps_Pve_Result <- model_gross_ps_Pve

model_gross_ps_Pve_Result

## S. Pistillata -----------------------------------------------------------

model_net_ps_Spi <- tab_psrate %>% 
  filter(Species=="Spi", sampling_point=="T1") %>% 
  lmer(log(net_ps)~treatment + (1|colony), data=.)

# Check Model
check_normality(model_net_ps_Spi)
check_heteroscedasticity(model_net_ps_Spi)

tidy(glht(model_net_ps_Spi, linfct = mcp(treatment = "Tukey"))) %>%
  add_significance("adj.p.value")

#Results for S. pistillata, respiration 

model_respiration_Spi <- tab_psrate %>% 
  filter(Species=="Spi", sampling_point=="T1") %>% 
  lmer(log(respiration)~treatment + (1|colony), data=.)

#Check Model

check_normality(model_respiration_Spi)
check_heteroscedasticity(model_respiration_Spi)

tidy(glht(model_respiration_Spi, linfct = mcp(treatment = 
    "Tukey"))) %>%  add_significance("adj.p.value")

# Results for S. pistillata, gross photosynthesis 

model_gross_ps_Spi <- tab_psrate %>% 
  filter(Species=="Spi", sampling_point=="T1") %>% 
  lmer(log(gross_ps)~treatment + (1|colony), data=.)

# Check Model

check_normality(model_gross_ps_Spi)
check_heteroscedasticity(model_gross_ps_Spi)

model_gross_ps_Spi <- tidy(glht(model_gross_ps_Spi, 
      linfct = mcp(treatment = "Tukey"))) %>%  
  add_significance("adj.p.value")

model_gross_ps_Spi_Result <- model_gross_ps_Spi


model_gross_ps_Spi_Result

# Null hypothesis - T0 Baseline -------------------------------------------
#All data analysis of the Baseline

## P. verrucosa ------------------------------------------------------------

#Results for Pve/ net photosynthesis T0(BASELINE)

tab_psrate <- na.omit(tab_psrate)

model_net_ps_Pve_T0 <- tab_psrate %>% 
  drop_na() %>% 
  filter(Species=="Pve", sampling_point=="T0") %>% 
  lmer((log(scale(net_ps)))~treatment + (1|colony), data=.)

# Check Model
check_normality(model_net_ps_Pve_T0)
check_heteroscedasticity(model_net_ps_Pve_T0)


model_net_ps_Pve_T0_result <- tidy(glht(model_net_ps_Pve_T0, linfct = mcp(treatment = "Tukey"))) %>%  
  add_significance("adj.p.value")

model_net_ps_Pve_T0_result


write.table(model_net_ps_Pve_T0_result, sep = ";", "Output/Tab_Result_Lmem_net_ps_Pve_T0.csv")

#Results for P. verrucosa, respiration / T0 baseline

tab_psrate <- na.omit(tab_psrate)

model_respiration_Pve_T0 <- tab_psrate %>% 
  filter(Species=="Pve", sampling_point=="T0") %>% 
  lmer((log(scale(respiration)))~treatment + (1|colony), data=.)

# Check Model
check_normality(model_respiration_Pve_T0)
check_heteroscedasticity(model_respiration_Pve_T0)


model_RES_Pve_T0_result <- tidy(glht(model_respiration_Pve_T0, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

model_RES_Pve_T0_result

write.table(model_RES_Pve_T0_result, sep = ";", "Output/Tab_Result_Lmem_res_Pve_T0.csv")

# Gross P / Pve / T0 Baseline

tab_psrate <- na.omit(tab_psrate)

model_gross_ps_Pve_T0 <- tab_psrate %>% 
  filter(Species=="Pve", sampling_point=="T0") %>% 
  lmer((log(scale(gross_ps)))~treatment + (1|colony), data=.)


# Check Model
check_normality(model_gross_ps_Pve_T0)
check_heteroscedasticity(model_gross_ps_Pve_T0)


model_gross_ps_Pve_T0 <- tidy(glht(model_gross_ps_Pve_T0, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

model_gross_ps_Pve_Result_T0 <- model_gross_ps_Pve_T0

model_gross_ps_Pve_Result_T0

write.table(model_gross_ps_Pve_Result_T0, sep = ";", "Output/Tab_Result_Lmem_gross_ps_Pve_T0.csv")

#Results for SPI/ net photosynthesis T0(BASELINE)

tab_psrate <- na.omit(tab_psrate)

model_net_ps_Spi_T0 <- tab_psrate %>% 
  filter(Species=="Spi", sampling_point=="T0") %>% 
  lmer((sqrt(scale(net_ps)))~treatment + (1|colony), data=.)


#Check Model
check_normality(model_net_ps_Spi_T0)
check_heteroscedasticity(model_net_ps_Spi_T0)


model_net_ps_Spi_T0_result <- tidy(glht(model_net_ps_Spi_T0, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

model_net_ps_Spi_T0_result


write.table(model_net_ps_Spi_T0_result, sep = ";", "Output/Tab_Result_Lmem_net_ps_Spi_T0.csv")

#Respiration Spi / T0 Baseline


model_respiration_Spi_T0 <- tab_psrate %>% 
  filter(Species=="Spi", sampling_point=="T0") %>% 
  lmer(log(respiration)~treatment + (1|colony), data=.)

# Check Model
check_normality(model_respiration_Spi_T0)
check_heteroscedasticity(model_respiration_Spi_T0)


model_RES_Spi_T0_result <- tidy(glht(model_respiration_Spi_T0, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

model_RES_Spi_T0_result

write.table(model_RES_Spi_T0_result, sep = ";", "Output/Tab_Result_Lmem_res_Spi_T0.csv")

# Gross P / Spi / T0 Baseline

model_gross_ps_Spi_T0 <- tab_psrate %>% 
  filter(Species=="Spi", sampling_point=="T0") %>% 
  lmer(log(gross_ps)~treatment + (1|colony), data=.)

 #not so good model fit, but best one found
check_normality(model_gross_ps_Spi_T0)
check_heteroscedasticity(model_gross_ps_Spi_T0)

#check_model(model_net_ps_Spi)

model_gross_ps_Spi_T0 <- tidy(glht(model_gross_ps_Spi_T0, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

model_gross_ps_Spi_Result_T0 <- model_gross_ps_Spi_T0


model_gross_ps_Spi_Result_T0

write.table(model_gross_ps_Spi_Result_T0, sep = ";", "Output/Tab_Result_Lmem_gross_ps_Spi_T0.csv")

## Data analysis PAM data ------------------------------------------------
#We analyse Effective quantum yield (Y(II)) from 
#light-adapted sampled (mean values from 3 
#measurements) and Maximum quantum yield (Fv/Fm)
#from dark-adapted samples.

## Descriptive statistics

### Check Outliers

outliers_YII <- tab_light_mean %>%
  drop_na() %>%
  group_by(Treatment, Species) %>%
  identify_outliers(YII)

outliers_YII
# no extreme outliers

tab_light_mean <- tab_light_mean[-79,]

###
outliers_FvFm <- tab_dark %>%
  drop_na() %>%
  group_by(Treatment, Species) %>%
  identify_outliers(FvFm)

outliers_FvFm
# one extreme outliers

## Overview mean and SD values ---------------------------------------------

skim_1 <- tab_light_mean %>% 
  skim()

summy_1 <- tab_light_mean %>%
  drop_na() %>%
  group_by(Species, Treatment, Date) %>%
  get_summary_stats(YII) 


#skim_1
datatable(summy_1, caption = "Table 5: Descriptive summary statistics Effective quantum yield (Y(II)).")

skim_1 <- tab_dark %>% 
  skim()

summy_1 <- tab_dark %>%
  drop_na() %>%
  group_by(Species, Treatment) %>%
  get_summary_stats(FvFm) 


#skim_1
datatable(summy_1, caption = "Table 6: 
Descriptive summary statistics Maximum quantum 
          yield (Fv/Fm).")

## Test null hypotheses -T1 -------------------------------------------------
#All analysis after 6 weeks 

# P. verrucosa ------------------------------------------------------------

model_YII_Pve <- tab_light_mean %>% 
  filter(Species=="Pve", Date=="T1") %>% 
  lmer((sqrt(scale(YII)))~Treatment + (1|Colony), data=.)

# Check Model
check_normality(model_YII_Pve)
check_heteroscedasticity(model_YII_Pve)

#check_model(model_YII_Pve)


model_FvFm_Pve <- tab_dark %>% 
  filter(Species=="Pve", Date=="T1") %>% 
  lmer(scale(FvFm)~Treatment + (1|Colony), data=.)

#Check Model
check_normality(model_FvFm_Pve)
check_heteroscedasticity(model_FvFm_Pve)


tidy(glht(model_YII_Pve, linfct = mcp(Treatment = 
"Tukey"))) %>%  add_significance("adj.p.value")

tidy(glht(model_FvFm_Pve, linfct = mcp(Treatment = 
"Tukey"))) %>%  add_significance("adj.p.value")

# S. pistillata -----------------------------------------------------------

model_YII_Spi <- tab_light_mean %>% 
  filter(Species=="Spi", Date=="T1") %>% 
  lmer(scale(YII)~Treatment + (1|Colony), data=.)

#Check Model
check_normality(model_YII_Spi)
check_heteroscedasticity(model_YII_Spi)

#check_model(model_YII_Spi)

model_FvFm_Spi <- tab_dark %>% 
  filter(Species=="Spi", Date=="T1") %>% 
  lmer(scale(FvFm)~Treatment + (1|Colony), data=.)

#Check Model
check_normality(model_FvFm_Spi)
check_heteroscedasticity(model_FvFm_Spi)

tidy(glht(model_YII_Spi, linfct = mcp(Treatment = 
    "Tukey"))) %>%  add_significance("adj.p.value")

tidy(glht(model_FvFm_Spi, linfct = mcp(Treatment = 
"Tukey"))) %>%  add_significance("adj.p.value")

## Null hypothesis T0 Baseline ---------------------------------------------

# P. verrucosa ------------------------------------------------------------

#Effective quantum yield / Pve / T0 Baseline

#Means

YII_PVE_Problem <- tab_light_mean %>% 
  filter(Species== "Spi", Date=="T0")%>%
  group_by(Treatment)%>%
  summarize(mean_II= mean(YII))

view(YII_PVE_Problem)

#Test

model_YII_Pve_T0 <- tab_light_mean %>% 
  filter(Species=="Pve", Date=="T0") %>% 
  lmer((log(scale(YII)))~Treatment + (1|Colony), data=.)

#not so good model fit, but best one found

check_normality(model_YII_Pve_T0)
check_heteroscedasticity(model_YII_Pve_T0)

#check_model(model_YII_Pve_T0)

model_YII_Pve_Result_T0 <- tidy(glht(model_YII_Pve_T0, linfct = mcp(Treatment = "Tukey"))) %>%  add_significance("adj.p.value")


model_YII_Pve_Result_T0

write.table(model_YII_Pve_Result_T0, sep = ";", "Output/Tab_Result_Lmem_YII_Pve_T0.csv")

#Maximum quantum yield / Pve / T0 Baseline

model_FvFm_Pve_T0 <- tab_dark %>% 
  filter(Species=="Pve", Date=="T0") %>% 
  lmer(scale(FvFm)~Treatment + (1|Colony), data=.)

# Check Modal
check_normality(model_FvFm_Pve_T0)
check_heteroscedasticity(model_FvFm_Pve_T0)

model_FvFm_Pve_Result_T0 <- tidy(glht(model_FvFm_Pve_T0, linfct = mcp(Treatment = 
                                                                        "Tukey"))) %>%  add_significance("adj.p.value")

model_FvFm_Pve_Result_T0

write.table(model_FvFm_Pve_Result_T0, sep = ";", "Output/Tab_Result_Lmem_FvFm_Pve_Result_T0.csv")

# S. pistillata -----------------------------------------------------------

#Effective quantum yield / Spi / T0 Baseline

model_YII_Spi_T0 <- tab_light_mean %>% 
  filter(Species=="Spi", Date=="T0") %>% 
  lmer(scale(YII)~Treatment + (1|Colony), data=.)

#Check Model
check_normality(model_YII_Spi_T0)
check_heteroscedasticity(model_YII_Spi_T0)

#check_model

model_YII_Spi_Result_T0 <- tidy(glht(model_YII_Spi_T0, linfct = mcp(Treatment = "Tukey"))) %>%  add_significance("adj.p.value")


model_YII_Spi_Result_T0

write.table(model_YII_Spi_Result_T0, sep = ";", "Output/Tab_Result_Lmem_YII_Spi_T0.csv")

#Maximum quantum yield / Spi / T0 Baseline

model_FvFm_Spi_T0 <- tab_dark %>% 
  filter(Species=="Spi", Date=="T0") %>% 
  lmer(scale(FvFm)~Treatment + (1|Colony), data=.)

#not so good model fit, but best one found
check_normality(model_FvFm_Spi_T0)
check_heteroscedasticity(model_FvFm_Spi_T0)


model_FvFm_Spi_Result_T0 <- tidy(glht(model_FvFm_Spi_T0, linfct = mcp(Treatment = "Tukey"))) %>%  add_significance("adj.p.value")

model_FvFm_Spi_Result_T0

write.table(model_FvFm_Spi_Result_T0, sep = ";", 
  "Output/Tab_Result_Lmem_FvFm_Spi_Result_T0.csv")

## Plot Big Graphic --------------------------------------------------------

species.labs <- c("P. verrucosa", "S. pistillata")
names(species.labs) <- c("Pve", "Spi")

levels(tab_psrate$treatment)
tab_psrate$treatment = factor(tab_psrate$treatment, levels=c("Control", "MPsallphases", "MPs+Food" ))
levels(tab_psrate$treatment)

levels(tab_light_mean$Treatment)
tab_light_mean$Treatment = factor(tab_light_mean$Treatment, levels=c("Control", "MPsallphases", "MPs+Food" ))
levels(tab_light_mean$Treatment)

levels(tab_psrate$treatment)
tab_psrate$treatment = factor(tab_psrate$treatment, levels=c("Control", "MPsallphases", "MPs+Food" ))
levels(tab_psrate$treatment)

levels(tab_dark$Treatment)
tab_dark$Treatment = factor( tab_dark$Treatment, levels=c("Control", "MPsallphases", "MPs+Food" ))
levels(tab_dark$Treatment)


#Plots 

#x labels 

a <- 11

#Y lebels 

b <- 9

c <- 9

#boxplot

d <- 0.8

#Photosynthesis_Graphic version A

#Net-Photosynthesis_Graphic version A

x = c(1,2.04,2)
xend = c(1.97,3,3)
y <- c(90,90,90)
yend <- c(90,90,90)
Species <- Species <- c("P. verrucosa", "P. verrucosa", "S. pistillata" )


df_segments <- data.frame(x, xend,y,yend, Species)


x= c(1.5,2.5,2.5)
y = c(90,90,90)
label = c("***", "***","**")
Species <- Species <- c("P. verrucosa", "P. verrucosa", "S. pistillata" )


df_segments_2 <- data.frame(x, y,label, Species)


tab_psrate$Species <- ifelse(tab_psrate$Species=="Pve", yes = "P. verrucosa", no = "S. pistillata")#Correr una vez


tab_psrate_T1 <- 
  tab_psrate %>% 
  filter(sampling_point == "T1")


plot_net_ps_Pve <- 
  tab_psrate_T1 %>% 
  filter(sampling_point == "T1") %>%
  ggplot()+
  geom_boxplot(data= tab_psrate_T1, aes(x=treatment, 
              y=net_ps, fill=treatment),
              width=d, outlier.shape = NA, 
              lwd=0.3, color="black")+
  geom_jitter(data= tab_psrate_T1, aes(x=treatment, 
                                   y=net_ps,
                                   fill=treatment),pch = 21,size= 1, stroke = 0.3)+
 geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 0.4, color = "black") + # Add a segment)
 geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 2.5)+
 facet_grid(~Species)+
 theme_bw()+
  theme(legend.position= "none",
        axis.ticks=element_line(size=0.1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        axis.text.x = element_blank(), #Quitar las etiquetas en x
        axis.title.x = element_blank(),#Quitar el titulo en x
        strip.text.x = element_text(size = a, face = "italic"),
        strip.background = element_blank(),
        strip.text.y = element_text(size=b),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        plot.margin=unit(c(-0.2,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(Net~photosynthesis, (µg~O[2]~cm^{-2}~h^{-1}))))+
  scale_y_continuous(limits=c(0, max(tab_psrate_T1$net_ps, na.rm = TRUE)*1.3))+
  scale_x_discrete(labels=c("", "", ""))

plot(plot_net_ps_Pve)
  
  
#Respiration_Graphic version A

x = c(2,1,2)
xend = c(2.93,3,3)
y <- c(46,50,46)
yend <- c(46,50,46)
Species <- Species <- c("P. verrucosa", "S. pistillata", "S. pistillata" )


df_segments <- data.frame(x, xend,y,yend, Species)


x= c(2.5,2,2.5)
y = c(46,50,46)
label = c("***", "***","**")
Species <- Species <- c("P. verrucosa", "S. pistillata", "S. pistillata" )


df_segments_2 <- data.frame(x, y,label, Species)



plot_respiration_Pve <- 
  tab_psrate %>% 
 filter(sampling_point=="T1") %>%
  ggplot()+
  geom_boxplot(data= tab_psrate, aes(x=treatment, 
                                     y=respiration, fill=treatment),
               width=d, outlier.shape = NA, 
               lwd=0.3, color="black")+
  geom_jitter(data= tab_psrate, aes(x=treatment, 
                                    y=respiration,
                                    fill=treatment),pch = 21,size= 1, stroke = 0.3)+
  geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 0.4, color = "black") + # Add a segment)
  geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 2.5)+
  facet_grid(~Species)+
  theme_bw()+
  theme(legend.position= "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_line(size=0.1),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=b),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        axis.title.x = element_blank(),
        plot.margin=unit(c(0.2,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(Respiration, (µg~O[2]~cm^{-2}~h^{-1}))))+
  scale_y_continuous(limits=c(0, max(tab_psrate$respiration, na.rm = TRUE)*1.3))+
  scale_x_discrete(labels=c("", "", ""))

plot(plot_respiration_Pve)



#Gross P_Graphic version A

x = c(1,1,2.1,1,2)
xend = c(3,1.98,3,3,3)
y <- c(145,138,138,145,138)
yend <- c(145,138,138,145,138)
Species <- Species <- c("P. verrucosa", "P. verrucosa", "P. verrucosa", "S. pistillata", "S. pistillata"   )


df_segments <- data.frame(x, xend,y,yend, Species)


x= c(2,1.5,2.5,2,2.5)
y = c(145,138,138,145,138)
label = c("**", "***","***", "*", "***" )
Species <- Species <- c("P. verrucosa", "P. verrucosa", "P. verrucosa", "S. pistillata" , "S. pistillata"   )

df_segments_2 <- data.frame(x, y,label, Species)


plot_gross_ps_Pve <- 
  tab_psrate %>% 
  filter(sampling_point=="T1") %>%
  ggplot()+
  geom_boxplot(data= tab_psrate, aes(x=treatment, 
                                     y=gross_ps, fill=treatment),
               width=d, outlier.shape = NA, 
               lwd=0.3, color="black")+
  geom_jitter(data= tab_psrate, aes(x=treatment, 
                                    y=gross_ps,
                                    fill=treatment),pch = 21,size= 1, stroke = 0.3)+
  geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 0.4, color = "black") + # Add a segment)
  geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 2.5)+
  facet_grid(~Species)+
  theme_bw()+
  theme(legend.position= "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_line(size=0.1),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=b),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        axis.title.x = element_blank(),
        plot.margin=unit(c(-0.2,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(Gross~photosynthesis, (µg~O[2]~cm^{-2}~h^{-1}))))+
  scale_y_continuous(limits=c(0, max(tab_psrate$gross_ps, na.rm = TRUE)*1.3))+
  scale_x_discrete(labels=c("", "", ""))

plot(plot_gross_ps_Pve)

#PAM_Graphic version A

tab_light_mean$Species <- ifelse(tab_light_mean$Species=="Pve", yes = "P. verrucosa", no = "S. pistillata")#Correr una vez



plot_YII_Pve <- 
  tab_light_mean %>% 
  filter(Date=="T1") %>%
  ggplot()+
  geom_boxplot(data=  tab_light_mean, aes(x=Treatment, 
                                     y=  YII, fill=Treatment),
               width=d, outlier.shape = NA, 
               lwd=0.3, color="black")+
  geom_jitter(data= tab_light_mean, aes(x=Treatment, 
                                    y=YII,
                                    fill=Treatment),pch = 21,size= 1, stroke = 0.3)+
  geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 1, color = "black") + # Add a segment)
  geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 5 )+
  facet_grid(~Species)+
  theme_bw()+
  theme(legend.position= "none", 
       panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_line(size=0.1),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=b),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        axis.title.x = element_blank(),
        plot.margin=unit(c(-0.2,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(Effective~quantum,yield~(Y(II)))))+
  scale_y_continuous(limits=c(min(tab_light_mean$YII, na.rm = TRUE)*1.2, max(tab_light_mean$YII, na.rm = TRUE)*1.1))+
  scale_x_discrete(labels=c("", "", ""))

plot(plot_YII_Pve)


#PAM_Graphic version B


species.labs <- c("P. verrucosa", "S. pistillata")
names(species.labs) <- c("Pve", "Spi")

plot_YII_Pve <- 
  tab_light_mean %>% 
  drop_na() %>% 
  filter(Date=="T1")%>%
  ggplot(aes(x=Treatment, 
             y=YII,
             fill=Treatment)) + 
  labs(title="",
       x="", y = "") +
  facet_grid(~Species, labeller = labeller(Species = species.labs)) +
  scale_fill_manual(values = c("#E7B800", "#008cFF", "#0027A8")) +
  theme_classic()+
  theme(legend.position= "none", 
        legend.justification = c(0, 0),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.ticks=element_line(size=0.1),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=b),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        plot.margin=unit(c(-0.5,0.5,-0.1,0.5), "cm") #top, right, bottom, left
  )+
  ylab(expression(atop(Effective~quantum,yield~(Y(II)))))+
  scale_y_continuous(limits=c(min(tab_light_mean$YII, na.rm = TRUE)*1.2, max(tab_light_mean$YII, na.rm = TRUE)*1.1))+
  geom_boxplot(width=d, outlier.shape = NA, lwd=0.3, color="black", aes(fill = Treatment)) +
  scale_x_discrete(labels=c("", "", ""))+
  geom_point(pch = 21, size= 1, stroke= 0.3,position = position_jitterdodge(), aes(fill = Treatment))+
  scale_color_manual(values=c("black", "black", "black", "black", "black"), name="")

plot(plot_YII_Pve)


#Maximun_Graphic version A


x = c(1,1)
xend = c(2,3)
y <- c(0.75,0.75)
yend <- c(0.75,0.75)
Species <- Species <- c("P. verrucosa", "S. pistillata")


df_segments <- data.frame(x, xend,y,yend, Species)


x= c(1.5,2)
y = c(0.75,0.75)
label = c("**", "*")
Species <- Species <- c("P. verrucosa", "S. pistillata")

df_segments_2 <- data.frame(x, y,label, Species)

tab_dark$Species <- ifelse(  tab_dark$Species=="Pve", yes = "P. verrucosa", no = "S. pistillata")#Correr una vez

plot_FvFm_Pve<- 
  tab_dark %>% 
  filter(Date=="T1") %>%
  ggplot()+
  geom_boxplot(data=    tab_dark, aes(x=Treatment, 
                                          y=  FvFm, fill=Treatment),
               width=d, outlier.shape = NA, 
               lwd=0.3, color="black")+
  geom_jitter(data= tab_dark, aes(x=Treatment, 
                                        y=FvFm,
                                        fill=Treatment),pch = 21,size= 1, stroke = 0.3)+
  geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 0.4, color = "black") + # Add a segment)
  geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 2.5)+
  facet_grid(~Species)+
  theme_bw()+
  theme(legend.position= "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_line(size=0.1),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=b),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        plot.margin=unit(c(-0.2,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(Maximum~quantum, yield~(Fv/Fm))))+
  labs(x='Treatment' )+
  scale_y_continuous(limits=c(min(tab_dark$FvFm, na.rm = TRUE)*1.2, 0.76))+
  scale_x_discrete(labels=c("Control", "MP", "MP+HF"))

plot(plot_FvFm_Pve)

#Join Plots

summaryplot_PSrate_full <- ggarrange (plot_net_ps_Pve, plot_respiration_Pve, plot_gross_ps_Pve, plot_YII_Pve, plot_FvFm_Pve,
                               
                                      align = "v", #common.legend = FALSE, 
                                      heights= c(1.3,1.3,1.3,1.3,1.4),
                                      ncol=1) 

summaryplot_PSrate_full



#save plot

ggsave(here ("Output", "summaryplot_PSrate_full_one_colum.png"), 
       plot = summaryplot_PSrate_full,  
       scale = 1, width = 11, height = 19, units = c("cm"),
       dpi = 600, limitsize = FALSE)

#The end