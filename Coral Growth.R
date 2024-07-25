
## Coral Growth -----------------------------------------------------------
# Entire Data analysis related to the Coral Growth

## Upload data tables -------------------------------------------------------

# Read table with infos to IDs, tanks, colonies and treatments

tab_info <- read.csv2("Data/info.csv",sep=";", dec=".",header=T) %>% 
  na.omit() %>% 

tab_info$treatment <- factor(tab_info$treatment, 
                             levels = c("Control", "MP+HF", "MP"))

# read buoyant weight data

tab_size <- read.csv2("Data/3D_Scanning_Data_Table.csv",sep=";", dec=".",header=T) %>% 
  na.omit() %>% 
  select(ID, Volume_t0, Surface_t0, Volume_t1, Surface_t1)

# Read buoyant weight data

tab_calcification <- read.csv2("Data/Weight_Data Table.csv",sep=";", 
                    dec=".",header=T) %>% 
  na.omit()

tab_calcification <-  merge(tab_calcification, tab_size, by="ID",
                            all = TRUE)

## Data processing for buoyant weight data ---------------------------------------------------------
# Here we standardize the difference in weight to the Surface_t1 area 
# of the coral colony. First we transform skeletal weight under water to 
# skeletal weight using the seacarb package with a the mean measured temperature 
# and salinity at the timepoints and a skeletal density of 2,93 g/cm³. 
# Then we use the photosynthetically active surface area at t1 (first measurement,
# actually corresponding to t0), after removing the necrotic areas from the 
# model (column surface_t1). The results are stated as mg/cm².

# Calculate and standardize

skeleton_weight <- function(buoyant_weight, S, T, P = 0, 
                            rho_aragonite = 2930){
  
  x <- seacarb::rho(S = S, T = T, P = P)
  y <- buoyant_weight / (1 - (x / rho_aragonite))
  attributes(y) <- NULL
  y
}

tab_calcification$weight_t0 <- skeleton_weight(tab_calcification$weight_.T0, S=34.3, T=23.8, P = 0, rho_aragonite = 2930)
tab_calcification$weight_t1 <- skeleton_weight(tab_calcification$weight_.T1, S=33, T=24.8, P = 0, rho_aragonite = 2930)


tab_calcification <- tab_calcification %>% 
  mutate(calcification = (((weight_.T1-weight_.T0)*1000)/(Surface_t1/100)))


write.table(tab_calcification, sep = ";", "Output/Tab_calcification_processed.csv")


## Data processing for growth in surface area and volume  ---------------
# Here we calculate the growth in surface area (%) and volume (cm³ per cm²) 
# of the coral colony.

# Calculate growth in surface area (%)

tab_size$growth_surf <- ((100/tab_size$Surface_t0)*tab_size$Surface_t1)-100

# Growth in volume (cm² per cm³)

tab_size$growth_vol <- (tab_size$Volume_t1 - tab_size$Volume_t0)/1000/(tab_size$Surface_t0/100)

str(tab_size)

tab_growth <- merge(tab_info, tab_size, by = "ID", all = TRUE)

write.table(tab_growth, sep = ";", "Output/Tab_growth_processed.csv")

## Means Buoyant weight ---------------------------------------------------

# Calcification mean

calcification_means <- tab_calcification %>% 
  
  drop_na() %>%
  
  group_by(Species, treatment) %>% summarize(calcification_means = mean(calcification))

calcification_means

# Calcification sd

calcification_sd <- tab_calcification %>% 
  
  drop_na() %>%
  
  group_by( Species) %>% summarize(calcification_sd = sd(calcification))

calcification_sd

## Data analysis Buoyant weight ------------------------------------------
# Here we analyze calcification rates of the corals. 
# These were recorded using buoyant weight at two 
# timepoints of the experiment.

# Descriptive statistics

### Check Outlier
### Check data for (extreme) outliers.

outliers_calcification <- tab_calcification %>%
  drop_na() %>%
  group_by(treatment, Species) %>%
  identify_outliers(calcification) 

outliers_calcification

### Results: 2 outliers; (Spi Control DT10) and (Spi MP CT4)

## Overview mean and SD values -----------------------------------------
# Here we get an impression on the data. 
# skimr is designed to provide summary statistics about variables 
# in data frames, tibbles, data tables and vectors. 
# It is opinionated in its defaults, but easy to modify.
# In base R, the most similar functions are summary() for 
# vectors and data frames and fivenum() for numeric vectors:

skim_1 <- tab_calcification %>% 
  skim()

summy_1 <- tab_calcification %>%
  drop_na() %>%
  group_by(Species, treatment) %>%
  get_summary_stats(calcification)  # package: rstatix

#skim_1

datatable(summy_1, caption = "Table 5: Descriptive summary 
          statistics of calcification (mg/cm²")

## Test null hypotheses of Calcification -------------------------------------------------
# Now we want to test if our treatments have a signifiant influence. 
# We use linear mixed effects models that include the coral colony as
# random factor. We construct the model and then test the distribution 
# of the model residuals for normality and heteroscedasticity. 
# We had to log transform the data to meet the model criteria. 
# This might change if the outliers are corrected or removed, 
# so go back here and check again if this is correct.

### Linear mixed effects models 

### Results for P. verrucosa

model_calcification_Pve <- tab_calcification %>% 
  filter(Species=="Pve") %>% 
  lmer(scale(calcification)~treatment + (1|colony), data=.)

### Check Model

check_normality(model_calcification_Pve)
check_heteroscedasticity(model_calcification_Pve)

#### Check_model(model_calcification_Pve)

result <- tidy(glht(model_calcification_Pve, 
linfct = mcp(treatment = "Tukey"))) %>%  add_significance
("adj.p.value")

result

write.table(result, sep = ";", 
            "Output/Tab_Result_Lmem_calcification_Pve.csv")

# Results for S. pistillata

model_calcification_Spi <- tab_calcification %>% 
  filter(Species=="Spi") %>% 
  lmer(sqrt(scale(calcification))~treatment + (1|colony), data=.)

### Check Model

check_normality(model_calcification_Spi)
check_heteroscedasticity(model_calcification_Spi)

result <- tidy(glht(model_calcification_Spi, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result

write.table(result, sep = ";", 
            "Output/Tab_Result_Lmem_calcification_spi.csv")


## Data analysis growth in surface area and volume  ----------------------
#Here we analyze growth rates of the corals. These were recorded 
#using 3D scanning at two timepoints (t0 and t1) of the experiment. 
#t0 corresponds to baseline  and  t1 corresponds to 6 weeks after 
#experiment

# Descriptive statistics

### Check Outlier

#Check data for (extreme) outliers. Surface growth is good, 
#no extreme outliers. Volume has still 2 extreme outliers, 
#A-T10-2 and A-T6-2 but they seem to be acceptable. 
#Fragment B-T1 was removed from all analyses.

outliers_growth_surf <- tab_growth %>%
  drop_na() %>%
  group_by(treatment, Species) %>%
  identify_outliers(growth_surf) 

outliers_growth_surf

# Result: surface is good, no extreme outliers

outliers_growth_vol <- tab_growth %>%
  drop_na() %>%
  group_by(treatment, Species) %>%
  identify_outliers(growth_vol) 

outliers_growth_vol

# Three extreme outliers: 
#(Pve MP+HF A-T6-2)
#(Pve MP C-T4-2)
#(Spi MP B-T14)


## Overview mean and SD values -----------------------------------------
#Here we get an impression on the data

skim_1 <- tab_growth %>% 
  skim()

summy_1 <- tab_growth %>%
  drop_na() %>%
  group_by(Species, treatment) %>%
  get_summary_stats(growth_surf)  # package: rstatix


#skim_2
datatable(summy_1, caption = "Table 5: Descriptive summary 
          statistics of surface growth (%)")

summy_2 <- tab_growth %>%
  drop_na() %>%
  group_by(Species, treatment) %>%
  get_summary_stats(growth_vol)  # package: rstatix


#skim_2

datatable(summy_2, caption = "Table 6: Descriptive summary 
          statistics of volume growth (cm³ per cm²)")


## Test null hypotheses of Surface and Volume  ------------------------------------------------

# Surface

# Linear mixed effects models 

###Results for P. verrucosa

model_surface_Pve <- tab_growth %>% 
  filter(Species=="Pve") %>% 
  lmer(scale(growth_surf)~treatment + (1|colony), data=.)

#Check Model

check_normality(model_surface_Pve)
check_heteroscedasticity(model_surface_Pve)


result <- tidy(glht(model_surface_Pve, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result

write.table(result, sep = ";", "Output/Tab_Result_Lmem_Surface_Pve.csv")


# Linear mixed effects models 
###Results for S. pistillata

model_surface_Spi <- tab_growth %>% 
  filter(Species=="Spi") %>% 
  lmer(scale(growth_surf)~treatment + (1|colony), data=.)

#Check Model

check_normality(model_surface_Spi)
check_heteroscedasticity(model_surface_Spi)

result <- tidy(glht(model_surface_Spi, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result

write.table(result, sep = ";", "Output/Tab_Result_Lmem_surface_spi.csv")


#Volume 

# Linear mixed effects models 

###Results for P. verrucosa

model_volume_Pve <- tab_growth %>% 
  filter(Species=="Pve") %>% 
  lmer(scale(growth_vol))~treatment + (1|colony), data=.)

# Check Model

check_normality(model_volume_Pve)
check_heteroscedasticity(model_volume_Pve)

result <- tidy(glht(model_volume_Pve, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result

write.table(result, sep = ";", "Output/Tab_Result_Lmem_volume_Pve.csv")

# Linear mixed effects models 

###Results for S. pistillata

model_volume_Spi <- tab_growth %>% 
  filter(Species=="Spi") %>% 
  lmer(scale(growth_vol)~treatment + (1|colony), data=.)

# Check Model

check_normality(model_volume_Spi)
check_heteroscedasticity(model_volume_Spi)

result <- tidy(glht(model_volume_Spi, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result

write.table(result, sep = ";", "Output/Tab_Result_Lmem_volume_spi.csv")

## Plot data Total - Big Graphic  ------------------------------------------------------------

# Coral Growth T1

species.labs <- c("P. verrucosa", "S. pistillata")
names(species.labs) <- c("Pve", "Spi")


#Plots 

#x labels 

a <- 11

#Y lebels 

b <- 9

c <- 9

#boxplot

d <- 0.8

### Plots

# Surface

plot_surface <- 
  tab_growth %>% 
  ggplot(aes(x=treatment, 
             y=growth_surf)) + 
  facet_grid(~Species, labeller = labeller(Species = species.labs)) +
  scale_fill_manual(values = c("#E7B800", "#008cFF", "#0027A8")) +
  theme_classic()+
  theme(legend.position= c(0, 3), 
        legend.justification = c(0, 0),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_line(size=0.1),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        strip.background = element_blank(),
        strip.text.x = element_text(size = a, face = "italic"),
        strip.text.y = element_text(size=b),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        plot.margin=unit(c(-0.1,0.3,-0.4,0.3), "cm"))+ #top, right, bottom, left
  ylab(expression(atop(Surface~area, growth~("%"))))+ 
  #scale_y_continuous(limits=c(-0.01,0.06))+
  geom_boxplot(width=d, outlier.shape = NA, lwd=0.3, color="black", aes(fill = treatment))+
  scale_x_discrete(labels=c("", "", ""))+
  geom_point(pch = 21, size= 1, stroke = 0.3,  position = position_jitterdodge(), aes(fill = treatment))+
  scale_color_manual(values=c("black", "black", "black", "black", "black"), name="")

plot_surface

# Volume

plot_volume <- 
  tab_growth %>% 
  ggplot(aes(x=treatment, 
             y=growth_vol)) + 
  facet_grid(~Species, labeller = labeller(Species = species.labs)) +
  scale_fill_manual(values = c("#E7B800", "#008cFF", "#0027A8")) +
  theme_classic()+
  theme(legend.position= c(0, 3), 
        legend.justification = c(0, 0),
        legend.direction = "horizontal",
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
        plot.margin=unit(c(-0.1,0.3,-0.4,0.3), "cm"))+
  ylab(expression(atop(Volume~growth, (cm^3~cm^-2))))+
  #scale_y_continuous(limits=c(-0.01,0.06))+
  geom_boxplot(width=d, outlier.shape = NA, lwd=0.3, color="black", aes(fill = treatment))+
  scale_x_discrete(labels=c("", "", ""))+
  geom_point(pch = 21, size= 1, stroke = 0.3,  position = position_jitterdodge(), aes(fill = treatment))+
  scale_color_manual(values=c("black", "black", "black", "black", "black"), name="")

plot_volume

# Calcification

plot_calcification <- 
  tab_calcification %>% 
  ggplot(aes(x=treatment, 
             y=calcification)) + 
  facet_grid(~Species) +
  scale_fill_manual(values = c("#E7B800", "#008cFF", "#0027A8")) +
  theme_classic()+
  theme(legend.position= "none", 
        legend.justification = c(0, 0),
        legend.direction = "horizontal",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        strip.text.y = element_text(size=b),
        strip.text.x = element_blank(),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        plot.margin=unit(c(-0.1,0.3,-0.07,0.3), "cm"))+ #top, right, bottom, left
  ylab(expression(atop(Calcification, (mg~CaCO^3~cm^-2))))+
  labs(x='Treatment')+
  #scale_y_continuous(limits=c(-0.01,0.06))+
  geom_boxplot(width=d, outlier.shape = NA, lwd=0.3, color="black", aes(fill = treatment))+
  scale_x_discrete(labels=c("Control", "MP", "MP+HF"))+
  geom_point(pch = 21, size= 1, stroke= 0.3, position = position_jitterdodge(), aes(fill = treatment))+
  scale_color_manual(values=c("black", "black", "black", "black", "black"), name="")

plot_calcification

### Summary

summaryplot <- ggarrange(plot_surface, plot_volume, plot_calcification,
                         #labels = c("a)", "b)", "c)", "d)"), 
                         align = "v", #common.legend = FALSE, 
                         heights= c(0.8,0.7,0.7),
                         ncol=1)

summaryplot

# Save plot
ggsave("Output/Summaryplot_growth_one_column.png", plot = summaryplot,
       scale = 1, width = 11, height = 13, units = c("cm"),
       dpi = 600, limitsize = FALSE)

### The End

