
## Energy bodged -----------------------------------------------------------
# Entire Data analysis related to the Coral Energy budged

## Upload data tables ---------------------------------------------------

tab_monaco <- read.csv("Data/Monaco_Samples.csv",sep=";", dec=".",header=T) 
  

## Data analysis -----------------------------------------------------------

### Check Outlier

##Check data for (extreme) outliers Carbohydrates

outliers_Carb_Total <- tab_monaco %>%
  drop_na() %>%
  group_by(treatment) %>%
  identify_outliers(Carb_Total) 

outliers_Carb_Total

##Check data for (extreme) outliers Lipids

outliers_Lipid_Total <- tab_monaco %>%
  drop_na() %>%
  group_by(treatment) %>%
  identify_outliers(Lip_total) 

outliers_Lipid_Total

##Check data for (extreme) outliers Proteins

outliers_Protein_Total <- tab_monaco %>%
  drop_na() %>%
  group_by(treatment) %>%
  identify_outliers(Prot_Total) 

outliers_Protein_Total


## Test null hypotheses ----------------------------------------------------

# Carbohydrates -----------------------------------------------------------

model_carb <- tab_monaco %>% 
  lmer(scale(Carb_Total)~treatment + (1|colony), data=.)


check_normality(model_carb)
check_heteroscedasticity(model_carb)


result_Carb <- tidy(glht(model_carb, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result_Carb

write.table(result_Carb, sep = ";", "Output/Tab_Result_Lmem_Monaco_Carb.csv")


## Lipids ------------------------------------------------------------------

model_lipids <- tab_monaco %>% 
  lmer(scale(Lip_total)~treatment + (1|colony), data=.)


check_normality(model_lipids)
check_heteroscedasticity(model_lipids)

#check_model(model_lipids)

result_Lipids <- tidy(glht(model_lipids, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result_Lipids

write.table(result_Lipids, sep = ";", "Output/Tab_Result_Lmem_Monaco_Lipids.csv")


## Proteins ----------------------------------------------------------------


model_proteins <- tab_monaco %>% 
  lmer(scale(Prot_Total)~treatment + (1|colony), data=.)


check_normality(model_proteins)
check_heteroscedasticity(model_proteins)


result_Proteins <- tidy(glht(model_proteins, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result_Proteins

write.table(result_Proteins, sep = ";", "Output/Tab_Result_Lmem_Monaco_Proteins.csv")


## Plots -------------------------------------------------------------------

levels(tab_monaco$treatment)
tab_monaco$treatment = factor(tab_monaco$treatment, 
        levels=c("Control", "MP", "MP + HF" ))
levels(tab_monaco$treatment)


#x labels 

a <- 7

#Y lebels 

b <- 6

c <- 6

#boxplot

d <- 0.7


## Carbohydrates -----------------------------------------------------------

#Carbohydrates_Graphic version A

x = c(1,1)
xend = c(3,2)
y <- c(0.27,0.25)
yend <- c(0.27,0.25)



df_segments <- data.frame(x, xend,y,yend)


x= c(2,1.5)
y = c(0.27,0.25)
label = c("**", "***")



df_segments_2 <- data.frame(x, y,label)


plot_Carb  <- 
  tab_monaco %>% 
  #filter(sampling_point == "T1") %>%
  ggplot()+
  geom_boxplot(data= tab_monaco, aes(x=treatment, 
                                        y=Carb_Total, fill=treatment),
               width=d, outlier.shape = NA, 
               lwd=0.3, color="black")+
  geom_jitter(data= tab_monaco, aes(x=treatment, 
                                       y=Carb_Total,
                                       fill=treatment),pch = 21,size= 1, stroke = 0.3)+
  geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 0.4, color = "black") + # Add a segment)
  geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 2.5)+
  #facet_grid(~Species)+
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
        plot.title = element_text(hjust = 0.5,size = a, face = "italic"),
        plot.margin=unit(c(-0.01,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
        labs(title="P. verrucosa", x="", y = "")+
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(Carbohydrates~total, (mg~cm^{-2}))))+
  scale_y_continuous(limits=c(min(tab_monaco$Carb_Total, na.rm = TRUE), max(tab_monaco$Carb_Total, na.rm = TRUE)))+
  scale_x_discrete(labels=c("", "", ""))

plot(plot_Carb)
  



## Lipids ------------------------------------------------------------------

#Lipids_Graphic version A

x = c(1)
xend = c(2)
y <- c(7)
yend <- c(7)


df_segments <- data.frame(x, xend,y,yend)


x= c(1.5)
y = c(7)
label = c("*")

df_segments_2 <- data.frame(x, y,label)


plot_Lip  <- 
  tab_monaco %>% 
  #filter(sampling_point == "T1") %>%
  ggplot()+
  geom_boxplot(data= tab_monaco, aes(x=treatment, 
                                     y=Lip_total, fill=treatment),
               width=d, outlier.shape = NA, 
               lwd=0.3, color="black")+
  geom_jitter(data= tab_monaco, aes(x=treatment, 
                                    y=Lip_total,
                                    fill=treatment),pch = 21,size= 1, stroke = 0.3)+
  geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 0.4, color = "black") + # Add a segment)
  geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 2.5)+
  #facet_grid(~Species)+
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
        plot.title = element_text(hjust = 0.5,size = a, face = "italic"),
        plot.margin=unit(c(-0.01,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
  #labs(title="P. verrucosa", x="", y = "")+
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(Lipids~total, (mg~cm^{-2}))))+
  scale_y_continuous(limits=c(min(tab_monaco$Lip_total, na.rm = TRUE), max(tab_monaco$Lip_total, na.rm = TRUE)))+
  scale_x_discrete(labels=c("", "", ""))

plot(plot_Lip)


# Proteins ----------------------------------------------------------------


#Protein_Graphic version A


x = c(1)
xend = c(2)
y <- c(1.55)
yend <- c(1.55)


df_segments <- data.frame(x, xend,y,yend)


x= c(1.5)
y = c(1.55)
label = c("*")

df_segments_2 <- data.frame(x, y,label)


plot_Prot  <- 
  tab_monaco %>% 
  #filter(sampling_point == "T1") %>%
  ggplot()+
  geom_boxplot(data= tab_monaco, aes(x=treatment, 
                                     y=Prot_Total, fill=treatment),
               width=d, outlier.shape = NA, 
               lwd=0.3, color="black")+
  geom_jitter(data= tab_monaco, aes(x=treatment, 
                                    y=Prot_Total,
                                    fill=treatment),pch = 21,size= 1, stroke = 0.3)+
  geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 0.4, color = "black") + # Add a segment)
  geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 2.5)+
  #facet_grid(~Species)+
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
        plot.title = element_text(hjust = 0.5,size = a, face = "italic"),
        plot.margin=unit(c(-0.01,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
  #labs(title="P. verrucosa", x="", y = "")+
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(Proteins~total, (mg~cm^{-2}))))+
  scale_y_continuous(limits=c(min(tab_monaco$Prot_Total, na.rm = TRUE), max(tab_monaco$Prot_Total, na.rm = TRUE)))+
  scale_x_discrete(labels=c("", "", ""))

plot(plot_Prot)


# Total Energy Reserves  --------------------------------------------------


# Upload data table -------------------------------------------------------


tab_monaco_2 <- read.csv("Data/Monaco_Samples_2.csv",sep=";", dec=".",header=T)  
 

# Data analysis -----------------------------------------------------------

### Check Outlier


outliers_Total_KJ <- tab_monaco_2 %>%
  drop_na() %>%
  group_by(treatment) %>%
  identify_outliers(Total_KJ) 

outliers_Total_KJ


# Test Null Hypothesys ----------------------------------------------------


model_Total_KJ <- tab_monaco_2 %>% 
  lmer(scale(Total_KJ)~treatment + (1|colony), data=.)


check_normality(model_Total_KJ)
check_heteroscedasticity(model_Total_KJ)


model_Total_KJ <- tidy(glht(model_Total_KJ, 
linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

model_Total_KJ

write.table(model_Total_KJ, sep = ";", "Output/Tab_Result_Lmem_Monaco_Total_KJ.csv")


# Plot --------------------------------------------------------------------

levels(tab_monaco_2$treatment)
tab_monaco_2$treatment = factor(tab_monaco_2$treatment, levels=c("Control", "MP", "MP + HF" ))
levels(tab_monaco_2$treatment)

#Total KJ_Graphic version A


x = c(1)
xend = c(2)
y<- c(299)
yend <- c(299)


df_segments <- data.frame(x, xend,y,yend)


x= c(1.5)
y = c(299)
label = c("*")

df_segments_2 <- data.frame(x, y,label)



plot_Total_KJ <- 
  tab_monaco_2 %>% 
  #filter(sampling_point == "T1") %>%
  ggplot()+
  geom_boxplot(data=  tab_monaco_2, aes(x=treatment, 
                                     y=Total_KJ, fill=treatment),
               width=d, outlier.shape = NA, 
               lwd=0.3, color="black")+
  geom_jitter(data=  tab_monaco_2, aes(x=treatment, 
                                    y=Total_KJ,
                                    fill=treatment),pch = 21,size= 1, stroke = 0.3)+
  geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 0.4, color = "black") + # Add a segment)
  geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 2.5)+
  #facet_grid(~Species)+
  theme_bw()+
  theme(legend.position= "none",
        axis.ticks=element_line(size=0.1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        axis.text.x = element_blank(), #Quitar las etiquetas en x
        axis.title.x = element_text(size = a),#Quitar el titulo en x
        strip.text.x = element_text(size = a, face = "italic"),
        strip.background = element_blank(),
        strip.text.y = element_text(size=b),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        plot.title = element_text(hjust = 0.5,size = a, face = "italic"),
        plot.margin=unit(c(-0.01,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
  #labs(title="P. verrucosa", x="", y = "")+
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(Total~energy~content, (J~cm^{-2}))))+
  labs(x='Treatment')+
  scale_y_continuous(limits=c(min(tab_monaco_2$Total_KJ, na.rm = TRUE), max(tab_monaco_2$Total_KJ, na.rm = TRUE)))+
  scale_x_discrete(labels=c("Control", "MP", "MP+HF"))

plot(plot_Total_KJ)



# Big Plot all components -------------------------------------------------

#Join Plots

summaryplot_EnergyKJ <- ggarrange (plot_Carb, plot_Lip, plot_Prot,
                                   plot_Total_KJ,
                                      #labels = c("A", "B", "C", "D", "E"),
                                      #font.label = list(size = 5),
                                      #hjust=-2,
                                      #vjust = 3, 
                                      align = "v", #common.legend = FALSE, 
                                      heights= c(1.2,0.92,1,1.20),
                                
                                      ncol=1)
                              

summaryplot_EnergyKJ

#save plot

ggsave(here ("Output", "summaryplot_EnergyKJ.png"), 
       plot = summaryplot_EnergyKJ,  
       scale = 1, width = 4.4, height = 9, units = c("cm"),
       dpi = 600, limitsize = FALSE)
