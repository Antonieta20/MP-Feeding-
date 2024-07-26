## Tissue brightness ---------------------------------
# Entire Data analysis related to the heat response

## Upload data tables -------------------------------------------------------------


tab_pictures <- read.csv2("Data/Data.csv",sep=";", dec=".",header=T)  %>%
  mutate(gray_sum= as.numeric(gray_sum))


tab_pictures$treatment <- factor(   tab_pictures$treatment, 
                                    levels = c("Control", "MP + HF", "MP"))

str(tab_pictures)


# Statistics analyses -----------------------------------------------------

## Pve -------------------------------------------------------------

model_pictures_Pve <- tab_pictures %>% 
  na.omit() %>% 
  filter(specie=="Pve") %>% 
  lmer(scale(gray_sum)~treatment + (1|colony), data=.)

#qqPlot(residuals(model_pictures_Pve))  #not so good model fit, but best one found
check_normality(model_pictures_Pve)
check_heteroscedasticity(model_pictures_Pve)

result_pictures_Pve <- tidy(glht(model_pictures_Pve, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result_pictures_Pve

write.table(result_pictures_Pve, sep = ";", "Output/Tab_Result_Lmem_pictures_Pve.csv")

## Spi -------------------------------------------------------------

model_pictures_Spi <- tab_pictures %>% 
  na.omit() %>% 
  filter(specie=="Spi") %>% 
  lmer(scale(gray_sum)~treatment + (1|colony), data=.)

#qqPlot(residuals(model_pictures_Pve))  #not so good model fit, but best one found
check_normality(model_pictures_Spi)
check_heteroscedasticity(model_pictures_Spi)

result_pictures_Spi <- tidy(glht(  model_pictures_Spi, linfct = mcp(treatment = "Tukey"))) %>%  add_significance("adj.p.value")

result_pictures_Spi

write.table(result_pictures_Spi, sep = ";", "Output/Tab_Result_Lmem_pictures_Spi.csv")  

## Plot -------------------------------------------------------------

tab_pictures <- read.csv2("Data/Data.csv",sep=";", dec=".",header=T) 

species.labs <- c("P. verrucosa", "S. pistillata")
names(species.labs) <- c("Pve", "Spi")

levels(tab_pictures$treatment)
tab_pictures$treatment = factor(tab_pictures$treatment, levels=c("Control", "MP", "MP + HF" ))
levels(tab_pictures$treatment)

#x labels 

a <- 5

#Y lebels 

b <- 5

c <- 5

#boxplot

d <- 0.6


#Plot_Version A


x = c(1,2.1)
xend = c(1.95,3)
y <- c(230,230)
yend <- c(230,230)
Species <- Species <- c("P. verrucosa", "P. verrucosa")


df_segments <- data.frame(x, xend,y,yend, Species)


x= c(1.5,2.5)
y = c(230,230)
label = c("*", "**")
Species <- Species <- c("P. verrucosa", "P. verrucosa")

df_segments_2 <- data.frame(x, y,label, Species)

tab_pictures$specie <- ifelse(tab_pictures$specie =="Pve", yes = "P. verrucosa", no = "S. pistillata")#Correr una vez

plot_pictures_heat <- 
  tab_pictures %>% 
  #filter(Date=="T1") %>%
  ggplot()+
  geom_boxplot(data=    tab_pictures, aes(x=treatment, 
                                      y=  gray_sum, fill=treatment),
               width=0.8, outlier.shape = NA, 
               lwd=0.3, color="black")+
  geom_jitter(data= tab_pictures, aes(x=treatment, 
                                  y=gray_sum,
                                  fill=treatment),pch = 21,size= 0.5, stroke = 0.3)+
  geom_point(data= tab_pictures, aes(x=treatment, 
                                     y=gray_sum,
                                     fill=treatment), pch = 21, size= 0.3, stroke = 0.2, position = position_jitterdodge())+
  geom_segment(data= df_segments, aes(x = x, xend = xend, y = y, yend = yend),size= 0.4, color = "black") + # Add a segment)
  geom_text(data= df_segments_2, aes(x = x, y=y, label=label), size= 2.5)+
  facet_grid(~specie)+
  theme_classic()+
  theme(legend.position= "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks=element_line(size=0.1),
        panel.background = element_rect(colour = "black", linewidth= 0.2),
        axis.line = element_line(colour = "black", linewidth= 0.2),
        strip.text.x = element_text(size = a, face = "italic"),
        strip.text.y = element_text(size=b),
        axis.text = element_text(size=c),
        axis.title.y = element_text(size=b),
        axis.title.x = element_text(size=b),
        plot.margin=unit(c(-0.1,0.5,-0,0.5), "cm"))+  #top, right, bottom, left
  scale_fill_manual("", values = c("#E7B800", "#008cFF", "#0027A8")) +
  ylab(expression(atop(~Bleaching~intensity)))+
  labs(x='Treatment')+
  scale_y_continuous(limits = c(0,255))+
  scale_x_discrete(labels=c("Control", "MP", "MP+HF"))

plot(plot_pictures_heat)

ggsave("Output/plot_pictures.png", plot =plot_pictures_heat,
       scale = 1, width = 6.3, height = 3, units = c("cm"),
       dpi = 600, limitsize = FALSE) 

