################################################################################
# TITLE: CODE FOR ANALYSIS OF 3D Scap data
# PURPOSE: FPA
# AUTHOR: JOE HEFNER & Micayla Spiros & Savannah Holcombe
# LAST UPDATED ON 2024/06/18 (SRH)
# UPDATES: 
# 1. 2019/06/14: CREATED INITIAL CODE (jth)
# 2. 2019/06/14: ADDED gpa, link file, morphologika reader, and comp code (jth)
# 3. 2020/10/16: fixed rgl plot issue and added new code for analysis (jth and srh)
# 4. 2020/10/21: fixed the tangent space plots and general discussion (jth,srh,mcs)
# 5. 2020/12/02: jth tackled the analysis...we shall see.
# 6. 2020/12/02: jth, mcs, sh meeting
# 7. 2020/12/17: PC issue addressed and outliers identified and subsequently removed 
#                (new file created, old file moved to legacy folder)
# 8. 2021/01/29: CREATED INITIAL CODE (MCS)
# 9. 2021/01/29: Create LDA (SEX~SCAPHT_L+SCAPBR_L)
# 10. 2021/01/29: Create LDA from principal components (GM)
# 11. 2021/01/29: MCC scores for each model
# 12. 2024/02/05: Linear summary stats (for pub)
# 13. 2024/02/13: scaling and centering added for PCA analysis to match LC results
# 14. 2024/05/07: consolidating scap code versions...
# 15. 2024/06/11: Pulling figures for publication (srh)
# 16. 2024/06/19: Updating classification rates (srh)
# NOTES: 
################################################################################
# START WITH A CLEAN WORKING DIRECTORY
rm(list=ls())
################################################################################
# OPEN USEFUL LIBRARIES 
library(geomorph)
library(Morpho)
library(abind)
library(curl)
library(knitr)
library(rgl)#Installed with geomorph
library(MASS) #lda for stepwise
library(ggplot2) #plots
library(mltools) #mcc
library(klaR) # for stepwise selection
library(psych)
library(gridExtra)
library(wesanderson)
library(scales)
################################################################################
getwd()

alinks <-read.table(file.choose(),sep=",") #we have 19 landmarks 
alinks<-as.matrix(alinks)

filelist <- list.files(pattern = ".txt")#This will search the working directory (default) for the files containing .tps
data<-read.morphologika(filelist)
classifier<- read.csv("groups.csv", header=T)
################################################################################
# initial analysis
plotAllSpecimens(data$coords,links=alinks, label=TRUE) 
# select subsets
group<-factor(paste(data$labels)) #levels(group) 

coords.raw<-coords.subset(data$coords,group) # sep male and female data

t<-coords.raw$Female #just female data
s<-coords.raw$Male #just male data

# now GPA on male and female (separate)
# generalized procrustes of female data to scale, translate, and rotate to same shape space
GPA.female<-gpagen(t)
GPA.female
summary(GPA.female)
plotAllSpecimens(GPA.female$coords, links=alinks)

female <- gm.prcomp(GPA.female$coords)
summary(female) # importance of principal components for females
loadings <- female$rotation
print(loadings) # extract the individual PC loadings (coefficients) for each variable (how much
                # each variable contributes to the individual PCs)

# generalized procrustes of male data to scale, translate, and rotate to same shape space
GPA.male<-gpagen(s)
GPA.male
summary(GPA.male)
gpagen(GPA.male$coords)

male <- gm.prcomp(GPA.male$coords) # importance of PCs
summary(male)
loadings <- male$rotation
print(loadings) # extract the individual PC loadings (coefficients) for each variable (how much
# each variable contributes to the individual PCs)

plotAllSpecimens(GPA.male$coords, links=alinks)

#now GPA for MEAN config (all)
GPA.all<-gpagen(data$coords)
GPA.all

summary(GPA.all)
gpagen(GPA.all$coords)


all <- gm.prcomp(GPA.all$coords)
summary(all) # importance of PCs
loadings <- all$rotation
print(loadings) # extract the individual PC loadings (coefficients) for each variable (how much
# each variable contributes to the individual PCs)
plot(all$d)



screeplot(all, npcs=50, main = "", xlab = "Scree Plot")


plotAllSpecimens(GPA.all$coords, links=alinks)

################################################################################
# mean shape for total sample
GPA<-gpagen(data$coords)
meanhead<-mshape(GPA$coords) # mean for all
plot(meanhead, links=alinks)

# mean shape for female
meanhead.female<-mshape(t) # mean for female
plot(meanhead.female, links=alinks)

# mean shape for male
meanhead.male<-mshape(s) # mean for male
plot<-plot(meanhead.male, links=alinks, color="red")

#rglwidget()

################################################################################
# comparisons of male to female

plotRefToTarget(M1=meanhead.male,M2=meanhead.female,mag=1.5,gridPars=GP1,
                label=F, method="points", axes=TRUE,links=alinks) 


rgl.snapshot(filename="comparison-plot_male_female.png")
#rglwidget()

################################################################################
#PCA analysis with plots (using the picknplot function)
# Multidimensional Trait Analysis
# basic analysis

Y.gpa <- gpagen(data$coords)
scap.pca <- gm.prcomp(Y.gpa$coords, scale=TRUE)
test<-scap.pca$x
write.csv(test,"pca.csv")
pca.plot<-cbind(classifier,test)

# Plot the the first 2 PCAs

plot(pca.plot$Comp1, pca.plot$Comp2, col = ifelse(pca.data$Sex == "Male","blue","red"), 
     pch = ifelse(pca.data$Sex == "Male",21,18), xlab = "PC1", ylab = "PC2") +
    legend("topleft", legend = c("Male", "Female"), col = c("blue", "red"), pch = c(21,18))

grand_budapest2_colors <- wes_palette("GrandBudapest2")

pca.figure<-ggplot(pca.plot,aes(pca.plot$Comp1,pca.plot$Comp2,color=Sex,shape=Sex))+
  geom_point(size=2)+
  xlab("PC1")+
  ylab("PC2")+
  theme_minimal()+ 
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = "black"),
        legend.background = element_rect(color = "black", fill = NA),
        axis.text = element_text(size = 13),  # Size for axis text (numbers)
        axis.title = element_text(size = 14),  # Size for axis labels (titles)
        legend.text = element_text(size = 12),  # Size for legend text
        legend.title = element_text(size = 14)) + # Size for legend title)
  scale_color_manual(values = c("steelblue4", "red3"))+
  stat_ellipse(type = "t", level = 0.90, linewidth=.75)  # 't' uses t-distribution, level for 95% confidence


pca.figure

#pick and plot for shape analysis
scap.pca.plot <- plot(scap.pca)
picknplot.shape(scap.pca.plot) # "Error in do.call(plot, x$plot.args) : second argument must be a list"

################################################################################
#possible plots
ggplot(pca.plot, aes(x=Comp2, y=Comp1, colour=Sex)) +
  geom_point(alpha=0.7)+
  geom_hline(yintercept=0, col="gray")+
  geom_vline(xintercept=0, col="gray")+
  theme_minimal()+
  theme(legend.position="bottom")+
  scale_color_manual(values=c('#E69F00','#56B4E9'))+
  stat_ellipse(type = "t", level = 0.99)  # 't' uses t-distribution, level for 95% confidence

################################################################################
#Open dataset
scap.data <- read.csv("Scap-Linear-new (FDB).csv")
scap_linear<-scap.data

#summary stats for linear data: Updated on 1/31/2024 by SRH
scap_linear$SEX<-as.factor(scap_linear$SEX)
scap_linear$ANCESTRY<-as.factor(scap_linear$ANCESTRY)

summary(scap_linear)

describeBy(scap_linear$SCAPHT_L, scap_linear$SEX)
describeBy(scap_linear$SCAPBR_L, scap_linear$SEX)
################################################################################
# for plots
# plots for GM data Age and Sex distribution
GM_Dem<-read.csv("GM_Dem.csv",head=TRUE)
# for plots
# Bar chart of average Age by Sex for GM data
GM_Dem$AgeGroup <- cut(GM_Dem$Age, breaks = seq(0, 100, by = 10), right = FALSE, 
                       labels = paste(seq(10, 100, by = 10) - 10, seq(10, 100, by = 10) - 1, sep = "-"))

bar_plot_GM <- ggplot(GM_Dem, aes(x = AgeGroup, fill = Sex)) +
  geom_bar(position="dodge", alpha=1) +
  scale_fill_manual(values = c("steelblue4", "red3")) +
  labs(title = NULL, x = "Age", y = "Count") +
  theme_minimal()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())
bar_plot_GM
# Violin plot of Age by PopAffin for GM data
violin_plot_GM <- ggplot(GM_Dem, aes(x = PopAffin, y = Age, fill = PopAffin)) +
  geom_violin(alpha=.8)+
  scale_fill_manual(values = c("red2", "steelblue3")) +
  labs(title = NULL, x = "", y = "Age") +
  theme_minimal()+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())
violin_plot_GM

# Bar chart of average age by sex for linear data
scap_linear$AgeGroup <- cut(scap_linear$AGE, breaks = seq(0, 100, by = 10), right = FALSE, 
                     labels = paste(seq(10, 100, by = 10) - 10, seq(10, 100, by = 10) - 1, sep = "-"))

bar_plot_linear <- ggplot(scap_linear, aes(x = AgeGroup, fill = SEX)) +
  geom_bar(position="dodge", alpha=1) +
  scale_fill_manual(values = c("steelblue4", "red3")) +
  labs(title = NULL, x = "", y = "Count") +
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())
bar_plot_linear
# Violin plot of age by population affinity for linear data
violin_plot_linear <- ggplot(scap_linear, aes(x = ANCESTRY, y = AGE, fill = ANCESTRY)) +
  geom_violin(alpha=.8)+
  scale_fill_manual(values = c("red2", "steelblue3")) +
  labs(title = NULL, x = "", y = "Age") +
  theme_minimal()+
  theme(legend.position = "none",
        panel.grid = element_blank())
violin_plot_linear
grid.arrange(bar_plot_linear, violin_plot_linear,bar_plot_GM, violin_plot_GM, nrow=2, ncol = 2)
#############################################################
#1/29/2020: Linear Discriminant Function Analysis- Linear Data
################################################################################

#Dependent variable: Sex
#Independent (predictor) variables: Scap Length, Scap Height

#Dependent (Sex) Groups-1 (i.e. male, female)

#check for missing data
missing<-is.na(scap.data)
summary(missing)

#save discriminant analysis as its own object
#dependent~then list variables of interest
scaplda<-lda(SEX~SCAPHT_L+SCAPBR_L, scap.data)


scaplda

p1<-predict(scaplda) # predict sex based on height and breadth
names(p1) #class (Sex), posterior (posterior probability), x is the linear discriminants
#for group membership (pp is the probability that the individual is correctly classified)

scaplda # priors, coefficients, group means

plot(scaplda) #separation of sexes based on the linear discriminants
  # most of the variability is on the first axis

qplot(p1$x[,1], color= SEX, data=scap.data) +
  labs(x = "Linear Discriminant 1", y = "Density", 
       title = NULL) +
  theme_minimal()

ggplot(scap.data, aes(x = LD1, fill = SEX)) +
  geom_density(alpha=.7) +
  scale_fill_manual(values = c("blue", "red"))+
  labs(x = "Linear Discriminant 1", y = "Density", 
       title = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

#Calculate Group Centroids per individual
LD1<- predict(scaplda)$x[,1]

LD1


#Classification tables
#table comparing actual classification vs. predicted classification
ct1<- table(scap.data$SEX, p1$class)
ct1

#    F    M
#F  390   44
#M  42   776


#misclasification of 44 females into male group
#misclasification of 42 males into female groups

#create proportion table for % CCR
diag(prop.table(ct1,1))

#F             M 
#0.8986175     0.9486553

#How well the overall MODEL is working (not for each group)
sum(diag(prop.table(ct1)))




#######################################################################
#1/29/2020: Linear Discriminant Function Analysis- Prin. Comp. GM Data
#######################################################################
#Open dataset
pca.data <- read.csv("pca_data.csv")
summary(pca.data)

#check for missing data
missing.1<-is.na(pca.data)
summary(missing.1)

#save discriminant analysis as its own object
#dependent~then list variables of interest


GMlda_step<-stepclass(Sex~PCA1+PCA2+PCA3+PCA4+PCA5+PCA6+PCA7+PCA8+PCA9+PCA10+PCA11+PCA12+PCA13+PCA14+PCA15+PCA16+
                        PCA17+PCA18+PCA19+PCA20+PCA21+PCA22+PCA23+PCA24+PCA25, pca.data, method="lda", direction="backward")

GMlda_step
  
GMlda<-lda(Sex ~ PCA1 + PCA2 + PCA3 + PCA4 + PCA5 + PCA6 + PCA7 + PCA9 + 
             PCA10 + PCA11 + PCA12 + PCA13 + PCA14 + PCA15 + PCA17 + PCA18 + 
             PCA19 + PCA20 + PCA21 + PCA22 + PCA24 + PCA25, pca.data, prior=c(1/2,1/2))


GMlda



p2<-predict(GMlda)
names(p2) #class (Sex), posterior (posterior probability)
#for group membership (pp is the probability that 
#the individual is correctly classified into the 
#group it was classified into)


plot(GMlda) #most of the variability is on the first axis





#Classification tables
#table comparing actual classification vs. predicted classification
ct2<- table(pca.data$Sex, p2$class)
ct2

#create proportion table for % CCR
diag(prop.table(ct2,1))

#F             M 
#0.9285714  0.8593750 

#How well the overall MODEL is working (not for each group)
sum(diag(prop.table(ct2)))

#################################################
#MCC Matthews correlation- Compare the  models
#################################################
#linear data
mcc.1<-mcc(confusionM = matrix(c(390,42,44,776), nrow = 2))
mcc.1
#GM data
mcc.2<-mcc(confusionM = matrix(c(39,9,3,55), nrow = 2))
mcc.2
################################################################################
# end of program
################################################################################
