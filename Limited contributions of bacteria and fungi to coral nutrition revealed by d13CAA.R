#####Code for Limited contributions of bacteria and fungi to ###
#####coral nutrition revealed by amino acid ¦Ä13C analysis####### 
#####(Wang et al., 2025)########################################
#####Here are code for LDA and Baysian mixing models############

###Package Required###
library(pacman)
library('knitr')
if (!require("pacman")) install.packages("pacman"); library(pacman) # for rapid install 
pacman::p_load(devtools, ellipse, ggbiplot, vqv, patchwork, graphics, plyr, effects, MASS, tidyverse, dplyr, plotrix, vegan, cowplot, caret, reshape, car, emmeans, lmerTest)
library(ggplot2)
library(ggpubr)
library(multcomp)
library(MixSIAR)
#### Please seperate each table form Supplement DATA 1 into single CSV before running the following code
#### LDA for overall samples (modified from Wall et al., 2021)
#### For Fig 3a 
#### Fraction: source and different host
EAA.raw.overall<-read.csv('OverallCSIA.csv')
EAA.raw.overall
#normalized
EAA.raw.overall$ID<-1:nrow(EAA.raw.overall)
for(i in 1:length(EAA.raw.overall$ID)){
  EAA.raw.overall$Ile.n[i] <- (EAA.raw.overall$Ile[i]-mean(as.numeric(EAA.raw.overall[i,3:8])))
  EAA.raw.overall$Phe.n[i] <- (EAA.raw.overall$Phe[i]-mean(as.numeric(EAA.raw.overall[i,3:8])))
  EAA.raw.overall$Leu.n[i] <- (EAA.raw.overall$Leu[i]-mean(as.numeric(EAA.raw.overall[i,3:8])))
  EAA.raw.overall$Lys.n[i] <- (EAA.raw.overall$Lys[i]-mean(as.numeric(EAA.raw.overall[i,3:8])))
  EAA.raw.overall$Thr.n[i] <- (EAA.raw.overall$Thr[i]-mean(as.numeric(EAA.raw.overall[i,3:8])))
  EAA.raw.overall$Val.n[i] <- (EAA.raw.overall$Val[i]-mean(as.numeric(EAA.raw.overall[i,3:8])))
}
EAA.raw.overall
#create a new table for normalized values
EAA.norm.overall<-EAA.raw.overall[,c(1:2,9:15)] 
EAA.norm.overall
#separate sources and animal
EAA.norm.overall.sources<-EAA.norm.overall[!(EAA.norm.overall$Fraction2=="Animal"),] 
EAA.norm.overall.animal<-EAA.norm.overall[( EAA.norm.overall$Fraction2=="Animal"),]
EAA.norm.overall.sources<-droplevels(EAA.norm.overall.sources)
EAA.norm.overall.animal<-droplevels(EAA.norm.overall.animal)
EAA.norm.overall.sources
EAA.norm.overall.animal
# remove additional information for now
EAA.norm.overall.sources.1<-EAA.norm.overall.sources[,c(-1,-3)]
EAA.norm.overall.animal.1<-EAA.norm.overall.animal[,c(-1,-3)]
EAA.norm.overall.sources.1
#LDA model
LDA.EAA.norm.overall.sources.1 <- lda(Fraction2 ~ Ile.n + Leu.n  +Phe.n +Thr.n+ Val.n+ Lys.n, data = EAA.norm.overall.sources.1, CV = TRUE)
# create a table which compares the classification of the LDA model to the actual spp
ct.prod.norm1 <- table(EAA.norm.overall.sources.1$Fraction2, LDA.EAA.norm.overall.sources.1$class)
ct.prod.norm1
sum(diag(prop.table(ct.prod.norm1)))
# what % of each species is being correctly classified
diag(prop.table(ct.prod.norm1, 1))
#training and predict
training.EAA.norm.samples <- lda(Fraction2 ~ Ile.n + Leu.n + Phe.n + Thr.n+ Val.n+ Lys.n, data = EAA.norm.overall.sources.1)
training.EAA.norm.samples$scaling
datPred.norm.sources <- data.frame(Fraction2=EAA.norm.overall.sources.1$Fraction2,
                                   predict(training.EAA.norm.samples)$x)
datPred.norm.sources$ID<-EAA.norm.overall.sources$ID
datPred.norm.sources$ID
host.norm.res <- predict(training.EAA.norm.samples, EAA.norm.overall.animal)
datPred2.norm.an <- data.frame(Fraction2='Animal', host.norm.res$x)
datPred2.norm.an$ID<-EAA.norm.overall.animal$ID
EAA.norm.overall.animal$class <- host.norm.res$class
EAA.norm.overall.animal$class
datPred.norm.sources
datPred2.norm.an
datPred3.norm.overall <- rbind(datPred.norm.sources, datPred2.norm.an)
datPred3.norm.overall       
#re-attach data by merging with ID
EAA.frac.LD.norm.meta1<-merge(datPred3.norm.overall,EAA.raw.overall, by="ID")
names(EAA.frac.LD.norm.meta1)[names(EAA.frac.LD.norm.meta1)=="Fraction.grouped.x"] <- "Fraction2"
EAA.frac.LD.norm.meta1 
EAA.frac.LD.norm.meta1<-EAA.frac.LD.norm.meta1[,-7]
EAA.frac.LD.norm.meta1
EAA.frac.LD.norm.meta1$Fraction.source2<-EAA.frac.LD.norm.meta1$Fraction2
EAA.frac.LD.norm.meta1
EAA.frac.LD.norm.meta1<-droplevels(EAA.frac.LD.norm.meta1)
EAA.frac.LD.norm.meta1
EAA.frac.LD.norm.meta1$Fraction.source2
EAA.frac.LD.norm.meta1$Fraction.grouped<-factor(EAA.frac.LD.norm.meta1$Fraction2,
                                                levels=c("Animal", "Symbiont", "Bacteria", "Fungi", "Particulate source"))
EAA.frac.LD.norm.meta1
#####shapes are different in the paper 
ggplot(EAA.frac.LD.norm.meta1, aes(LD1, LD2)) +
  geom_point(aes(color = Fraction.grouped,shape=Fraction), size=7)+
  scale_shape_manual(values=c(1,3,15,5,4,0,16,7))+
  scale_color_manual(values = c("#EF476F", "#FFD166", "#AD91BB", "#99BC29", "#26547C", "black"))+
  stat_ellipse(aes(color=Fraction.grouped),type = "norm",level = 0.95)+
  theme_minimal()

##### Overall mixing models using mean-centred ¦Ä13CEAA
##### For Fig.3b
##### Load mix data
# Seperate MN data in SUPPLEMENTARY DATA 1
##### Species code as follow
##### Pdcool and Pdwarm means P. damicornis in cool or warm season
##### Pm= Pocillopora meandrina; Mc=Montipora capitata; Pmns=P.meandrina without paried symbiont
mix <- load_mix_data(filename="overall_consumer_mn.csv",
                     iso_names=c("Leu","Ile","Thr","Val","Lys","Phe"),
                     factors="Species",
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)
# Load source data
source <- load_source_data(filename="overall_sources_mn.csv",
                           source_factors="Species",
                           conc_dep=FALSE,
                           data_type="means",
                           mix)
# Load discrimination/TDF data
discr <- load_discr_data(filename ="overall_discrimination_1.csv",mix)
# Plot prior
plot_prior(alpha.prior=1,source)
# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)
# Run the JAGS model ("test" first, then "long")
jags.1 <- run_model(run="long", mix, source, discr, model_filename,alpha.prior=1)
output_JAGS(jags.1, mix, source)
source$source_names
original <- combine_sources(jags.1 , mix, source, 1,
                            groups=list(Particulate="Particulate source",Symbiont="Symbiont",
                                        Fungi="Fungi",Bacteria="Bacteria"))
plot_intervals(original,toplot="fac1")
####For sponges and bamboo corals
mix <- load_mix_data(filename="BC_consumer_mn.csv",
                     iso_names=c("Leu","Ile","Thr","Val","Lys","Phe"),
                     factors="Species",
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)
# Load source data
source <- load_source_data(filename="BC_sources_mn.csv",
                           source_factors="Species",
                           conc_dep=FALSE,
                           data_type="means",
                           mix)
# Load discrimination/TDF data
discr <- load_discr_data(filename ="BC_discrimination_1.csv",mix)
# Plot prior
plot_prior(alpha.prior=1,source)
# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)
# Run the JAGS model ("test" first, then "normal")
jags.1 <- run_model(run="long", mix,  source, discr, model_filename,alpha.prior=1)
output_JAGS(jags.1, mix, source)
source$source_names
original <- combine_sources(jags.1 , mix, source, 1,
                            groups=list(Particulate="Particulate source",
                                        Fungi="Fungi",Bacteria="Bacteria"))
plot_intervals(original,toplot="fac1")

# Estimate regional sources contribution for P. damicornis
# For Fig4b
# Load mix data
# Seperate raw ¦Ä13CEAA data in SUPPLEMENTARY DATA 1
mix <- load_mix_data(filename="coral_consumer_raw.csv",
                     iso_names=c("Leu","Ile","Thr","Val","Lys","Phe"),
                     factors="Season",
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)
# Load source data
source <- load_source_data(filename="coral_sources.csv",
                           source_factors="Season",
                           conc_dep=FALSE,
                           data_type="means",
                           mix)
# Load discrimination/TDF data
discr <- load_discr_data(filename ="coral_discrimination_1.csv",mix)
# Plot prior
plot_prior(alpha.prior=1,source)
# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)
# Run the JAGS model ("test" first, then "normal"or "Long")
jags.1 <- run_model(run="Long", mix, source, discr, model_filename,alpha.prior=1)
output_JAGS(jags.1, mix, source)
source$source_names
original <- combine_sources(jags.1 , mix, source, 1,
                            groups=list(Heter="Heter",
                                        Symbiont="Symbiont"))

plot_intervals(original,toplot="fac1")


