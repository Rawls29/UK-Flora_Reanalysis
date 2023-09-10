rm(list=ls())
setwd('C:/Users/samra/OneDrive/Documents/Uni/Edinburgh Honours Project Write-Up')

####Libraries####
library(car)
library(beeswarm)
library(BIFloraExplorer)
library(ggplot2)
library(MASS)

####Loading Datasets####
flora_all <- readRDS('new_data_new_fert_mode.rds')
flora_all <- merge(flora_all, BI_main, by.x="main_species_name", by.y="taxon_name_binom")

flora <- flora_all[complete.cases(flora_all$main_species_name),]
flora <- flora[complete.cases(flora$myFert3),]

flora$total_hectads_2010_2019 <- flora$CI_hectads_2010_2019+flora$Ire_hectads_2010_2019+flora$GB_Man_hectads_2010_2019
flora$total_hectads_1987_1999 <- flora$CI_hectads_1987_1999+flora$Ire_hectads_1987_1999+flora$GB_Man_hectads_1987_1999
flora$total_hectads_2000_2009 <- flora$CI_hectads_2000_2009+flora$Ire_hectads_2000_2009+flora$GB_Man_hectads_2000_2009
flora$total_hectads_post2000 <- flora$CI_hectads_post2000+flora$Ire_hectads_post2000+flora$GB_Man_hectads_post2000
flora$total_hectads_increase <- ((flora$total_hectads_2010_2019-flora$total_hectads_1987_1999)/flora$total_hectads_2010_2019)*100

genome <- flora[complete.cases(flora$'2C DNA content (pg)'),]

for(i in 1:nrow(genome)){
value <-genome$`2C DNA content (pg)`[i]
value2 <- strsplit(value, ", ")
value2 <- as.numeric(unlist(value2))
mean2C <- mean(value2)
genome$'2C_DNA_numeric'[i]<-mean2C}

genome$`2C_DNA_numeric` <- as.numeric(genome$`2C_DNA_numeric`)

range <- flora[complete.cases(flora$GB_Man_hectads_2010_2019),]
range <- flora[which(flora$GB_Man_hectads_1987_1999!=0),]
range <- flora[which(flora$GB_Man_hectads_2000_2009!=0),]
range <- flora[which(flora$GB_Man_hectads_2010_2019!=0),]
range$range_increase <- ((range$GB_Man_hectads_2010_2019-range$GB_Man_hectads_1987_1999)/range$GB_Man_hectads_2010_2019)*100

####Genome Complete Analysis####
#Analysis
shapiro.test(genome$'2C_DNA_numeric') #p<0.05 therefore non-normal
shapiro.test(log(genome$'2C_DNA_numeric')) #p<0.05 therefore non-normal
leveneTest(genome$'2C_DNA_numeric'~myFert3, data=genome) #p>0.05 so variances are homogenous
kruskal.test(genome$'2C_DNA_numeric'~myFert3, data=genome) #p>0.05 therefore median genome size not sig. different
#Beeswarm Plot
beeswarm(log(genome$'2C_DNA_numeric')~myFert3, data=genome)
boxplot(log(genome$'2C_DNA_numeric')~myFert3, data=genome, add=T, 
        col="#0000ff22")

####Genome Family Specific Analysis####
genome_family <- genome[complete.cases(genome$family),]
family_freq <- as.data.frame(table(genome_family$family))
large_fams <- family_freq[which(family_freq$Freq>=15),]

for(i in 1:nrow(large_fams)){
  fam <- genome_family[which(genome_family$family==large_fams[i,1]),]
  tab <- as.data.frame(table(fam$myFert3))
  if(nrow(tab)==3){
    if(tab[(which(tab$Var1=="outcrossing")),2]>3&tab[(which(tab$Var1=="selfing")),2]>3){
    if(shapiro.test(fam$`2C_DNA_numeric`)$p.value>0.05){
      fam_anova <- aov(fam$`2C_DNA_numeric`~myFert3, data=fam)
      print(unique(fam$family))
      print("non-transformed")
      print(summary(fam_anova))
      beeswarm(fam$`2C_DNA_numeric`~myFert3, data=fam)
      boxplot(fam$`2C_DNA_numeric`~myFert3, data=fam, add=T, 
              col="#0000ff22")} else{
        if(shapiro.test(log(fam$`2C_DNA_numeric`))$p.value>0.05){
          fam_anova <- aov(log(fam$`2C_DNA_numeric`)~myFert3, data=fam)
          print(unique(fam$family))
          print("log transformed")
          print(summary(fam_anova))
          beeswarm(log(fam$`2C_DNA_numeric`)~myFert3, data=fam)
          boxplot(log(fam$`2C_DNA_numeric`)~myFert3, data=fam, add=T, 
                  col="#0000ff22")} else {print(unique(fam$family))
            print("non-normal")
            if(leveneTest(fam$`2C_DNA_numeric`~myFert3, data=fam)$`Pr(>F)`[1]>0.05){print(kruskal.test(fam$`2C_DNA_numeric`~myFert3, data=fam))}else{
              if(leveneTest(log(fam$`2C_DNA_numeric`)~myFert3, data=fam)$`Pr(>F)`[1]>0.05){print(kruskal.test(log(fam$`2C_DNA_numeric`)~myFert3, data=fam))}else{print("Non-Homogenous Variances")}
            }
            beeswarm(log(fam$`2C_DNA_numeric`)~myFert3, data=fam)
            boxplot(log(fam$`2C_DNA_numeric`)~myFert3, data=fam, add=T, 
                    col="#0000ff22")}}
    }else{print(unique(fam$family))
      print("too few selfers/outcrossers")}
}}

####Genome Pair Analysis####
self_cross_pairs_genome <- data.frame(self_species=c("Briza minor",
                                                     "Callitriche brutia",
                                                     "Cephalanthera damasonium",
                                                     "Juncus squarrosus",
                                                     "Luzula campestris",
                                                     #"Papaver dubium",
                                                     #"Senecio vulgaris",
                                                     #"Trifolium dubium",
                                                     #"Trifolium fragiferum",
                                                     #"Trifolium glomeratum",
                                                     "Veronica polita"),
                                      cross_species=c("Briza media",
                                                      "Callitriche stagnalis",
                                                      "Cephalanthera longifolia",
                                                      "Juncus bufonius",
                                                      "Luzula forsteri",
                                                      #"Papaver rhoeas",
                                                      #"Senecio squalidus",
                                                      #"Trifolium campestre",
                                                      #"Trifolium hybridum",
                                                      #"Trifolium repens",
                                                      "Veronia chamaedrys"),
                                      genus=c("Briza",
                                              "Callitriche",
                                              "Cephalanthera",
                                              "Juncus",
                                              "Luzula",
                                              #"Papaver",
                                              #"Senecio",
                                              #"Trifolium",
                                              #"Trifolium",
                                              #"Trifolium",
                                              "Veronica"),
                                      self_2c=c(mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Briza minor")]),
                                                mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Callitriche brutia")]),
                                                mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Cephalanthera damasonium")]),
                                                mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Juncus squarrosus")]),
                                                mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Luzula campestris")]),
                                                #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Papaver dubium")]),
                                                #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Senecio vulgaris")]),
                                                #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Trifolium dubium")]),
                                                #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Trifolium fragiferum")]),
                                                #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Trifolium glomeratum")]),
                                                mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Veronica polita")])),
                                      cross_2c=c(mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Briza media")]),
                                                 mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Callitriche stagnalis")]),
                                                 mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Cephalanthera longifolia")]),
                                                 mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Juncus bufonius")]),
                                                 mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Luzula forsteri")]),
                                                 #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Papaver rhoeas")]),
                                                 #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Senecio squalidus")]),
                                                 #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Trifolium campestre")]),
                                                 #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Trifolium hybridum")]),
                                                 #mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Trifolium repens")]),
                                                 mean(genome$'2C_DNA_numeric'[which(genome$main_species_name=="Veronica chamaedrys")])),
                                      pair=c("B. minor - B. media",
                                             "C. brutia - C. stagnalis",
                                             "C. damasonium - C. longifolia",
                                             "J. squarrosus - J. bufonius",
                                             "L. campestris - L. forsteri",
                                             #"P. dubium - P. rhoeas",
                                             #"S. vulgaris - S. squalidus",
                                             #"T. dubium - T. campestre",
                                             #"T. fragigerum - T. hybridum",
                                             #"T. glomeratum - T. repens",
                                             "V. polita - V. chamaedrys"))
self_cross_pairs_genome$difference_2c <- self_cross_pairs_genome$self_2c-self_cross_pairs_genome$cross_2c

library(ggplot2)
ggplot(data= self_cross_pairs_genome, aes(x=pair, y=difference_2c, fill=genus))+
  geom_bar(stat="identity")+
  xlab("")+
  ylab("Difference in 2C DNA Content (pg)")+
  labs(fill="Genus")+
  theme(axis.text.x=element_text(angle=90, face="italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.text=element_text(face="italic"))

t.test(self_cross_pairs_genome$difference_2c, mu=0)
#p>0.05, so mean difference between pairs doesn't differ sig. from 0

####Genome Family Difference Plot####
fam_nums <- as.data.frame(table(genome$family))

mean_fam_2c_cs <- data.frame(matrix(NA,
                                    nrow=0,
                                    ncol=3))
colnames(mean_fam_2c_cs) <- c("family", "cross_av_2c", "self_av_2c")
for(i in 1:nrow(fam_nums)){
  if(fam_nums$Freq[i]>1){
    fam <- fam_nums$Var1[i]
    temp <- as.data.frame(table(genome$myFert3[which(genome$family==fam_nums$Var1[i])]))
    if(nrow(as.data.frame(temp$Freq[which(temp$Var1=="outcrossing")]))!=0&
       nrow(as.data.frame(temp$Freq[which(temp$Var1=="selfing")]))!=0){
      if(temp$Freq[which(temp$Var1=="outcrossing")]>=1&
         temp$Freq[which(temp$Var1=="selfing")]>=1){
        family <- unique(genome$family[which(genome$family==fam_nums$Var1[i])])
        cross_av_2c <- mean(genome$`2C_DNA_numeric`[which(genome$family==fam&genome$myFert3=="outcrossing")])
        self_av_2c <- mean(genome$`2C_DNA_numeric`[which(genome$family==fam&genome$myFert3=="selfing")])
        temp2 <- data.frame(t(as.data.frame(c(family, cross_av_2c, self_av_2c))))
        colnames(temp2) <- c("family", "cross_av_2c", "self_av_2c")
        mean_fam_2c_cs <- rbind(mean_fam_2c_cs, temp2)
      }}
  }
}

mean_fam_2c_cs$self_av_2c <- as.numeric(mean_fam_2c_cs$self_av_2c)
mean_fam_2c_cs$cross_av_2c <- as.numeric(mean_fam_2c_cs$cross_av_2c)
mean_fam_2c_cs$diff <- NA
for(i in 1:nrow(mean_fam_2c_cs)){
  mean_fam_2c_cs$diff <- mean_fam_2c_cs$self_av_2c - mean_fam_2c_cs$cross_av_2c
}

dev.off()
ggplot(data=mean_fam_2c_cs, aes(x=family, y=diff, 
                                fill=family)) +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.text=element_text(face="italic")) +
  geom_bar(stat = "identity") +
  xlab("") + ylab("Difference in Genome Size")+
  scale_fill_hue(name="Family")

####Range Outliers####
range_endemic_sorbus <- range[which(range$main_species_name=="Sorbus arranensis"|
                                      range$main_species_name=="Sorbus porrigentiformis"|
                                      range$main_species_name=="Sorbus devonensis")]
range_montanes <- range[which(range$Biome=="Boreal montane"|range$Biome=="Arctic Montane"|range$Biome=="Boreo-Arctic Montane"),]
range_coastals
range_aquatics <- range[which(range$'Life-form'=="hydrophyte, therophyte"|
                                range$'Life-form'=="hydrophyte"|
                                range$'Life-form'=="helophyte, hydrophyte"|
                                range$'Life-form'=="helophyte, hemicryptophyte, hydrophyte"|
                                range$'Life-form'=="helophyte, hydrophyte, therophyte" |
                                range$'Life-form'=="hemicryptophyte, hydrophyte" ),]

range_reduced <- range[-c(which(range$main_species_name%in%range_endemic_sorbus$main_species_name|
                                  #range$main_species_name%in%range_coastals$main_species_names|
                                  range$main_species_name%in%range_montanes$main_species_name|
                                  range$main_species_name%in%range_aquatics$main_species_name)),]

#Montanes
shapiro.test(range_montanes$total_hectads_2010_2019) #p<0.05 therefore non-normal
shapiro.test(log(range_montanes$total_hectads_2010_2019)) #p<0.05 therefore non-normal
leveneTest(total_hectads_2010_2019~myFert3, data=range_montanes) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_2019~myFert3, data=range_montanes) #p>0.05 therefore median genome size not different
#pairwise.wilcox.test(range$total_hectads_2010_2019, range$myFert3, p.adjust.method = "BH")
#Beeswarm Plot
beeswarm(total_hectads_2010_2019~myFert3, data=range_montanes)
boxplot(total_hectads_2010_2019~myFert3, data=range_montanes, add=T, 
        col="#0000ff22")

#Coastals

#Aquatics
shapiro.test(range_aquatics$total_hectads_2010_2019) #p<0.05 therefore non-normal
shapiro.test(log(range_aquatics$total_hectads_2010_2019)) #p<0.05 therefore non-normal
leveneTest(total_hectads_2010_2019~myFert3, data=range_aquatics) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_2019~myFert3, data=range_aquatics) #p>0.05 therefore median genome size not different
#pairwise.wilcox.test(range$total_hectads_2010_2019, range$myFert3, p.adjust.method = "BH")
#Beeswarm Plot
beeswarm(total_hectads_2010_2019~myFert3, data=range_aquatics)
boxplot(total_hectads_2010_2019~myFert3, data=range_aquatics, add=T, 
        col="#0000ff22")

#Others

#Beeswarm Plot
beeswarm(total_hectads_2010_2019~myFert3, data=range_reduced)
boxplot(total_hectads_2010_2019~myFert3, data=range_reduced, add=T, 
        col="#0000ff22")

#Annual vs Perennial

#Diploids vs Polyploids

#Natives and Non-Natives
range_natives <- range_reduced[which(range_reduced$StaceIV_nativity=="N"|
                                       range_reduced$StaceIV_nativity=="Arch-colonist"|
                                       range_reduced$StaceIV_nativity=="Arch-denizen"),]
shapiro.test(range_natives$total_hectads_2010_2019) #p<0.05 therefore non-normal
shapiro.test(log(range_natives$total_hectads_2010_2019)) #p<0.05 therefore non-normal
leveneTest(total_hectads_2010_2019~myFert3, data=range_natives) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_2019~myFert3, data=range_natives) #p>0.05 therefore median genome size not different
#pairwise.wilcox.test(range$total_hectads_2010_2019, range$myFert3, p.adjust.method = "BH")
#Beeswarm Plot
beeswarm(total_hectads_2010_2019~myFert3, data=range_natives)
boxplot(total_hectads_2010_2019~myFert3, data=range_natives, add=T, 
        col="#0000ff22")

range_non_natives <- range_reduced[which(range_reduced$StaceIV_nativity=="Neo-casual"|
                                       range_reduced$StaceIV_nativity=="Neo-natd"|
                                       range_reduced$StaceIV_nativity=="Neo-surv"|
                                       range_reduced$StaceIV_nativity=="Neonative"),]
shapiro.test(range_non_natives$total_hectads_2010_2019) #p<0.05 therefore non-normal
shapiro.test(log(range_non_natives$total_hectads_2010_2019)) #p<0.05 therefore non-normal
leveneTest(total_hectads_2010_2019~myFert3, data=range_non_natives) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_2019~myFert3, data=range_non_natives) #p>0.05 therefore median genome size not different
#pairwise.wilcox.test(range$total_hectads_2010_2019, range$myFert3, p.adjust.method = "BH")
#Beeswarm Plot
beeswarm(total_hectads_2010_2019~myFert3, data=range_non_natives)
boxplot(total_hectads_2010_2019~myFert3, data=range_non_natives, add=T, 
        col="#0000ff22")

####Range Complete Analysis####
shapiro.test(range$total_hectads_2010_2019) #p<0.05 therefore non-normal
shapiro.test(log(range$total_hectads_2010_2019)) #p<0.05 therefore non-normal
leveneTest(total_hectads_2010_2019~myFert3, data=range) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_2010_2019~myFert3, data=range) #p>0.05 therefore median genome size not different
#pairwise.wilcox.test(range$total_hectads_2010_2019, range$myFert3, p.adjust.method = "BH")
#Beeswarm Plot
beeswarm(total_hectads_2010_2019~myFert3, data=range)
boxplot(total_hectads_2010_2019~myFert3, data=range, add=T, 
        col="#0000ff22")

####Range Family Specific Analysis####
range_family <- range[complete.cases(range$family),]
family_freq_r <- as.data.frame(table(range_family$family))
large_fams_r <- family_freq_r[which(family_freq_r$Freq>=15),]

for(i in 1:nrow(large_fams_r)){
  fam <- range_family[which(range_family$family==large_fams_r[i,1]),]
  tab <- as.data.frame(table(fam$myFert3))
  if(nrow(tab)==3){
    if(tab[(which(tab$Var1=="outcrossing")),2]>3&tab[(which(tab$Var1=="selfing")),2]>3){
      if(shapiro.test(fam$total_hectads_2010_2019)$p.value>0.05){
        fam_anova <- aov(fam$total_hectads_2010_2019~myFert3, data=fam)
        print(unique(fam$family))
        print("non-transformed")
        print(summary(fam_anova))
        beeswarm(fam$total_hectads_2010_2019~myFert3, data=fam)
        boxplot(fam$total_hectads_2010_2019~myFert3, data=fam, add=T, 
                col="#0000ff22")} else{
                  if(shapiro.test(log(fam$total_hectads_2010_2019))$p.value>0.05){
                    fam_anova <- aov(log(fam$total_hectads_2010_2019)~myFert3, data=fam)
                    print(unique(fam$family))
                    print("log transformed")
                    print(summary(fam_anova))
                    beeswarm(log(fam$total_hectads_2010_2019)~myFert3, data=fam)
                    boxplot(log(fam$total_hectads_2010_2019)~myFert3, data=fam, add=T, 
                            col="#0000ff22")} else {print(unique(fam$family))
                              print("non-normal")
                              if(leveneTest(fam$total_hectads_2010_2019~myFert3, data=fam)$`Pr(>F)`[1]>0.05){print(kruskal.test(fam$GB_Man_hectads_2010_2019~myFert3, data=fam))}else{
                                if(leveneTest(log(fam$total_hectads_2010_2019)~myFert3, data=fam)$`Pr(>F)`[1]>0.05){print(kruskal.test(log(fam$GB_Man_hectads_2010_2019)~myFert3, data=fam))}else{print("Non-Homogenous Variances")}
                              }
                              beeswarm(log(fam$total_hectads_2010_2019)~myFert3, data=fam)
                              boxplot(log(fam$total_hectads_2010_2019)~myFert3, data=fam, add=T, 
                                      col="#0000ff22")}}
    }else{print(unique(fam$family))
      print("too few selfers/outcrossers")}
  }}


####Range PGLS Analysis####

###Range Change Complete Analysis####
shapiro.test(range$total_hectads_increase) #p<0.05 therefore non-normal
leveneTest(total_hectads_increase~myFert3, data=range) #p>0.05 so variances are homogenous
kruskal.test(total_hectads_increase~myFert3, data=range) #p>0.05 therefore median genome size not different
#pairwise.wilcox.test(range$GB_Man_hectads_2010_2019, range$myFert3, p.adjust.method = "BH")
#Beeswarm Plot
beeswarm(total_hectads_increase~myFert3, data=range)
boxplot(total_hectads_increase~myFert3, data=range, add=T, 
        col="#0000ff22")

###Range Change Family Specific Analysis####
range_family <- range[complete.cases(range$family),]
family_freq_r <- as.data.frame(table(range_family$family))
large_fams_r <- family_freq_r[which(family_freq_r$Freq>=15),]

for(i in 1:nrow(large_fams_r)){
  fam <- range_family[which(range_family$family==large_fams_r[i,1]),]
  tab <- as.data.frame(table(fam$myFert3))
  if(nrow(tab)==3){
    if(tab[(which(tab$Var1=="outcrossing")),2]>3&tab[(which(tab$Var1=="selfing")),2]>3){
      if(shapiro.test(fam$total_hectads_increase)$p.value>0.05){
        fam_anova <- aov(fam$total_hectads_increase~myFert3, data=fam)
        print(unique(fam$family))
        print("non-transformed")
        print(summary(fam_anova))
        beeswarm(fam$total_hectads_increase~myFert3, data=fam)
        boxplot(fam$total_hectads_increase~myFert3, data=fam, add=T, 
                col="#0000ff22")} else {print(unique(fam$family))
                              print("non-normal")
                              if(leveneTest(fam$total_hectads_increase~myFert3, data=fam)$`Pr(>F)`[1]>0.05){print(kruskal.test(fam$GB_Man_hectads_2010_2019~myFert3, data=fam))}else{print("Non-Homogenous Variances")}
                              beeswarm(fam$total_hectads_increase~myFert3, data=fam)
                              boxplot(fam$total_hectads_increase~myFert3, data=fam, add=T, 
                                      col="#0000ff22")}}
    }else{print(unique(fam$family))
      print("too few selfers/outcrossers")}
  }

###Range Change Outliers####

###Continent Presence####
for(i in 1:nrow(flora)){
  if(is.na(flora$`Range: 3. continents where native`[i])==FALSE){
    continents <- strsplit(flora$`Range: 3. continents where native`[i], ", ")
    flora$regions_native[i] <- length(continents[[1]])
    }
}

bc<-boxcox(flora$total_hectads_2010_2019~flora$regions_native)
(lambda <- bc$x[which.max(bc$y)])
model <- lm(((flora$total_hectads_2010_2019^lambda-1)/lambda)~flora$regions_native)
par(mfrow=c(2,2))
plot(model)
par(mfrow=c(1,1))
summary(model)

beeswarm(((flora$total_hectads_2010_2019^lambda-1)/lambda)~regions_native, data=flora)
boxplot(((flora$total_hectads_2010_2019^lambda-1)/lambda)~regions_native, data=flora, add=T, 
        col="#0000ff22")
abline(lm(((flora$total_hectads_2010_2019^lambda-1)/lambda)~flora$regions_native), data=flora, add=T)
