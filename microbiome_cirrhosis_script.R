library(readxl)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(vegan)

library(tibble)
library(ggpubr)
library(wesanderson)
library(tidyverse)
library(reshape2)
library(gplots)
library(broom)

install.packages("broom") #install 
setwd("~/Google Drive/MSc/2 - Microbiome Analysis/Report/LCReport")

main=read.table('merged_abundance_table.txt',header=TRUE,sep="\t")
# t = separate(main, ID, into = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species","Type"), sep="\\|")
t = separate(main, clade_name, into = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species","Type"), sep="\\|")

t=t[!is.na(t$Species), ]
t=t[!duplicated(t$Species),]

x=1:dim(t)[1]
OTUs= paste("OTU",x)
# t=add_column(t, OTUs = OTUs, .before = "s1")
t=add_column(t, OTUs = OTUs, .after = "Type")


# OTU =as.data.frame(t[,-which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Type"))])
OTU =as.data.frame(t[,-which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Type","NCBI_tax_id"))])
row.names(OTU) = OTU$OTUs
OTU =as.data.frame(OTU[,-which(names(OTU) %in% c("OTUs","Species"))])
CountMatrix = OTU %>% as.matrix()
mode(CountMatrix) <- 'integer'

# TAX =as.data.frame(t[,which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species","OTUs"))])
TAX =as.data.frame(t[,which(names(t) %in% c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species","OTUs","NCBI_tax_id"))])

row.names(TAX) = TAX$OTUs
TAX =as.data.frame(TAX[,-which(names(TAX) %in% c("OTUs"))])
TaxaMatrix <- as.matrix(TAX)

Metadata <- read.table('metadata_reports.txt',header=TRUE,sep="\t")
rownames(Metadata) <- Metadata$SampleID

otuTABLE <- otu_table(CountMatrix, taxa_are_rows = TRUE)
taxTABLE <- tax_table(TaxaMatrix)
sampleDATA <- sample_data(Metadata)

phylo_obj <- phyloseq(otuTABLE, taxTABLE, sampleDATA)

#Extracting data from phyloseq object 
OTUdata <- abundances(phylo_obj)
TAXAData <- as.data.frame(tax_table(phylo_obj)@.Data)
SampleData <- meta(phylo_obj)

save(phylo_obj, file = "phylo_obj.Rdata")

########################## DATA VISUALISATION ############################
names(wes_palettes)
colors <- wes_palette(n=45, name="Darjeeling1", type = "continuous")
colors2 <- wes_palette(n=5, name="Darjeeling1", type ="continuous")
my_palette <- colorRampPalette(c("#e30019", "black", "#00A08A"))(n = 299)

####### 1. Heatmap based on abundance ###
rownames(CountMatrix) = TAX$Species
mean_value<-rowMeans(CountMatrix)
drop_index<-which(mean_value<1)
Countmatrix_new<-CountMatrix[-drop_index,]
#pdf ("Countmatrix6.pdf",width = 20, height = 20 )
#what does heatmap show?
png(filename = "abundance_heatmap.png", width = 10, height = 8, units = "in", res = 300)
heatmap.2(as.matrix (Countmatrix_new), margins = c(10,12), col=my_palette,
          key.title = "Colour Key",
          keysize = 1.3,
          main = "Heat Map of Genus Abundances by Sample")

dev.off()


################
#Visualisation of abundances
phyobj_wrkshp = phylo_obj
gp = subset_taxa(phyobj_wrkshp)
pd <- psmelt(gp) 
# I don't know what 'psmelt' does
pd_filter <- pd %>% filter(Abundance > 0)

######## 2a. Stack bar plot: % of each genus in all samples #######
fig2a <- ggplot(data=pd_filter, aes(x=Sample, y=Abundance, fill = Genus)) + # define which columns you will display in plot
  geom_bar(stat="identity") + # plotting bar plot
  scale_fill_manual(values = colors) +
  # name of x and y axis and title
  labs(x = "Sample", y = "Relative mean abundance (%)",
       title = "Bar Plot of Genus Abundance by Sample", fill = "Genus") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.1), axis.text.x = element_text(angle = 90)) # define the style of plot title. 
fig2a


ggsave("genus_abundance.pdf", plot = fig2a, width =6, height=6, dpi=300) 

###### 2b. Stack bar plot: % of each family in all samples #######
fig2b <- ggplot(data=pd_filter, aes(x=Sample, y=Abundance, fill = Family)) + # define which columns you will display in plot
  geom_bar(stat="identity") + # plotting bar plot
  scale_fill_manual(values = colors) +
  labs(x = "Sample", y = "Relative mean abundance (%)",
       title = "Bar Plot of Family Abundance by Sample", fill = "Family") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.1), axis.text.x = element_text(angle = 90)) # define the style of plot title. 
fig2b
ggsave("family_abundance.pdf", plot = fig2b, width =6, height=6, dpi=300) 

fig2c <- ggplot(data=pd_filter, aes(x=Sample, y=Abundance, fill = Phylum)) + # define which columns you will display in plot
  geom_bar(stat="identity") + # plotting bar plot
  scale_fill_manual(values = colors2) +
  # name of x and y axis and title
  labs(x = "Sample", y = "Relative mean abundance (%)",
       title = "Bar Plot of Phylum Abundance by Sample", fill = "Phylum") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.1), axis.text.x = element_text(angle = 90)) # define the style of plot title. 
fig2c


######## 3. T-Test for each bacteria 
#T-test to test if significant difference between bacteria in healthy and cirrhosis groups
disease<-cbind(t$ERR528726,t$ERR528725,t$ERR528324)
health <- cbind(t$ERR527028, t$ERR527031, t$ERR527081)
t_list<-NULL
p_list<-NULL

for (i in 1:263){
  result<-t.test(disease[i,],health[i,])
  t_abundance<-result$statistic
  p_abundance<-result$p.value
  t_list<-rbind(t_list,t_abundance)
  p_list<-rbind(p_list,p_abundance)
}
#make table to find species with p-values smaller than 0.05 (less abundant species in disease)
index<-which(p_list<0.05)
result<-cbind(t$Species[index],p_list[index],t_list[index])
colnames(result) = c("Species", "P-Value", "T-Value")
#make table to transfer to word
write.table(result, file = "species_p_t_table", sep = ",", quote = FALSE, row.names = F)

#genus table
resultg<-cbind(t$Genus[index],p_list[index],t_list[index])
colnames(resultg) = c("Genus", "P-Value", "T-Value")
#make table to transfer to word
write.table(resultg, file = "genus_p_t_table", sep = ",", quote = FALSE, row.names = F)

#family table 
#2 families coming up more than once
resultf<-cbind(t$Family[index],p_list[index],t_list[index])
colnames(resultf) = c("Family", "P-Value", "T-Value")
#make table to transfer to word
write.table(resultf, file = "family_p_t_table", sep = ",", quote = FALSE, row.names = F)

#phylum table 
#One phylum is coming up more than once
resultp<-cbind(t$Phylum[index],p_list[index],t_list[index])
colnames(resultp) = c("Phylum", "P-Value", "T-Value")
#make table to transfer to word
write.table(resultp, file = "phylum_p_t_table", sep = ",", quote = FALSE, row.names = F)

int_matrix<-rbind(t(disease[index,]),t(health[index,]))
Condition<-as.matrix(c("disease","disease","disease","health","health","health"))
int_matrix<-cbind(Condition,int_matrix)
colnames(int_matrix)<-c("Condition",result[,1])
int_matrix<-data.frame(int_matrix)

int_matrix_2<-cbind(int_matrix$Condition,int_matrix$s__Barnesiella_intestinihominis)
int_matrix_2 <- data.frame(int_matrix_2)
#Make boxplot for the abundances
#ggboxplot(p_list_ttest, x= "Condition", y ="p_list")
labelsch <- c("Cirrhosis", "Health")
box1 <- ggplot(data=int_matrix_2, aes(x=X1, y=as.numeric(X2), fill = X1))+
  geom_boxplot() +
  theme(strip.text = element_text(size=8),
        axis.text.x = element_text(face="bold", size=10, angle = 45, hjust = 1), 
        axis.text.y = element_text(colour = "grey30", size = 10, face = "italic"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.position = "none") +
  scale_x_discrete(labels = labelsch)
box1
ggsave("boxplot_Barnesiella_intestinihominis.pdf", plot = box1, width =6, height=6, dpi=300)
#scale_fill_manual(values = colors) +

##### 4. Alpha diversity
#calculate alpha diversity estimators: chao1 and shannon
Adiv <- estimate_richness(phylo_obj, measures=c("Observed", "Chao1",  "Shannon", "InvSimpson"))
#retrieve the meta data + bind to Adiv
RE_MetaData <- meta(phylo_obj)
Adiv$SampleID <- rownames(Adiv)
Adiv_DF <- full_join(RE_MetaData,Adiv)

#normal distribution test of Chao1 
stChao1<-shapiro.test(Adiv_DF$Chao1)
stChao1
p_value=stChao1$p.value

resultt<-t.test(disease[i,],health[i,])
t.test(c(12.750,18.000,83.000),c(23.000,42.125,26.000))
shannonttest <- t.test(c(1.868929, 0.911557, 1.964093), c(2.764006, 3.057923, 2.561824))
#significant
#Boxplot

ggplot(data=Adiv_DF, aes(x= Condition, y=as.numeric(Shannon), fill = Condition))+
  geom_boxplot()

#Normal distribution of shannon
stShannon <- shapiro.test(Adiv$Shannon)
stShannon

#ABUNDANCE OF SIGNIFICANT SPECIES BY CONDITION MASS BOXPLOT
pd_sig <- pd %>%
  
  filter(OTU == 'OTU 51' | OTU == 'OTU 154' | OTU == 'OTU 165' | OTU == 'OTU 179' | OTU == 'OTU 183' | OTU == 'OTU 185' | OTU == 'OTU 206' | OTU == 'OTU 236' )

speciesboxplot <- ggboxplot(pd_sig, x = "Species", y = "Abundance",
                            color = "Condition",palette = c("#E30019", "#00A08A"), add = "jitter")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = "Species", y = "Mean Abundance",
       title = "Abundance of Significant Species by Condition", fill = "Condition") +
  theme(axis.text = element_text(size = 8))  
speciesboxplot

dev.off()

########### 5. T-Test for alpha diversity estimators - whether estimator levels are significantly different 

index_l<-which(Adiv_DF$Condition=="LC samples")
index_h<-which(Adiv_DF$Condition=="health")
shannon_1<-Adiv_DF$Shannon[index_l]
shannon_2 <- Adiv_DF$Shannon[index_h]
#t test
ttest_result<-t.test(shannon_1,shannon_2)
ttest_result
#p-value = 0.2513 > 0.05 so accept null hypothesis, no difference

chao1_1<-Adiv_DF$Chao1[index_l]
chao1_2 <- Adiv_DF$Chao1[index_h]
#t test
ttest_result2<-t.test(chao1_1,chao1_2)
ttest_result2

t_stat<-ttest_result$statistic
p_ttest<-ttest_result$p.value

###### 6. Linear model between two continuous metadata/ between metadata with alpha diversity estimators

#CRP c reactive protein 
glmShannonCrea <- glm(Shannon ~ Crea..umol.L., data=Adiv_DF, family="quasipoisson")
summary(glmShannonCrea) #Not significant

glmChao1Crea <- glm(Chao1 ~ Crea..umol.L., data=Adiv_DF, family="quasipoisson")
summary(glmChao1Crea) #Significant
#p-value = 0.00774 **

glmShannonAlb <- glm(Shannon ~ Alb..g.L., data=Adiv_DF, family="quasipoisson")
summary(glmShannonAlb) 
#significant p-value = 0.0386

glmChao1Alb <- glm(Chao1 ~ Alb..g.L., data=Adiv_DF, family="quasipoisson")
summary(glmChao1Alb) # Not significant

glmChao1TB <- glm(Chao1 ~ TB..umol.L., data=Adiv_DF, family="quasipoisson")
summary(glmChao1TB) # Not significant

glmShannonTB <- glm(Shannon ~ TB..umol.L., data=Adiv_DF, family="quasipoisson")
summary(glmShannonTB) # Not significant

plot(Adiv_DF$Crea..umol.L., Adiv_DF$Chao1, pch = 19)
abline(glmShannonCRPLevel$coefficients, col = "red", lty = 2)

p3 <- ggplot(Adiv_DF, aes(x=Alb..g.L., y=Shannon, color=Condition)) +
  geom_point() +
  geom_smooth(method=lm , color="#E30019", fill="#00A08A", se=TRUE) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size =6))
p3

####### 7. Beta Diversity 
aovShannonCondition = aov(Shannon ~ Condition, data=Adiv_DF)
summary(aovShannonCondition)

#comparing patients from one patient group with e.g. obese group 
#have to transform data 
phyobj_wrkshp_shift <- microbiome::transform(phyobj_wrkshp, transform="shift", shift=1)
phylo_BrayDis <- phyloseq::distance(phyobj_wrkshp_shift, "bray")
wrkshp__MetaData <- meta(phyobj_wrkshp_shift)

#From this should be able to plot PCoA 
#if you have time found out how and put it in but it's not crucial
head(phylo_BrayDis)

#permANOVA analysis for beta diversity (Condition)
adonis(phylo_BrayDis ~ Condition, data=wrkshp__MetaData, permutations=999)

#permANOVA analysis for beta diversity (Gender)
adonis(phylo_BrayDis ~ Gender, data=wrkshp__MetaData, permutations=999)

#permANOVA analysis for beta diversity (Alcohol)
adonis(phylo_BrayDis ~ Alcohol.related..Y.or.N., data=wrkshp__MetaData, permutations=999)

#permANOVA analysis for beta diversity (Bilirubin)
adonis(phylo_BrayDis ~ TB..umol.L., data=wrkshp__MetaData, permutations=999)
###Significant
