#Task_4

#Calculate the square root of 10 in R
sqrt(10)
#Calculate the logarithm of 32 to the base 2 in R
log(32, base = 2)
#Calculate the sum of the numbers from 1 to 1000 in R
sum(1:1000)
#Calculate the sum of all even numbers from 2 to 1000 in R
sum(seq(2,1000,by=2))
#Claculate pairwise comparisions number for 100 genes
choose(100,2)
#Claculate 100 genes arranged in triples
factorial(100)/(factorial(100-3))

#Task_5

#read data CO2 to Rstudio
data(CO2)
#Open help documents for CO2 dataset
?CO2
#The CO2 data frame has 84 rows and 5 columns of data from an experiment on the 
#cold tolerance of the grass species Echinochloa crus-galli.
#By head(CO2) function, the data is clean data
#Noticed the Type only contains Quebec and Mississippi regarding to 
#the question 3
#Show the aeverage and median CO2 uptake of the plants from Quebec and Mississippi#
library(tidyverse)
CO2 %>% 
  group_by(Type) %>% 
  summarise(
    mean_uptake=mean(uptake),
    median_uptake=median(uptake))
#Don't know why my default pipe operator |> not work in VScode
#Summary how many genes are expressed in each sample in "airway" example data
#and how many genes are not expressed in and sample
BiocManager::install("airway")
library(airway)
data(airway)
#library a S4 container for the airway data
library(SummarizedExperiment)
#create a data.frame to store the col and row with ENSMBLE id listed in the table
airway_counts<-assay(airway)%>%
  as.data.frame()%>%
  rownames_to_column(var="ENSMBLE_id")
#tidy the data and count for how many genes are expressed in each sample
gene_expressed<-airway_counts %>%
  pivot_longer(-ENSMBLE_id,names_to="sample_SRR",values_to="raw_reads")%>%
  group_by(sample_SRR) %>%
  summarise(sum(raw_reads>0))
gene_expressed
#The responses from R:
#  sample_SRR `sum(raw_reads > 0)`
#    <chr>                     <int>
#  1 SRR1039508                24633
#  2 SRR1039509                24527
#  3 SRR1039512                25699
#  4 SRR1039513                23124
#  5 SRR1039516                25508
#  6 SRR1039517                25998
#  7 SRR1039520                24662
#  8 SRR1039521                23991
#For genes not expressed, filter and select the gene id that in any sample raw counts are always 0.
#Functions applied to each row (ENSMBLE_id).
genes_not_expressed<-airway_counts %>%
  rowwise() %>%
  filter(all(c_across(-ENSMBLE_id) == 0)) %>%
  select(ENSMBLE_id) 
genes_not_expressed
#The responses from R:
# A tibble: 30,633 × 1
# Rowwise: 
#ENSMBLE_id     
#<chr>          
#1 ENSG00000000005
#2 ENSG00000004809
#3 ENSG00000004848
#4 ENSG00000004939
#5 ENSG00000004948
#6 ENSG00000005001
#7 ENSG00000005981
#8 ENSG00000006059
#9 ENSG00000006071
#10 ENSG00000006075

#Task 6

#Function that calculates the ratio of the mean and the median of a given vector
mean_median_ratio<-function(x){
  mean_vector<-mean(x)
  median_vector<-median(x)
  mean_vector/median_vector
}

#Write a function that ignores the lowest and the highest value from a given vector and calculate the mean
mean_filtered<-function(x){
  vector_max=max(x)
  vector_min=min(x)
  mean(x[x!=vector_max&x!=vector_min])
}

#Usage of pipes

#The pipes operator allows the pass of augment from one function to another.  
#This usage also results in some scenarios where the pipe should be avoided during R coding.  
#First of all, the intermediate results are important and should be stored in a new date frame.  
#As pipe ignores all the intermediate results, you can not reuse the intermediate results when you apply pipe for your coding.  
#Therefore, to some extent, pipe decreases the reusability of your coding.

#Second, many of R functions do not have the first argument as a data frame for you to pass.  
#Such as tapply, in this scenario, using intermediate augments is easier during coding.

#Moreover, using lots of pipes makes your coding less readable and difficult to debug.  
#Whenever you are going to debug a script with lots of pipes, you always need to read the code from the very beginning through the end of the coding.  
#It is very difficult to find a bug from the “pipe forests”.  When you miss the bug, you need to read from the very beginning and debug again, decreasing the possibility of your coding working properly.

#Although it is less common in the bioinformatics area, the usage of pipes will increase the usage of memory as it creates intermediate variables.  
#Therefore, system memory is a concern, it is better not to use pipes (or maybe not use R coding in this scenario).  

#Apply-family functions
#The apply family of functions allows you to code with R concisely and efficiently.  
#By using the apply family of functions, one can avoid writing difficult loops to apply functions to columns, rows, and matrics in one command line.  
#Although without knowledge of apply-family functions, the coding is very difficult to understand.  
#However, after learning it, it is the most understandable function, even better than lots of functions in dplyr package (for example the group_by() function vs the tapply()).  
#The apply-family functions will be more powerful when apply complicated functions to the datasheet, for example when doing group-wise gene expression analysis during scRNAseq analysis with the Seurat package.

#Task 7

#Compare the distributions of the body heights of the two species by using hist, ggplot with geom_histogram plot
library("ggplot2")
#read data
magic_guys<-read.csv("magic_guys.csv")
#extract targeted parameters
jedi_heights <- magic_guys$length[magic_guys$species == "jedi"]
sith_heights <- magic_guys$length[magic_guys$species == "sith"]
#Show the figure in a same figure, room for 2 figures
par(mfrow = c(1, 2))
#Generate figures
hist(jedi_heights,col="blue",xlab="length",main="The distribution of length for jedi")
hist(sith_heights,col="red",xlab="length",main="The distribution of length for sith")
#use ggplot
ggplot(magic_guys, aes(x=length, fill=species)) +
  geom_histogram(alpha=0.5,position="stack",bins=15) +
  labs(title="Species:Jedi vs Sith", x="Hights", y="Count")+
  theme_classic()
#use ggplot                      
ggplot(magic_guys, aes(x=species, y=length, fill=species))+
  geom_boxplot(alpha = 0.5)+
  labs(title = "Species:Jedi vs Sith",
       x="Species",
       y="Height") +
  theme_classic()
#save as png format due to the background is transparent.

#Load the gene expression data matrix from 'microarray_data.tab' and do the visualization as instructed.
microarray<-read.delim2("microarray_data.tab",header=T,sep='\t',na.strings = "NA")
dim(microarray)
library(tidyverse)
#Only filter the genes with missing_counts
missing_counts <- microarray %>%
  summarise_all(~ sum(is.na(.))) %>% 
  gather(key="Gene", value="Missing") %>% 
  filter(Missing>0)
#plot the reults
ggplot(missing_counts, aes(x = reorder(Gene, -Missing), y=Missing)) +
  geom_bar(stat="identity", fill="blue") +
  theme_classic() +
  labs(title="missing values of genes", y="Count of missing values", x = "Genes")  #This task is too difficult for me, I give up, the figure looks weird.
#find 10%,20%,50% missing value
missing_counts_0.1 <- microarray %>%
  summarise_all(~ sum(is.na(.))) %>% 
  gather(key="Gene", value="Missing") %>% 
  filter(Missing>0.1)
missing_counts_0.1                 #missing value>10%

missing_counts_0.2 <- microarray %>%
  summarise_all(~ sum(is.na(.))) %>% 
  gather(key="Gene", value="Missing") %>% 
  filter(Missing>0.2)
missing_counts_0.2                 #missing value>20%

missing_counts_0.5 <- microarray %>%
  summarise_all(~ sum(is.na(.))) %>% 
  gather(key="Gene", value="Missing") %>% 
  filter(Missing>0.5)
missing_counts_0.5                 #missing value>50%

#Replace the missing value by the average expression value of particular gene.
library(tidyverse)
microarray<-read.delim2("microarray_data.tab",header=T,sep='\t',na.strings = "NA")
microarray_copy<-microarray
microarray_copy%>%
  mutate(across(where(is.numeric),~ ifelse(is.na(.),mean(., na.rm = TRUE),.))) #mutate from tidyverse to forming new dataframe
write.table(microarray_copy,"gene expression with NA replaced", sep='\t')      #write into new table rather than replace old data

#visualize the CO2 dataset
data(CO2)
ggplot(CO2, aes(x=conc, y=uptake, color=Type)) +
  geom_point(size=3, alpha=0.5) +  
  geom_smooth(method="glm", se=FALSE, linetype="dashed") +
  theme_classic() +
  labs(title = "CO2 Uptake vs. Concentration",
       x = "CO2 Concentration", 
       y = "CO2 Uptake")
#looks the plant in Quebec can uptake more CO2 when CO2 conc. is same in Mississippi.

#Task 8

devtools::install_github("hirscheylab/tidybiology")
library(tidybiology)
library(tidyverse)
library(patchwork)
data(chromosome)
data(proteins)
summarise(h)
#summarise the data
chromosome_summary<-chromosome %>%
  summarise(mean_variations=mean(variations),
            median_variations=median(variations),
            max_variations=max(variations),
            mean_proteincodingenes=mean(protein_codinggenes),
            median_proteincodinggenes=median(protein_codinggenes),
            max_proteincodinggenes=max(protein_codinggenes),
            mean_mi_rna=mean(mi_rna),
            median_mi_rna=median(mi_rna),
            max_mi_rna=max(mi_rna)
  )

#plot the chromosome distribution
ggplot(chromosome, aes(x=basepairs)) +
  geom_density(fill="blue", alpha=0.5) +
  theme_classic() +
  labs(title="Chromosome distribution",
       x="Chromosome Size (bp)",
       y="Density")

#correlation between coding genes vs length & miRNA vs length.
plot1<-ggplot(chromosome, aes(x=length_mm, y=protein_codinggenes))+
  geom_point(color="blue", alpha=0.5, size=3)+
  geom_smooth(method="lm", color="red")+
  theme_minimal()
plot2<-ggplot(chromosome, aes(x=length_mm, y=mi_rna))+
  geom_point(color="green", alpha=0.5, size=3)+
  geom_smooth(method="lm", color="red")+
  theme_minimal()                             #creative plots one by one
plot=plot1+plot2                              #combine with patchwork
plot

#Summary statistics for proteins
summary_proteins<-proteins %>%
  summarise(mean_length=mean(length),
            median_length=median(length),
            max_length=max(length),
            mean_mass=mean(mass),
            median_mass=median(mass),
            max_mass=max(mass)
  )
summary_proteins
#plot of correlation plots
ggplot(proteins, aes(x=length, y=mass))+
  geom_point(color="blue", alpha=0.5, size=3)+
  geom_smooth(method="lm", color="red",linetype=dashed)+
  theme_classic()+
  labs(titles="length vs mass",
       x="protein length",
       y="protein mass")