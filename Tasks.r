#Task_4
#Calculate the square root of 10 in R#
sqrt(10)
#Calculate the logarithm of 32 to the base 2 in R#
log(32, base = 2)
#Calculate the sum of the numbers from 1 to 1000 in R#
sum(1:1000)
#Calculate the sum of all even numbers from 2 to 1000 in R#
sum(seq(2,1000,by=2))
#Claculate pairwise comparisions number for 100 genes#
choose(100,2)
#Claculate 100 genes arranged in triples#
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