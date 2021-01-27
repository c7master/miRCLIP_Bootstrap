##################################################################################################
# Yifei Yan
# McGill University
# University of Montreal
#
# Dec 12, 2020
#
# This script takes Lian Mignacca's RNA-seq identified miR-146a targets and performs bootstrapping 
# upon the lncRNA sequences predicted by different miRNA target prediction software tools, 
# and test the significance of the fraction of lncRNAs that were detected. 

# 7500 lnRNAs, 

# This script requies the following input files:
# 1. List of RNAseq-identified lncRNAs (about 7500): lncRNA_in_Seq.csv
# 2. The number of lncRNA genes selected by 2-fold change cut-off from the identified lncRNAs: 1425
# 3. List of all lncRNA targets predicted by EACH software tool: predicted_mir146_targets.csv
# 4. Number of lncRNAs selected by 2-fold-change cut-off and identified by at least one of the prediction tools (1121)

###################################################################################################


# read in all lncRNAs identified by clip-seq as a dataframe
all_lncRNAs <- read.csv("lncRNA_in_Seq.csv")
summary(all_lncRNAs)
head(all_lncRNAs)

# Parse the $Gene column into a vector to save the gene names. 
# Use vector to look up the list of predicted targets.
mir146a_lncRNAs_by_name <- data.frame(all_lncRNAs$Gene)
head(mir146a_lncRNAs_by_name)
summary(mir146a_lncRNAs_by_name)
##########################
# Parse the lncRNA sequences predited as miR-146a target by all 
# software tools into a dataframe. Retrieve results of individual tools from 
# the columns.
all_predicted_targets <- read.csv(file="predicted_mir146_targets.csv")
summary(all_predicted_targets)
##########################


########################### 
# make a function that counts the number of lncRNAs identified by 
# at least one software tool 
# count_identified()
############################

count_identified <- function (sample1_mat) {  
# Next, check membership of entries in sample1 with the tool-predicted
# list of miR-146. Use intersect function from matrix.

number_of_identified <- 0
################################### miRWalk
# Parse the prediction of individual tools into a vector
miRWalk_targets <- all_predicted_targets$miRWalk
summary(miRWalk_targets)
# parse prediction tool results into matrices (vectors)
miRWalk_tar_mat <- as.matrix(miRWalk_targets)
summary(miRWalk_tar_mat)
head(miRWalk_tar_mat[,1])

#### compare intersection: bootstrap sample with miRWalk prediction
intersection_s1_miRWalk <- intersect(sample1_mat[,1], miRWalk_tar_mat[,1])
head(intersection_s1_miRWalk)
summary(intersection_s1_miRWalk)
length(intersection_s1_miRWalk)

### update the union_of_targets
union_of_targets <- union(union_of_targets,intersection_s1_miRWalk)
number_of_identified <- length(union_of_targets)

################################### Targetscan
Targetscan_targets <- all_predicted_targets$Targetscan
summary(Targetscan_targets)
Targetscan_tar_mat <- as.matrix(Targetscan_targets)
summary(Targetscan_tar_mat)
head(Targetscan_tar_mat[,1])

#### compare intersection: bootstrap sample with Targetscan prediction
intersection_s1_Targetscan <- intersect(sample1_mat[,1], Targetscan_tar_mat[,1])
head(intersection_s1_Targetscan)
summary(intersection_s1_Targetscan)

### update the union_of_targets
union_of_targets <- union(union_of_targets,intersection_s1_Targetscan)
number_of_identified <-length(union_of_targets)

################################### Miranda
Miranda_targets <- all_predicted_targets$Miranda
summary(Miranda_targets)
Miranda_tar_mat <- as.matrix(Miranda_targets)
summary(Miranda_tar_mat)
head(Miranda_tar_mat[,1])

#### compare intersection: bootstrap sample with Miranda prediction
intersection_s1_Miranda <- intersect(sample1_mat[,1], Miranda_tar_mat[,1])
head(intersection_s1_Miranda)
summary(intersection_s1_Miranda)

### update the union_of_targets
union_of_targets <- union(union_of_targets,intersection_s1_Miranda)
number_of_identified <-length(union_of_targets)

################################### RNAhybrid
RNAhybrid_targets <- all_predicted_targets$RNAhybrid
summary(RNAhybrid_targets)
RNAhybrid_tar_mat <- as.matrix(RNAhybrid_targets)
summary(RNAhybrid_tar_mat)
head(RNAhybrid_tar_mat[,1])

#### compare intersection: bootstrap sample with RNAhybrid prediction
intersection_s1_RNAhybrid <- intersect(sample1_mat[,1], RNAhybrid_tar_mat[,1])
head(intersection_s1_RNAhybrid)
summary(intersection_s1_RNAhybrid)

### update the union_of_targets
union_of_targets <- union(union_of_targets,intersection_s1_RNAhybrid)
number_of_identified <-length(union_of_targets)

################################### DianalncBase
DianalncBase_targets <- all_predicted_targets$Diana.lncBase
summary(DianalncBase_targets)
DianalncBase_tar_mat <- as.matrix(DianalncBase_targets)
summary(DianalncBase_tar_mat)
head(DianalncBase_tar_mat[,1])

#### compare intersection: bootstrap sample with DianalncBase prediction
intersection_s1_DianalncBase <- intersect(sample1_mat[,1], DianalncBase_tar_mat[,1])
head(intersection_s1_DianalncBase)
summary(intersection_s1_DianalncBase)

### update the union_of_targets
union_of_targets <- union(union_of_targets,intersection_s1_DianalncBase)
number_of_identified <-length(union_of_targets)

################################### RNA22
RNA22_targets <- all_predicted_targets$RNA22
summary(RNA22_targets)
RNA22_tar_mat <- as.matrix(RNA22_targets)
summary(RNA22_tar_mat)
head(RNA22_tar_mat[,1])

#### compare intersection: bootstrap sample with RNA22 prediction
intersection_s1_RNA22 <- intersect(sample1_mat[,1], RNA22_tar_mat[,1])
head(intersection_s1_RNA22)
summary(intersection_s1_RNA22)

### update the union_of_targets
union_of_targets <- union(union_of_targets,intersection_s1_RNA22)
number_of_identified<-length(union_of_targets)
number_of_identified
##################

return(number_of_identified)
} # end of count_identified function


##########################  The main program starts here ###################################
# generate bootstrapped samples n times and save all of them in a 
# dataframe/table of size  1425 X n. Perform checking on each column
# by a loop scripted.

n<-1000
X<-1425

nubmer_of_bootstrap_Subsets <- n
number_of_genes_per_set <-X

bootstrap_lncRNA_sets <- replicate(n,sample(all_lncRNAs$Gene, X, replace = TRUE))
summary(bootstrap_lncRNA_sets)

bootstrap_counts <- c(1:n)
# initialize a vector to hold all the counts of identified lncRNAs from all boostrapped subsets
identity_counts <- c()

# initialize the set of union of targets identified by all different tools
union_of_targets <- c()

# for loop go through every column of the dataframe of bootstrapped subsets, one subset per column
for (col_index in bootstrap_counts) {
  #sample <- sample(all_lncRNAs$Gene, 1425, replace = TRUE)
  subset_sample <- bootstrap_lncRNA_sets[,col_index]
  # head(subset_sample)
  summary(subset_sample)

  # Save gene names in a matrix and use all rows of the first column for membership check. 
  subset_sample_mat <- as.matrix(subset_sample)
  # head(subset_sample_mat)
  summary(subset_sample_mat)
  head(subset_sample_mat[,1]) # The matrix has only one column, so it's always [,1] in any loop

  software_detected_count <- count_identified(subset_sample_mat) # pass the 1-col matrix of 1425 lncRNAs to the counting function
  software_detected_count
  # call the count function and store results in a vector by appending new results to the vector
  identity_counts <- append(identity_counts, software_detected_count)

}# end of for loop

identity_counts
summary(identity_counts)
write.csv(identity_counts, "IdentityCounts_bootstrap.csv")

############ analysis and plot of the result ############
hist(identity_counts)
bootstrap_pop_mean<-mean(identity_counts)
bootstrap_pop_sd<-sd(identity_counts)

###### test stats functions ######
real_count <- 1121
z_val <- (real_count - bootstrap_pop_mean)/bootstrap_pop_sd
z_val
p_value <- round(pnorm(z_val, bootstrap_pop_mean, bootstrap_pop_sd), 10)
p_value

#### output: ### 
# > z_val
# [1] 6.840985
# > p_value
# [1] 0
##############