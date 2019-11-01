library(tidyverse)
library(stringr)
library(ggfortify)
library(cluster)
library(factoextra)
library(magrittr)
library(tibble)
library(NbClust)
library(viridis)
library(dendextend)
library(heatmaply)
library(superheat)
# Sometimes the xlsx doesn't load well on different computers. 
# It's just used for exporting an Excel file at the end, so don't worry if it's not working for you
library(xlsx)

############ Step 1: Importing Data and Calculating Fold Changes ####################
#####################################################################################
# Before you begin, you want to make sure your Skyline report is in the right format
# There is a Skyline report file included with this script which can be imported into Skyline
# The peptide note is the name of each peptide. 
# This workflow will be easiest if your names are similar to the ones here, but you can also modify this workflow to reflect your peptide names

# BEFORE YOU BEGIN  - open the Skyline report in Excel and replace all #N/A with 0 or NA
# R doesn't recognize #N/A as NA and makes everything factors and it's THE WORST

# Reading in all the data
data <- read.csv("Histone_Report.csv")

# Turning peptide notes and column names into something that R can read
data <-  as.data.frame(data)
data$Peptide.Note <- unlist(lapply(data$Peptide.Note, function(x){unlist(str_replace_all(x, "[ :().,]", ""))}))
data$Peptide.Note <- unlist(lapply(data$Peptide.Note, function(x){unlist(str_replace_all(x, "K18ac/", "K18"))}))
data$Peptide.Note <- unlist(lapply(data$Peptide.Note, function(x){unlist(str_replace_all(x, "H4", "H4__"))}))
colnames(data) <- sub('[.]', '_', make.names(colnames(data), unique = TRUE))

# Spreading the table so that each column is a different sample with areas for each mod
data <- data %>% 
  spread(File_Name, Total_Area.MS1)

# Summing mods with multiple precursors
data[is.na(data)] <- 0
data <- data %>% 
  select(-Precursor_Charge) %>% 
  group_by(Peptide_Note) %>% 
    summarise_all(sum)

# NOTE: I didn't include the data check step in this script (like there is in EpiProfile)
# If you've processed your data through Skyline, you should know what the data quality is already
# So just delete any bad quality samples before you export your data
# With great power comes great responsibility

# Filter out peptides that weren't idenitifed in many samples
# This command calculates the mean value for each peptide and the number of times it wasn't identified
# You can check which peptides these are by opening the data dataframe and sorting by the Count column
data <- data %>% 
  mutate(Mean = rowMeans(data[,-1]), 
         Count = rowSums(data == 0))

# We'll filter out any peptides that were at very low abundance, or were present in less than 70% of the samples
# You can change this filter by modifying the number after "(ncol(data)-4)*"
data <- data %>% 
  filter(Mean > 1e+6) %>% 
  filter(Count <= (ncol(data)-4)*.3) %>% 
  select(-Mean, -Count)

# Caculating sums for each peptide family
# To do this, first make a label for each peptide family
data$peptide_id <-unlist(lapply(data$Peptide_Note, function(x){
  unlist(str_sub(x, 1, 4))
}))

# Now normalize the data
normalized <- data %>%
  ungroup() %>%
  mutate_if(is.double, list(~./sum(., na.rm = TRUE))) %>%
  select(-peptide_id)

# Or you if you want to normalize by peptide
# normalized <- data %>%
#   group_by(peptide_id) %>%
#   mutate_if(is.double, list(~./sum(., na.rm = TRUE))) %>%
#   ungroup() %>%
#   select(-peptide_id)

# Transpose the data
normalized <- normalized %>% 
  gather(Sample, value, 2:ncol(normalized)) %>% 
  spread(Peptide_Note, value)

# Now create a column that contains your sample groups
# Basically, you want to split all your biological/technical replicates into seperate groups so you can average and compare them
# This is done by finding a unique ID that matches the names of all the samples in each group and then assigning a number to the samples it matches
# For example, if I wanted to choose all the B_CutC samples from B_CutC_1, B_CutC_2, C_CutC_3, and B_WT_1, I would use "B_CutC" as the identifier
# You can use more than one set of characters for your identifier if you want, just input grepl("ID1|ID2", Sample) (| means or). See regex for more details
# The last number at the end of the command will be assigned to any samples which aren't found by any of your IDs (useful for checking that everything worked)
# Delete or add as many columns as you need to this command, just make sure there are enough parantheses at the end to close it off
groups <- normalized %>% 
  mutate(Sample_Type = ifelse(grepl("ID_1", Sample), "set01",
                       ifelse(grepl("ID_2", Sample), "set02", 
                       ifelse(grepl("ID_3", Sample), "set03", 
                       ifelse(grepl("ID_4", Sample), "set04",
                       ifelse(grepl("ID_5", Sample), "set05",
                       ifelse(grepl("ID_6", Sample), "set06",
                       ifelse(grepl("ID_7", Sample), "set07", 
                       ifelse(grepl("ID_8", Sample), "set08", 
                       ifelse(grepl("ID_9", Sample), "set09", 
                       ifelse(grepl("ID_10", Sample), "set10",
                       ifelse(grepl("ID_11", Sample), "set11", 
                       ifelse(grepl("ID_12", Sample), "set12", 
                       ifelse(grepl("ID_13", Sample), "set13",
                       ifelse(grepl("ID_14", Sample), "set14", 
                       ifelse(grepl("ID_15", Sample), "set15",
                       ifelse(grepl("ID_16", Sample), "set16", "set17")))))))))))))))))

# It's useful to check that each sample is in its appropriate group, you can use this command to pull out just the sample names and groups
peptide_check <- groups %>% 
  select(Sample, Sample_Type) %>% 
  arrange(Sample_Type)
rm(peptide_check)

# It's often useful to take this normalized data table and plot a correlation plot using the command below
# This helps you see if your biological replicates are similar to each other and whether any samples are significant outliers
distance <- as.data.frame(normalized)

# This will create a matrix where the most similar samples will have a 0 value and the value increases with increasing difference
distance <- column_to_rownames(distance, var = "Sample")
res.dist <- get_dist(distance, method = "spearman")
fviz_dist(res.dist, lab_size = 5, order = FALSE)

# After you've looked at the data above, you can look at which samples have high disimilarity from the rest of your samples using this command.
dist_check <- as.matrix(res.dist)
dist_check <- as.data.frame(dist_check) %>% 
  rownames_to_column()
dist_check <- as_tibble(dist_check)

distance <- dist_check %>% 
  ungroup() %>% 
  mutate(Mean = select(dist_check, -rowname) %>%
           rowMeans()) %>% 
  select(rowname, Mean) %>% 
  arrange(desc(Mean))  

# In general, I like to say that a sample with a disimilarity above 0.05 may be an outlier, but this number is arbitrary and will depend on your samples
# However, it may just be that one set of samples is different from the rest (all your brain samples are different than colon samples, etc)
# Before you throw something out, it is useful to see if it is more similar to samples of the same tissue, treatment, etc. 
# You'll need to use an ID that matches all the samples you want to compare
distance_group <- dist_check %>% 
  mutate(Mean1 = select(dist_check, contains("ID")) %>% 
           rowMeans()) %>% 
  select(rowname, Mean1) %>% 
  arrange(desc(Mean1)) %>%
  filter(str_detect(rowname, 'ID')) 

# You can then join these two matrices to find samples aren't similar to their biological replicates or the rest of the data
poss_outliers <- distance %>%
  inner_join(distance_group, by = "rowname")

rm(poss_outliers)
rm(distance_group)
rm(distance)
rm(dist_check)
rm(res.dist)

# If any samples are outliers (they don't correlate well with any of your other data), it might be useful to remove them using the following command
# Change "row_name_1", etc. to the name of the samples you want to remove
# Obviously, you should also try to figure out why any sample is an outlier before discarding it :)
# normalized <- normalized %>% 
#   filter(!Sample %in% c("row_name_1", "row_name_2","etc"))

# Now you can start doing calculations.  
average <- groups %>% 
  group_by(Sample_Type) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  select(-c(Sample))
# You'll get a warning on this command, it just means that R couldn't average your sample names, which is probably good since they're letters and not numbers

# Next calculate standard deviations for each set of replicates
stdev <- groups %>% 
  group_by(Sample_Type) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  select(-c(Sample))
stdev$Sample_Type <- unlist(lapply(stdev$Sample_Type, function(x){unlist(str_replace_all(x, "$", "_Sd"))}))

# Now you'll caluculate fold changes. It's easiest to transpose the data before you do this.
average <- average %>% 
  gather(Peptide, value, 2:ncol(average)) %>% 
  spread(Sample_Type, value)

# This will calculate all possible fold changes in your dataset
# Obviously, we don't need all of these (eg, we know dividing a sample by itself will always equal 1)
# We'll filter out the unecessary fold changes later
fold_change <- average
fold_change <- fold_change %>% 
  select_if(is.numeric) %>% 
  map(~./fold_change %>% select_if(is.numeric)) %>% 
  imap(~set_names(.x, paste0(.y, "-", names(.x)))) %>% 
  bind_cols()
fold_change <- bind_cols(average[1], fold_change)

# Calculate p values for the fold changes you just calculated above
# This requires you to change the format of the table
long <- groups %>% 
  select(-c(Sample)) %>% 
  select(Sample_Type, everything())
long <- long %>%   
  gather(Peptide, value, 2:ncol(long)) 
colnames(long) <- c("Sample_Type", "Peptide", "value")

# Now perform an Anova
# There are two types of anovas: one-way and two-way. 
# Use one-way if you're just comparing factors to one control (eg, three cell treatments to untreated), use two-way for everything else
# One-way Anova ####################################################################
# Convert to factor
long_2 <- long
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

# Make into dataframe of lists
long_3 <- split(long_2, paste(long_2$Peptide))

# Temp function for running anova
tmpfn <- function(dat) {
	 fit <- lm(value ~ Sample_Type, dat)
	 rsq <- summary(fit)$r.squared
	 if(rsq < 0.99)
	    pval <- anova(fit)[,"Pr(>F)"][1]
	 else
	    pval <- NA
	 coefs <- t(coef(fit))
	 data.frame(coefs,rsq, pval)
}

long_stat <- long_3 %>%
    map(tmpfn) %>%
    bind_rows(.id = "Peptide")

long_stat <- long_stat %>%
  filter(pval < 0.05)
list <- c(long_stat[,1])
long <- filter(long, Peptide %in% list)

# Perform Tukey's HSD to get adjusted pvalues for all possible comparisons
long_2 <- long
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

long_3 <- split(long_2, paste(long_2$Peptide))

tmpfn <- function(dat) {
  fit <- aov(value ~ Sample_Type, dat)
  Tukey <- TukeyHSD(fit)$Sample_Type
}

long_stat <- long_3 %>%
  map(tmpfn)

p_adj <- data.frame(long_stat)
p_adj <- p_adj %>% 
  select(contains("p.adj"))
names(p_adj) <- sub(".p.adj", "", names(p_adj))
p_adj <- rownames_to_column(p_adj, "Sample")

# You can use this command to write results of Anova to a .csv if you wish
# write.csv(long_stat, file = "filename_anova2.csv")
### End One-way Anova #################################################################################

# Two-way Anova ######################################################################
# Separate sample type column into seperate condition columns
# This command seperates by a "_", so if you've named each sample set as Factor1_Factor2, this will seperate everything correctly
long_2 <- long
long_2 <-separate(long_2, Sample_Type, c("Factor1", "Factor2"), sep = "_")

# Convert to factors
long_2$Factor1 <- as.factor(long_2$Factor1)
long_2$Factor2 <- as.factor(long_2$Factor2)

# Make into dataframe of lists
long_3 <- split(long_2, paste(long_2$Peptide))

# Create temp function for anova
# Output is pvalues from F test
tmpfn <- function(dat) {
  fit <- aov(value ~ Factor2 + Factor1  + Factor2:Factor1, dat)
  sum_aov <- summary.aov(fit)
  pval_Factor2 <- sum_aov[[1]][["Pr(>F)"]][1]
  pval_Factor1 <- sum_aov[[1]][["Pr(>F)"]][2]
  pval_pair <- sum_aov[[1]][["Pr(>F)"]][3]
  res <- (fit)$res
  data.frame(pval_Factor2,pval_Factor1, pval_pair, res)
}

long_stat <- long_3 %>%
    map(tmpfn) %>%
    bind_rows(.id = "Peptide")

# Check that residuals are evenly distributed
hist(long_stat$res,main="Histogram of
     residuals",xlab="Residuals")

long_stat <- long_stat %>%
  select(-res) %>%
  unique()

# Two-way anovas can have main effects and interaction
# Significant interaction effect indicates relationship between dependent variable and main effect (Factor 1) depends on other main effect (Factor 2)
# We can filter out those peptides
long_interaction <- long_stat %>%
  filter(pval_pair <= 0.05)
list <- c(long_interaction[,1])

# Plot the interaction
# This script will make an interaction plot for each peptide that showed significant interaction between Factor1 and Factor2
long_2_plot <- long_2
long_2_plot <- filter(long_2, Peptide %in% list)
long_2_plot <- long_2_plot %>%
  group_by(Factor1, Factor2, Peptide) %>%
  summarise_at(c("value"), mean)
# If there's a lot of Peptides in this category (>15), the plot can get a little crowded
# You can use the commands below to separate out peptides belongs in Histone H3, H4, etc.
long_2_plot_H3 <- long_2_plot %>%
  filter(str_detect(Peptide, "H3"))
long_2_plot_H4 <- long_2_plot %>%
  filter(str_detect(Peptide, "H4"))
long_2_plot_H2 <- long_2_plot %>%
  filter(str_detect(Peptide, "H2"))
# Make sure to reference these filtered plots in the first line of the ggplot command below

ggplot(long_2_plot) +
  aes(x = Factor1, y = value, color = Factor2) +
  geom_line(aes(group = Factor2)) +
  geom_point() +
  facet_grid(vars(), vars(Peptide))

# Now you can look at peptides with significantly changing main interactions (eg, changing in either Factor1 or Factor2)
long_main <- long_stat %>%
  mutate(pval = apply(long_stat[,2:3], 1, FUN=min))

long_main <- long_main %>%
  filter(pval <= 0.05)

list <- c(long_main[,1])
long <- filter(long, Peptide %in% list)

# Perform Tukey's HSD to get adjusted pvalues for all possible comparisons
long_2 <- long
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

long_3 <- split(long_2, paste(long_2$Peptide))

tmpfn <- function(dat) {
  fit <- aov(value ~ Sample_Type, dat)
  Tukey <- TukeyHSD(fit)$Sample_Type
}

long_stat <- long_3 %>%
  map(tmpfn)

p_adj <- data.frame(long_stat)
p_adj <- p_adj %>%
  select(contains("p.adj"))
names(p_adj) <- sub(".p.adj", "", names(p_adj))
p_adj <- rownames_to_column(p_adj, "Sample")

# You can use this command to write results of Anova to a .csv if you wish
# write.csv(long_stat, file = "filename_anova2.csv")
### End Two-way Anova #################################################################################
rm(list)
rm(tmpfn)
rm(long_2)
rm(long_3)
rm(long_stat)

# Combine average, standard deviation, fold change, and p-values into one table
# First, make sure all your data is transposed
stdev <- stdev %>% 
  gather(Peptide, value, 2:ncol(stdev)) %>% 
  spread(Sample_Type, value)

# Filter out the extra fold change comparisons we don't need
fold_change <- fold_change %>% 
  gather(Sample, value, 2:ncol(fold_change)) %>% 
  spread(Peptide, value)
list <- c(p_adj[,1])
fold_change <- filter(fold_change, Sample %in% list)
fold_change$Sample <- unlist(lapply(fold_change$Sample, function(x){unlist(str_replace_all(x, "-", "/"))}))
fold_change <- fold_change %>%
  gather(Peptide, value, 2:ncol(fold_change)) %>% 
  spread(Sample, value)
rm(list)
p_adj <- p_adj %>%
  gather(Peptide, value, 2:ncol(p_adj)) %>% 
  spread(Sample, value)
p_adj$Peptide <- unlist(lapply(p_adj$Peptide, function(x){unlist(str_replace_all(x, "\\.Sample_Type", ""))}))

# Then combine
final <- average %>% 
  full_join(stdev, by = "Peptide") %>% 
  full_join(fold_change, by = "Peptide") %>%
  full_join(p_adj, by = "Peptide")

# You can write this write to .csv using the command below if you'd like (we'll create an excel file later)
# write.csv(final, file = "filname.csv")

# You can now use this dataset to make graphs in MatLab, Excel or R.

########## Step 2: Deconvoluting peptides ###########################################
#####################################################################################

deconvoluted <- groups %>%  
  ungroup() %>% 
  mutate(H3_K9me1 = select(groups, contains("K9me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_K9me2 = select(groups, contains("K9me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_K9me3 = select(groups, contains("K9me3")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_K27me1 = select(groups, contains("H3K27me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_K27me2 = select(groups, contains("H3K27me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_K27me3 = select(groups, contains("H3K27me3")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_K36me1 = select(groups, matches("H3K.*K36me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_K36me2 = select(groups, matches('H3K.*K36me2')) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_K36me3 = select(groups, matches("H3K.*K36me3")) %>%
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_K27K36un = select(groups, matches("H3K27unK36un")) %>%
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_K27me1 = select(groups, contains("H33K27me1")) %>% 
           rowSums(na.rm = TRUE)) %>% 
  mutate(H33_K27me2 = select(groups, contains("H33K27me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_K27me3 = select(groups, contains("H33K27me3")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_K36me1 = select(groups, matches("H33.*K36me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_K36me2 = select(groups, matches("H33.*K36me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_K36me3 = select(groups, matches("H33.*K36me3")) %>% 
           rowSums(na.rm = TRUE)) %>% 
  mutate(H33_K27K36un = select(groups, matches("H33K27unK36un")) %>%
           rowSums(na.rm = TRUE))

deconvoluted <- deconvoluted %>% 
  select(-matches("H3K9me(1|2|3)")) %>% 
  select(-matches("H3K(27|36)(un|me1|me2|me3)")) %>% 
  select(-matches("H33K(27|36)(un|me1|me2|me3)"))

# This is also a useful file to look at in other formats
# write.csv(deconvoluted, "filename.csv")

# Now that the data is deconvoluted, you can calculate average, standard dev, etc for your sample sets
# This will be very similar to the averaging you did above
average_deconv <- deconvoluted %>% 
  group_by(Sample_Type) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  select(-c(Sample))
# You'll get a warning on this command, it just means that R couldn't average your sample names, which is probably good since they're letter and not numbers

# Next calculate standard deviations for each set of replicates
stdev_deconv <- deconvoluted %>% 
  group_by(Sample_Type) %>% 
  summarise_all(sd, na.rm = TRUE) %>% 
  select(-c(Sample))
stdev_deconv$Sample_Type <- unlist(lapply(stdev_deconv$Sample_Type, function(x){unlist(str_replace_all(x, "$", "_Sd"))}))

# Now you'll caluculate fold changes. It's easiest to transpose the data before you do this.
average_deconv <- average_deconv %>% 
  gather(Peptide, value, 2:ncol(average_deconv)) %>% 
  spread(Sample_Type, value)

# Calculate fold changes
fold_change_deconv <- average_deconv
fold_change_deconv <- fold_change_deconv %>% 
  select_if(is.numeric) %>% 
  map(~./fold_change_deconv %>% select_if(is.numeric)) %>% 
  imap(~set_names(.x, paste0(.y, "-", names(.x)))) %>% 
  bind_cols()
fold_change_deconv <- bind_cols(average_deconv[1], fold_change_deconv)

# Calculate p values for the fold changes you just calculated above
# This requires you to change the format of the table
long_deconv <- deconvoluted %>% 
  select(-c(Sample)) %>% 
  select(Sample_Type, everything())
long_deconv <- long_deconv %>%   
  gather(Peptide, value, 2:ncol(long_deconv)) 
colnames(long_deconv) <- c("Sample_Type", "Peptide", "value")

# Now perform an Anova
# One-way Anova ####################################################################
# Convert to factor
long_2 <- long_deconv
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

# Make into dataframe of lists
long_3 <- split(long_2, paste(long_2$Peptide))

# Temp function for running anova
tmpfn <- function(dat) {
	 fit <- lm(value ~ Sample_Type, dat)
	 rsq <- summary(fit)$r.squared
	 if(rsq < 0.99)
	    pval <- anova(fit)[,"Pr(>F)"][1]
	 else
	    pval <- NA
	 coefs <- t(coef(fit))
	 data.frame(coefs,rsq, pval)
}

long_stat <- long_3 %>%
    map(tmpfn) %>%
    bind_rows(.id = "Peptide")

long_stat <- long_stat %>%
  filter(pval < 0.05)
list <- c(long_stat[,1])
long_deconv <- filter(long_deconv, Peptide %in% list)

# Perform Tukey's HSD to get adjusted pvalues for all possible comparisons
long_2 <- long_deconv
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

long_3 <- split(long_2, paste(long_2$Peptide))

tmpfn <- function(dat) {
  fit <- aov(value ~ Sample_Type, dat)
  Tukey <- TukeyHSD(fit)$Sample_Type
}

long_stat <- long_3 %>%
  map(tmpfn)

p_adj_deconv <- data.frame(long_stat)
p_adj_deconv <- p_adj_deconv %>% 
  select(contains("p.adj"))
names(p_adj_deconv) <- sub(".p.adj", "", names(p_adj_deconv))
p_adj_deconv <- rownames_to_column(p_adj_deconv, "Sample")

# You can use this command to write results of Anova to a .csv if you wish
# write.csv(long_stat, file = "filename_anova2.csv")
### End One-way Anova #################################################################################

# Two-way Anova ######################################################################
# Separate sample type column into seperate condition columns
# This command seperates by a "_", so if you've named each sample set as Factor1_Factor2, this will seperate everything correctly

long_2 <- long_deconv
long_2 <-separate(long_2, Sample_Type, c("Factor1", "Factor2"), sep = "_")

# Convert to factors
long_2$Factor1 <- as.factor(long_2$Factor1)
long_2$Factor2 <- as.factor(long_2$Factor2)

# Make into dataframe of lists
long_3 <- split(long_2, paste(long_2$Peptide))

# Create temp function for anova
# Output is pvalues from F test
tmpfn <- function(dat) {
  fit <- aov(value ~ Factor2 + Factor1  + Factor2:Factor1, dat)
  sum_aov <- summary.aov(fit)
  pval_Factor2 <- sum_aov[[1]][["Pr(>F)"]][1]
  pval_Factor1 <- sum_aov[[1]][["Pr(>F)"]][2]
  pval_pair <- sum_aov[[1]][["Pr(>F)"]][3]
  res <- (fit)$res
  data.frame(pval_Factor2,pval_Factor1, pval_pair, res)
}

long_stat <- long_3 %>%
    map(tmpfn) %>%
    bind_rows(.id = "Peptide")

# Check that residuals are evenly distributed
hist(long_stat$res,main="Histogram of
     residuals",xlab="Residuals")
long_stat <- long_stat %>%
  select(-res) %>%
  unique()

# Two-way anovas can have main effects and interaction
# Significant interaction effect indicates relationship between dependent variable and main effect (Factor 1) depends on other main effect (Factor 2)
# We can filter out those peptides
long_interaction <- long_stat %>%
  filter(pval_pair <= 0.05)
list <- c(long_interaction[,1])

# Plot the interaction
# This script will make an interaction plot for each peptide that showed significant interaction between Factor1 and Factor2
long_2_plot <- long_2
long_2_plot <- filter(long_2, Peptide %in% list)
long_2_plot <- long_2_plot %>%
  group_by(Factor1, Factor2, Peptide) %>%
  summarise_at(c("value"), mean)
# If there's a lot of Peptides in this category (>15), the plot can get a little crowded
# You can use the commands below to separate out peptides belongs in Histone H3, H4, etc.
long_2_plot_H3 <- long_2_plot %>%
  filter(str_detect(Peptide, "H3"))
long_2_plot_H4 <- long_2_plot %>%
  filter(str_detect(Peptide, "H4"))
long_2_plot_H2 <- long_2_plot %>%
  filter(str_detect(Peptide, "H2"))
# Make sure to reference these filtered plots in the ggplot command below

ggplot(long_2_plot) +
  aes(x = Factor1, y = value, color = Factor2) +
  geom_line(aes(group = Factor2)) +
  geom_point() +
  facet_grid(vars(), vars(Peptide))

# Now you can look at peptides with significantly changing main interactions (eg, changing in either Factor1 or Factor2)
long_main <- long_stat %>%
  mutate(pval = apply(long_stat[,2:3], 1, FUN=min))

long_main <- long_main %>%
  filter(pval <= 0.05)

list <- c(long_main[,1])
long_deconv <- filter(long_deconv, Peptide %in% list)

# Perform Tukey's HSD to get adjusted pvalues for all possible comparisons
long_2 <- long_deconv
long_2$Sample_Type <- as.factor(long_2$Sample_Type)

long_3 <- split(long_2, paste(long_2$Peptide))

tmpfn <- function(dat) {
  fit <- aov(value ~ Sample_Type, dat)
  Tukey <- TukeyHSD(fit)$Sample_Type
}

long_stat <- long_3 %>%
  map(tmpfn)

p_adj_deconv <- data.frame(long_stat)
p_adj_deconv <- p_adj_deconv %>%
  select(contains("p.adj"))
names(p_adj_deconv) <- sub(".p.adj", "", names(p_adj_deconv))
p_adj_deconv <- rownames_to_column(p_adj_deconv, "Sample")

# You can use this command to write results of Anova to a .csv if you wish
# write.csv(long_stat, file = "filename_anova2.csv")
### End Two-way Anova #################################################################################
rm(list)
rm(tmpfn)
rm(long_2)
rm(long_3)
rm(long_stat)

# Combine average, standard deviation, fold change, and p-values into one table
# First, make sure all your data is transposed
stdev_deconv <- stdev_deconv %>% 
  gather(Peptide, value, 2:ncol(stdev_deconv)) %>% 
  spread(Sample_Type, value)

# Filter out the extra fold change comparisons we don't need
fold_change_deconv <- fold_change_deconv %>% 
  gather(Sample, value, 2:ncol(fold_change_deconv)) %>% 
  spread(Peptide, value)
list <- c(p_adj_deconv[,1])
fold_change_deconv <- filter(fold_change_deconv, Sample %in% list)
fold_change_deconv$Sample <- unlist(lapply(fold_change_deconv$Sample, function(x){unlist(str_replace_all(x, "-", "/"))}))
fold_change_deconv <- fold_change_deconv %>%
  gather(Peptide, value, 2:ncol(fold_change_deconv)) %>% 
  spread(Sample, value)
rm(list)
p_adj_deconv <- p_adj_deconv %>%
  gather(Peptide, value, 2:ncol(p_adj_deconv)) %>% 
  spread(Sample, value)
p_adj_deconv$Peptide <- unlist(lapply(p_adj_deconv$Peptide, function(x){unlist(str_replace_all(x, "\\.Sample_Type", ""))}))

# Then combine
final_deconv <- average_deconv %>% 
  full_join(stdev_deconv, by = "Peptide") %>% 
  full_join(fold_change_deconv, by = "Peptide") %>%
  full_join(p_adj_deconv, by = "Peptide")

# Write everything to an excel file
write.xlsx(final, file = "Skyline PIC-Pr.xlsx", sheetName = "Final")
write.xlsx(final_deconv, file = "Skyline PIC-Pr.xlsx", sheetName = "Final_deconv", append = TRUE)
write.xlsx(normalized, file = "Skyline PIC-Pr.xlsx", sheetName = "Normalized", append = TRUE)

################### Step 3: Visualizing your results! ###############################################################################################
#####################################################################################################################################################
# You can use either the average or average_deconv dataframe here 
fig <- average %>% 
  gather(Sample, value, 2:ncol(average)) %>% 
  spread(Peptide, value)
fig <- column_to_rownames(fig, var = "Sample")
# It will give you a warning here, it's okay to ignore
fig <- as.data.frame(fig)

# Plot distance matrix. This will show how similar each of your sample sets are to each other. 0 = most similar
res.dist <- get_dist(fig[, -1], method = "spearman")
fviz_dist(res.dist, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"), lab_size = 5, order = FALSE)

# Assess clustering tendency
# A $hopkins_stat result closer to 0 means the data is clusterable, a value over 0.5 says the data is random. 
# The plot should show patterns of similar coloring as well.
fig %>%
  get_clust_tendency(n = (nrow(fig)-1))

# Principal Components Analysis
# First, find out if PCA is even worth it. You want each PC (principal component) to explain as much variance as possible
# Check the "Cumulative Proportion" row after the summary command to see how much of the variance you see with just two PCs. 
# A high cumulative proportion (closer to 1) means that those first two PCs expain most of the variance in your data
pca <- prcomp(fig)
summary(pca)

# Now you can plot the first two principle components
autoplot(prcomp(fig), data = fig, label = TRUE, label.hjust = 1, label.vjust =1)

# Now that you've seen how your sample sets relate to each other, you can figure out which histone PTMs are most important for separating groups
# This is done by looking at the contribution of each PTM to the the principle components, the higher the contribution, the more important those PTMs are to grouping your samples 
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/ has a lot more suggestions on how to look at PCA data, give it a try!
pca_var <- get_pca_var(pca)
# Shows the top 20 PTMs that contribute to your first PC (or Dim.1)
fviz_contrib(pca, choice = "var", axes = 1, top = 20)
# Shows the top 20 PTMs that contribute to your first PC (or Dim.2)
fviz_contrib(pca, choice = "var", axes = 2, top = 20)
# Shows the top 20 PTMs that contribute to both the first and second PCs
fviz_contrib(pca, choice = "var", axes = 1:2, top = 20)
# You can also just look at the contributions for each PC using this command. The values are the percent each PTM contributes to that PC
contributions <-  as.data.frame(pca_var$contrib)

# Now you can cluster the data. First, determine optimal amount of clusters
# Optimal number of clusters is shown by the dotted line.
clust <- average[,2:ncol(average)]
fviz_nbclust(clust, pam, method = "gap_stat")
rm(clust)

# You can also use this method. The optimal number of clusters is shown by the red point. 
d <- dist(fig)
dend <- as.dendrogram(hclust(d))
dend_k <- find_k(dend)
plot(dend_k)

# Once you've determined your optimal number of clusters, write it here
k <- 3

# Plotting pam (more robust kmeans). 
autoplot(pam(fig, k), data = fig, label = TRUE)

# Plotting hierarchical clusters
res.hc <- fig %>% 
  dist(method = "euclidean") %>% 
  hclust(method = "ward.D2")
fviz_dend(res.hc, k = k, cex = 0.5, color_labels_by_k = TRUE, rect = TRUE)

# Plot heatmap with clustering
dend <- fig %>% 
  dist %>% 
  hclust %>% 
  as.dendrogram() %>% 
  color_branches(k = k)
superheat(fig, left.label.text.size = 2, bottom.label.text.size = 1.5, pretty.order.rows = TRUE,
          bottom.label.text.angle = 90, grid.hline.col = "white", force.grid.vline = TRUE, grid.vline.col = "white", 
          row.dendrogram = TRUE)

# Plot an INTERACTIVE heatmap (so fancy!). You can change the margins by changing l and b.
heatmaply(fig, Rowv = dend, color = viridis(200), k_row = k, Colv = NULL, xlab = "Peptide", ylab = "Conditions", grid_gap = 0.5, fontsize_col = 6) %>% 
  layout(margin = list(l = 170, b = 150))

# It can also be useful to plot heatmaps using fold change data
# To do so, tranpose the fold_change or fold_change_deconv dataframe
fig_fold <- fold_change %>% 
  gather(Sample, value, 2:ncol(fold_change)) %>% 
  spread(Peptide, value)
fig_fold <- column_to_rownames(fig_fold, var = "Sample")
fig_fold <- as.data.frame(fig_fold)
# Log transform
fig_fold <- log2(fig_fold)

# Plot heatmap (without clustering this time)
superheat(fig_fold, left.label.text.size = 2, bottom.label.text.size = 1.5, pretty.order.rows = TRUE,
          bottom.label.text.angle = 90, grid.hline.col = "white", force.grid.vline = TRUE, grid.vline.col = "white")

# Or interactive heatmap
heatmaply(fig_fold, Rowv = NULL, color = viridis(200), Colv = NULL, xlab = "Peptide", ylab = "Conditions", grid_gap = 0.5, key.title = "log2(Fold Change)", fontsize_col = 6) %>% 
  layout(margin = list(l = 170, b = 150))

# If you have a bunch of conditions, these heatmaps can become complicated
# If you only want to look at a subset of conditions, you can select specific rows using this command
fig_fold <- rownames_to_column(fig_fold, "Sample")
fig_fold <- fig_fold %>% 
  filter(str_detect(Sample, "condition1/condition2|condition3/condition1|etc"))
fig_fold <- column_to_rownames(fig_fold, var = "Sample")
fig_fold <- as.data.frame(fig_fold)
# You can now rerun the commands above
# Have fun analyzing your data!
