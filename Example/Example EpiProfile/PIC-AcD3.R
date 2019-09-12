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
# Before you begin, you'll need to save the histone_ratios file from EpiProfile as a .csv file

# Reading in all the data. You'll need to change the filename to whatever matches your file 
# You'll get a few warnings doing this, they're okay to ignore
data <- read.csv("histone_ratios_cell_line_basic.csv", skip = 1)

# This spreadsheet is in a weird format, and we really only want the total area in the middle columns 
# So we filter out the rest of it 
data <- data %>% 
  select(1, contains("Area")) %>% 
  filter(!is.na(Area))

# Now add your sample names
# This script will pull the file names from the EpiProfiler Results
labels <- read.csv("histone_ratios_cell_line_basic.csv", header = FALSE)
labels <- labels[1,2:ncol(data)]
# Now you'll chose a subset of these labels to use as your sample IDs. The command below uses the last parts of each word to create a label
# You can choose how long you want this label to be by changing the first number (-5).
# -5 keeps the last 5 letters, -11 would keep the last 11 letters, etc. 
labels <- unlist(sapply(labels, function(x){unlist(str_sub(x, -6, -1))}))
colnames(data)[-1] <- labels
rm(labels)
colnames(data) <- sub('[.]', '_', make.names(colnames(data), unique = TRUE))

# Make the peptide names something R can read
data <-  as.data.frame(data)
data$Peptide <- unlist(lapply(data$Peptide, function(x){unlist(str_replace_all(x, "[ :().]", "_"))}))
data <- as_tibble(data)

########################### You can skip the next section if you're only using PIC-Pr labeling ###########################################
# But you should do Ac-D3 labeling because it's cool
# Save your basic Ac-D3 dataset
data_basic <- data

# Reading in the Ac-D3 multi dataset. You'll need to change the filename to whatever matches your file 
# You'll also get a few warnings doing this, they're okay to ignore
data <- read.csv("histone_ratios_cell_line_multi.csv", skip = 1)

# This excel spreadsheet is in a weird format, and we really only want the total area in the middle columns 
# So we filter out the rest of it 
data <- data %>% 
  select(1, contains("Area")) %>% 
  filter(!is.na(Area))

# Now add your sample names
# This script will pull the file names from the EpiProfiler Results
labels <- read.csv("histone_ratios_cell_line_multi.csv", header = FALSE)
labels <- labels[1,2:ncol(data)]
# Now you'll chose a subset of these labels to use as your sample IDs. The command below uses the last parts of each word to create a label
# You can choose how long you want this label to be by changing the first number (-5).
# -5 keeps the last 5 letters, -11 would keep the last 11 letters, etc. 
labels <- unlist(sapply(labels, function(x){unlist(str_sub(x, -6, -1))}))
colnames(data)[-1] <- labels
rm(labels)
colnames(data) <- sub('[.]', '_', make.names(colnames(data), unique = TRUE))

# Make the peptide names something R can read
data <-  as.data.frame(data)
data$Peptide <- unlist(lapply(data$Peptide, function(x){unlist(str_replace_all(x, "[ :().]", "_"))}))
data <- as_tibble(data)

# Save the results seperately
data_multi <- data

# Now join them. In general, I use the basic dataset values for any duplicate mods
# I've checked and they both call the same peaks, so using one or the other is arbitrary
data <- anti_join(data_multi, data_basic, by = "Peptide")
data_basic <- data_basic %>% 
  gather(Sample, value, 2:ncol(data_basic)) %>% 
  spread(Peptide, value)
data <- data %>% 
  gather(Sample, value, 2:ncol(data)) %>% 
  spread(Peptide, value)
data <- inner_join(data, data_basic, by = "Sample")
data <- data %>% 
  gather(Peptide, value, 2:ncol(data)) %>% 
  spread(Sample, value)

# Remove a couple duplicate rows that didn't get filtered out correctly
data <- data %>% 
  filter(!str_detect(Peptide, "4_17_(1|2)ac")) %>% 
  filter(!str_detect(Peptide, "9_17_K(9|14)hi"))
rm(data_multi)
rm(data_basic)

######################### You can start here again if you're just doing PIC-Pr (boring) #############################################
######################## Data check to make sure your data is good enough to run with EpiProfile ##################################
# Filter out representative mods that should be consistently identified in your samples
data_check <- data %>% 
  filter(str_detect(Peptide, "H3_18_26_K(18|23)ac|H3_27_40_K27me|H33_27_40_K27me|H3_9_17_K(9|14)ac|H4_4_17_unmod|H3_3_8_unmod|H3_73_83_K79me1"))

# Transpose
data_check <- data_check %>% 
  gather(Sample, value, 2:ncol(data_check)) %>% 
  spread(Peptide, value)

# Check how many of these were identified in each sample
data_check <- data_check %>%
  mutate(zeroes = rowSums(data_check==0))

# All of these mods should be identified in good quality samples, although it's okay if one or two are missing
# This dataframe filters all the samples that have more than two mods unidentified
bad_samples <- data_check %>% 
  filter(zeroes > 2) %>% 
  select(Sample, zeroes)

# If the number of bad samples is small (less than 5% of your samples) 
# AND you have enough biological replicates to remove them without compromising your statistics, you can proceed with this script
# Otherwise, you may need to re-rerun your samples or analyze them in Skyline

# NOTE: H33 is difficult to identify, especially in tissue samples.
# If all of the unidentified mods are on H33, it may be possible to continue analysis, but all H33 mods should be deleted. 
# You can use the following command to remove H33 from your dataset:
# data <- data %>% 
#   filter(!str_detect(Peptide, "H33"))
# After removing these mods, rerun the data-check steps 

# To remove bad samples from your analyses, first transpose the data 
data <- data %>% 
  gather(Sample, value, 2:ncol(data)) %>% 
  spread(Peptide, value)

# Then remove bad samples
data <- data %>% 
  filter(!data$Sample %in% bad_samples$Sample)
rm(data_check)
rm(bad_samples)

# Now you can start filtering your data ######################################################################################################################
# Combine mods that are difficult to disambiguate
# If you just did propionyl labeling, use this command
# data <- data %>% 
#   mutate(H3_9_17_K9K14ac = apply(data[c("H3_9_17_K9ac", "H3_9_17_K14ac")], 1, sum, na.rm = TRUE),
#          H3_18_26_K18K23ac = apply(data[c("H3_18_26_K18ac", "H3_18_26_K23ac")], 1, sum, na.rm = TRUE),
#          H2A1_4_11_K5K9ac = apply(data[c("H2A1_4_11_K5ac", "H2A1_4_11_K9ac")], 1, sum, na.rm = TRUE),
#          H2AX_4_11_K5K9ac = apply(data[c("H2AX_4_11_K5ac", "H2AX_4_11_K9ac")], 1, sum, na.rm = TRUE),
#          H2AJ_4_11_K5K9ac = apply(data[c("H2AJ_4_11_K5ac", "H2AJ_4_11_K9ac")], 1, sum, na.rm = TRUE)) %>% 
#   select(-matches("H3_18_26_K18ac$")) %>% 
#   select(-matches("H3_18_26_K23ac$")) %>%
#   select(-matches("H3_9_17_K9ac$")) %>%
#   select(-matches("H3_9_17_K14ac$")) %>% 
#   select(-matches("H2A1_4_11_K5ac$")) %>% 
#   select(-matches("H2A1_4_11_K9ac$")) %>% 
#   select(-matches("H2AX_4_11_K5ac$")) %>% 
#   select(-matches("H2AX_4_11_K9ac$")) %>% 
#   select(-matches("H2AJ_4_11_K5ac$")) %>% 
#   select(-matches("H2AJ_4_11_K9ac$"))

# If you're doing Ac-D3 labeling, use this one
data <- data %>% 
  mutate(H3_9_17_K9K14ma = apply(data[c("H3_9_17_K9ma", "H3_9_17_K14ma")], 1, sum, na.rm = TRUE),
         H3_9_17_K9K14ac = apply(data[c("H3_9_17_K9ac", "H3_9_17_K14ac")], 1, sum, na.rm = TRUE),
         H3_9_17_K9K14pr = apply(data[c("H3_9_17_K9pr", "H3_9_17_K14pr")], 1, sum, na.rm = TRUE),
         H3_9_17_K9K14gl = apply(data[c("H3_9_17_K9gl", "H3_9_17_K14gl")], 1, sum, na.rm = TRUE),
         H3_9_17_K9K14cr = apply(data[c("H3_9_17_K9cr", "H3_9_17_K14cr")], 1, sum, na.rm = TRUE),
         H3_9_17_K9K14su = apply(data[c("H3_9_17_K9su", "H3_9_17_K14su")], 1, sum, na.rm = TRUE),
         H3_9_17_K9K14bu = apply(data[c("H3_9_17_K9bu", "H3_9_17_K14bu")], 1, sum, na.rm = TRUE),
         H3_18_26_K18K23ac = apply(data[c("H3_18_26_K18ac", "H3_18_26_K23ac")], 1, sum, na.rm = TRUE),
         H3_18_26_K18K23pr = apply(data[c("H3_18_26_K18pr", "H3_18_26_K23pr")], 1, sum, na.rm = TRUE),
         H3_18_26_K18K23cr = apply(data[c("H3_18_26_K18cr", "H3_18_26_K23cr")], 1, sum, na.rm = TRUE),
         H3_18_26_K18K23gl = apply(data[c("H3_18_26_K18gl", "H3_18_26_K23gl")], 1, sum, na.rm = TRUE),
         H3_18_26_K18K23hi = apply(data[c("H3_18_26_K18hi", "H3_18_26_K23hi")], 1, sum, na.rm = TRUE),
         H3_18_26_K18K23ma = apply(data[c("H3_18_26_K18ma", "H3_18_26_K23ma")], 1, sum, na.rm = TRUE),
         H3_18_26_K18K23su = apply(data[c("H3_18_26_K18su", "H3_18_26_K23su")], 1, sum, na.rm = TRUE),
         H3_18_26_K18K23bu = apply(data[c("H3_18_26_K18bu", "H3_18_26_K23bu")], 1, sum, na.rm = TRUE),
         H3_27_40_K27K36hi = apply(data[c("H3_27_40_K27hi", "H3_27_40_K36hi")], 1, sum, na.rm = TRUE),
         H3_27_40_K27K36ma = apply(data[c("H3_27_40_K27ma", "H3_27_40_K36ma")], 1, sum, na.rm = TRUE),
         H3_27_40_K27K36cr = apply(data[c("H3_27_40_K27cr", "H3_27_40_K36cr")], 1, sum, na.rm = TRUE),
         H3_27_40_K27K36bu = apply(data[c("H3_27_40_K27bu", "H3_27_40_K36bu")], 1, sum, na.rm = TRUE),
         H3_27_40_K27K36pr = apply(data[c("H3_27_40_K27pr", "H3_27_40_K36pr")], 1, sum, na.rm = TRUE),
         H2A1_4_11_K5K9ac = apply(data[c("H2A1_4_11_K5ac", "H2A1_4_11_K9ac")], 1, sum, na.rm = TRUE),
         H2AX_4_11_K5K9ac = apply(data[c("H2AX_4_11_K5ac", "H2AX_4_11_K9ac")], 1, sum, na.rm = TRUE),
         H2AJ_4_11_K5K9ac = apply(data[c("H2AJ_4_11_K5ac", "H2AJ_4_11_K9ac")], 1, sum, na.rm = TRUE)) %>% 
  select(-matches("H3_9_17_K(9|14)(ma|pr|gl|cr|su|bu)")) %>% 
  select(-matches("H3_18_26_K(18|23)(pr|cr|gl|hi|ma|su|bu)")) %>% 
  select(-matches("H3_27_40_K(27|36)(ac|hi|ma|cr|bu|pr)")) %>% 
  select(-matches("H3_18_26_K18ac$")) %>% 
  select(-matches("H3_18_26_K23ac$")) %>%
  select(-matches("H3_9_17_K9ac$")) %>%
  select(-matches("H3_9_17_K14ac$")) %>% 
  select(-matches("H2A1_4_11_K5ac$")) %>% 
  select(-matches("H2A1_4_11_K9ac$")) %>% 
  select(-matches("H2AX_4_11_K5ac$")) %>% 
  select(-matches("H2AX_4_11_K9ac$")) %>% 
  select(-matches("H2AJ_4_11_K5ac$")) %>% 
  select(-matches("H2AJ_4_11_K9ac$"))

# Tranpose the data back 
data <- data %>% 
  gather(Peptide, value, 2:ncol(data)) %>% 
  spread(Sample, value)

# Average H4 peptides
H4 <- data %>% 
  filter(str_detect(Peptide, "H4.4.17.*ac"))
H4$Acetyl_count <- str_count(H4$Peptide, "ac")
H4 <- H4 %>%
  select(-Peptide) %>% 
  group_by(Acetyl_count) %>% 
  summarise_all(sum)
H4$Peptide = (c("H4_4_17_1ac", "H4_4_17_2ac", "H4_4_17_3ac", "H4_4_17_4ac"))
H4 <- H4 %>% 
  ungroup() %>% 
  select(-Acetyl_count) %>% 
  select(Peptide, everything())
data <- rbind(H4, data) %>% 
  filter(!str_detect(Peptide, "H4.4.17.K.*ac"))
rm(H4)

# Filter out peptides that weren't idenitifed in many samples
# This command calculates the mean value for each peptide and the number of times it wasn't identified
data <- data %>% 
  mutate(Mean = rowMeans(data[,-1]), 
         Count = rowSums(data == 0))

# We'll filter out any peptides that were at very low abundance, or were present in less than 70% of the samples
data <- data %>% 
  filter(Mean > 1e+6) %>% 
  filter(Count <= (ncol(data)-4)*.3) %>% 
  select(-Mean, -Count)

# Caculating sums for each peptide family
# To do this, first make a label for each peptide family
data$peptide_id <-unlist(lapply(data$Peptide, function(x){
  unlist(str_sub(x, 1, 7))
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
  spread(Peptide, value)

# Now create a column that contains your sample groups
# Basically, you want to split all your biological/technical replicates into seperate groups so you can average and compare them
# This is done by finding a unique ID that matches the names of all the samples in each group and then assigning a number to the samples it matches
# For example, if I wanted to choose all the B_CutC samples from B_CutC_1, B_CutC_2, C_CutC_3, and B_WT_1, I would use "B_CutC" as the identifier
# You can use more than one set of characters for your identifier if you want, just input grepl("ID1|ID2", Sample) (| means or). See regex for more details
# The last number at the end of the command will be assigned to any samples which aren't found by any of your IDs (useful for checking that everything worked)
# Delete or add as many columns as you need to this command, just make sure there are enough parantheses at the end to close it off
groups <- normalized %>% 
  mutate(Sample_Type = ifelse(grepl("HCT116", Sample), "HCT116",
                       ifelse(grepl("HEK293", Sample), "HEK293", 
                       ifelse(grepl("HepG2", Sample), "HepG2", 
                       ifelse(grepl("PANC1", Sample), "PANC1",
                       ifelse(grepl("MCF7", Sample), "MCF7", "set17"))))))

# It's useful to check that each sample is in its appropriate group, you can use this command to pull out just the sample names and groups
peptide_check <- groups %>% 
  select(Sample, Sample_Type) %>% 
  arrange(Sample_Type)
rm(peptide_check)

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

# Identify mods with high %CV within each sample group
# It's important to do this by sample group since you don't want to throw out mods that change between treatments
CV <- stdev[2:ncol(stdev)]/average[2:ncol(average)]
CV <- cbind(average[,1], CV)
CV <- CV %>%
  summarise_all(median, na.rm = TRUE)
# You'll get a warning message here about an argument being not numeric or logical, you can ignore it. 
CV <- CV %>%
  gather(Peptide, value, 2:ncol(CV)) %>% 
  spread(Sample_Type, value)
colnames(CV) <- c("Peptide", "Median")
# This command filters out any mods whose median standard deviation is more than half of its average values in each sample set
# You can make this more or less stringent by changing the value after "Median >" 
CV <- CV %>% 
  filter(Median > 0.5)

# Filter out those mods in the original dataframe
data <- data %>% 
  filter(!Peptide %in% CV$Peptide)
rm(CV)

# Filter out an NA's and replace with small value
data[is.na(data)] <- 0
data[data == 0] <- 1E5

# Then remove any peptides with just one mod identified
data <- data %>%
  group_by(peptide_id) %>% 
  mutate(ids = n())
data <- data %>% 
  filter(ids >= 2) %>% 
  select(-ids)

# Do the normalization again
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
  spread(Peptide, value)

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
normalized <- normalized %>% 
  filter(!Sample %in% c("row_name_1", "row_name_2","etc"))

# Now create a column that contains your sample groups again
# Here, you need to decide if you're going to do a one-way or a two-way anova
# Use one-way if you're just comparing factors to one control (eg, three cell treatments to untreated), use two-way for everything else
# You can change the names ("set01", "set02", etc) to whatever you want if you're running a one-way anova
# If you're running a two-way anova, you'll need to name each sample group by its factors in the format "Factor1_Factor2
# eg. if you're comparing two cell lines and two treatments, the name would be "Line1_Treatment1", then "Line2_Treatment1", etc.
# The following commands will reorder things in alphabetical order, so if you want them in a particular order, make sure the names are alphabetical
groups <- normalized %>% 
  mutate(Sample_Type = ifelse(grepl("HCT116", Sample), "HCT116",
                       ifelse(grepl("HEK293", Sample), "HEK293", 
                       ifelse(grepl("HepG2", Sample), "HepG2", 
                       ifelse(grepl("PANC1", Sample), "PANC1",
                       ifelse(grepl("MCF7", Sample), "MCF7", "set17"))))))

# Now you can start doing calculations. First, average all the replicates in each sample set 
average <- groups %>% 
  group_by(Sample_Type) %>% 
  summarise_all(mean, na.rm = TRUE) %>% 
  select(-c(Sample))
# You'll get a warning on this command, it just means that R couldn't average your sample names, which is probably good since they're letter and not numbers

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
# long_2 <- long
# long_2 <-separate(long_2, Sample_Type, c("Factor1", "Factor2"), sep = "_")
# 
# # Convert to factors
# long_2$Factor1 <- as.factor(long_2$Factor1)
# long_2$Factor2 <- as.factor(long_2$Factor2)
# 
# # Make into dataframe of lists
# long_3 <- split(long_2, paste(long_2$Peptide))
# 
# # Create temp function for anova
# # Output is pvalues from F test
# tmpfn <- function(dat) {
#   fit <- aov(value ~ Factor2 + Factor1  + Factor2:Factor1, dat)
#   sum_aov <- summary.aov(fit)
#   pval_Factor2 <- sum_aov[[1]][["Pr(>F)"]][1]
#   pval_Factor1 <- sum_aov[[1]][["Pr(>F)"]][2]
#   pval_pair <- sum_aov[[1]][["Pr(>F)"]][3]
#   res <- (fit)$res
#   data.frame(pval_Factor2,pval_Factor1, pval_pair, res)
# }
# 
# long_stat <- long_3 %>%
#     map(tmpfn) %>%
#     bind_rows(.id = "Peptide")
# 
# # Check that residuals are evenly distributed
# hist(long_stat$res,main="Histogram of
#      residuals",xlab="Residuals")
# 
# long_stat <- long_stat %>%
#   select(-res) %>%
#   unique()
# 
# # Two-way anovas can have main effects and interaction
# # Significant interaction effect indicates relationship between dependent variable and main effect (Factor 1) depends on other main effect (Factor 2)
# # We can filter out those peptides
# long_interaction <- long_stat %>%
#   filter(pval_pair <= 0.05)
# list <- c(long_interaction[,1])
# 
# # Plot the interaction 
# # This script will make an interaction plot for each peptide that showed significant interaction between Factor1 and Factor2
# long_2_plot <- long_2
# long_2_plot <- filter(long_2, Peptide %in% list)
# long_2_plot <- long_2_plot %>% 
#   group_by(Factor1, Factor2, Peptide) %>% 
#   summarise_at(c("value"), mean)
# # If there's a lot of Peptides in this category (>15), the plot can get a little crowded
# # You can use the commands below to separate out peptides belongs in Histone H3, H4, etc.
# long_2_plot_H3 <- long_2_plot %>% 
#   filter(str_detect(Peptide, "H3"))
# long_2_plot_H4 <- long_2_plot %>% 
#   filter(str_detect(Peptide, "H4"))
# long_2_plot_H2 <- long_2_plot %>% 
#   filter(str_detect(Peptide, "H2"))
# # Make sure to reference these filtered plots in the first line of the ggplot command below
# 
# ggplot(long_2_plot) +
#   aes(x = Factor1, y = value, color = Factor2) + 
#   geom_line(aes(group = Factor2)) +
#   geom_point() +
#   facet_grid(vars(), vars(Peptide))
# 
# # Now you can look at peptides with significantly changing main interactions (eg, changing in either Factor1 or Factor2)
# long_main <- long_stat %>%
#   mutate(pval = apply(long_stat[,2:3], 1, FUN=min))
# 
# long_main <- long_main %>%
#   filter(pval <= 0.05)
# 
# list <- c(long_main[,1])
# long <- filter(long, Peptide %in% list)
# 
# # Perform Tukey's HSD to get adjusted pvalues for all possible comparisons
# long_2 <- long
# long_2$Sample_Type <- as.factor(long_2$Sample_Type)
# 
# long_3 <- split(long_2, paste(long_2$Peptide))
# 
# tmpfn <- function(dat) {
#   fit <- aov(value ~ Sample_Type, dat)
#   Tukey <- TukeyHSD(fit)$Sample_Type
# }
# 
# long_stat <- long_3 %>%
#   map(tmpfn)
# 
# p_adj <- data.frame(long_stat)
# p_adj <- p_adj %>% 
#   select(contains("p.adj"))
# names(p_adj) <- sub(".p.adj", "", names(p_adj))
# p_adj <- rownames_to_column(p_adj, "Sample")

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
  mutate(H3_9_18_K9me1 = select(groups, contains("K9me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_9_18_K9me2 = select(groups, contains("K9me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_9_18_K9me3 = select(groups, contains("K9me3")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K27me1 = select(groups, contains("H3_27_40_K27me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K27me2 = select(groups, contains("H3_27_40_K27me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K27me3 = select(groups, contains("H3_27_40_K27me3")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K36me1 = select(groups, matches("H3_.*K36me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K36me2 = select(groups, matches('H3_.*K36me2')) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H3_27_41_K36me3 = select(groups, matches("H3_.*K36me3")) %>%
           rowSums(na.rm = TRUE)) %>% 
  mutate(H33_27_41_K27me1 = select(groups, contains("H33_27_40_K27me1")) %>% 
           rowSums(na.rm = TRUE)) %>% 
  mutate(H33_27_41_K27me2 = select(groups, contains("H33_27_40_K27me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_27_41_K27me3 = select(groups, contains("H33_27_40_K27me3")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_27_41_K36me1 = select(groups, matches("H33_.*K36me1")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_27_41_K36me2 = select(groups, matches("H33_.*K36me2")) %>% 
           rowSums(na.rm = TRUE)) %>%
  mutate(H33_27_41_K36me3 = select(groups, matches("H33_.*K36me3")) %>% 
           rowSums(na.rm = TRUE))

deconvoluted <- deconvoluted %>% 
  select(-matches("H3_9_17_K9me(1|2|3)")) %>% 
  select(-matches("H3_27_40_K(27|36)(me1|me2|me3)")) %>% 
  select(-matches("H33_27_40_K(27|36)(me1|me2|me3)"))

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

# long_2 <- long_deconv
# long_2 <-separate(long_2, Sample_Type, c("Factor1", "Factor2"), sep = "_")
# 
# # Convert to factors
# long_2$Factor1 <- as.factor(long_2$Factor1)
# long_2$Factor2 <- as.factor(long_2$Factor2)
# 
# # Make into dataframe of lists
# long_3 <- split(long_2, paste(long_2$Peptide))
# 
# # Create temp function for anova
# # Output is pvalues from F test
# tmpfn <- function(dat) {
#   fit <- aov(value ~ Factor2 + Factor1  + Factor2:Factor1, dat)
#   sum_aov <- summary.aov(fit)
#   pval_Factor2 <- sum_aov[[1]][["Pr(>F)"]][1]
#   pval_Factor1 <- sum_aov[[1]][["Pr(>F)"]][2]
#   pval_pair <- sum_aov[[1]][["Pr(>F)"]][3]
#   res <- (fit)$res
#   data.frame(pval_Factor2,pval_Factor1, pval_pair, res)
# }
# 
# long_stat <- long_3 %>%
#     map(tmpfn) %>%
#     bind_rows(.id = "Peptide")
# 
# # Check that residuals are evenly distributed
# hist(long_stat$res,main="Histogram of
#      residuals",xlab="Residuals")
# long_stat <- long_stat %>%
#   select(-res) %>%
#   unique()
# 
# # Two-way anovas can have main effects and interaction
# # Significant interaction effect indicates relationship between dependent variable and main effect (Factor 1) depends on other main effect (Factor 2)
# # We can filter out those peptides
# long_interaction <- long_stat %>%
#   filter(pval_pair <= 0.05)
# list <- c(long_interaction[,1])
# 
# # Plot the interaction 
# # This script will make an interaction plot for each peptide that showed significant interaction between Factor1 and Factor2
# long_2_plot <- long_2
# long_2_plot <- filter(long_2, Peptide %in% list)
# long_2_plot <- long_2_plot %>% 
#   group_by(Factor1, Factor2, Peptide) %>% 
#   summarise_at(c("value"), mean)
# # If there's a lot of Peptides in this category (>15), the plot can get a little crowded
# # You can use the commands below to separate out peptides belongs in Histone H3, H4, etc.
# long_2_plot_H3 <- long_2_plot %>% 
#   filter(str_detect(Peptide, "H3"))
# long_2_plot_H4 <- long_2_plot %>% 
#   filter(str_detect(Peptide, "H4"))
# long_2_plot_H2 <- long_2_plot %>% 
#   filter(str_detect(Peptide, "H2"))
# # Make sure to reference these filtered plots in the ggplot command below
# 
# ggplot(long_2_plot) +
#   aes(x = Factor1, y = value, color = Factor2) + 
#   geom_line(aes(group = Factor2)) +
#   geom_point() +
#   facet_grid(vars(), vars(Peptide))
# 
# # Now you can look at peptides with significantly changing main interactions (eg, changing in either Factor1 or Factor2)
# long_main <- long_stat %>%
#   mutate(pval = apply(long_stat[,2:3], 1, FUN=min))
# 
# long_main <- long_main %>%
#   filter(pval <= 0.05)
# 
# list <- c(long_main[,1])
# long_deconv <- filter(long_deconv, Peptide %in% list)
# 
# # Perform Tukey's HSD to get adjusted pvalues for all possible comparisons
# long_2 <- long_deconv
# long_2$Sample_Type <- as.factor(long_2$Sample_Type)
# 
# long_3 <- split(long_2, paste(long_2$Peptide))
# 
# tmpfn <- function(dat) {
#   fit <- aov(value ~ Sample_Type, dat)
#   Tukey <- TukeyHSD(fit)$Sample_Type
# }
# 
# long_stat <- long_3 %>%
#   map(tmpfn)
# 
# p_adj_deconv <- data.frame(long_stat)
# p_adj_deconv <- p_adj_deconv %>% 
#   select(contains("p.adj"))
# names(p_adj_deconv) <- sub(".p.adj", "", names(p_adj_deconv))
# p_adj_deconv <- rownames_to_column(p_adj_deconv, "Sample")

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
write.xlsx(final, file = "PIC-AcD3.xlsx", sheetName = "Final")
write.xlsx(final_deconv, file = "PIC-AcD3.xlsx", sheetName = "Final_deconv", append = TRUE)
write.xlsx(normalized, file = "PIC-AcD3.xlsx", sheetName = "Normalized", append = TRUE)

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
k <- 2

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
