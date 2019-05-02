# R script to perform analyze stoperator data for Mavrich & Hatfull, mBio, 2019.
# The scripts computes position weight matrices for stoperators and
# compares similarities between PWMs.
# Distinct analyses and code blocks separated by "###".

### 1. Prepare environment.
### 2. 
### 3. 
### 4. 
### 5. 
### 6. 


### 1. Prepare environment.

# Install dependencies
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# biocLite("TFBSTools")
library("Biostrings")
library("TFBSTools")




# Set working directory variables specific to local directory structure.
DIR_INPUT_MISC = "~/scratch/process_stoperator_data/input_misc/"
DIR_INPUT_GENOMES = "~/scratch/process_stoperator_data/input_genomes/"
DIR_OUTPUT = "~/scratch/process_stoperator_data/output/"





#TODO fix variables
#Set paths for all input files.

PHAGE_METADATA_FILENAME = 
  paste(DIR_INPUT_MISC,
        "phage_metadata.csv",
        sep="")

STOPERATOR_DATA_FILENAME = 
  paste(DIR_INPUT_MISC,
        "stoperator_data.csv",
        sep="")








### Import datasets


# Import all 327 Cluster A genomes from the Actino1321 database.
# Read individual fasta files by generating a list of files in the directory,
# and passing that list to the readDNAStringSet function.


setwd(DIR_INPUT_GENOMES)
actino1321_genomes <- readDNAStringSet(list.files())
setwd(DIR_OUTPUT)





# Import list of MEME-identified stoperators.
# Data structure:
# 1. phage name (imported as factor)
# 2. stoperator_id (imported as factor)
# 3. sequence (imported as factor)
meme_sites <- read.csv(STOPERATOR_DATA_FILENAME,sep=",",header=TRUE)



# Tally # of stoperators identified by MEME as input.
meme_freq <- as.data.frame(table(meme_sites$phage))
names(meme_freq) <- c("phage","meme_frequency")





# Import phage metadata.
# This data is derived from the same metadata input used for the
# analyze_immunity_data.R. It contains data for all phages in the Actino1321
# database.
# Data structure:
# 1. phageid (imported as factor)
# 2. size (imported as int)
# 3. coordinate_genome_center (imported as int)
phage_metadata <- read.csv(PHAGE_METADATA_FILENAME,sep=",",header=TRUE)


#Phages not in Cluster A do not have a specific genome center coordinate.
phage_metadata[phage_metadata == "Unspecified"] <- NA
phage_metadata$coordinate_genome_center <- 
  as.numeric(as.character(phage_metadata$coordinate_genome_center))


# QC: Confirm presence of sites in each genomes in case genome sequence has
# changed since original MEME search.

# Create an empty dataframe to store all the search data
bios_sites <- data.frame(phage="empty",
                         forward_seq="empty",
                         reverse_seq="empty",
                         start=0,
                         end=0,
                         strand="empty")


# For each phage, 
# 1. create list of unique 13bp stoperators.
# 2. search the corresponding parent genome, both strands.
# 3. save list of sites and coordinates.
for (phageid in levels(meme_sites$phage)){
  
  
  reduced_table_seqs <- subset(meme_sites,
                               meme_sites$phage == phageid,
                               select = c("sequence"))
  
  reduced_table_seqs$sequence <- factor(reduced_table_seqs$sequence)
  
  for (stoperator in levels(reduced_table_seqs$sequence)){
    
    stop_forward <- DNAString(stoperator)
    stop_reverse <- reverseComplement(stop_forward)
    genome_to_search <- actino1321_genomes[[phageid]]
    
    hits_plus <- matchPattern(stop_forward,genome_to_search)
    hits_minus <- matchPattern(stop_reverse,genome_to_search)
    
    if (length(hits_plus) > 0){
      
      hits_plus_df <- 
        data.frame(phage=phageid,
                   forward_seq = as.character(stop_forward),
                   reverse_seq = as.character(stop_reverse),
                   start=start(hits_plus),
                   end=end(hits_plus),
                   strand="forward")
      bios_sites <- rbind(bios_sites,hits_plus_df)
    }

    if (length(hits_minus) > 0){
      
      hits_minus_df <- 
        data.frame(phage=phageid,
                   forward_seq=as.character(stop_forward),
                   reverse_seq=as.character(stop_reverse),
                   start=start(hits_minus),
                   end=end(hits_minus),
                   strand="reverse")
      bios_sites <- rbind(bios_sites,hits_minus_df)
    }
  }
}


# Remove the row of empty data used to initiate the table.
bios_sites <- bios_sites[bios_sites$phage != "empty",]

bios_sites$phage <- factor(bios_sites$phage)

bios_sites$strand <- factor(bios_sites$strand)

bios_sites$forward_seq <- factor(bios_sites$forward_seq)

bios_sites$reverse_seq <- factor(bios_sites$reverse_seq)



# Create unique identifier for each site.
bios_sites$stop_site_id <- paste(bios_sites$phage,
                                 bios_sites$strand,
                                 bios_sites$start,
                                 sep="_")

bios_sites$stop_site_id <- factor(bios_sites$stop_site_id)


# Add metadata.
bios_sites <- merge(bios_sites,phage_metadata,by.x="phage",by.y="phageid")

bios_sites$dist_from_center <- 
  bios_sites$start - bios_sites$coordinate_genome_center


# Tally # of stoperators confirmed.
bios_freq <- as.data.frame(table(bios_sites$phage))

names(bios_freq) <- c("phage","biostrings_freq")


# Tally # of stoperators on each side of genome center.
bios_sites_left <- subset(bios_sites,bios_sites$dist_from_center <= 0)
bios_sites_right <- subset(bios_sites,bios_sites$dist_from_center > 0)

bios_freq_left <- as.data.frame(table(bios_sites_left$phage))
bios_freq_right <- as.data.frame(table(bios_sites_right$phage))

names(bios_freq_left) <- c("phage","biostrings_freq_left")
names(bios_freq_right) <- c("phage","biostrings_freq_right")


# QC: should equal 0.
nrow(bios_sites) - nrow(bios_sites_left) - nrow(bios_sites_right)




# Tally # of stoperators on each strand of each side of genome center.
bios_sites_left_for <- subset(bios_sites_left,
                              bios_sites_left$strand == "forward")
bios_sites_left_rev <- subset(bios_sites_left,
                              bios_sites_left$strand == "reverse")
bios_sites_right_for <- subset(bios_sites_right,
                               bios_sites_right$strand == "forward")
bios_sites_right_rev <- subset(bios_sites_right,
                               bios_sites_right$strand == "reverse")


bios_freq_left_for <- as.data.frame(table(bios_sites_left_for$phage))
bios_freq_left_rev <- as.data.frame(table(bios_sites_left_rev$phage))
bios_freq_right_for <- as.data.frame(table(bios_sites_right_for$phage))
bios_freq_right_rev <- as.data.frame(table(bios_sites_right_rev$phage))


names(bios_freq_left_for) <- c("phage","biostrings_freq_left_forward")
names(bios_freq_left_rev) <- c("phage","biostrings_freq_left_reverse")
names(bios_freq_right_for) <- c("phage","biostrings_freq_right_forward")
names(bios_freq_right_rev) <- c("phage","biostrings_freq_right_reverse")




# Combine MEME and Biostrings data.
freq_table <- merge(meme_freq,bios_freq,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_left,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_right,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_left_for,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_left_rev,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_right_for,by.x="phage",by.y="phage")
freq_table <- merge(freq_table,bios_freq_right_rev,by.x="phage",by.y="phage")



freq_table$meme_biostrings_diff <- 
  freq_table$meme_frequency - freq_table$biostrings_freq

freq_table$biostrings_percent_left <- 
  freq_table$biostrings_freq_left / freq_table$biostrings_freq

freq_table$biostrings_percent_right <- 
  freq_table$biostrings_freq_right / freq_table$biostrings_freq

freq_table$biostrings_percent_left_forward <- 
  freq_table$biostrings_freq_left_forward / freq_table$biostrings_freq_left

freq_table$biostrings_percent_left_reverse <- 
  freq_table$biostrings_freq_left_reverse / freq_table$biostrings_freq_left

freq_table$biostrings_percent_right_forward <- 
  freq_table$biostrings_freq_right_forward / freq_table$biostrings_freq_right

freq_table$biostrings_percent_right_reverse <- 
  freq_table$biostrings_freq_right_reverse / freq_table$biostrings_freq_right


# QC: Tally checks should equal 0.
freq_table$biostrings_check1 <- 
  freq_table$biostrings_freq - 
  freq_table$biostrings_freq_left - 
  freq_table$biostrings_freq_right

freq_table$biostrings_check2 <- 
  freq_table$biostrings_freq_left - 
  freq_table$biostrings_freq_left_forward - 
  freq_table$biostrings_freq_left_reverse


freq_table$biostrings_check3 <- 
  freq_table$biostrings_freq_right - 
  freq_table$biostrings_freq_right_forward - 
  freq_table$biostrings_freq_right_reverse

summary(freq_table$biostrings_check1)
summary(freq_table$biostrings_check2)
summary(freq_table$biostrings_check3)



# QC: Percent checks should equal 1.
freq_table$biostrings_check4 <- 
  freq_table$biostrings_percent_left +
  freq_table$biostrings_percent_right

freq_table$biostrings_check5 <- 
  freq_table$biostrings_percent_left_forward +
  freq_table$biostrings_percent_left_reverse

freq_table$biostrings_check6 <- 
  freq_table$biostrings_percent_right_forward +
  freq_table$biostrings_percent_right_reverse


summary(freq_table$biostrings_check4)
summary(freq_table$biostrings_check5)
summary(freq_table$biostrings_check6)








#HERE
# QC: 

nrow(meme_sites)
nrow(bios_sites)

# If this does not equal 0, then there is a difference between the number
# of stoperators present in the two sets of stoperators.
nrow(meme_sites) - nrow(bios_sites)



#TODO delete or refine
hist(freq_table$meme_biostrings_diff)
setdiff(levels(meme_freq$phage),
        levels(bios_freq$phage))
setdiff(levels(bios_freq$phage),
        levels(meme_freq$phage))



# MEME sites tally = 9266. Biostrings sites tally = 9273.
# Biostrings contains 7 more stoperators. Overall,  only a few genomes contain
# small differences in stoperators than expected. Seven phages contain 1 fewer
# stoperator from MEME than Biostrings:
# BeesKnees (A1)
# Changeling (A2)
# Doom (A1)
# KBG (A1)
# Petruchio (A1)
# Saintus (A8)
# SarFire (A1)


# Compare # MEME sites versus # Biostrings sites per genome.
plot(freq_table$meme_frequency,
     freq_table$biostrings_freq)



# Plot distribution of number of stoperators.

par(mar=c(4,8,8,4))
hist(table(bios_sites$phage),col="black",
     xlim=c(0,50),ylim=c(0,80),
     breaks=10,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,paste(DIR_OUTPUT,
                   'stoperator327_number_stops_per_genome.pdf',
                   sep=""))
dev.off()

# Output biostrings data for reference
write.table(bios_sites,
            paste(DIR_OUTPUT,"stoperator_biostrings_coordinates.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)









###
###
###
###
###
###
###
###
###
### Create PWMs
# Use the sites that were identified in each genome by Biostrings search.
# For each phage in the table, create a pwm.
# Create two lists of PWMs.

# Log2ProbRatio PWM.
stoperator_pwm_list <- vector("list",nlevels(bios_sites$phage))


# Prob PWM.
stoperator_pwm2_list <- vector("list",nlevels(bios_sites$phage))


# For each phage genome, the list of stoperators must contain at least 
# one of each of the 4 nucleotides A,C,T,G. Otherwise, when PFMatrix is called,
# an error is encountered. The for loop will probably crash if one or more
# PWMs does not meet this criteria, but a QC step has been used to confirm it.
stoperator_letters_check_sum <- 0



count <- 1
for (stoperator_phage in levels(bios_sites$phage)){
  
  reduced_table_seqs <- subset(bios_sites,
                               bios_sites$phage == 
                                 stoperator_phage,
                               select = c("forward_seq"))
  
  reduced_table_seqs <- as.character(reduced_table_seqs$forward_seq)
  
  consensus_matrix <- consensusMatrix(reduced_table_seqs)
  stoperator_pfm <- PFMatrix(ID=stoperator_phage,
                                  name=stoperator_phage,
                                  profileMatrix = consensus_matrix)
  
  stoperator_pwm <- toPWM(stoperator_pfm,type="log2probratio") 
  stoperator_pwm2 <- toPWM(stoperator_pfm,type="prob") 
  
  
  
  # Append the pwm to the list of pwms.
  stoperator_pwm_list[[count]] <- stoperator_pwm
  stoperator_pwm2_list[[count]] <- stoperator_pwm2  
  
  
  # Append DNA letter frequency check.
  reduced_table_seqs_dnastringset <- DNAStringSet(reduced_table_seqs)
  all_4_nucleotides <- 
    ifelse(length(uniqueLetters(reduced_table_seqs_dnastringset)) == 4,1,0)
  
  stoperator_letters_check_sum <- 
    stoperator_letters_check_sum + all_4_nucleotides
  
  count <- count + 1
  
}  


# QC: should equal 0.
nlevels(bios_sites$phage) - stoperator_letters_check_sum









### Compute Log2ProbRatio PWM distances
# This log2-based score seems to be commonly used instead of a linear score.
stoperator_pwm_list1 <- stoperator_pwm_list
stoperator_pwm_list2 <- stoperator_pwm_list


pwm_output_list <- 
  vector("list",length(stoperator_pwm_list1)*length(stoperator_pwm_list2))

count <- 1
for (pwm1 in stoperator_pwm_list1){
  for (pwm2 in stoperator_pwm_list2){
    
    dist_pearson <- 1 - PWMSimilarity(pwm1,
                                      pwm2,
                                      method="Pearson")
    
    dist_euc <- PWMSimilarity(pwm1,
                              pwm2,
                              method="Euclidean")
    
    # Append the distance data to the output list.
    pwm_output_list[[count]] <- c(ID(pwm1),
                                  ID(pwm2),
                                  round(dist_pearson,2),
                                  round(dist_euc,2))
    count <- count + 1
  }
}



pwm_output_df <- as.data.frame(t(as.data.frame(pwm_output_list)))
rownames(pwm_output_df) <- NULL
names(pwm_output_df) <- c("phage1","phage2","dist_pearson","dist_euc")


# Export data
setwd("~/scratch/stoperator_distance/output/")
write.table(pwm_output_df,
            paste(DIR_OUTPUT,"stoperator_pwm_distances.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)








# QC: Compare the euclidean distance of PWM probabilities and the
# log2probratio probabilities.

pwm_output_df$phage1_phage2 <- paste(pwm_output_df$phage1,
                                     "_",
                                     pwm_output_df$phage2,
                                     sep="")

pwm_output_df$dist_pearson <- as.numeric(as.character(pwm_output_df$dist_pearson))
pwm_output_df$dist_euc <- as.numeric(as.character(pwm_output_df$dist_euc))







### Compute Prob PWM distances.
stoperator_pwm2_list1 <- stoperator_pwm2_list
stoperator_pwm2_list2 <- stoperator_pwm2_list


pwm2_output_list <- 
  vector("list",length(stoperator_pwm2_list1)*length(stoperator_pwm2_list2))

count <- 1
for (pwm1 in stoperator_pwm2_list1){
  for (pwm2 in stoperator_pwm2_list2){
    
    dist_pearson <- 1 - PWMSimilarity(pwm1,
                                      pwm2,
                                      method="Pearson")
    
    dist_euc <- PWMSimilarity(pwm1,
                              pwm2,
                              method="Euclidean")

    # Append the distance data to the output list.
    pwm2_output_list[[count]] <- c(ID(pwm1),
                                   ID(pwm2),
                                   round(dist_pearson,2),
                                   round(dist_euc,2))
    count <- count + 1
  }
}



pwm2_output_df <- as.data.frame(t(as.data.frame(pwm2_output_list)))
rownames(pwm2_output_df) <- NULL
names(pwm2_output_df) <- 
  c("prob_phage1","prob_phage2","prob_dist_pearson","prob_dist_euc")

pwm2_output_df$phage1_phage2 <- paste(pwm2_output_df$prob_phage1,
                                      "_",
                                      pwm2_output_df$prob_phage2,
                                      sep="")

pwm2_output_df$prob_dist_pearson <- 
  as.numeric(as.character(pwm2_output_df$prob_dist_pearson))

pwm2_output_df$prob_dist_euc <- 
  as.numeric(as.character(pwm2_output_df$prob_dist_euc))






# Merge both PWM lists.

combined_pwm_output_df <- merge(pwm_output_df,
                                pwm2_output_df,
                                by.x="phage1_phage2",
                                by.y="phage1_phage2")


plot(combined_pwm_output_df$dist_euc,
     combined_pwm_output_df$prob_dist_euc)

plot(combined_pwm_output_df$dist_euc,
     combined_pwm_output_df$dist_pearson)

plot(combined_pwm_output_df$prob_dist_euc,
     combined_pwm_output_df$prob_dist_pearson)
abline(0,1)

#TODO
#Summary of comparisons


###
###
###
###
###
###
###
###
###
### After creating all stoperator PWMs, iterate through all genomes for each
# PWM and predict sites.

# Convert object type of stoperator PWMs from standard list to PWMatrixList.
stoperator_pwm_matrixlist <- do.call(PWMatrixList,stoperator_pwm_list)


stop_site_set80_list <- searchSeq(stoperator_pwm_matrixlist,actino1321_genomes)
stop_site_set80_df <- as(stop_site_set80_list,"data.frame")
# Note: it takes about an hour to run the two commands above using stop327 data.




# Dataframe output: each row is a predicted site in a genome from one PWM
# Column names:
# 1. seqnames = factor, 327 levels are individual phage genomes
# 2. source = factor, 1 level ("TFBS")
# 3. feature = factor, 1 level ("TFBS")
# 4. start = int
# 5. end = int
# 6. absScore = number (from ~4 to ~25)
# 7. relScore = number (from 0.8 to 1). It reflects the range of possible
#               scores that could be obtained with the motif and is
#               normalized between 0 and 1.
# 8. strand = factor, 2 levels ("+","-")
# 9. ID = factor, 327 levels are stoperator PWM names. Identical to "TF"
# 10. TF = factor, 327 levels are stoperator PWM names. Identical to "ID"
# 11. class = factor, 1 level ("Unknown")
# 12. siteSeqs = chr, 13bp stoperator DNA sequence




#TODO 
# QC
# stop_site_set95_list <- searchSeq(stoperator_pwm_matrixlist,actino1321_genomes,min.score = "95%")
# stop_site_set95_df <- as(stop_site_set95_list,"data.frame")


# actino1319_analysis:
# For 80%: 1797889 rows. 
# For 95%: 350763 rows.


# Reduce data for analysis with immunity data
stoperator_site_predictions80 <- subset(stop_site_set80_df,
                                     select=c("seqnames",
                                              "start",
                                              "end",
                                              "strand",
                                              "siteSeqs",
                                              "absScore",
                                              "relScore",
                                              "ID"))

names(stoperator_site_predictions80) <- c("stoperator_target",
                                       "site_start",
                                       "site_end",
                                       "site_strand",
                                       "site_seq",
                                       "site_abs_score",
                                       "site_rel_score",
                                       "stoperator_motif")
  

stoperator_site_predictions80$motif_target <- 
  paste(stoperator_site_predictions80$stoperator_motif,
        "_",
        stoperator_site_predictions80$stoperator_target,
        sep="")

stoperator_site_predictions80$motif_target <- 
  as.factor(stoperator_site_predictions80$motif_target)


#TODO
# QC: Check number of unique comparisons.
# This list has 106,300 levels, but in the Actino1321 database there are
# 106,929 unique possible pairwise combinations (327 PWMs * 327 target genomes).


all_levels <- 
  expand.grid(stoperator_target = 
                as.character(levels(stoperator_site_predictions80$stoperator_target)),
              stoperator_motif = 
                as.character(levels(stoperator_site_predictions80$stoperator_motif)))

all_levels$motif_target <- paste(all_levels$stoperator_motif,
                                 "_",
                                 all_levels$stoperator_target,
                                 sep="")

all_levels$motif_target <- as.factor(all_levels$motif_target)


# Re-factor motif_target using all possible pairwise combinations.
# This will enable quantification for all combinations, 
# even those that are missing.
stoperator_site_predictions80$motif_target <- 
  factor(stoperator_site_predictions80$motif_target,
         levels(all_levels$motif_target))
# Now there are 106,929 levels.




# Create unique site ids
stoperator_site_predictions80$site_strand2 <- 
  ifelse(stoperator_site_predictions80$site_strand == "+",
         "forward",
         "reverse")
stoperator_site_predictions80$site_strand2 <- 
  factor(stoperator_site_predictions80$site_strand2)

stoperator_site_predictions80$stop_site_id <- 
  paste(stoperator_site_predictions80$stoperator_target,
        stoperator_site_predictions80$site_strand2,
        stoperator_site_predictions80$site_start,
        sep="_")

stoperator_site_predictions80$stop_site_id <- 
  factor(stoperator_site_predictions80$stop_site_id)



#QC
#temp_l5 <- subset(stoperator_site_predictions80,stoperator_site_predictions80$stoperator_target == "l5" & stoperator_site_predictions80$stoperator_motif == "l5")






# Table creates a two column table = the levels and the frequency (int)
stoperator_site_predictions80_self <- 
  subset(stoperator_site_predictions80,
         as.character(stoperator_site_predictions80$stoperator_target) == 
           as.character(stoperator_site_predictions80$stoperator_motif))


stoperator_site_predictions80_self$stoperator_target <- 
  factor(stoperator_site_predictions80_self$stoperator_target)
stoperator_site_predictions80_self$stoperator_motif <- 
  factor(stoperator_site_predictions80_self$stoperator_motif)
stoperator_site_predictions80_self$motif_target <- 
  factor(stoperator_site_predictions80_self$motif_target)
stoperator_site_predictions80_self$site_strand2 <- 
  factor(stoperator_site_predictions80_self$site_strand2)

stoperator_site_predictions85_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.85)
stoperator_site_predictions86_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.86)
stoperator_site_predictions87_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.87)
stoperator_site_predictions88_self <- subset(stoperator_site_predictions80_self,
                                            stoperator_site_predictions80_self$site_rel_score >= 0.88)
stoperator_site_predictions89_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.89)
stoperator_site_predictions90_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.90)
stoperator_site_predictions91_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.91)
stoperator_site_predictions92_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.92)
stoperator_site_predictions93_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.93)
stoperator_site_predictions94_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.94)
stoperator_site_predictions95_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 0.95)
stoperator_site_predictions100_self <- subset(stoperator_site_predictions80_self,
                                             stoperator_site_predictions80_self$site_rel_score >= 1)




stoperator_site_predictions80_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions80_self$stoperator_target))
names(stoperator_site_predictions80_self_tfbstools_frequency) <- c("phage","tfbstools80_frequency")

stoperator_site_predictions85_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions85_self$stoperator_target))
names(stoperator_site_predictions85_self_tfbstools_frequency) <- c("phage","tfbstools85_frequency")

stoperator_site_predictions86_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions86_self$stoperator_target))
names(stoperator_site_predictions86_self_tfbstools_frequency) <- c("phage","tfbstools86_frequency")

stoperator_site_predictions87_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions87_self$stoperator_target))
names(stoperator_site_predictions87_self_tfbstools_frequency) <- c("phage","tfbstools87_frequency")

stoperator_site_predictions88_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions88_self$stoperator_target))
names(stoperator_site_predictions88_self_tfbstools_frequency) <- c("phage","tfbstools88_frequency")

stoperator_site_predictions89_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions89_self$stoperator_target))
names(stoperator_site_predictions89_self_tfbstools_frequency) <- c("phage","tfbstools89_frequency")

stoperator_site_predictions90_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions90_self$stoperator_target))
names(stoperator_site_predictions90_self_tfbstools_frequency) <- c("phage","tfbstools90_frequency")

stoperator_site_predictions91_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions91_self$stoperator_target))
names(stoperator_site_predictions91_self_tfbstools_frequency) <- c("phage","tfbstools91_frequency")

stoperator_site_predictions92_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions92_self$stoperator_target))
names(stoperator_site_predictions92_self_tfbstools_frequency) <- c("phage","tfbstools92_frequency")

stoperator_site_predictions93_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions93_self$stoperator_target))
names(stoperator_site_predictions93_self_tfbstools_frequency) <- c("phage","tfbstools93_frequency")

stoperator_site_predictions94_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions94_self$stoperator_target))
names(stoperator_site_predictions94_self_tfbstools_frequency) <- c("phage","tfbstools94_frequency")

stoperator_site_predictions95_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions95_self$stoperator_target))
names(stoperator_site_predictions95_self_tfbstools_frequency) <- c("phage","tfbstools95_frequency")

stoperator_site_predictions100_self_tfbstools_frequency <- as.data.frame(table(stoperator_site_predictions100_self$stoperator_target))
names(stoperator_site_predictions100_self_tfbstools_frequency) <- c("phage","tfbstools100_frequency")





#Combine TFBSTools data with MEME and Biostrings data
freq_table <- merge(freq_table,
                              stoperator_site_predictions80_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")


freq_table <- merge(freq_table,
                              stoperator_site_predictions85_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions86_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions87_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions88_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions89_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions90_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions91_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions92_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions93_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions94_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions95_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

freq_table <- merge(freq_table,
                              stoperator_site_predictions100_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")





# Check how well TFBSTools matches with Biostrings.

plot_scatter <- function(table,field1,field2,x_range,y_range,title){
  
  plot(table$field1,
       table$field2,
       xlim = x_range,ylim = y_range,main= title)
  abline(0,1)  
  
  correlation <- lm(field1 ~
                      field2,
                    data = table)
  summary(correlation)
  
  print(paste("Total number of sites in Field 1: ",sum(table$field1)))
  print(paste("Total number of sites in Field 2: ",sum(table$field2)))
}


x_coords <- c(0,100)
y_coords <- c(0,100)

plot_scatter("freq_table","biostrings_freq","tfbstools80_frequency",x_coords,y_coords,"80")
plot_scatter("freq_table","biostrings_freq","tfbstools85_frequency",x_coords,y_coords,"85")
plot_scatter("freq_table","biostrings_freq","tfbstools86_frequency",x_coords,y_coords,"86")
plot_scatter("freq_table","biostrings_freq","tfbstools87_frequency",x_coords,y_coords,"87")
plot_scatter("freq_table","biostrings_freq","tfbstools88_frequency",x_coords,y_coords,"88")
plot_scatter("freq_table","biostrings_freq","tfbstools89_frequency",x_coords,y_coords,"89")
plot_scatter("freq_table","biostrings_freq","tfbstools90_frequency",x_coords,y_coords,"90")
plot_scatter("freq_table","biostrings_freq","tfbstools91_frequency",x_coords,y_coords,"91")
plot_scatter("freq_table","biostrings_freq","tfbstools92_frequency",x_coords,y_coords,"92")
plot_scatter("freq_table","biostrings_freq","tfbstools93_frequency",x_coords,y_coords,"93")
plot_scatter("freq_table","biostrings_freq","tfbstools94_frequency",x_coords,y_coords,"94")
plot_scatter("freq_table","biostrings_freq","tfbstools95_frequency",x_coords,y_coords,"95")
plot_scatter("freq_table","biostrings_freq","tfbstools100_frequency",x_coords,y_coords,"100")



# QC:
# biostrings:              ; # sites = 9273
# 80 = Multiple R2 = 0.1943; # sites = 13193
# 85 = Multiple R2 = 0.8067; # sites = 9686
# 86 = Multiple R2 = 0.876; # sites = 9536
# 87 = Multiple R2 = 0.9224; # sites = 9399
# 88 = Multiple R2 = 0.9532; # sites = 9202
# 89 = Multiple R2 = 0.8966; # sites = 8965
# 90 = Multiple R2 = 0.7866; # sites = 8720
# 91 = Multiple R2 = 0.7791; # sites = 8535
# 92 = Multiple R2 = 0.7465; # sites = 8373
# 93 = Multiple R2 = 0.6517; # sites = 8124
# 94 = Multiple R2 = 0.5796; # sites = 7869
# 95 = Multiple R2 = 0.5526; # sites = 7493
# 100 = Multiple R2 = 3.081e-05; # sites = 2361
# Using TFBSTools cutoff of 88% identifies nearly the same number of
# stoperators in each genome (9202 compared to 9273), and it produces the most
# correlated data, with an R2 of 0.9532, So using 88% is probably the best 
# cutoff to use.






#TODO delete
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools80_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="80")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools85_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="85")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools86_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="86")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools87_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="87")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools88_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="88")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools89_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="89")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools90_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="90")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools91_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="91")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools92_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="92")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools93_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="93")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools94_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="94")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools95_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="95")
# abline(0,1)  
# plot(freq_table$biostrings_freq,
#      freq_table$tfbstools100_frequency,
#      xlim=c(0,100),ylim=c(0,100),main="100")
# abline(0,1)  
# 
# 
# 
# lm_biostrings_tfbs80_cor <- lm(biostrings_freq ~
#                                  tfbstools80_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs80_cor)
# #Stop264 = Multiple R2 = 0.1551
# #Stop327 = Multiple R2 = 0.1943
# 
# lm_biostrings_tfbs85_cor <- lm(biostrings_freq ~
#                                  tfbstools85_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs85_cor)
# #Stop264 = Multiple R2 = 0.7822
# #Stop327 = Multiple R2 = 0.8067
# 
# lm_biostrings_tfbs86_cor <- lm(biostrings_freq ~
#                                  tfbstools86_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs86_cor)
# #Stop264 = Multiple R2 = 0.8615
# #Stop327 = Multiple R2 = 0.876
# 
# lm_biostrings_tfbs87_cor <- lm(biostrings_freq ~
#                                  tfbstools87_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs87_cor)
# #Stop264 = Multiple R2 = 0.9162
# #Stop327 = Multiple R2 = 0.9224
# 
# lm_biostrings_tfbs88_cor <- lm(biostrings_freq ~
#                                  tfbstools88_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs88_cor)
# #Stop264 = Multiple R2 = 0.954
# #Stop327 = Multiple R2 = 0.9532
# 
# lm_biostrings_tfbs89_cor <- lm(biostrings_freq ~
#                                  tfbstools89_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs89_cor)
# #Stop264 = Multiple R2 = 0.9054
# #Stop327 = Multiple R2 = 0.8966
# 
# 
# lm_biostrings_tfbs90_cor <- lm(biostrings_freq ~
#                                  tfbstools90_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs90_cor)
# #Stop264 = Multiple R2 = 0.8098
# #Stop327 = Multiple R2 = 0.7866
# 
# lm_biostrings_tfbs91_cor <- lm(biostrings_freq ~
#                                  tfbstools91_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs91_cor)
# #Stop264 = Multiple R2 = 0.8005
# #Stop327 = Multiple R2 = 0.7791
# 
# lm_biostrings_tfbs92_cor <- lm(biostrings_freq ~
#                                  tfbstools92_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs92_cor)
# #Stop264 = Multiple R2 = 0.7619
# #Stop327 = Multiple R2 = 0.7465
# 
# lm_biostrings_tfbs93_cor <- lm(biostrings_freq ~
#                                  tfbstools93_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs93_cor)
# #Stop264 = Multiple R2 = 0.6733
# #Stop327 = Multiple R2 = 0.6517
# 
# lm_biostrings_tfbs94_cor <- lm(biostrings_freq ~
#                                  tfbstools94_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs94_cor)
# #Stop264 = Multiple R2 = 0.6078
# #Stop327 = Multiple R2 = 0.5796
# 
# 
# 
# lm_biostrings_tfbs95_cor <- lm(biostrings_freq ~
#                                  tfbstools95_frequency,
#                                data = freq_table)
# summary(lm_biostrings_tfbs95_cor)
# #Stop264 = Multiple R2 = 0.5812
# #Stop327 = Multiple R2 = 0.5526
# 
# 
# lm_biostrings_tfbs100_cor <- lm(biostrings_freq ~
#                                   tfbstools100_frequency,
#                                 data = freq_table)
# summary(lm_biostrings_tfbs100_cor)
# #Stop264 = Multiple R2 = 2.752e-06
# #Stop327 = Multiple R2 = 3.081e-05
# 
# 
# 
# sum(freq_table$biostrings_freq) #Stop264 = 7500 sites; Stop327 = 9273 sites
# sum(freq_table$tfbstools80_frequency) #Stop264 = 10705 sites; Stop327 = 13193 sites
# sum(freq_table$tfbstools85_frequency) #Stop264 = 7835 sites; Stop327 = 9686 sites
# sum(freq_table$tfbstools86_frequency) #Stop264 = 7715 sites; Stop327 = 9536 sites
# sum(freq_table$tfbstools87_frequency) #Stop264 = 7602 sites; Stop327 = 9399 sites
# sum(freq_table$tfbstools88_frequency) #Stop264 = 7445 sites; Stop327 = 9202 sites
# sum(freq_table$tfbstools89_frequency) #Stop264 = 7249 sites; Stop327 = 8965 sites
# sum(freq_table$tfbstools90_frequency) #Stop264 = 7055 sites; Stop327 = 8720 sites
# sum(freq_table$tfbstools91_frequency) #Stop264 = 6903 sites; Stop327 = 8535 sites
# sum(freq_table$tfbstools92_frequency) #Stop264 = 6768 sites; Stop327 = 8373 sites
# sum(freq_table$tfbstools93_frequency) #Stop264 = 6569 sites; Stop327 = 8124 sites
# sum(freq_table$tfbstools94_frequency) #Stop264 = 6367 sites; Stop327 = 7869 sites
# sum(freq_table$tfbstools95_frequency) #Stop264 = 6063 sites; Stop327 = 7493 sites
# sum(freq_table$tfbstools100_frequency) #Stop264 = 1856 sites; Stop327 = 2361 sites
# 













names(stoperator_site_predictions88_self) <- 
  paste("tfbs88_",names(stoperator_site_predictions88_self),sep="")


#TODO
# QC
# temp_tfbs_l5 <- subset(stoperator_site_predictions88_self,stoperator_site_predictions88_self$stoperator_target == "l5")
# temp_biostrings_l5 <- subset(bios_sites,bios_sites$phage == "l5")







# Compute the exact position of sites present in biostrings and 
# tfbstools88 datasets.


stoperator_biostrings_tfbs88 <- merge(bios_sites,
                                      stoperator_site_predictions88_self,
                                      by.x="stop_site_id",
                                      by.y="tfbs88_stop_site_id",
                                      all.x=TRUE,
                                      all.y=TRUE)


stoperator_biostrings_tfbs88$source_biostrings <- 
  ifelse(is.na(stoperator_biostrings_tfbs88$phage),
         FALSE,TRUE)

stoperator_biostrings_tfbs88$source_tfbs88 <- 
  ifelse(is.na(stoperator_biostrings_tfbs88$tfbs88_stoperator_target),
         FALSE,TRUE)

stoperator_biostrings_tfbs88$source_both <- 
  ifelse(stoperator_biostrings_tfbs88$source_biostrings & 
           stoperator_biostrings_tfbs88$source_tfbs88,
         TRUE,FALSE)

summary(stoperator_biostrings_tfbs88$source_biostrings)
summary(stoperator_biostrings_tfbs88$source_tfbs88)
summary(stoperator_biostrings_tfbs88$source_both)


# QC:
# There are 9123 sites in both Biostrings and TFBS88 (98%).
# There are 150 sites in Biostrings that are not present in TFBS88. 
# There are 79 sites in TFBS88 that are not present in Biostrings.



#TODO
#QC
#Stop264:
#temp_tfbs_sites <- subset(stoperator_biostrings_tfbs88,stoperator_biostrings_tfbs88$source_tfbs88 == TRUE & stoperator_biostrings_tfbs88$source_biostrings == FALSE)
#temp_biostrings_sites <- subset(stoperator_biostrings_tfbs88,stoperator_biostrings_tfbs88$source_tfbs88 == FALSE & stoperator_biostrings_tfbs88$source_biostrings == TRUE)
#temp_tfbs88_l5_sites <- subset(stoperator_site_predictions88_self,stoperator_site_predictions88_self$tfbs88_stoperator_target == 'l5')
#Of the 30 L5 site originally reported in Brown et al. 1997, the TFBS88 dataset does not include Site #3 (EMSA bound) or Site #30 (EMSA unbound)







# Now that the best TFBS cutoff has been determined (88%), search all genomes
# with each PWM and determine how many sites of each motif are in each genome.
# Re-factor columns in the self-site dataset.
stoperator_site_predictions88_self$tfbs88_stoperator_target <- 
  factor(stoperator_site_predictions88_self$tfbs88_stoperator_target)

stoperator_site_predictions88_self$tfbs88_stoperator_motif <- 
  factor(stoperator_site_predictions88_self$tfbs88_stoperator_motif)

stoperator_site_predictions88_self$tfbs88_motif_target <- 
  factor(stoperator_site_predictions88_self$tfbs88_motif_target)

stoperator_site_predictions88_self$tfbs88_site_strand2 <- 
  factor(stoperator_site_predictions88_self$tfbs88_site_strand2)



# Subset all sites predicted in each genome from each PWM above 88% threshold.
stoperator_site_predictions80$site_self <- 
  ifelse(is.element(stoperator_site_predictions80$stop_site_id,
                    stoperator_site_predictions88_self$tfbs88_stop_site_id),
         TRUE,FALSE)

stoperator_site_predictions88 <- 
  subset(stoperator_site_predictions80,
         stoperator_site_predictions80$site_rel_score >= 0.88)



stoperator_site_predictions88$stoperator_target <- 
  factor(stoperator_site_predictions88$stoperator_target)

stoperator_site_predictions88$stoperator_motif <- 
  factor(stoperator_site_predictions88$stoperator_motif)

stoperator_site_predictions88$motif_target <- 
  factor(stoperator_site_predictions88$motif_target)

stoperator_site_predictions88$site_strand2 <- 
  factor(stoperator_site_predictions88$site_strand2)

names(stoperator_site_predictions88) <- 
  paste("tfbs88_",names(stoperator_site_predictions88),sep="")



stoperator_site_tfbs88_reduced <- 
  subset(stoperator_site_predictions88,select=c("tfbs88_stop_site_id",
                                                "tfbs88_motif_target",
                                                "tfbs88_stoperator_target",
                                                "tfbs88_stoperator_motif",
                                                "tfbs88_site_start",
                                                "tfbs88_site_end",
                                                "tfbs88_site_strand2",
                                                "tfbs88_site_seq",
                                                "tfbs88_site_abs_score",
                                                "tfbs88_site_rel_score",
                                                "tfbs88_site_self"))


stoperator_site_tfbs88_reduced$tfbs88_site_abs_score <- 
  round(stoperator_site_tfbs88_reduced$tfbs88_site_abs_score,2)

stoperator_site_tfbs88_reduced$tfbs88_site_rel_score <- 
  round(stoperator_site_tfbs88_reduced$tfbs88_site_rel_score,2)



# Export data for analysis with immunity data.
write.table(stoperator_site_tfbs88_reduced,
            paste(DIR_OUTPUT,"stoperator_site_tfbs88_reduced.csv",sep=""),
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)




###
