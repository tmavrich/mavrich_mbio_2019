#To compute position weight matrices for stoperators and compare similarities between PWMs

###Load packages

#Ensure packages are installed
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")
# biocLite("TFBSTools")


library("Biostrings")
library("TFBSTools")






###Import datasets


#Import all Cluster A genomes from actino1321 database
#Read individual fasta files by generating a list of files in the directory,
#and passing that list to the readDNAStringSet function
setwd("~/scratch/stoperator_distance/input_genomes/")
actino1321_genomes <- readDNAStringSet(list.files())





#Import list of MEME-identified stoperators
setwd("~/scratch/stoperator_distance/input/")

#stoperator data columns:
# "phage" (imported as factor)
# "stoperator_id" (imported as factor)
# "sequence" (imported as factor)

stoperator_data <- read.csv("stoperator_data.csv",sep=",",header=TRUE)



#Tally # of stoperators identified by MEME as input
stoperator_MEME_frequency <- as.data.frame(table(stoperator_data$phage))
names(stoperator_MEME_frequency) <- c("phage","meme_frequency")






#phage metadata columns:
# "phageid" (imported as factor)
# "size" (imported as int)
# "genome_center_alignment_reference" (imported as int)

phage_metadata <- read.csv("phage_metadata.csv",sep=",",header=TRUE)
















#Confirm presence of sites in each genomes in case genome sequence has changed since original MEME search

#For each phage, 
#1. create list of unique 13bp stoperators
#2. search the corresponding parent genome, both strands
#3. save list of sites and coordinates

#Create an empty dataframe to store all the search data
stoperator_biostrings_df <- data.frame(phage="empty",
                                       stoperator_forward_seq="empty",
                                       stoperator_reverse_seq="empty",
                                       start=0,
                                       end=0,
                                       strand="empty")




for (phageid in levels(stoperator_data$phage)){
  
  
  reduced_table_seqs <- subset(stoperator_data,
                               stoperator_data$phage == phageid,
                               select = c("sequence"))
  
  reduced_table_seqs$sequence <- factor(reduced_table_seqs$sequence)
  
  
  for (stoperator in levels(reduced_table_seqs$sequence)){
    
    
    stop_forward <- DNAString(stoperator)
    stop_reverse <- reverseComplement(stop_forward)
    genome_to_search <- actino1321_genomes[[phageid]]

    
    hits_plus <- matchPattern(stop_forward,genome_to_search)
    hits_minus <- matchPattern(stop_reverse,genome_to_search)
    
    
    if (length(hits_plus) > 0){
      
      hits_plus_df <- data.frame(phage=phageid,
                                 stoperator_forward_seq=as.character(stop_forward),
                                 stoperator_reverse_seq=as.character(stop_reverse),
                                 start=start(hits_plus),
                                 end=end(hits_plus),
                                 strand="forward")
      stoperator_biostrings_df <- rbind(stoperator_biostrings_df,hits_plus_df)
    }
    
    
    
    if (length(hits_minus) > 0){
      
      hits_minus_df <- data.frame(phage=phageid,
                                  stoperator_forward_seq=as.character(stop_forward),
                                  stoperator_reverse_seq=as.character(stop_reverse),
                                  start=start(hits_minus),
                                  end=end(hits_minus),
                                  strand="reverse")
      stoperator_biostrings_df <- rbind(stoperator_biostrings_df,hits_minus_df)
    }
    
    
    
    
    
  }
  
}


#Remove the row of empty data used to initiate the table
stoperator_biostrings_df <- stoperator_biostrings_df[stoperator_biostrings_df$phage != "empty",]
stoperator_biostrings_df$phage <- factor(stoperator_biostrings_df$phage)
stoperator_biostrings_df$strand <- factor(stoperator_biostrings_df$strand)
stoperator_biostrings_df$stoperator_forward_seq <- factor(stoperator_biostrings_df$stoperator_forward_seq)
stoperator_biostrings_df$stoperator_reverse_seq <- factor(stoperator_biostrings_df$stoperator_reverse_seq)



#Create unique identifier for each site
stoperator_biostrings_df$stop_site_id <- paste(stoperator_biostrings_df$phage,
                                               stoperator_biostrings_df$strand,
                                               stoperator_biostrings_df$start,
                                               sep="_")

stoperator_biostrings_df$stop_site_id <- factor(stoperator_biostrings_df$stop_site_id)


#add metadata
stoperator_biostrings_df <- merge(stoperator_biostrings_df,phage_metadata,by.x="phage",by.y="phageid")

stoperator_biostrings_df$site_dist_from_center <- stoperator_biostrings_df$start - stoperator_biostrings_df$genome_center_alignment_reference


                                               
#Tally # of stoperators confirmed
stoperator_biostrings_frequency <- as.data.frame(table(stoperator_biostrings_df$phage))
names(stoperator_biostrings_frequency) <- c("phage","biostrings_frequency")


#Tally # of stoperators on each side of genome center
stoperator_biostrings_left_sites <- subset(stoperator_biostrings_df,stoperator_biostrings_df$site_dist_from_center <= 0)
stoperator_biostrings_right_sites <- subset(stoperator_biostrings_df,stoperator_biostrings_df$site_dist_from_center > 0)



stoperator_biostrings_left_sites_frequency <- as.data.frame(table(stoperator_biostrings_left_sites$phage))
names(stoperator_biostrings_left_sites_frequency) <- c("phage","biostrings_left_sites_frequency")

stoperator_biostrings_right_sites_frequency <- as.data.frame(table(stoperator_biostrings_right_sites$phage))
names(stoperator_biostrings_right_sites_frequency) <- c("phage","biostrings_right_sites_frequency")

#QC: should equal 0
nrow(stoperator_biostrings_df) - nrow(stoperator_biostrings_left_sites) - nrow(stoperator_biostrings_right_sites)








#Tally # of stoperators on each strand of each side of genome center
stoperator_biostrings_left_sites_forward <- subset(stoperator_biostrings_left_sites,stoperator_biostrings_left_sites$strand == "forward")
stoperator_biostrings_left_sites_reverse <- subset(stoperator_biostrings_left_sites,stoperator_biostrings_left_sites$strand == "reverse")


stoperator_biostrings_right_sites_forward <- subset(stoperator_biostrings_right_sites,stoperator_biostrings_right_sites$strand == "forward")
stoperator_biostrings_right_sites_reverse <- subset(stoperator_biostrings_right_sites,stoperator_biostrings_right_sites$strand == "reverse")



stoperator_biostrings_left_sites_forward_frequency <- as.data.frame(table(stoperator_biostrings_left_sites_forward$phage))
names(stoperator_biostrings_left_sites_forward_frequency) <- c("phage","biostrings_left_sites_forward_frequency")

stoperator_biostrings_left_sites_reverse_frequency <- as.data.frame(table(stoperator_biostrings_left_sites_reverse$phage))
names(stoperator_biostrings_left_sites_reverse_frequency) <- c("phage","biostrings_left_sites_reverse_frequency")


stoperator_biostrings_right_sites_forward_frequency <- as.data.frame(table(stoperator_biostrings_right_sites_forward$phage))
names(stoperator_biostrings_right_sites_forward_frequency) <- c("phage","biostrings_right_sites_forward_frequency")

stoperator_biostrings_right_sites_reverse_frequency <- as.data.frame(table(stoperator_biostrings_right_sites_reverse$phage))
names(stoperator_biostrings_right_sites_reverse_frequency) <- c("phage","biostrings_right_sites_reverse_frequency")











#Combine MEME and Biostrings data
stoperator_frequency <- merge(stoperator_MEME_frequency,
                              stoperator_biostrings_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_biostrings_left_sites_frequency,
                              by.x="phage",
                              by.y="phage")
stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_biostrings_right_sites_frequency,
                              by.x="phage",
                              by.y="phage")
stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_biostrings_left_sites_forward_frequency,
                              by.x="phage",
                              by.y="phage")
stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_biostrings_left_sites_reverse_frequency,
                              by.x="phage",
                              by.y="phage")
stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_biostrings_right_sites_forward_frequency,
                              by.x="phage",
                              by.y="phage")
stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_biostrings_right_sites_reverse_frequency,
                              by.x="phage",
                              by.y="phage")





stoperator_frequency$meme_biostrings_diff <- stoperator_frequency$meme_frequency - stoperator_frequency$biostrings_frequency
stoperator_frequency$biostrings_left_sites_percent <- stoperator_frequency$biostrings_left_sites_frequency/stoperator_frequency$biostrings_frequency
stoperator_frequency$biostrings_right_sites_percent <- stoperator_frequency$biostrings_right_sites_frequency/stoperator_frequency$biostrings_frequency


stoperator_frequency$biostrings_left_sites_forward_percent <- stoperator_frequency$biostrings_left_sites_forward_frequency/stoperator_frequency$biostrings_left_sites_frequency
stoperator_frequency$biostrings_left_sites_reverse_percent <- stoperator_frequency$biostrings_left_sites_reverse_frequency/stoperator_frequency$biostrings_left_sites_frequency

stoperator_frequency$biostrings_right_sites_forward_percent <- stoperator_frequency$biostrings_right_sites_forward_frequency/stoperator_frequency$biostrings_right_sites_frequency
stoperator_frequency$biostrings_right_sites_reverse_percent <- stoperator_frequency$biostrings_right_sites_reverse_frequency/stoperator_frequency$biostrings_right_sites_frequency


#QC
#Site tally checks should equal 0
stoperator_frequency$site_tally_check <- stoperator_frequency$biostrings_frequency - stoperator_frequency$biostrings_left_sites_frequency - stoperator_frequency$biostrings_right_sites_frequency
summary(stoperator_frequency$site_tally_check)

stoperator_frequency$site_left_tally_check <- stoperator_frequency$biostrings_left_sites_frequency - stoperator_frequency$biostrings_left_sites_forward_frequency - stoperator_frequency$biostrings_left_sites_reverse_frequency
summary(stoperator_frequency$site_left_tally_check)

stoperator_frequency$site_right_tally_check <- stoperator_frequency$biostrings_right_sites_frequency - stoperator_frequency$biostrings_right_sites_forward_frequency - stoperator_frequency$biostrings_right_sites_reverse_frequency
summary(stoperator_frequency$site_right_tally_check)



#Site percent checks should equal 1
stoperator_frequency$site_percent_check <- stoperator_frequency$biostrings_left_sites_percent + stoperator_frequency$biostrings_right_sites_percent
summary(stoperator_frequency$site_percent_check)


stoperator_frequency$site_left_percent_check <- stoperator_frequency$biostrings_left_sites_forward_percent + stoperator_frequency$biostrings_left_sites_reverse_percent
summary(stoperator_frequency$site_left_percent_check)

stoperator_frequency$site_right_percent_check <- stoperator_frequency$biostrings_right_sites_forward_percent + stoperator_frequency$biostrings_right_sites_reverse_percent
summary(stoperator_frequency$site_right_percent_check)



#QC: if this does not equal 0, then there is a difference between the number of stoperators present in the two sets of stoperators.
nrow(stoperator_data) - nrow(stoperator_biostrings_df)

hist(stoperator_frequency$meme_biostrings_diff)
setdiff(levels(stoperator_MEME_frequency$phage),levels(stoperator_biostrings_frequency$phage))
setdiff(levels(stoperator_biostrings_frequency$phage),levels(stoperator_MEME_frequency$phage))


#Stoperator327 data: MEME sites tally = 9266. Biostrings sites tally = 9273.
#biostrings contains 7 more stoperators. 
#Seven phages contain 1 fewer stoperator from MEME than Biostrings
#BeesKnees (A1)
#Changeling (A2)
#Doom (A1)
#KBG (A1)
#Petruchio (A1)
#Saintus (A8)
#SarFire (A1)
#Overall,  only a few genomes contain small differences in stoperators than expected

#Compare # MEME sites versus # Biostrings sites per genome
plot(stoperator_frequency$meme_frequency,stoperator_frequency$biostrings_frequency)



#Plot distribution of number of stoperators
setwd("~/scratch/stoperator_distance/output/")
par(mar=c(4,8,8,4))
hist(table(stoperator_biostrings_df$phage),col="black",xlim=c(0,50),ylim=c(0,80),breaks=10,cex.axis=2,ann=FALSE,main=NULL,las=1)
dev.copy(pdf,'stoperator327_number_stops_per_genome.pdf')
dev.off()







#Plot distribution of percent of stoperators on each side of genome
hist(stoperator_frequency$biostrings_left_sites_forward_percent,xlim=c(0,1),col="black",breaks=100)
hist(stoperator_frequency$biostrings_left_sites_reverse_percent,xlim=c(0,1),col="black",breaks=100)


hist(stoperator_frequency$biostrings_right_sites_forward_percent,xlim=c(0,1),col="black",breaks=100)
hist(stoperator_frequency$biostrings_right_sites_reverse_percent,xlim=c(0,1),col="black",breaks=100)

par(mar=c(4,8,8,4))
plot(stoperator_frequency$biostrings_right_sites_reverse_percent,
     stoperator_frequency$biostrings_left_sites_forward_percent,
     xlim=c(0,1),ylim=c(0,1),
     cex.axis=2,ann=FALSE,main=NULL,las=1,
     col="black",pch=16,cex=2)
dev.copy(pdf,'stoperator327_percent_txn_oriented_sites_per_genome.pdf')
dev.off()



#Output biostrings data for reference
write.table(stoperator_biostrings_df,
            "stoperator_biostrings_coordinates.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)

























###Create PWMs
#Use the sites that were identified in each genome by Biostrings search
#For each phage in the table, create a pwm
#Create two lists of PWMs.

#Log2ProbRatio PWM
stoperator_pwm_list <- vector("list",nlevels(stoperator_biostrings_df$phage))


#Prob PWM
stoperator_pwm2_list <- vector("list",nlevels(stoperator_biostrings_df$phage))


#For each phage genome, the list of stoperators must contain at least one of each
#of the 4 nucleotides A,C,T,G. Otherwise, when PFMatrix is called, it will throw an error.
#The for loop will probably crash if one or more PWMs does not meet this criteria,
#but a QC step has been used to confirm it.
stoperator_letters_check_sum <- 0



count <- 1
for (stoperator_phage in levels(stoperator_biostrings_df$phage)){
  
  reduced_table_seqs <- subset(stoperator_biostrings_df,
                               stoperator_biostrings_df$phage == stoperator_phage,
                               select = c("stoperator_forward_seq"))
  
  reduced_table_seqs <- as.character(reduced_table_seqs$stoperator_forward_seq)
  
  consensus_matrix <- consensusMatrix(reduced_table_seqs)
  stoperator_pfm <- PFMatrix(ID=stoperator_phage,
                                  name=stoperator_phage,
                                  profileMatrix = consensus_matrix)
  
  stoperator_pwm <- toPWM(stoperator_pfm,type="log2probratio") 
  stoperator_pwm2 <- toPWM(stoperator_pfm,type="prob") 
  
  
  
  #Append the pwm to the list of pwms
  stoperator_pwm_list[[count]] <- stoperator_pwm
  stoperator_pwm2_list[[count]] <- stoperator_pwm2  
  
  
  #Append DNA letter frequency check
  reduced_table_seqs_dnastringset <- DNAStringSet(reduced_table_seqs)
  all_4_nucleotides <- ifelse(length(uniqueLetters(reduced_table_seqs_dnastringset)) == 4,1,0)
  
  stoperator_letters_check_sum <- stoperator_letters_check_sum + all_4_nucleotides
  
  
  count <- count + 1
  
}  


#QC: should equal 0
nlevels(stoperator_biostrings_df$phage) - stoperator_letters_check_sum














###Compute Log2ProbRatio PWM distances
#This log2-based score seems to be commonly used instead of the straight up probability.
stoperator_pwm_list1 <- stoperator_pwm_list
stoperator_pwm_list2 <- stoperator_pwm_list


pwm_output_list <- vector("list",length(stoperator_pwm_list1)*length(stoperator_pwm_list2))

count <- 1
for (pwm1 in stoperator_pwm_list1){
  for (pwm2 in stoperator_pwm_list2){
    
    dist_pearson <- 1 - PWMSimilarity(pwm1,
                                      pwm2,
                                      method="Pearson")
    
    dist_euc <- PWMSimilarity(pwm1,
                              pwm2,
                              method="Euclidean")
    
    
    
    
    #print(c(ID(pwm1),ID(pwm2),dist_pearson,dist_euc))

    #Append the distance data to the output list
    pwm_output_list[[count]] <- c(ID(pwm1),ID(pwm2),round(dist_pearson,2),round(dist_euc,2))
    count <- count + 1

  }
}



pwm_output_df <- as.data.frame(t(as.data.frame(pwm_output_list)))
rownames(pwm_output_df) <- NULL
names(pwm_output_df) <- c("phage1","phage2","dist_pearson","dist_euc")


#Output data
setwd("~/scratch/stoperator_distance/output/")
write.table(pwm_output_df,
            "stoperator_pwm_distances.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)








#QC
#Compare the euclidean distance of PWM probabilities and the log2probratio probabilities

pwm_output_df$phage1_phage2 <- paste(pwm_output_df$phage1,
                                     "_",
                                     pwm_output_df$phage2,
                                     sep="")

pwm_output_df$dist_pearson <- as.numeric(as.character(pwm_output_df$dist_pearson))
pwm_output_df$dist_euc <- as.numeric(as.character(pwm_output_df$dist_euc))









###Compute Prob PWM distances
stoperator_pwm2_list1 <- stoperator_pwm2_list
stoperator_pwm2_list2 <- stoperator_pwm2_list


pwm2_output_list <- vector("list",length(stoperator_pwm2_list1)*length(stoperator_pwm2_list2))

count <- 1
for (pwm1 in stoperator_pwm2_list1){
  for (pwm2 in stoperator_pwm2_list2){
    
    dist_pearson <- 1 - PWMSimilarity(pwm1,
                                      pwm2,
                                      method="Pearson")
    
    dist_euc <- PWMSimilarity(pwm1,
                              pwm2,
                              method="Euclidean")
    

    #Append the distance data to the output list
    pwm2_output_list[[count]] <- c(ID(pwm1),ID(pwm2),round(dist_pearson,2),round(dist_euc,2))
    count <- count + 1
    
  }
}



pwm2_output_df <- as.data.frame(t(as.data.frame(pwm2_output_list)))
rownames(pwm2_output_df) <- NULL
names(pwm2_output_df) <- c("prob_phage1","prob_phage2","prob_dist_pearson","prob_dist_euc")
pwm2_output_df$phage1_phage2 <- paste(pwm2_output_df$prob_phage1,
                                      "_",
                                      pwm2_output_df$prob_phage2,
                                      sep="")

pwm2_output_df$prob_dist_pearson <- as.numeric(as.character(pwm2_output_df$prob_dist_pearson))
pwm2_output_df$prob_dist_euc <- as.numeric(as.character(pwm2_output_df$prob_dist_euc))






#Merge both PWM lists

combined_pwm_output_df <- merge(pwm_output_df,
                                pwm2_output_df,
                                by.x="phage1_phage2",
                                by.y="phage1_phage2")


plot(combined_pwm_output_df$dist_euc,combined_pwm_output_df$prob_dist_euc)

plot(combined_pwm_output_df$dist_euc,combined_pwm_output_df$dist_pearson)

plot(combined_pwm_output_df$prob_dist_euc,combined_pwm_output_df$prob_dist_pearson)
abline(0,1)
























###After creating all stoperator PWMs, iterate through all genomes for each PWM and predict sites

#Convert object type of stoperator PWMs from standard list to PWMatrixList
stoperator_pwm_matrixlist <- do.call(PWMatrixList,stoperator_pwm_list)


stop_site_set80_list <- searchSeq(stoperator_pwm_matrixlist,actino1321_genomes)
stop_site_set80_df <- as(stop_site_set80_list,"data.frame")
#Note: it takes about an hour to run the two commands above using stop327 data




#dataframe output: each row is a predicted site in a genome from one PWM
#Column names in dataframe
# "seqnames" = factor, 327 levels are individual phage genomes
# "source" = factor, 1 level ("TFBS")
# "feature" = factor, 1 level ("TFBS")
# "start" = int
# "end" = int
# "absScore" = number (from ~4 to ~25)
# "relScore" = number (from 0.8 to 1). Seems like it reflects the range of possible scores that could be
#obtained with the motif and is normalized between 0 and 1.
# "strand" = factor, 2 levels ("+","-")
# "ID" = factor, 327 levels are stoperator PWM names. Identical to "TF"
# "TF" = factor, 327 levels are stoperator PWM names. Identical to "ID"
# "class" = factor, 1 level ("Unknown")
# "siteSeqs" = chr, 13bp stoperator DNA sequence





#QC analysis
#stop_site_set95_list <- searchSeq(stoperator_pwm_matrixlist,actino1321_genomes,min.score = "95%")
#stop_site_set95_df <- as(stop_site_set95_list,"data.frame")


#actino1319_analysis:
#For 80%: 1797889 rows. 
#For 95%: 350763 rows.
#






#reduce data for analysis with immunity data
# "seqnames" = factor, 327 levels are individual phage genomes
# "start" = int
# "end" = int
# "strand" = factor, 2 levels ("+","-")
# "siteSeqs" = chr, 13bp stoperator DNA sequence
# "absScore" = number (from ~4 to ~25)
# "relScore" = number (from 0.8 to 1). Seems like it reflects the range of possible scores that could be
# "ID" = factor, 327 levels are stoperator PWM names. Identical to "TF"
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
  
  

stoperator_site_predictions80$motif_target <- paste(stoperator_site_predictions80$stoperator_motif,
                                                    "_",
                                                    stoperator_site_predictions80$stoperator_target,
                                                    sep="")

stoperator_site_predictions80$motif_target <- as.factor(stoperator_site_predictions80$motif_target)

#actino1319 = 248 PWMs * 325 target genomes = 80,600 combinations. This list only contains 80,165 levels.
#This likely reflects the fact that target genomes may not contain any predicted sites.

#actino1321 = 264 PWMs * 327 target genomes = 86,328 combinations. 
#This list only has 85,892 levels though.

#actino1321, stop327 data = 327 PWMs * 327 target genomes = 106929 combinations.
#This list has 106300 levels.

all_levels <- expand.grid(stoperator_target = as.character(levels(stoperator_site_predictions80$stoperator_target)),
                          stoperator_motif = as.character(levels(stoperator_site_predictions80$stoperator_motif)))

all_levels$motif_target <- paste(all_levels$stoperator_motif,
                                 "_",
                                 all_levels$stoperator_target,
                                 sep="")

all_levels$motif_target <- as.factor(all_levels$motif_target)


#Re-factor motif_target using all possible pairwise combinations.
#This will enable quantification for all combinations, even those that are missing.
stoperator_site_predictions80$motif_target <- factor(stoperator_site_predictions80$motif_target,levels(all_levels$motif_target))
#Now there are 106929 levels.




#Create unique site ids
stoperator_site_predictions80$site_strand2 <- ifelse(stoperator_site_predictions80$site_strand == "+",
                                                          "forward",
                                                          "reverse")
stoperator_site_predictions80$site_strand2 <- factor(stoperator_site_predictions80$site_strand2)

stoperator_site_predictions80$stop_site_id <- paste(stoperator_site_predictions80$stoperator_target,
                                                    stoperator_site_predictions80$site_strand2,
                                                    stoperator_site_predictions80$site_start,
                                                         sep="_")

stoperator_site_predictions80$stop_site_id <- factor(stoperator_site_predictions80$stop_site_id)



#QC
#temp_l5 <- subset(stoperator_site_predictions80,stoperator_site_predictions80$stoperator_target == "l5" & stoperator_site_predictions80$stoperator_motif == "l5")






#Table creates a two column table = the levels and the frequency (int)
stoperator_site_predictions80_self <- subset(stoperator_site_predictions80,
                                             as.character(stoperator_site_predictions80$stoperator_target) == as.character(stoperator_site_predictions80$stoperator_motif))


stoperator_site_predictions80_self$stoperator_target <- factor(stoperator_site_predictions80_self$stoperator_target)
stoperator_site_predictions80_self$stoperator_motif <- factor(stoperator_site_predictions80_self$stoperator_motif)
stoperator_site_predictions80_self$motif_target <- factor(stoperator_site_predictions80_self$motif_target)
stoperator_site_predictions80_self$site_strand2 <- factor(stoperator_site_predictions80_self$site_strand2)

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
stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions80_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")


stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions85_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions86_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions87_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions88_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions89_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions90_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions91_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions92_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions93_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions94_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions95_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")

stoperator_frequency <- merge(stoperator_frequency,
                              stoperator_site_predictions100_self_tfbstools_frequency,
                              by.x="phage",
                              by.y="phage")





#Check how well TFBSTools matches with Biostrings



plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools80_frequency,xlim=c(0,100),ylim=c(0,100),main="80")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools85_frequency,xlim=c(0,100),ylim=c(0,100),main="85")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools86_frequency,xlim=c(0,100),ylim=c(0,100),main="86")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools87_frequency,xlim=c(0,100),ylim=c(0,100),main="87")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools88_frequency,xlim=c(0,100),ylim=c(0,100),main="88")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools89_frequency,xlim=c(0,100),ylim=c(0,100),main="89")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools90_frequency,xlim=c(0,100),ylim=c(0,100),main="90")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools91_frequency,xlim=c(0,100),ylim=c(0,100),main="91")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools92_frequency,xlim=c(0,100),ylim=c(0,100),main="92")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools93_frequency,xlim=c(0,100),ylim=c(0,100),main="93")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools94_frequency,xlim=c(0,100),ylim=c(0,100),main="94")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools95_frequency,xlim=c(0,100),ylim=c(0,100),main="95")
abline(0,1)  
plot(stoperator_frequency$biostrings_frequency,stoperator_frequency$tfbstools100_frequency,xlim=c(0,100),ylim=c(0,100),main="100")
abline(0,1)  



lm_biostrings_tfbs80_cor <- lm(biostrings_frequency ~
                                 tfbstools80_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs80_cor)
#Stop264 = Multiple R2 = 0.1551
#Stop327 = Multiple R2 = 0.1943

lm_biostrings_tfbs85_cor <- lm(biostrings_frequency ~
                                 tfbstools85_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs85_cor)
#Stop264 = Multiple R2 = 0.7822
#Stop327 = Multiple R2 = 0.8067

lm_biostrings_tfbs86_cor <- lm(biostrings_frequency ~
                                 tfbstools86_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs86_cor)
#Stop264 = Multiple R2 = 0.8615
#Stop327 = Multiple R2 = 0.876

lm_biostrings_tfbs87_cor <- lm(biostrings_frequency ~
                                 tfbstools87_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs87_cor)
#Stop264 = Multiple R2 = 0.9162
#Stop327 = Multiple R2 = 0.9224

lm_biostrings_tfbs88_cor <- lm(biostrings_frequency ~
                                 tfbstools88_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs88_cor)
#Stop264 = Multiple R2 = 0.954
#Stop327 = Multiple R2 = 0.9532

lm_biostrings_tfbs89_cor <- lm(biostrings_frequency ~
                                 tfbstools89_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs89_cor)
#Stop264 = Multiple R2 = 0.9054
#Stop327 = Multiple R2 = 0.8966


lm_biostrings_tfbs90_cor <- lm(biostrings_frequency ~
                                 tfbstools90_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs90_cor)
#Stop264 = Multiple R2 = 0.8098
#Stop327 = Multiple R2 = 0.7866

lm_biostrings_tfbs91_cor <- lm(biostrings_frequency ~
                                 tfbstools91_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs91_cor)
#Stop264 = Multiple R2 = 0.8005
#Stop327 = Multiple R2 = 0.7791

lm_biostrings_tfbs92_cor <- lm(biostrings_frequency ~
                                 tfbstools92_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs92_cor)
#Stop264 = Multiple R2 = 0.7619
#Stop327 = Multiple R2 = 0.7465

lm_biostrings_tfbs93_cor <- lm(biostrings_frequency ~
                                 tfbstools93_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs93_cor)
#Stop264 = Multiple R2 = 0.6733
#Stop327 = Multiple R2 = 0.6517

lm_biostrings_tfbs94_cor <- lm(biostrings_frequency ~
                                 tfbstools94_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs94_cor)
#Stop264 = Multiple R2 = 0.6078
#Stop327 = Multiple R2 = 0.5796



lm_biostrings_tfbs95_cor <- lm(biostrings_frequency ~
                                 tfbstools95_frequency,
                               data = stoperator_frequency)
summary(lm_biostrings_tfbs95_cor)
#Stop264 = Multiple R2 = 0.5812
#Stop327 = Multiple R2 = 0.5526


lm_biostrings_tfbs100_cor <- lm(biostrings_frequency ~
                                  tfbstools100_frequency,
                                data = stoperator_frequency)
summary(lm_biostrings_tfbs100_cor)
#Stop264 = Multiple R2 = 2.752e-06
#Stop327 = Multiple R2 = 3.081e-05



sum(stoperator_frequency$biostrings_frequency) #Stop264 = 7500 sites; Stop327 = 9273 sites
sum(stoperator_frequency$tfbstools80_frequency) #Stop264 = 10705 sites; Stop327 = 13193 sites
sum(stoperator_frequency$tfbstools85_frequency) #Stop264 = 7835 sites; Stop327 = 9686 sites
sum(stoperator_frequency$tfbstools86_frequency) #Stop264 = 7715 sites; Stop327 = 9536 sites
sum(stoperator_frequency$tfbstools87_frequency) #Stop264 = 7602 sites; Stop327 = 9399 sites
sum(stoperator_frequency$tfbstools88_frequency) #Stop264 = 7445 sites; Stop327 = 9202 sites
sum(stoperator_frequency$tfbstools89_frequency) #Stop264 = 7249 sites; Stop327 = 8965 sites
sum(stoperator_frequency$tfbstools90_frequency) #Stop264 = 7055 sites; Stop327 = 8720 sites
sum(stoperator_frequency$tfbstools91_frequency) #Stop264 = 6903 sites; Stop327 = 8535 sites
sum(stoperator_frequency$tfbstools92_frequency) #Stop264 = 6768 sites; Stop327 = 8373 sites
sum(stoperator_frequency$tfbstools93_frequency) #Stop264 = 6569 sites; Stop327 = 8124 sites
sum(stoperator_frequency$tfbstools94_frequency) #Stop264 = 6367 sites; Stop327 = 7869 sites
sum(stoperator_frequency$tfbstools95_frequency) #Stop264 = 6063 sites; Stop327 = 7493 sites
sum(stoperator_frequency$tfbstools100_frequency) #Stop264 = 1856 sites; Stop327 = 2361 sites







#QC Analysis: 
#Stop264: using TFBSTools cutoff of 88% identifies nearly the same number of stoperators in each genome (7,445 compared to 7,500)
#and it produces the most correlated data, with an R2 of 0.954. So using 88% is probably the best cutoff to use.

#Stop327: using TFBSTools cutoff of 88% identifies nearly the same number of stoperators in each genome (9202 compared to 9273)
#and it produces the most correlated data, with an R2 of 0.9532, So using 88% is probably the best cutoff to use.







names(stoperator_site_predictions88_self) <- paste("tfbs88_",names(stoperator_site_predictions88_self),sep="")

#QC
#temp_tfbs_l5 <- subset(stoperator_site_predictions88_self,stoperator_site_predictions88_self$stoperator_target == "l5")
#temp_biostrings_l5 <- subset(stoperator_biostrings_df,stoperator_biostrings_df$phage == "l5")







#Compute the exact position of sites present in biostrings and tfbstools88 datasets


stoperator_biostrings_tfbs88 <- merge(stoperator_biostrings_df,
                                      stoperator_site_predictions88_self,
                                      by.x="stop_site_id",
                                      by.y="tfbs88_stop_site_id",
                                      all.x=TRUE,
                                      all.y=TRUE)


stoperator_biostrings_tfbs88$source_biostrings <- ifelse(is.na(stoperator_biostrings_tfbs88$phage),FALSE,TRUE)
stoperator_biostrings_tfbs88$source_tfbs88 <- ifelse(is.na(stoperator_biostrings_tfbs88$tfbs88_stoperator_target),FALSE,TRUE)
stoperator_biostrings_tfbs88$source_both <- ifelse(stoperator_biostrings_tfbs88$source_biostrings & stoperator_biostrings_tfbs88$source_tfbs88,TRUE,FALSE)

summary(stoperator_biostrings_tfbs88$source_biostrings)
summary(stoperator_biostrings_tfbs88$source_tfbs88)
summary(stoperator_biostrings_tfbs88$source_both)


#Stop264:
#There are 7378 sites in both Biostrings and TFBS88 (98%).
#There are 122 sites in Biostrings that are not present in TFBS88. 
#There are 67 sites in TFBS88 that are not present in Biostrings.


#Stop327:
#There are 9123 sites in both Biostrings and TFBS88 (98%).
#There are 150 sites in Biostrings that are not present in TFBS88. 
#There are 79 sites in TFBS88 that are not present in Biostrings.







#QC
#Stop264:
#temp_tfbs_sites <- subset(stoperator_biostrings_tfbs88,stoperator_biostrings_tfbs88$source_tfbs88 == TRUE & stoperator_biostrings_tfbs88$source_biostrings == FALSE)
#temp_biostrings_sites <- subset(stoperator_biostrings_tfbs88,stoperator_biostrings_tfbs88$source_tfbs88 == FALSE & stoperator_biostrings_tfbs88$source_biostrings == TRUE)
#temp_tfbs88_l5_sites <- subset(stoperator_site_predictions88_self,stoperator_site_predictions88_self$tfbs88_stoperator_target == 'l5')
#Of the 30 L5 site originally reported in Brown et al. 1997, the TFBS88 dataset does not include Site #3 (EMSA bound) or Site #30 (EMSA unbound)







#Now that the best TFBS cutoff has been determined (88%), 
#search all genomes with each PWM and determine how many sites of each motif are in each genome.

#Re-factor columns in the self-site dataset
stoperator_site_predictions88_self$tfbs88_stoperator_target <- factor(stoperator_site_predictions88_self$tfbs88_stoperator_target)
stoperator_site_predictions88_self$tfbs88_stoperator_motif <- factor(stoperator_site_predictions88_self$tfbs88_stoperator_motif)
stoperator_site_predictions88_self$tfbs88_motif_target <- factor(stoperator_site_predictions88_self$tfbs88_motif_target)
stoperator_site_predictions88_self$tfbs88_site_strand2 <- factor(stoperator_site_predictions88_self$tfbs88_site_strand2)



#Subset all sites predicted in each genome from each PWM above 88% threshold
stoperator_site_predictions80$site_self <- ifelse(is.element(stoperator_site_predictions80$stop_site_id,stoperator_site_predictions88_self$tfbs88_stop_site_id),TRUE,FALSE)

stoperator_site_predictions88 <- subset(stoperator_site_predictions80,
                                             stoperator_site_predictions80$site_rel_score >= 0.88)



stoperator_site_predictions88$stoperator_target <- factor(stoperator_site_predictions88$stoperator_target)
stoperator_site_predictions88$stoperator_motif <- factor(stoperator_site_predictions88$stoperator_motif)
stoperator_site_predictions88$motif_target <- factor(stoperator_site_predictions88$motif_target)
stoperator_site_predictions88$site_strand2 <- factor(stoperator_site_predictions88$site_strand2)

names(stoperator_site_predictions88) <- paste("tfbs88_",names(stoperator_site_predictions88),sep="")



stoperator_site_tfbs88_reduced <- subset(stoperator_site_predictions88,select=c("tfbs88_stop_site_id",
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


stoperator_site_tfbs88_reduced$tfbs88_site_abs_score <- round(stoperator_site_tfbs88_reduced$tfbs88_site_abs_score,2)
stoperator_site_tfbs88_reduced$tfbs88_site_rel_score <- round(stoperator_site_tfbs88_reduced$tfbs88_site_rel_score,2)



#output data for analysis with immunity data
setwd("~/scratch/stoperator_distance/output/")
write.table(stoperator_site_tfbs88_reduced,
            "stoperator_site_tfbs88_reduced.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)



# #output data for analysis with immunity data
# setwd("~/scratch/stoperator_distance/output/")
# write.table(stop_site_set80_df_reduced,
#             "stop_site_set80_df_reduced.csv",
#             sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)




































###misc code for QC analysis



stop_site_set80_df_self <- subset(stop_site_set80_df,
                                  as.character(stop_site_set80_df$seqnames) == as.character(stop_site_set80_df$ID))






# stop_site_set80_df$check1 <- as.character(stop_site_set80_df$ID) == as.character(stop_site_set80_df$TF)
# stop_site_set80_df$check2 <- as.character(stop_site_set80_df$seqnames) == as.character(stop_site_set80_df$ID)
# stop_site_set80_df$check3 <- as.character(stop_site_set80_df$seqnames) == as.character(stop_site_set80_df$TF)




stop_site_set80_df_self$seqnames <- factor(stop_site_set80_df_self$seqnames)
stop_site_set80_df_self$ID <- factor(stop_site_set80_df_self$ID)


stop_site_set80_df_bxb1_l5 <- subset(stop_site_set80_df,
                                  stop_site_set80_df$seqnames == 'bxb1' &
                                  stop_site_set80_df$ID == 'l5')

stop_site_set80_df_l5_bxb1 <- subset(stop_site_set80_df,
                                     stop_site_set80_df$seqnames == 'l5' &
                                       stop_site_set80_df$ID == 'bxb1')


setwd("~/scratch/stoperator_distance/output/")

write.table(stop_site_set80_df_self,
            "stop_site_set80_df_self.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)







stop_site_set95_df_self <- subset(stop_site_set95_df,
                                  as.character(stop_site_set95_df$seqnames) == as.character(stop_site_set95_df$ID))

stop_site_set95_df_self$seqnames <- factor(stop_site_set95_df_self$seqnames)
stop_site_set95_df_self$ID <- factor(stop_site_set95_df_self$ID)
write.table(stop_site_set95_df_self,
            "stop_site_set95_df_self.csv",
            sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)







