The repository 'mavrich_mbio_2019' contain the scripts
used for the analyses in the publication Mavrich & Hatfull, mBio, 2019. 
Below is a brief description of each script in the repository,
the requisite input files, and the output files.
The requisite input and output data files are available upon request.






1. Script: process_genome_distances.R

Description: compile whole genome nucleotide distance and whole genome
gene content dissimilarity metrics.



Input files:


    1. phage_metadata.csv
    
    Description: information about all phages in the database. Derived from 
    the same metadata file used for the analyze_immunity_data.R script.
    
    Data structure:
        1. Phage name
        2. Host
        3. Cluster
        4. Subcluster
        
        
    2. processed_mash_output.csv
    
    Description: pairwise whole genome nucleotide distance data generated by 
    Mash and processed using scripts published in Mavrich & Hatfull 2017.
    
    Data structure:
        1. Reference phage name
        2. Query phage name
        3. Mash distance
        4. Mash p-value
        5. Mash kmer count
        6. Reference-Query identifier
    
    
    3. pairwise_pham_proportions.csv
    
    Description: pairwise whole genome averaged pham proportion data generated
    by scripts published in Mavrich & Hatfull 2017.
    
    Data structure:
        1. Phage1 name
        2. Phage1 number of unshared phams
        3. Phage1 shared proportion
        4. Phage2 name
        5. Phage2 number of unshared phams
        6. Phage2 shared proportion
        7. Number of shared phams
        8. Average shared proportion
        9. Jaccard similarity
        10. Shared pham distribution mean
        11. Shared pham distribution median
        12. Shared pham distribution max
        13. Unshared pham distribution mean
        14. Unshared pham distribution median
        15. Unshared pham distribution max
        16. Unshared orpham count



Output files:


    1. genomic_distance_data.csv

    Description: pairwise comparisons that are used in the 
    analyze_immunity_data.R script. Dataset contains self comparisons
    (i.e. L5 compared to L5) as well as reciprocal comparisons (i.e.
    L5 compared to Trixie and Trixie compared to L5).
    
    Data structure:
        1. Phage1-Phage2 ID
        2. Modified Mash distance
        3. Gene content dissimilarity


    2. data_for_mode_prediction.csv

    Description: pairwise comparisons within intra-cluster boundaries
    (as defined by Mavrich & Hatfull 2017) so that evolutionary mode can be
    identified, as published in Mavrich & Hatfull 2017 
    (but not used in this publication).
    
    Data structure:
        1. Reference phage name
        2. Query phage name
        3. Modified Mash distance
        4. Gene content dissimilarity


    3. data_for_maxgcdgap_all.csv

    Description: all pairwise comparisons, which can be used to compute 
    MaxGCDGap as published in Pope et al. 2017
    (but not used in this publication).

    Data structure:
        1. Reference phage name
        2. Query phage name
        3. Modified Mash distance
        4. Gene content dissimilarity


    4. data_for_mash_network_nuc005.csv
    
    Description: all pairwise comparisons below Mash distance of 0.05. 
    This can be used to create phage networks
    (but not used in this publication).
    
    Data structure:
        1. Reference phage name
        2. Query phage name
        3. Modified Mash distance
        4. Mash p-value
        5. Mash kmer count
        6. Reference_Query ID






2. Script: process_stoperator_data.R

Description: using previously identified potential binding sites, create 
binding site position weight matrices (PWMs), compute pairwise distances
between PWMs, and scan all genomes for potential binding sites for each
PWM.



Input files:


    1. meme_stoperators.csv
    
    Description: table of stoperators within each genome, predicted by MEME.
    
    Data structure:
        1. Phage Name
        2. Unique stoperator ID 
        3. Stoperator sequence (13bp)


    2. folder of fasta-formatted genome sequences
    
    Description: all genomes that need to be analyzed. One file per genome,
    fasta-formatted.



Output files:


    1. biostrings_stoperators.csv
    
    Description: table of stoperator binding sites for each genome, derived
    from the input MEME stoperator data, and containing confirmed coordinate and
    sequence data using Biostrings tools.
    (This serves as a reference dataset and is not used in downstream analyses).
    
    
    Data structure:
        1. Phage name
        2. Forward sequence
        3. Reverse sequence
        4. Start coordinate
        5. End coordinate
        6. Strand
        7. Unique stoperator ID


    2. stoperator_pwm_distances.csv
    
    Description: pairwise genetic distances between PWMs that were generated
    using TFBSTools log2probratio.
    
    Data structure:
        1. Phage name #1
        2. Phage name #2
        3. Pearson distance
        4. Euclidean distance


    3. stoperator_site_predictions.csv
    
    Description: table of predicted stoperator sites in every genome that match
    every PWM, using TFBSTools and an 88% similarity score cutoff.
    
    Data structure:
        1. Unique stoperator ID
        2. Motif-Target ID
        3. Stoperator target genome
        4. Stoperator motif genome
        5. Start coordinate
        6. End coordinate
        7. Strand
        8. Sequence
        9. Absolute score
        10. Relative score
        11. Endogenous stoperator






3. Script: analyze_immunity_data.R

Description: compare superinfection immunity phenotypes with genomic 
characteristics to identify how immunity changes as phages evolve.



Input files:


    1. immunity_data.csv
    
    Description: table containing all immunity assay data used for this
    analysis.
    
    Data structure:
        1. Unique immunity assay identifier
        2. Unique immunity assay set identifier
        3. Date of immunity assay
        4. Lab notebook of immunity assay
        5. Lab notebook page of immunity assay
        6. Unique bacterial strain identifier
        7. Prophage (if applicable)
        8. Cloned Repressor (if applicable)
        9. Strain type (lysogen or cloned repressor)
        10. Defending phage
        11. Infecting phage
        12. Assay type (multiple titer or single titer)
        13. Notes regarding bacterial lawn
        14. Lawn reliability score (0 = unreliable; 3 = reliable)
        15. Tested phage titer
        16. Phage spot reliability score (0 = unreliable; 3 = reliable)
        17. Observed infection strength
        18. Observed turbidity
        19. Observed plaque/spot size
        20. Observed plaques (yes or no)
        21. Infection score (systematically scored infection phenotype)


    2. genomic_distance_data.csv

    Description: pairwise Mash-based nucleotide distances and Pham-based gene
    content dissimilarities for all phages in the Actino1321 database,
    computed as previously described in Mavrich & Hatfull 2017.
    This file is generated by the process_genome_distances.R script.


    3. phage_metadata.csv

    Description: manual compilation of misc. data for each phage in the 
    Actino1321 database, generated from multiple tools.

    Data structure:
        1. Phage name
        2. Host
        3. Cluster
        4. Subcluster
        5. Genome size
        6. Lysogen type (extrachromosomal, integration, none)
        7. Integrase pham number
        8. ParA pham number 
        9. Phage source (environment or lab)
        10. Parent phage name
        11. Repressor is functional (yes, no, NA)
        12. Lysogen can be generated (no, unknown, yes, NA)
        13. Repressor helix-turn-helix domain sequence
        14. Length of repressor full sequence
        15. Length of repressor N-terminal sequence
        16. Length of repressor C-terminal sequence
        17. ParB pham number
        18. Gene content clade (clade2 = the "L5 clade")
        19. Coordinate of Pleft locus for manual alignment
        20. Coordinate of repressor locus for manual alignment
        21. Coordinate of genome center for manual alignment


    4. repressor_336_distance_data.csv

    Description: pairwise genetic distances of 336 Immunity Repressor homologs.

    Data structure:
        1. Phage1-Phage2 ID
        2. PhyML distance of full sequence from MAFFT alignment
        3. PhyML distance of N-terminal sequence from MAFFT alignment
        4. PhyML distance of C-terminal sequence from MAFFT alignment
        5. EMBOSS distance of full sequence from MAFFT alignment
        6. EMBOSS distance of N-terminal sequence from MAFFT alignment
        7. EMBOSS distance of C-terminal sequence from MAFFT alignment


    5. cas4_311_distance_data.csv

    Description: pairwise genetic distances of 311 Cas4-family homologs.

    Data structure:
        1. Phage1-Phage2 ID
        2. EMBOSS distance of full sequence from MAFFT alignment


    6. endovii_306_distance_data.csv

    Description: pairwise genetic distances of 306 Endonuclease VII homologs.

    Data structure:
        1. Phage1-Phage2 ID
        2. EMBOSS distance of full sequence from MAFFT alignment


    7. dnapol_311_distance_data.csv

    Description: pairwise genetic distances of 311 DNA Polymerase homologs.

    Data structure:
        1. Phage1-Phage2 ID
        2. EMBOSS distance of full sequence from MAFFT alignment


    8. portal_311_distance_data.csv

    Description: pairwise genetic distances of 311 Portal homologs.

    Data structure:
        1. Phage1-Phage2 ID
        2. EMBOSS distance of full sequence from MAFFT alignment


    9. stoperator_pwm_distances.csv

    Description: pairwise distances between 327 position weight matrices 
    generated from predicted stoperators within each of the 327 Cluster A 
    phage genomes. This file is generated by the process_stoperator_data.R
    script.


    10. stoperator_site_predictions.csv

    Description: table of predicted stoperator sites from all 327 Cluster A 
    phage genomes using all 327 position weight matrices. This file is 
    generated by the process_stoperator_data.R script.


    11. reciprocal_infection_matrix.csv

    Description: table of averaged reciprocal infection scores used to 
    generate the heatmap in Figure 5a. 

    Data structure:
        1. Rows = infecting phage
        2. Columns = defending phage



Output files:


    1. Table_S1_averaged_immunity_data.csv
    
    Description: Processed infection data, with low-confidence data removed
    and all replicate assay results averaged. Data from this table was used
    to generate heatmaps in several figure panels. This file is also used to
    generate the reciprocal_infection_matrix.csv file.

    Data structure:
        1. Strain type
        2. Defending phage
        3. Infecting phage
        4. Whole genome nucleotide distance
        5. Whole genome gene content dissimilarity
        6. Repressor distance
        7. Stoperator motif distance
        8. Immunity assay replicates
        9. Minimum infection score
        10. Maximum infection score
        11. Average infection score
        
    2. Misc. data analysis plots used to create figures for the publication.





