The repository 'mavrich_mbio_2019' contain the scripts
used for the analyses in the publication Mavrich & Hatfull, mBio, 2019. 
The requisite input and output data files are available upon request.
Below is a brief description of each file in the repository. 

Notes: 





Input files:


1. genomic_distance_data.csv

Description: Pairwise Mash-based nucleotide distances or Pham-based gene
content dissimilarities for all phages in the Actinobacteriophage_1321
database, computed as previously described in Mavrich & Hatfull 2017.


2. phage_metadata.csv

Compilation of misc. data for each phage in the Actinobacteriophage_1321
database, generated from multiple tools.



3. repressor_distance_data.csv

4. portal_distance_data.csv

5. recb_distance_data.csv

6. repressor336_distance_data.csv

Description: Pairwise genetic distances of 336 Immunity Repressor homologs.

7. cas4311_distance_data.csv

Description: Pairwise genetic distances of 311 Cas4-family homologs.


8. endovii306_distance_data.csv

Description: Pairwise genetic distances of 306 Endonuclease VII homologs.


9. dnapol311_distance_data.csv

Description: Pairwise genetic distances of 311 DNA Polymerase homologs.


10. portal311_distance_data.csv

Description: Pairwise genetic distances of 311 Portal homologs.


11. stoperator_pwm_data.csv

Description: Pairwise distances between 327 position weight matrices generated
from predicted stoperators within each of the 327 Cluster A phage genomes.


12. stoperator_site_predictions.csv

Description: Table of predicted stoperator sites from all 327 Cluster A phage
genomes using the position weight matrices listed above.






Data analysis scripts:

1. analyze_immunity_data.R

Description: RStudio script to ...





    Output files:
    
    1. conf_assay_strain_def_chal_average.csv
    
    Description: Processed infection data, with low-confidence data removed
    and all replicate assay results averaged.
    
    
    2. reciprocal_data_alpha_ordered.csv
    
    Description:
    
    
    3. mutant_analysis.csv
    
    
    Description:








