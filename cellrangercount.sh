
# Get Mouse reference (mm10) dataset required for Cell Ranger.
# Download - 9.6 GB - md5sum: 8ce6bc561e2554701fc43871301042e6

curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-3.0.0.tar.gz
# or 
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-3.0.0.tar.gz

# make a Cell Ranger compatible "pre-mRNA" reference package according to the instructions:

# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/3.0/advanced/references#premrna

### bash script for cellranger count 

#!/bin/bash

module add bioinfo-tools
module add cellranger/3.0.1

# Define variables
REFERENCE="/path/to/refdata-mm10_premrna"
FASTQ_DIR="/path/to/fastq_files"

# Get list of sample names from their paths
SAMPLES=($(ls -d "$FASTQ_DIR"/* | awk -F'/' '{print $NF}'))

# Loop through sample names
for SAMPLE_ID in "${SAMPLES[@]}"; do
    echo "Processing sample: $SAMPLE_ID"
    
    # Run cellranger count command
    cellranger count \
        --id="$SAMPLE_ID" \
        --fastqs="$FASTQ_DIR/$SAMPLE_ID" \
        --transcriptome="$REFERENCE" 
    
    echo "Sample $SAMPLE_ID processing complete."
done

echo "All samples processed."

