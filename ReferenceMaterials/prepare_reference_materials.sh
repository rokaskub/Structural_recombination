# Needs ncbi-acc-download --> pip install ncbi-acc-download
# NCBI TAIR10.1 --> https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000001735.4/

# Refseq ids 
CHLORO_REFSEQ=NC_000932.1
MITO_REFSEQ=NC_037304.1

# Wanted filenames
CHLORO=ChrC.$CHLORO_REFSEQ.fa
MITO=ChrM.$MITO_REFSEQ.fa

# Download and saved them
echo "Downloading $CHLORO_REFSEQ & $MITO_REFSEQ from NBCI..."
ncbi-acc-download --format fasta $CHLORO_REFSEQ --out $CHLORO
ncbi-acc-download --format fasta $MITO_REFSEQ --out $MITO

# Requires module load, samtools, gatk, hisat2, bwa
createIndexes (){
    echo "Creating indexes for $1..." 
    module load samtools hisat2 bwa gatk  
    samtools faidx $1
    gatk CreateSequenceDictionary --REFERENCE $1
    hisat2-build $1 $1
    bwa index $1
    echo "✅  Indexes done for $1"
}

createIndexes $CHLORO
createIndexes $MITO

# Mask the MITO
REPEATS=repeats.sorted.bed

MITO_MASKED=${MITO/.fa/.repeat_masked.fa}

module load bedtools
bedtools maskfasta -fi $MITO -bed $REPEATS -fo $MITO_MASKED 

createIndexes $MITO_MASKED









