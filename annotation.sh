#! /bin/bash

INPUT_GENOME_FASTA="..."

# RepeatMasker
#
# This process should be done prior to pipeline execution, the reference repeat database can be reused.
#
# FamDB can be downloaded from https://www.dfam.org/releases/current/families/FamDB/
#
# famdb.py -i Libraries/RepeatMaskerLib.h5 families -f embl  -a -d Insecta  > Insecta_ad.embl
# util/buildRMLibFromEMBL.pl Insecta_ad.embl > $FASTA_REPEAT_SEQUENCE_MASK
#
FASTA_REPEAT_SEQUENCE_MASK="..."

BUSCO_LINAGE="..."
PERFORM_MASKING=true

ROOT_WORKING_DIR=$(pwd)

busco --cpu $NCPU \
	-l /gpfs/home/meiyang/opt/insecta_odb10 \
	-m genome --force -o $ROOT_WORKING_DIR/BUSCO \
	-i $INPUT_GENOME_FASTA \
	--offline

if [ $PERFORM_MASKING = true ] ; then
    # 3. Repeat annotation and genome mask

    ## 3.1 RepeatModeler v2 & RepeatMasker
    #  Input:
    #     - genome.fa
    #  Output:
    #     - genome.fa.masked

    # 3.1.1 Building reference repeat database
    mkdir $ROOT_WORKING_DIR/01_RepeatModeler
    cd $ROOT_WORKING_DIR/01_RepeatModeler
    # BuildDatabase - Format FASTA files for use with RepeatModeler
    REPEAT_DB_NAME="GDB" # Given by the user, can be any name, has nothing to do with external databaase.
    BuildDatabase -name $REPEAT_DB_NAME -engine ncbi $INPUT_GENOME_FASTA | tee BuildDatabase.log
    # OUTPUT:
    #    - $REPEAT_DB_NAME.nhr
    #    - $REPEAT_DB_NAME.nin
    #    - $REPEAT_DB_NAME.nnd
    #    - $REPEAT_DB_NAME.nni
    #    - $REPEAT_DB_NAME.nog
    #    - $REPEAT_DB_NAME.nsq


    # RepeatModeler
    RepeatModeler -engine ncbi -threads $NCPU -database $REPEAT_DB_NAME -LTRStruct | tee RepeatModele.log
    # OUTPUT: Produce soft-masked genome.
    #    - $REPEAT_DB_NAME-families.fa - Consensus sequences for each family identified.
    #    - $REPEAT_DB_NAME-families.fa - Seed alignments for each family identified.
    #    - $REPEAT_DB_NAME-rmod.log - Execution log.  Useful for reproducing results.

    mkdir $ROOT_WORKING_DIR/02_RepeatMasker
    cd $ROOT_WORKING_DIR/02_RepeatMasker

    if [ $FASTA_REPEAT_SEQUENCE_MASK != "" ]; then
        cat $ROOT_WORKING_DIR/01_RepeatModeler/$REPEAT_DB_NAME-families.fa $FASTA_REPEAT_SEQUENCE_MASK > repeat_db.fa
    else
        cp $ROOT_WORKING_DIR/01_RepeatModeler/$REPEAT_DB_NAME-families.fa $ROOT_WORKING_DIR/02_RepeatMasker/repeat_db.fa
    fi

    # 3.1.2 Running RepeatMasker for genome masking
    # RepeatMasker (output: genome.fa.masked)
    RepeatMasker -xsmall -gff -html -lib repeat_db.fa -pa $NCPU -dir $ROOT_WORKING_DIR/02_RepeatMasker $INPUT_GENOME_FASTA | tee RepeatMasker.log

    MASKED_GENOME="$ROOT_WORKING_DIR/02_RepeatMasker/$(basename ${INPUT_GENOME_FASTA}).masked"

    # 3.2 EDTA - Identify and annotate transposable elements (TEs) and repetitive sequences. Produce soft-masked genome.
    # EDTA.pl --genome $INPUT_GENOME_FASTA --species others --sensitive 1 --anno 1 --evaluate 1 --threads 30
fi


# 4. Gene prediction

## 4.1 RNA-seq based gene prediction

### 4.1.1 HISAT2 & StringTie
#         Input:
#            - masked geome (genome.fa.masked)
#            - transcriptome (reads.fq)
#         Output:
#            - stringtie.gtf

# Build genome index
if [ $PERFORM_MASKING = true ] && [ -e $MASKED_GENOME ] ; then
    hisat2-build -p $NCPU $MASKED_GENOME genome
else
    hisat2-build -p $NCPU $INPUT_GENOME_FASTA genome
fi

# Mapping to genome
# single end
hisat2 --threads $NCPU -x genome --summary-file hisat2_summary.txt --quiet --dta -U reads.fq | samtools sort -@ $NCPU > reads.bam 

## OR ##

# paired end
hisat2 --threads $NCPU -x genome --summary-file hisat2_summary.txt --quiet --dta -1 reads_1.fq -2 reads_2.fq | samtools sort -@ $NCPU > reads.bam

# GTF merging
samtools merge -@ $NCPU merged.bam `ls *bam`
stringtie -p $NCPU -o stringtie.gtf merged.bam # OUTPUT: stringtie.gtf

### 4.1.2 TransDecoder
#         Input:
#            - masked geome (genome.fa.masked)
#            - stringtie.gtf
#         Output:
#            - transcript_alignments.gff3

gtf_genome_to_cdna_fasta.pl stringtie.gtf genome.fa > transcripts.fasta
gtf_to_alignment_gff3.pl stringtie.gtf > transcripts.gff3
TransDecoder.LongOrfs -t transcripts.fasta

# homology search
blastp -query transdecoder_dir/longest_orfs.pep -db uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 28 > blastp.outfmt6
hmmscan --cpu 28 --domtblout pfam.domtblout Pfam-A.hmm transdecoder_dir/longest_orfs.pep

TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcript_alignments.gff3

## OUTPUT: transcript_alignments.gff3

## 4.2 Ab initio gene prediction
### BRAKER v3
#   Input:
#   
#   - masked geome (genome.fa.masked)
#   - homology protein (OrthoDB), Arthropoda.fa
#   
#   Parameters
#   
#   ``` bash
#   --species=<species_name>
#   --min_contig, less than genome N50
#   ```

braker.pl --genome=genome.fa \
	--species=Sfru \
	--prot_seq=Arthropoda.fa \
	--bam=merged.bam \
	--threads 30 \
	--gff3 --workingdir=out

python gff_rename.py braker.gff3 sfru > gene_predictions.gff3

## 4.3 Homology-based gene prediction

### miniprot
#   Input:
#   
#   - masked geome (genome.fa.masked)
#   - homology protein (OrthoDB), Arthropoda.fa

miniprot -t28 -d genome.mpi genome.fa.masked 
miniprot -It28 --gff genome.mpi protein.fasta > miniprot.gff3 ## OUTPUT: miniprot.gff3

# 5. EVidenceModeler (EVM)

# Input:
#    - masked geome (genome.fa.masked)
#    - weights.txt
#    - GFF3 file
#        - gene_predictions.gff3 (from step 4.2)
#        - protein_alignments.gff3 (from EVidenceModeler)
#        - transcript_alignments.gff3 (from step 4.1.2)

# Example of weights.txt

# ```
# PROTEIN	miniprot	5
# ABINITIO_PREDICTION	BRAKER3	10
# OTHER_PREDICTION	transdecoder	10
# ```

# The gff3 file of miniprot should be reformated.
python ~/software/EVidenceModeler-v2.1.0/EvmUtils/misc/miniprot_GFF_2_EVM_alignment_GFF3.py miniprot.gff3 > protein_alignmentss.gff3

EVidenceModeler \
	--sample_id speceis \
	--genome genome.fa \
	--weights weights.txt  \
	--gene_predictions gff/gene_predictions.gff3 \
	--protein_alignments gff/protein_alignments.gff3 \
	--transcript_alignments gff/transcript_alignments.gff3 \
	--segmentSize 100000 --overlapSize 10000 --CPU 20

# 6. OGS annotation

## 6.1 OGS annotation updates
### 6.1.1 PASApipeline
#           - masked geome (genome.fa)
#           - species.evm.gff3
#           - stringtie.gtf

# PASA alignment Assembly
util/gtf_genome_to_cdna_fasta.pl stringtie.gtf genome.fa > transcripts.fasta
bin/seqclean transcripts.fasta

Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g genome.fa -t transcripts.fasta.clean -T -u transcripts.fasta --ALIGNERS blat --CPU 16

# two cycles !!! of annotation loading, annotation comparison, and annotation updates
# check gff3
misc_utilities/pasa_gff3_validator.pl species.evm.gff3

# load annotation
scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g genome.fa -P species.evm.gff3

# update
# annotCompare.config, set up the mysql database name same as alignAssembly.config
Launch_PASA_pipeline.pl -c annotCompare.config -A -g genome.fa -t transcripts.fasta.clean

peaks2utr -p 20 species.evm.gff3 merged.bam

# rename gff3
# species name, Sfru
python gff_rename.py gene_structures_post_PASA_updates.gff3 Sfru > Sfru.gff3


## 6.2 Collect GFF, cds, PEP

# Input:
#    - gene_structures_post_PASA_updates.gff3
# Output:
#     - Sfru.gff3 (Sfru_no_alt.gff3)
#     - cds.fa (cds_no_alt.fa)
#     - pep.fa (pep_no_alt.fa)

# collect
gffread Sfru.gff3 -g genome.fa -x Sfru_cds.fa -y Sfru_pep.fa

# collect no alt gff, cds, pep 
python Collect_no_alt.py pep.fa cds.fa Sfru.gff3
# no_alt.gff3, cds_no_alt.fa, pep_no_alt.fa

# rename gff3
# species name, Sfru
python gff_rename.py gene_structures_post_PASA_updates.gff3 Sfru > Sfru.gff3

# collect
gffread Sfru.gff3 -g genome.fa -x Sfru_cds.fa -y Sfru_pep.fa

# collect no alt gff, cds, pep 
python Collect_no_alt.py pep.fa cds.fa Sfru.gff3
# no_alt.gff3, cds_no_alt.fa, pep_no_alt.fa
