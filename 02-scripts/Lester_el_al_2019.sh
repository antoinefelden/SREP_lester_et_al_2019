###########################################################################################################################################
################ Step 1: Generate a sorted sorted file with all the reads mapped on the reference genome

# copy reference genome
cp /srv/global/scratch/groups/sbs/hymenoptera_genomes/argentine_ant/RefSeq/GCF_000217595.1_Lhum_UMD_V04_genomic.fna .

# copy raw reads
raw_dir=/home/scifachpc-fs01/feldenan/raw_data/northland_viral_loads/Cleandata
find $raw_dir -iname '*.fq.gz' -exec cp {} . \;

hisat2-build GCF_000217595.1_Lhum_UMD_V04_genomic.fna GCF_000217595.1_Lhum_UMD_V04_genomic.fna

for SAMPLE in CS297 CS298 CS300 CS301 CS302 CS303 CS305 CS306 CS308 CS320 CS321 CS323;

do
R1=${SAMPLE}_R1.fq.gz
R2=${SAMPLE}_R2.fq.gz
BAM=nth_samples_${SAMPLE}.bam
SAM=nth_samples_${SAMPLE}.sam
UMPD=nth_samples_unmapped_${SAMPLE}.fq.gz

# align reads with Hisat2
hisat2 -1 $R1 -2 $R2 -p $NSLOTS GCF_000217595.1_Lhum_UMD_V04_genomic.fna | samtools sort > $BAM
hisat2 -1 $R1 -2 $R2 -p $NSLOTS GCF_000217595.1_Lhum_UMD_V04_genomic.fna -S $SAM --un-conc-gz $UMPD
samtools sort $SAM > $BAM
rm $SAM
done

mkdir ./unmapped_reads
mv nth_samples_unmapped_* ./unmapped_reads
rm *.sam
rm *.fq.gz

###########################################################################################################################################
################ Step 2:  StringTie -e to get read coverage table (simple step, no novel transcripts)

input=/srv/global/scratch/feldenan/Northland_viral_loads/xo10_northland_viral_loads_alignemnt_Hisat2_752372
cp $input/*.bam .

cp /srv/global/scratch/groups/sbs/hymenoptera_genomes/argentine_ant/RefSeq/GCF_000217595.1_Lhum_UMD_V04_genomic.gff .

for SAMPLE in CS297 CS298 CS300 CS301 CS302 CS303 CS305 CS306 CS308 CS320 CS321 CS323;
do
    mkdir ./nth_${SAMPLE}
    BAM=nth_samples_${SAMPLE}.bam
    GTF=./nth_${SAMPLE}/nth_${SAMPLE}.gtf
    stringtie $BAM -p $NSLOTS -o $GTF -G GCF_000217595.1_Lhum_UMD_V04_genomic.gff -e
done

cp /srv/global/scratch/feldenan/ZZ_programs/StringTie/prepDE.py .
chmod +x prepDE.py

./prepDE.py -i ./ -l 150 -p "nth_"


###########################################################################################################################################
################ Step 3:  Assemble unmapped transcripts with Trinity

cp /srv/global/scratch/feldenan/Northland_viral_loads/samples_viral_loads.txt .

 $TRINITY_HOME/Trinity \
 --samples_file samples_viral_loads.txt \
 --CPU=$NSLOTS \
 --max_memory=256G \
 --seqType=fq

###########################################################################################################################################
################ Step 4:  BLAST de-novo assembled transcripts against the NCBI virus database

cp /home/scifachpc-fs01/feldenan/references/viral/viral.1.protein.faa.gz .
cp /home/scifachpc-fs01/feldenan/references/viral/viral.2.protein.faa.gz .
cat viral.1.protein.faa.gz viral.2.protein.faa.gz > viral.all.protein.faa.gz
gunzip viral.all.protein.faa.gz
makeblastdb -in viral.all.protein.faa -dbtype prot

input=/srv/global/scratch/feldenan/Northland_viral_loads/xo21_trinity_assembly_752416
cp $input/Trinity.fasta .

blastx -query Trinity.fasta -db viral.all.protein.faa -num_threads $NSLOTS -max_target_seqs 10 -outfmt 6 -evalue 1e-5 > viruses_nthland_refseq.blastx.outfmt6

rm *.gz

###########################################################################################################################################
################ Step 4':  BLAST de-novo assembled transcripts against the Viljakainen virus database

#input files
input=/srv/global/scratch/feldenan/Northland_viral_loads/xo21_trinity_assembly_752416
cp $input/Trinity.fasta .

cp /home/scifachpc-fs01/feldenan/references/viral/Viljakainen/* .

# prot databae
cat *_prots.fasta > viljakainen_all_viral_prots.fasta
makeblastdb -in viljakainen_all_viral_prots.fasta -dbtype prot
blastx -query Trinity.fasta -db viljakainen_all_viral_prots.fasta -num_threads $NSLOTS -max_target_seqs 10 -outfmt 6 -evalue 1e-5 > viljakainen_viruses_nthland.blastx.outfmt6

###########################################################################################################################################
################ Step 5:  Transcript Abundance

# Copy Trinity assembly
cp /srv/global/scratch/feldenan/Northland_viral_loads/xo21_trinity_assembly_752416/Trinity.fasta .

# copy samples_file
cp /srv/global/scratch/feldenan/Northland_viral_loads/samples_viral_loads.txt .

perl $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq \
--samples_file samples_viral_loads.txt \
--est_method eXpress --aln_method bowtie2 --trinity_mode --prep_reference \
--thread_count $NSLOTS

################ Stage output back for collection

cp -p -r . $output_dir

history

mv /srv/global/scratch/feldenan/*$JOB_ID $output_dir
