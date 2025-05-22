### **TASK 1**
Tuber[Organism]  -->  118562 seqs
Tuber[Organism] AND (ribosomal RNA) --> 5746 seqs
Tuber[Organism] AND (ribosomal RNA) AND ( "500"[SLEN] : "5000"[SLEN] )   --> 5087 seqs
Tuber[Organism] AND (ribosomal RNA) AND ( "500"[SLEN] : "5000"[SLEN] ) NOT environmental  --> 4773 seqs
(Sequence number generated on Apr 15 2025)

### **TASK 2**
- Select blastn
- lead the sequence or sequecne fasta
- Enter in the Organism window the scientific name (or taxa id from NCBI Taxonomy) to limit the results to that species, try a blast, if you get no results try to progressively relax blast parameters . Click on **Algorithm parameters** and progressively relax "word size" and "match/mismatch" 
- Select the best hit. If you have problem finding a trustable sequence for yeast you may think is because is too distantly related to animals, well it is quite distant, but this is not the case. The reality is COX1 is not very often sequences for yeasts, so, where would you search instead of  the nucleotide collection?
- Download the best hit in fasta: tick on the sequence -> Download -> FASTA
### **Optional**
concatenate fastas
```bash
cat *.fas > cat_besthit_COX1.fasta
```
align the multifasta file
```bash
mafft --localpair --maxiterate 1000 cat_besthit_COX1.fasta > cat_besthit_COX1_mafft_aligned.fasta
```
check if the alignment makes sense (aliview, jalview, MEGA... ) and delete unaligned blocks, in this case with TrimAL heuristic method
```bash
trimal -in cat_besthit_COX1_mafft_aligned.fasta -out cat_besthit_COX1_mafft_aligned_trimmed.fas -automated1
```
phylogenetic inference based on 
```bash
iqtree2 -s cat_besthit_COX1_mafft_aligned_trimmed.fas -m MFP -bb 1000 -nt AUTO
```
Load the .treefile in ITOL to visualize the phylogenetic tree

### **TASK 3**
First unzip the databases
```shell
for i in ./data/NCBI_BLAST_databases/*.gz ; do tar -xvzf $i -C ./data/NCBI_BLAST_databases; done
```
Check the fasta
```bash
less ./data/sequence.fasta  
```
BLAST on the databases the query sequence
```bash
blastn -num_threads 4 -query ./data/sequence.fasta -db ./data/NCBI_BLAST_databases/16S_ribosomal_RNA -out ./results/task3/blastout_sequence_vs_16S.txt

```
and so on...
Or you can loop through databases!
```bash
# Define the names of the databases in a list
databases=("16S_ribosomal_RNA" "ITS_eukaryote_sequences" "LSU_eukaryote_rRNA" "LSU_prokaryote_rRNA" "SSU_eukaryote_rRNA" "18S_fungal_sequences" "28S_fungal_sequences" "ITS_RefSeq_Fungi")

# Loop over them
for db in "${databases[@]}"; do
    # Construct your blastn command
    blastn \
        -query ./data/sequence.fasta \
        -db "./data/NCBI_BLAST_databases/${db}" \
        -num_threads 4 \
        -out "./results/task3/blastout_sequence_vs_${db}.txt"
done
```
Using default setting you should not get any result, try to relax search parameters, as you did on BLAST webserver, for example with:
```bash
# Define the names of the databases in a list
databases=("16S_ribosomal_RNA" "ITS_eukaryote_sequences" "LSU_eukaryote_rRNA" "LSU_prokaryote_rRNA" "SSU_eukaryote_rRNA" "18S_fungal_sequences" "28S_fungal_sequences" "ITS_RefSeq_Fungi")

# Loop over them
for db in "${databases[@]}"; do
    # Construct your blastn command
    blastn \
        -query ./data/sequence.fasta \
        -db "./data/NCBI_BLAST_databases/${db}" \
        -num_threads 4 \
        -out "./results/task3/blastout_sequence_vs_${db}_relaxed.txt" \
        -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -word_size 7
done
```


Get the first n hits of each blast result from the blast output file:
```bash
$ find ./results/ -name blastout*_relaxed.txt -exec  grep -A 2 --with-filename "significant alignments" {} \;
```


















