### **TASK 1**
Tuber[Organism]  ->  118562 seqs

Tuber[Organism] AND (ribosomal RNA) -> 5746 seqs

Tuber[Organism] AND (ribosomal RNA) AND ( "500"[SLEN] : "5000"[SLEN] )   -> 5087 seqs

Tuber[Organism] AND (ribosomal RNA) AND ( "500"[SLEN] : "5000"[SLEN] ) NOT environmental  -> 4773 seqs

(Sequence number generated on Apr 15 2025)
Download the sequences ordered by length (useful later) in fasta 
### **TASK 2**
```bash
$ mkdir results/task2
```
- Select blastn
- load or copy/paste into the window sequence.fasta
- Enter in the Organism window the scientific name (or taxa id from NCBI Taxonomy) to limit the results to that species, try a blast, if you get no results try to progressively relax blast parameters . Click on **Algorithm parameters** and progressively relax "word size" and "match/mismatch" 
- Select the best hit. If you have problem finding a trustable sequence for yeast you may think is because is too distantly related to animals, well it is quite distant, but this is not the case. The reality is COX1 is not very often sequences for yeasts, but you get hist from complete mitochondrial yests's genomes. Go to "Alignments" and download as fasta only the portion of the mitochondrial genomes that aligns to your query.
-  For the other species just download the best hit in fasta: tick on the sequence -> Download -> FASTA
### **Optional**
Cat fastas
```bash
 $ cat ./results/task2/*.fasta > results/task2/cat_besthit_COX1.fasta
```

Align the multifasta file
```bash
$ singularity pull https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_1

$ singularity exec mafft\:7.525--h031d066_1 mafft --localpair --maxiterate 1000 --thread 4 results/task2/cat_besthit_COX1.fasta > results/task2/cat_besthit_COX1_mafft_aligned.fasta
```

Check if the alignment makes sense (aliview, jalview, MEGA... ) and delete unaligned blocks, in this case with TrimAl heuristic method
```bash
$ singularity pull https://depot.galaxyproject.org/singularity/trimal:1.5

$ singularity exec trimal:1.5 trimal -in results/task2/cat_besthit_COX1_mafft_aligned.fasta -out results/task2/cat_besthit_COX1_mafft_aligned_trimmed.fas -automated1 -keepheader
```

Phylogenetic inference
```bash
$ sed -i 's/ /_/g' results/task2/cat_besthit_COX1_mafft_aligned_trimmed.fas #substitute spaces with _ or name will be truncated at first space

$ singularity pull https://depot.galaxyproject.org/singularity/iqtree:2.3.6--h503566f_1
$ singularity exec iqtree\:2.3.6--h503566f_1 iqtree2 -s results/task2/cat_besthit_COX1_mafft_aligned_trimmed.fas -m MFP -bb 1000  -nt 4 -af fasta
```

Load the .treefile in ITOL to visualize the phylogenetic tree
```bash
$ less ./results/task2/cat_besthit_COX1_mafft_aligned_trimmed.fas.treefile 
```
### **TASK 3**
```bash
$ mkdir results/task3
```
First check file integrity and unzip the databases 
```shell
$ cd results/NCBI_databases/
$ md5sum -c *.md5
$ cd -

$ for i in ./results/NCBI_databases/*.gz ; do tar -xvzf $i -C ./results/NCBI_databases; done

```
Check the fasta:
```bash
$ cat ./data/sequence.fasta  
```
- Is a DNA sequence of lenght 166 nt
```bash
$ grep -v "^>" data/sequence.fasta | wc -c
```

- Blastn, we do not know if it is coding DNA, so we should stick to blastn for now
BLAST on the databases the query sequence
```bash
singularity exec blast:2.16.0--h66d330f_4 blastn -num_threads 4 -query ./data/sequence.fasta -db ./data/NCBI_BLAST_databases/16S_ribosomal_RNA -out ./results/task3/blastout_sequence_vs_16S.txt

```
and so on...
or you can loop through databases!
```bash
# Names of the databases in a list
databases=("16S_ribosomal_RNA" "ITS_eukaryote_sequences" "LSU_eukaryote_rRNA" "LSU_prokaryote_rRNA" "SSU_eukaryote_rRNA" "18S_fungal_sequences" "28S_fungal_sequences" "ITS_RefSeq_Fungi")

# Loop over them
for db in "${databases[@]}"; do
    #blastn command
    singularity exec blast:2.16.0--h66d330f_4 blastn \
        -query ./data/sequence.fasta \
        -db "./data/NCBI_databases/${db}" \
        -num_threads 4 \
        -out "./results/task3/blastout_sequence_vs_${db}.txt"
done
```
check results:
```bash 
$ cat results/task3/blastout_sequence_vs_*
```
no hits!

- spacer are very variable (purifying selection does not select negatively substitutions), relaxed parameter are necessary to find hits in these smaller rDNA databases we are using (exactly as we already did on BLAST webserver):
```bash
# Names of the databases in a list
databases=("16S_ribosomal_RNA" "ITS_eukaryote_sequences" "LSU_eukaryote_rRNA" "LSU_prokaryote_rRNA" "SSU_eukaryote_rRNA" "18S_fungal_sequences" "28S_fungal_sequences" "ITS_RefSeq_Fungi")

# Loop over them
for db in "${databases[@]}"; do
    # blastn command
    singularity exec blast:2.16.0--h66d330f_4 blastn \
        -query ./data/sequence.fasta \
        -db "./data/NCBI_databases/${db}" \
        -num_threads 4 \
        -out "./results/task3/blastout_relaxed_sequence_vs_${db}.txt" \
        -reward 2 -penalty -3 -gapopen 5 -gapextend 2 -word_size 7
done
```
check results:
```bash 
$ cat results/task3/blastout_relaxed_sequence_vs_* 
```
we got hits! quite short but we have some.
Parse the output (but we do not know from which database the hits come from, except for the name of hits):
``` bash
$ cat results/task3/blastout_relaxed_sequence_vs_* | grep -A 2 --with-filename "significant alignments"
```
Better luck with:
```bash
$ find ./results/task3 -name blastout_relaxed_sequence_vs_* -exec  grep -A 2 --with-filename "significant alignments" {} \;
```
decent e-value from ITS eukaryote and refseq fungi, same bitscore (same sequence aligns) but slightly different e-value (e-value depends on the length of the database).
It is an ITS from Fungi, putatively *Tuber cryptobrumale*


To clearly delimit which ITS we can build a MSA that more clearly shows the boundaries of each locus in rDNA, for example using task1 tuber_rDNA.fas file:
```bash
# get the first (longest) seqs from file
$ head -n 3000 ./results/task1/tuber_rDNA.fas > results/task3/tuber_rDNA_longest_plus_sequence.fas
# add the query sequence
$ cat ./data/sequence.fasta >> results/task3/tuber_rDNA_longest_plus_sequence.fas

$ singularity exec mafft\:7.525--h031d066_1 mafft --thread 4 results/task3/tuber_rDNA_longest_plus_sequence.fas > results/task3/tuber_rDNA_longest_plus_sequence_mafft_aligned.fas
```
download and visualize in aliview (or similar software). 
The file in results/backup is manually edited to keep only meaningful sequences, as downloaded sequence are often non-overlapping portions of rDNA. It is an ITS1.

### **TASK 4**

Show blast results:
```bash
cat results/blastout_sequence_vs_GCA_000151645.1_tabfmt.txt
```
A total of 5 hits are present on three scaffolds. 
Grepping the headers of the hits:
```
$ grep  "FN430142.1" results/GCA_000151645.1_ASM15164v1_genomic.fna
$ grep  "FN430206.1" results/GCA_000151645.1_ASM15164v1_genomic.fna
$ grep  "FN430270.1" results/GCA_000151645.1_ASM15164v1_genomic.fna 
```
One hit on scaffold 24 , two on scaffolds 298 and 355, hits on smaller scaffolds possibly presenting two repeated rDNA sequences (two ITS1 hits ~8000 bp apart)
Multiple hits are common, as rDNA is usually in tandem‐repeat array (two or more copies of the same sequence that sit right next to one another).
It is not possible to say if they are originally on the same chromosome, as hit are on different scaffolds that can originally be on different, or on the same chromosome. 
A chromosome level assembly would be needed. Also, being rDNA characterized by repeated copies, mis-assemblies are likely, if using only short reads.










