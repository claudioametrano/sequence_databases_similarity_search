### **TASK 1**
Tuber[Organism]  ->  118562 seqs

Tuber[Organism] AND (ribosomal RNA) -> 5746 seqs

Tuber[Organism] AND (ribosomal RNA) AND ( "500"[SLEN] : "5000"[SLEN] )   -> 5087 seqs

Tuber[Organism] AND (ribosomal RNA) AND ( "500"[SLEN] : "5000"[SLEN] ) NOT environmental  -> 4773 seqs

(Sequence number generated on Apr 15 2025)
Download the sequences ordered by length in fasta (send to:)

### **TASK 2**
```bash
$ mkdir results/task2
```
- Select blastn
- load or copy/paste into the window sequence.fasta (cat or less can be useful)
- Enter in the Organism window the scientific name (or taxa id from NCBI Taxonomy) to limit the results to that species, try a blast, if you get no results try to progressively relax blast parameters . Click on **Algorithm parameters** and progressively relax "word size" and "match/mismatch" 
- Select the best hit. If you have problem finding a trustable sequence for yeast you may think is because is too distantly related to animals, well it is quite distant (relax BLAST algorithm!), but this is not the case. The reality is COX1 is not very often sequences for yeasts, but you get hist from complete mitochondrial yeasts's genomes. Go to "Alignments" and download as fasta only the portion of the mitochondrial genomes that aligns to your query.
-  For the other species just download the best hit in fasta: tick on the sequence -> Download -> FASTA
### **Optional**
Concatenate fastas in a multi-fasta file
```bash
 $ cat ./results/task2/*.fasta > results/task2/cat_besthit_COX1.fas
```

Align the multifasta file
```bash
$ singularity pull https://depot.galaxyproject.org/singularity/mafft:7.525--h031d066_1

$ singularity exec mafft\:7.525--h031d066_1 mafft --localpair --maxiterate 1000 --thread 4 results/task2/cat_besthit_COX1.fas > results/task2/cat_besthit_COX1_mafft_aligned.fas
```

Check if the alignment makes sense (aliview, jalview, MEGA... ) and delete unaligned blocks, in this case with TrimAl heuristic method
```bash
$ singularity pull https://depot.galaxyproject.org/singularity/trimal:1.5

$ singularity exec trimal:1.5 trimal -in results/task2/cat_besthit_COX1_mafft_aligned.fas -out results/task2/cat_besthit_COX1_mafft_aligned_trimmed.fas -automated1 -keepheader
```

Phylogenetic inference
```bash
$ sed -i 's/ /_/g' results/task2/cat_besthit_COX1_mafft_aligned_trimmed.fas #substitute spaces with _ or name will be truncated at first space

$ singularity pull https://depot.galaxyproject.org/singularity/iqtree:2.3.6--h503566f_1
$ singularity exec iqtree\:2.3.6--h503566f_1 iqtree2 -s results/task2/cat_besthit_COX1_mafft_aligned_trimmed.fas -m MFP -bb 1000  -T AUTO
```

Load the .treefile in ITOL to visualize the phylogenetic tree
```bash
$ cat ./results/task2/cat_besthit_COX1_mafft_aligned_trimmed.fas.treefile 
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
$ singularity exec blast:2.16.0--h66d330f_4 blastn -num_threads 4 -query ./data/sequence.fasta -db ./results/NCBI_databases/SSU_eukaryote_rRNA -out ./results/task3/blastout_sequence_vs_SSUeuk_rRNA.txt

```
and so on...
or you can loop through databases!
```bash
# Names of the databases in a list
databases=("ITS_eukaryote_sequences" "LSU_eukaryote_rRNA" "SSU_eukaryote_rRNA")

# Loop over them
for db in "${databases[@]}"; do
    #blastn command
    singularity exec blast:2.16.0--h66d330f_4 blastn \
        -query ./data/sequence.fasta \
        -db "./results/NCBI_databases/${db}" \
        -num_threads 4 \
        -out "./results/task3/blastout_sequence_vs_${db}.txt"
done
```
check results:
```bash 
$ cat results/task3/blastout_sequence_vs_*
```
no hits!

- spacer are very variable (purifying selection does not select negatively substitutions), relaxed parameter are necessary to find hits in these smaller rDNA databases we are using (exactly as we already did on the BLAST webserver):
```bash
# Names of the databases in a list
databases=("ITS_eukaryote_sequences" "LSU_eukaryote_rRNA" "SSU_eukaryote_rRNA")

# Loop over them
for db in "${databases[@]}"; do
    # blastn command
    singularity exec blast:2.16.0--h66d330f_4 blastn \
        -query ./data/sequence.fasta \
        -db "./results/NCBI_databases/${db}" \
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
Parse the output  using grep:
``` bash
$ grep  -A 2 --with-filename "significant alignments" results/task3/blastout_relaxed_sequence_vs_*
```
Or also...
```bash
$ find ./results/task3 -name blastout_relaxed_sequence_vs_* -exec  grep -A 2 --with-filename "significant alignments" {} \;
```
good e-value from ITS eukaryote. It looks like it is  an ITS sequence from Fungi, putatively *Tuber cryptobrumale*

### **TASK 4**
Show blast results:
```bash
cat results/task4/blastout_sequence_vs_GCA_000151645.1_tabfmt.txt
```
A total of 5 hits are present on three scaffolds. 
Searching the headers of the hits:
```
$ grep  "NW_003299036.1" results/task4/GCF_000151645.1_ASM15164v1_genomic.fna
$ grep  "NW_003298972.1" results/task4/GCF_000151645.1_ASM15164v1_genomic.fna
$ grep  "NW_003298908.1" results/task4/GCF_000151645.1_ASM15164v1_genomic.fna 
```
One hit on scaffold 24 , two on scaffolds 298 and 355, hits on smaller scaffolds possibly presenting two repeated rDNA sequences (two ITS1 hits ~8000 bp apart)
Multiple hits are common, as rDNA is usually in tandem‐repeat array (two or more copies of the same sequence that sit right next to one another).
It is not possible to say if they are originally on the same chromosome, as hit are on different scaffolds that can originally be on different, or on the same chromosome. 
A chromosome level assembly would be needed. Also, being rDNA characterized by repeated copies, mis-assemblies are likely, if  only short reads assembly is available.

#### **TASK 5**
The lowest value for the match state is 0.39301, in correspondence to a M (Methionine).
Yes, it makes sense, start codon produce an M.








