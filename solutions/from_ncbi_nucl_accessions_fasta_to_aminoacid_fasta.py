#!/usr/bin/env python

import sys
from Bio import Entrez, SeqIO

def fetch_protein_sequences(accession):
    """
    Fetches a GenBank record for 'accession' from NCBI,
    extracts all protein translations (CDS features),
    and returns a list of protein sequences (strings).
    """
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    gb_record = SeqIO.read(handle, "genbank")
    handle.close()

    protein_seqs = []
    for feature in gb_record.features:
        # Only look at CDS features with a translation qualifier
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            protein_seqs.append(feature.qualifiers["translation"][0])

    return protein_seqs

def main(input_fasta, output_fasta):
    # Set your email to something valid (NCBI requirement)
    Entrez.email = "your_email@example.com"

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        # Parse the input FASTA
        for record in SeqIO.parse(infile, "fasta"):
            accession = record.id  # e.g., "XM_062772314.1"
            
            # Fetch and parse protein translations
            proteins = fetch_protein_sequences(accession)
            
            # Write each translated CDS to the output FASTA
            for idx, prot_seq in enumerate(proteins, start=1):
                # Example header: >XM_062772314.1_CDS1
                outfile.write(f">{accession}_CDS{idx}\n{prot_seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_fasta> <output_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    main(input_fasta, output_fasta)
