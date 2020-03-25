import pandas as pd
import numpy as np
import os
import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import translate, complement, reverse_complement
from Bio import SeqIO

## load config and input files
configfile:"config.yaml"
workdir: config['workdir']

df = pd.read_table(config["query"], header=None, names=["ID", "path"])
id_list = df["ID"].to_list()
path_list = df["path"].to_list()
query_dict = dict(zip(id_list, path_list))
#print(query_dict)

## prepare queries
os.makedirs("database", exist_ok=True)
for id in id_list:
    link="database/{id}.fa".format(id=id)
    if not os.path.exists(link):
        shell("ln -s {target} {link_name}".format(target=query_dict[id], link_name=link))

## target rules
rule all:
    input:
        "result/protein_alignment.clustal",
        "result/protein_alignment_top10.clustal"

### makeblastdb
rule makedb:
    input:
        protein=config["db_file"]
    output:
        "db/db.dmnd"
    params:
        db="db/db"
    shell:
        "diamond makedb --in {input.protein} --db {params.db}"

### run blast
rule blast:
    input:
        db="db/db.dmnd",
        query="database/{sample}.fa"
    output:
        "blast/{sample}.tsv"
    threads: 1
    run:
        if config["querytype"] == "nucl":
            shell("diamond blastx --threads {threads} --db {input.db} --query {input.query} --out {output} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen qseq --header")
        
        #for proteomic, use blastp, get full_qseq instead of qseq.
        elif config["querytype"] == "protein":
            shell("diamond blastp --threads {threads} --db {input.db} --query {input.query} --out {output} --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen full_qseq --header")      

### extract seq
rule extract_seq:
    input:
        tsv="blast/{sample}.tsv"
    output:
        fa="seq/{sample}.fa"
    run:
        blast_out = pd.read_table(input.tsv, header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'slen', 'qlen', 'qseq'], comment="#")
        
        # If there's no blast result, create an empty .fa file
        if blast_out.dropna().empty:
            with open(output.fa, 'w') as fp:
                pass

        # iterate over each blast result
        sequences_to_save = []
        for index, data in blast_out.iterrows():
            id = data['sseqid']
            description = "{file}|contig: {contig}|{start}..{end}".format(file=wildcards.sample, contig=data['qseqid'], start=data['qstart'], end=data['qend'])


            if config["querytype"] == "protein":
                seq = Seq(data['qseq'], IUPAC.extended_protein)

            elif config["querytype"] == "nucl":
                seq_dna = Seq(data['qseq'], IUPAC.extended_dna)

                if data['qstart'] > data['qend']:
                    seq_dna = seq_dna.reverse_complement()

                seq=seq_dna.translate(table=config["translate_table"])
            
            sequence = SeqRecord(seq, id=id, description=description)
            sequences_to_save.append(sequence)

        if sequences_to_save != []:
            SeqIO.write(sequences_to_save, output.fa, "fasta")

### summarize

rule summarize:
    input:
        expand("seq/{sample}.fa", sample=id_list)
    output:
        aln="result/protein_alleles.aln",
        tsv="result/summary.tsv"
    run:
        aln_to_save = []
        summary = []
        for fasta_file in input:
            file_name = re.search("seq/(.+).fa", fasta_file).group(1)
            
            protein_types = []

            # if the fasta file is empty
            if os.path.getsize(fasta_file) == 0:
                protein_types.append('None')

            for seq_record in SeqIO.parse(fasta_file, "fasta"):
                seq_string = str(seq_record.seq)
                seq_is_new = True

                # find whether seq is in aln
                for i in range(len(aln_to_save)):
                    protein_string = str(aln_to_save[i].seq)
                    protein_id = aln_to_save[i].id
                    if seq_string == protein_string:
                        seq_is_new = False
                        protein_types.append(protein_id)
                        break

                    elif seq_string in protein_string:
                        seq_is_new = False
                        start = protein_string.find(seq_string)
                        end = start + len(seq_string)
                        protein_types.append("{0}({1}..{2})".format(protein_id, start, end))
                        break
                    
                if seq_is_new:
                    seq_id = seq_record.id+"_{}".format(len(aln_to_save)+1)
                    seq_record.description = seq_record.description.split(' ',1)[1]
                    seq_record.id = seq_id
                    aln_to_save.append(seq_record)
                    protein_types.append(seq_id)
                    
            summary.append("{0}\t{1}\n".format(file_name, ",".join(protein_types)))

        
        SeqIO.write(aln_to_save, output.aln, "fasta")
        with open(output.tsv, 'w') as f:
            f.writelines(summary)
                
### alignment
rule alignment:
    input:
        "result/protein_alleles.aln"
    output:
        "result/protein_alignment.clustal"
    shell:
        "clustalw -INFILE={input} -TYPE=protein -OUTFILE={output}"
    
### alignment top 10
rule get_top10:
    input:
        ref=config['db_file'],
        aln="result/protein_alleles.aln",
        tsv="result/summary.tsv"
    output:
        aln="result/protein_top10.aln"
    run:
        summary = pd.read_table(input.tsv, header=None, names=['id', 'protein_type'])
        values,counts = np.unique(summary['protein_type'], return_counts=True)
        values_dict = dict(zip(values, counts))
        values_sort_by_counts = [k for k,v in sorted(values_dict.items(), key=lambda item: item[1], reverse=True)]
        top10 = values_sort_by_counts[0:10]
        
        if 'None' in top10:
            top10.remove('None')
        
        top10_aln = []
        with open(input.ref, 'r') as ref:
            for seq in SeqIO.parse(ref, "fasta"):
                seq.id = "reference"
                top10_aln.append(seq)

        with open(input.aln, 'r') as fi:
            for seq in SeqIO.parse(fi, "fasta"):
                for id in top10:
                    if seq.id == id:
                        top10_aln.append(seq)

        SeqIO.write(top10_aln, output.aln, "fasta")

rule align_top10:
    input:
        "result/protein_top10.aln"
    output:
        "result/protein_alignment_top10.clustal"
    shell:
        "clustalw -INFILE={input} -TYPE=protein -OUTFILE={output}"

        


