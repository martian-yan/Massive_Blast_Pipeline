# Massive Blast Pipeline

This is a pipeline based on [`SnakeMake`](https://snakemake.readthedocs.io/en/stable/). The subject of this pipeline is to search a protein sequence in thousands of bacterial genomes/proteomes. Making use of the parallel processing property of `SnakeMake`, this pipeline can run >8000 `diamond` search in ~30 minutes with 16 cores.

## Usage

- Dependencies
    - SnakeMake
    - Python 3
        - biopython
        - pandas
    - [`Diamond`](http://www.diamondsearch.org/index.php)
    - Clustalw

    Recommend install dependencies with `Conda`
- Download this pipeline

    ```sh
    git clone git@github.com:martian-yan/Massive_Blast_Pipeline.git/
    ```

- Config the inputs: edit `config.yaml`, recommend use a customized config file for each task.

    ```yaml
    # Working directory: use absolute path, or relative path based from "Snakefile" folder.
    workdir: ./test

    # input files, use absolute path, or relative path based from "workdir"
    query: test.tsv # a table give the sequences ID and path
    reference: SopF.fa # a protein sequence

    # querytype: protein or nucl
    querytype: nucl

    # translate_table if querytype is nucl
    # it's a paramater of Bio.Seq.translate module 
    translate_table: Bacterial
    ```

- Edit the "query" table, use a "tab" separated 2-column table, first column is the ID of each sequence, second column is path, for example in my "test.tsv":
    
    ```plain
    FD01851500      /pub28/yan/04_Genome_collections/10K_data/V3_merged_analysis/8_fasta_to_send/fasta/FD01851500.fa
    FD01851501      /pub28/yan/04_Genome_collections/10K_data/V3_merged_analysis/8_fasta_to_send/fasta/FD01851501.fa
    ...
    ```

    **Notice**: when using absolute path, start from "root folder" (/), avoid using "home folder" (\~) cause `SnakeMake` can't recognise "\~".

- Run snakemake: go to the folder with `Snakefile`, run

    ```bash
    # simplest with 1 thread
    snakemake
    ```

    or

    ```bash
    # Use 16 threads and customized config file
    snakemake --cores 16 --configfile /path/to/your/config.yaml
    ```

## Result

The final result is in the `result` folder, including:

  - `protein_alleles.aln`: all the protein alleles,
  - `protein_alignment.clustal`: an `clustalw` alignment of all the protein alleles,
  - `protein_top10.aln`: top 10 protein alleles with most appearance,
  - `protein_alignment_top10.clustal`: an `clustalw` alignment of top 10 protein alleles,
  - `summary.tsv`: a summary of the protein type of each sequence.

## Todo

The planned updates in the future:

  - [ ] Use `--conda` feature of `SnakeMake` to automatically install dependencies,
  - [ ] Support get all queries from a single folder rather than config the path one by one,
  - [ ] Support use [`Last`](http://last.cbrc.jp/) rather than `Diamond`, `Last` is another sequence comparison program which is faster and can recognise frameshift mutations. 