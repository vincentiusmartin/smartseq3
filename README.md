

# Smart-seq pipeline

This repository contains alignment and downstream analysis pipeline for Smart-seq
RNA seq.

<details><summary><b>Table of content</b></summary>

- [Directory structure](#directory)
- [Running the pipeline](#running-the-pipeline)
    - [Step-1: Generate the config file](#step1)
    - [Step-2: Running Nextflow pipeline](#step2)
    - [Step-3: Downstream analysis](#step3)
    
</details>

## Directory structure <a name="directory"></a>
    .
    ├── samplesheets       # example of samplesheets
    ├── modules            # script for nextflow module
    ├── analysis           # Currently contain template for downstream analysis
    ├── bin                # contains executable for the nextflow pipeline
    └── README.md

## Running the pipeline <a name="running-the-pipeline"></a>

### Step-1: Generate the config file <a name="step1"></a>

The pipeline accepts a samplesheet containing paths to the input files. An example samplesheet is shown below. 
An example script for generating the samplesheet is provided in `createsamplesheet.R`.

| sample | fastq1 | fastq2 |
| --- | --- | --- |
| cell1 | /path/to/cell1_L001_R1.fastq.gz | /path/to/cell1_L001_R2.fastq.gz |
| cell1 | /path/to/cell1_L002_R1.fastq.gz | /path/to/cell1_L002_R2.fastq.gz |
| cell2 | /path/to/cell2_L001_R1.fastq.gz | /path/to/cell2_L001_R2.fastq.gz |
| cell2 | /path/to/cell2_L002_R1.fastq.gz | /path/to/cell2_L002_R2.fastq.gz |


### Step-2: Running Smart-seq pipeline <a name="step2"></a>

1. Edit `nextflow.config` file. Parameter definition for now is defined in the config file
to simplify running the pipeline.
2. Run the nextflow pipeline: `nextflow run main.nf`

### Step-3: Interpreting results <a name="step3"></a>

An example of script to perform downstream analysis is in the `analysis` folder.
Example of resulting analysis is also provided in the same folder.



