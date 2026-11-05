

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

