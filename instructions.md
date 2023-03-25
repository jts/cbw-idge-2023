---
layout: tutorial_page
title: CBW Infectious Disease Epidemiology 2023
header1: Workshop Pages for Students
header2: Informatics for High-throughput Sequencing Data Analysis 2019 Module 6 Lab
home: https://bioinformaticsdotca.github.io/
description: CBW IDE 2023 Module 3 - Virus Sequencing
author: Jared Simpson
modified: March 21, 2023
---

# Virus Short-Read Genome Assembly and Variant Analysis

by [Jared Simpson](https://simpsonlab.github.io)

## Introduction

In this lab we will perform reference based assembly of Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2) short-read sequencing data. You will be guided through the genome assembly and subsequent analysis of output visualizations for assessing variants and lineages. At the end of the lab you will know:

1. How to run a short-read assembly using Illumina data
2. How to assess the quality of an assembly
3. How to assess single nucleotide polymorphisms (SNPs) for SARS-CoV-2 variant calls
4. How to view the evolutionary changes of SARS-CoV-2 and interpret differences between lineages belonging to variants of concern (i.e., Alpha vs. Delta variants)

## Data Sets

In this lab we will use subset of data from the COVID-19 Genomics UK COnsortium (COG-UK) to demonstrate genome assembly of SARS-CoV-2. This data set consists of short-read sequencing reads collected as part of the COVID-19 Pandemic response and will be a mix of different variants of concern (VOCs) up to the end of summer 2021. This would include predominantly Alpha (also known as PANGO lineage B.1.1.7; first designated December 2020) and Delta (also known as PANGO lineage B.1.617.2; first designated in May 2021) variants. We have provided accession IDs to Illumina short-reads and the assemblies will be done using the SARS-CoV-2 Illumina GeNome Assembly Line (`SIGNAL`) with additional quality control and assessment using `ncov-tools`. We will perform some additional steps to downsample our reads to speed up the overall run time; however, running the overall pipelines use the exact same set of commands and the results you will obtain are comparable to assembling using the full dataset.

## Data Download and Preparation

First, lets create and move to a directory that we'll use to organize our raw data:

```
mkdir -p Module4
cd Module4
```

From within your `Module4` directory, you can view the list of accessions that will be downloaded, you can view the following file (found within `etc/accessions.txt`):

```
more ../etc/accessions.txt
```

Note that the `..` means to go one directory above, taking us back to our starting position. You can also use `less` to view the file and use the space bar to scroll down, but use the `q` key to exit the program when you cannot scroll any further.

To download our seqeucning data, we'll use a provided script with our list of accession IDs as input:

```
bash ../bin/get_raw_data.sh ../etc/accessions.txt
```

If you run `ls` you should now be able to see a directory has been created called `raw_data`. If you run `ls raw_data` you should be able to see a series of FASTQ files, including "R1" and matching "R2" files.

We are going to add a step to downsample our reads prior to running through the assembly pipeline. This will results in fewer overall reads provided as input; however, the reads put forward will be sufficient for our analysis. This assumes that the `covid-19-signal` directory is located one level above, where we first started this lab practical.

```
bash ../bin/downsample.sh -a ../etc/accessions.txt -r ../covid-19-signal
```

This script will take our list of accessions (-a) and the location of the SIGNAL repository (-r) and produce another directory called `raw_reads_downsampled`, which you can verify with `ls`. Also, if you run `ls raw_reads_downsampled`, you should be able to see a series of FASTQ files, including "R1" and matching "R2" files.

Navigate to the `covid-19-signal` directory where our results will be located (subject to change given filesystem...) and activate the conda environment required to run both `SIGNAL` and `ncov-tools`:

```
cd ../covid-19-signal
conda activate signalcovtools
```

In order to run SIGNAL, we first need to prepare two files: a configuration file, where all of our assembly parameters will be assigned, and a sample table, which will list the indivdual samples and the location of corresponding R1 and R2 FASTQs. Remember that our sequencing data is located one directory level up (i.e., `../Module4/raw_reads_downsampled`). Generating the required files can all be done using doing the following:

```
python signalexe.py --directory ../Module4/raw_reads_downsampled --config-only
```

If you run `ls` you should see `raw_reads_downsampled_config.yaml` and `raw_reads_downsamples_sample_table.csv` files have been created. You can use `more` or `less` to examine the input files.

(Need confirmation of `data` directory for SIGNAL to work)

## Reference-based assembly using SIGNAL

Using our configuatrion file as input, we can begin our assembly of SARS-CoV-2 sequencing reads. Run the following:

```
python signalexe.py --configfile raw_reads_downsampled_config.yaml --cores 4 all postprocess
```

We can now start assessing the quality of our assembly. We typically measure the quality of an assembly using three factors:

- Contiguity: Long contigs are better than short contigs as long contigs give more information about the structure of the genome (for example, the order of genes)
- Completeness: Most of the genome should be assembled into contigs with few regions missing from the assembly (remember for this exercise we are only assembling one megabase of the genome, not the entire genome)
- Accuracy: The assembly should have few large-scale _misassemblies_ and _consensus errors_ (mismatches or insertions/deletions)

## Additional quality control and assessments using ncov-tools

Similarly to how our configuatrion file was used as input, we can similarly run `ncov-tools` through SIGNAL. Run the following:

```
python signalexe.py --configfile raw_reads_downsampled_config.yaml --cores 4 ncov_tools
```

## Interpretation of the data
