---
ospool:
    path: software_examples/bioinformatics/tutorial-minimap2/README.md
---

# Minimap2 Tutorial 🧬  
**Mapping Long Reads at Scale using Minimap2, HTCondor, and the Open Science Pool**

This repository contains a Jupyter Notebook tutorial, `Minimap2_Tutorial.ipynb`, that walks through a complete **long-read sequencing read mapping workflow** using [Minimap2](https://github.com/lh3/minimap2).  
It demonstrates how to run mapping tasks locally and then scale them out to hundreds or thousands of jobs on distributed computing systems like the **Open Science Pool (OSPool)** using **HTCondor** and **Apptainer** containers.

---

## 📘 Tutorial Overview

You will learn how to:
- Map long Oxford Nanopore reads to a reference genome using **Minimap2**
- Break down large sequencing datasets into many smaller sub-jobs
- Submit hundreds to thousands of mapping jobs efficiently using **HTCondor**
- Manage software environments with **Apptainer**
- Use the **Open Science Data Federation (OSDF)** to stage and transfer data efficiently across distributed sites

This exercise uses realistic genomic data and highlights principles of **performance, reproducibility, and scalability** in bioinformatics workflows.

---

## 🧠 Background

**Minimap2** is a versatile sequence alignment program for mapping long and short sequencing reads to reference genomes.  
It supports ONT, PacBio, and Illumina data and can be run efficiently in parallel across many compute nodes.

The dataset used here includes simulated Oxford Nanopore reads from the **humpback whale (_Megaptera novaeangliae_)** genome, mapped against its reference assembly.

---

## 🧩 Tutorial Structure

The notebook includes the following sections:

1. **Setup & Environment**
   - Ensuring you are working in your `tutorial-minimap2` directory
   - Using provided setup scripts (`download_data.sh`) to fetch data and software

2. **Data Preparation**
   - Example data derived from the _Megaptera novaeangliae_ genome - [GenBank: GCA_041834305.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_041834305.1/)
   - Simulated ONT reads generated with [pbsims](https://github.com/yukiteruono/pbsim3)
   - FASTQ files are split into manageable subsets:
     ```bash
     split -l 40000 inputs/m_novaeangliae.fastq --additional-suffix=_m_novaeangliae.fastq subset_
     ```

3. **Software Environment**
   - Minimap2 is executed within an **Apptainer container** (`software/minimap2.def`)
   - Demonstrates how to use reproducible environments across nodes

4. **Parallelization Strategy**
   - Subset the full read dataset into smaller chunks for parallel processing
   - Lists all input subsets:
     ```bash
     ls inputs/subset_* | xargs -n 1 basename > listOfReads.txt
     ```
   - Generates an HTCondor submit file to launch mapping jobs in parallel

5. **Job Execution with HTCondor**
   - Overview of the HTCondor submission system
   - How job lists are created and managed
   - Monitoring job progress and viewing outputs

6. **OSDF Integration**
   - Discussion of distributed data access through the **Open Science Data Federation**
    - Efficient data staging and transfer across multiple sites

---

## 🧬 Requirements

- Access to OSPool via the [Guest Notebook Service]()
- Provided with the Tutorial:
  - Simulated ONT read datasets
  - Reference genome files
  - HTCondor submit files and scripts
  - Apptainer definition files for Minimap2
- Basic familiarity with command-line operations, Jupyter Notebooks, and bioinformatics concepts

---

## 📊 Expected Outputs

Each job produces:
- `.sam` alignment files
- Optional `.bam` and `.bai` outputs after conversion with Samtools
- Summary logs (`.out`, `.err`, `.log`) for each job

You’ll also learn to:
- Validate successful job completion
- Aggregate outputs
- Scale the workflow to larger datasets
- Leverage distributed computing effectively
- Run reproducible bioinformatics analyses on the OSPool

---

## 🤝 Acknowledgments

This tutorial is adapted for training in high-throughput computing and bioinformatics within the **Open Science Grid (OSG)** ecosystems.  
It draws from realistic genomics workflows and data from public resources such as **Data Carpentry** and **NCBI**.

---

## 🧾 License

This tutorial is released under the [MIT License](LICENSE).

---

## 📬 Contact

For questions or contributions:
- **Author:** Daniel Morales  
- **Institution:** Open Science Pool + Center for High Throughput Computing (CHTC), University of Wisconsin–Madison  
- **Email:** damorales4@wisc.edu  
- **Website:** [https://chtc.cs.wisc.edu](https://chtc.cs.wisc.edu)


