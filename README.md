# ISMB Tutorial: 3D Genome Data Processing, Analysis, and Visualization Tutorial

Student Requirements: June 27th, 2017 
* Install Docker
  * https://www.docker.com/community-edition
* Install Miniconda
  * https://conda.io/miniconda.html  
* Windows users
  * Putty (for ssh)
  * https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html

## Files in this repository

* **Tutorial Part 1 (Introduction to HiC)**: ismb-workshop.pptx
* **Tutorial Part 2 (Data Processing and Visualization)**: index.html or [https://hms-dbmi.github.io/3d-genome-processing-tutorial/](https://hms-dbmi.github.io/3d-genome-processing-tutorial/)
* **Tutorial Part 3 (Data Analysis for Nuclear Compartmentalization)**: Ma-ISMB-tutorial-2017.pdf
  

## Agenda

**10:00 - 10:25 - Introduction and Overview**

* Who are we? [5 minutes]
* What we’ll cover [5 minutes]
* AWS accounts and ssh in [10 minutes]

**10:25 - 11:30 Hi-C Analysis**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Intro to Hi-C [Nezar]**

* 15 minutes:
  * Introduction
    * Protocol
    * What Hi-C sees
    * Levels of genome organization inferred from patterns
      * Checkering: compartments
      * Enriched squares without checkering: TADs
      * Corner peaks: loops
    * How Hi-C is processed
    * How features are assessed
  * Hi-C processing steps (informational) [Nezar]
    * Mapping
    * Filtering
    * Creating a list of contacts
    * Binning
    * Normalization
    * Feature analysis
      * TADs
      * Loops
    * QC

* 45min
  * Practical - processing pipeline [Soo]: 
    * Mapping
      * bwa mem -SP5M index fastq1 fastq2 | samtools view -bhS - > output.bam
      * samtools view output.bam | head -n 5
      * output: bam file
    * Filtering / sorting / Creating a list of contacts
      * pairsamtools
      * outputs: pairs, bam
      * Pairs / pairix (indexes the pairs file)
    * Binning
      * Cooler / cool
    * Normalization
      * Cooler / cool
    * map a small Hi-C dataset using distiller (https://github.com/mirnylab/distiller) and generate contact matrices using cooler (https://github.com/mirnylab/cooler)
  * Practical - Feature analysis [Nezar]: 
    * Jupyter notebook
    * Cis vs trans and scaling (contact probability vs genomic distance)
    * Feature analysis
      * Compartments, saddle plots
      * Insulation, TADs
      * Loops
      * Pileups


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**11:30 - 11:45 Coffee Break**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**11:45 - 12:45 Visualization**

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Existing tools for contact matrix exploration**

* 15 minutes [Nils]: 
  * 3D genome browser	
  * WashU epigenome browser
  * Juicebox
  * HiGlass 

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**Using HiGlass (http://higlass.io) to display contact maps [Peter]**

* 25 minutes: 
  * Overview of common operations such as adding tracks, removing tracks, adding views, removing view, linking views by zoom and location
  * Practical: 
    * Create an interactive version of a figure
* 20 minutes: Installing HiGlass
  * Overview of the HiGlass architecture and description of the infrastructure used to run it
  * Practical: 
    * Create a local HiGlass instance
    * Convert a contact map to multi-resolution format and import it
    * Convert a bigWig file to hitile format and import it
    * Open both files in the client and navigate to an interesting location

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**12:45 - 1:30 - Data Analysis for Nuclear Compartmentalization [Jian Ma]**

* Introduction
* DamID analysis
* Repli-seq analysis
* Data from emerging technologies



