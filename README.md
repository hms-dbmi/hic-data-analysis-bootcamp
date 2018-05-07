# Hi-C Data Analysis Bootcamp

> Assessing, analyzing, and visualizing the 3D genome.

![Funky Colromaps](/teaser.jpg?raw=true "Some funky colormaps")

## Files in this repository

* **Tutorial Part 1 (Hi-C Protocol)**: ???
* **Tutorial Part 2 (From fastqs to contact matrices)**: [Introduction](https://github.com/hms-dbmi/3d-genome-processing-tutorial/blob/master/ISMB-Tutorial-Visualization-Intro.pdf)
* **Tutorial Part 3 (From contact matrices to biology)**: ???
* **Tutorial Part 4 (Hi-C Data Visualization - HiGlass)**: [https://hms-dbmi.github.io/hic-data-analysis-bootcamp/](https://hms-dbmi.github.io/hic-data-analysis-bootcamp/)
* **Tutorial Part 5 (Hi-C Data Visualization - HiPiler)**: hipiler.pdf

## Presenters

* Johan Gibcus, Research Instructor, Universy Massachusetts Medical School
* [Nezar Abdennur](http://nvictus.me/), PhD student, MIT
* [Soo Lee](https://compbio.hms.harvard.edu/people/soohyun-lee), Senior Bioinformatics Scientist, Harvard Medical School
* [Peter Kerpedjiev](http://emptypipes.org/about), Postdoctoral Research Fellow, Harvard Medical School
* [Fritz Lekschas](https://lekschas.de/) PhD Student, Harvard University
* [Leonid Mirny](http://mirnylab.mit.edu/) Professor, MIT

## Organizers

* Burak Alver, Scientific Project Manager, Harvard Medical School
* [Nils Gehlenborg](http://gehlenborglab.org/), Assistant Professor, Harvard Medical School
* [Peter Park](https://compbio.hms.harvard.edu/), Professor, Harvard Medical School

## Motivation

Due in large part to the explanatory power of chromosome organization in gene regulation, its association with disease and disorder as well as the unanswered questions regarding the mechanisms behind its maintenance and function, the 3D structure and function of the genome are becoming increasingly target of scientific scrutiny. With efforts such as the 4D Nucleome Project and ENCODE 4 already beginning to generate large amounts of data, the ability to analyze and visualize it will be a valuable asset to any computational biologist tasked with interpretation of experimental results.

## Objectives

The objectives of this tutorial are

* To introduce the theoretical concepts related to 3D genome data analysis
* To familiarize participants with the data types, analysis pipeline, and common tools for analysis and visualization of 3D genome data
* To provide a hands on experience in data analysis by walking through some common use cases of existing tools for data analysis and visualization.

## Goal

After the workshop participants should be able to obtain, process, analyze, and visualize 3D genome data on their own as well as to understand some of the logic, motivation and pitfalls associated with common operations such as matrix balancing and multi-resolution visualization.

## Intended audience and level

The subject matter and practical exercises presented in this tutorial will be accessible to a broad audience. Prior experience with next generation sequencing and the data it produces will be helpful for understanding the subsequent processing steps used to derive contact maps as well as some of the artifacts that can arise during data processing.

The material will be most useful to computational biologists and biologists working on genomics-related topics.

## Student Requirements:

* Install Docker
  * https://www.docker.com/community-edition
* Install Miniconda
  * https://conda.io/miniconda.html
* Windows users
  * Putty (for ssh)
  * https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html

## Agenda

**09:00 - 09:10** - Introduction and Overview (Peter Park and Burak Alver, Harvard)

**09:10 - 09:30** - Hi-C Protocol (Johan Gibcus, UMass)

**10:30 - 10:45** - _Break_

**10:45 - 12:15** - From fastqs to contact matrices (Soohyun Lee, Harvard)

**12:15 - 13:00** - _Lunch_

**13:00 - 14:00** - From contact matrices to biology (Nezar Abdennur, MIT)

**14:00 - 15:00** - Hi-C Data Visualization - HiGlass (Peter Kerpedjiev, Harvard)

**15:00 - 15:15** - _Break_

**15:15 - 16:00** - Hi-C Data Visualization - HiPiler (Fritz Lekschas, Harvard)

**16:00 - 17:00** - Keynote Speaker - Leonid Mirny, MIT


## Instructor Bios

### Johan Gibcus

Johan Gibcus is a Research Instructor at the University of Massachussetts Medical School. He has not only used but also refined the Hi-C protocol to answer important biological questions about chromosome organization and replication. Web: [http://www.dekkerlab.org/](http://www.dekkerlab.org/)

### Soo Lee

Soo Lee is a Senior Bioinformatics Scientist in the Department of Biomedical Informatics at Harvard Medical School. She is creating cloud-based pipelines for Hi-C and other genomic data and developing infrastructure for automation of such pipelines as part of the 4D Nucleome Data Coordination and Integration Center. Web: [compbio.hms.harvard.edu/people/soohyun-lee](https://compbio.hms.harvard.edu/people/soohyun-lee)

### Nezar Abdennur

Nezar Abdennur is a PhD candidate in Computational and Systems Biology at MIT. His research focuses on the determinants of 3D genome organization and the development of tools for dealing with large Hi-C datasets. Twitter: [@nv1ctus](https://twitter.com/nv1ctus) Web: [nvictus.me](http://nvictus.me)

### Peter Kerpedjiev

Peter Kerpedjiev is a postdoctoral researcher working on creating tools (such as HiGlass) for visualizing large genomic data sets. Twitter: [@pkerpedjiev](https://twitter.com/pkerpedjiev) Web: [emptypipes.org](http://emptypipes.org)

### Fritz Lekschas

Fritz is a PhD student working on biomedical information visualization with focus on large multiscale genomic data sets. Twitter: [@flekschas](https://twitter.com/flekschas) Web: [lekschas.de](https://lekschas.de)

### Leonid Mirny

Leonid Mirny is a professor at MIT's Institute for Medical Engineering & Science. His lab studies the three dimensional organization of chromosomes using a combination of computational analysis and simulation. Twitter: [@leonidmirny](https://twitter.com/leonidmirny) Web: [mirnylab.mit.edu](http://mirnylab.mit.edu/)


## Resources

**Software:**

* [bwa](https://github.com/lh3/bwa) and [SAM spec](https://samtools.github.io/hts-specs/SAMv1.pdf)
* [pairsamtools](https://github.com/mirnylab/pairsamtools)
* [pairix](https://github.com/4dn-dcic/pairix)
* [cooler](https://github.com/mirnylab/cooler) and [docs](http://cooler.readthedocs.io/en/latest/)
* [HiGlass](http://cooler.readthedocs.io/en/latest/) and [wiki](https://github.com/hms-dbmi/higlass/wiki)
* [HiPiler](http://hipiler.lekschas.de) and [source code](https://github.com/flekschas/hipiler) and [docs](https://github.com/flekschas/hipiler/wiki)


**Bioinformatics workflow managers:**

* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [nextflow](https://www.nextflow.io/)


**Package and environment management:**

* [conda](https://conda.io/miniconda.html)
* [bioconda](https://bioconda.github.io/)


**Papers:**

* Imakaev, Maxim, et al. "Iterative correction of Hi-C data reveals hallmarks of chromosome organization." Nature methods 9.10 (2012): 999-1003. doi:[10.1038/nmeth.2148](https://doi.org/10.1038/nmeth.2148)
* Lajoie, Bryan R., Job Dekker, and Noam Kaplan. "The Hitchhiker’s guide to Hi-C analysis: practical guidelines." Methods 72 (2015): 65-75. doi:[10.1016/j.ymeth.2014.10.031](https://doi.org/10.1016/j.ymeth.2014.10.031)
* Kerpedjiev, Peter, et al. "HiGlass: Web-based Visual Comparison And Exploration Of Genome Interaction Maps" bioRxiv. doi:[10.1101/121889](https://doi.org/10.1101/121889)
* Lekschas, Fritz et al. "HiPiler: Visual Exploration Of Large Genome Interaction Matrices With Interactive Small Multiples" IEEE Transactions on Visualization and Computer Graphics, 24(1), 522-531. [doi:10.1109/TVCG.2017.2745978](https://doi.org/10.1109/TVCG.2017.2745978)
