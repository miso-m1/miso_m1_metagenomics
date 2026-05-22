# M1 project

# Statistical Analysis of Shotgun Metagenomic Data

This repository contains a complete statistical analysis pipeline for shotgun metagenomic datasets of human gut microbiota, carried out as part of the Master’s program MISO in Bioinformatics at Université de Lille.

The project is based on publicly available human metagenomics datasets from Kaggle Human Metagenomics Dataset.

Our work combines statistical analysis, alpha and beta diversity exploration, and network analysis through an interactive blog built with Quarto.

# Project Overview

The objective of this project was to perform a complete analysis of metagenomic data from intestinal microbiomes. The workflow includes:

Data import, cleaning, and preprocessing,
Alpha-diversity analysis (within-sample diversity), statistical testing using ANOVA,
Beta-diversity analysis (between-sample diversity), statistical testing using PERMANOVA, 
Construction and interpretation of co-occurrence networks.

The analyses were conducted using R and various bioinformatics/statistical packages.

#Blog Website

The project is presented as an interactive blog available here : https://miso-m1.github.io/miso_m1_metagenomics/ .

To properly understand the project workflow, we strongly recommend beginning with the Getting Started section of the blog.

The posts are organized to follow the logical progression of a metagenomic analysis workflow. For the best reading experience, we recommend consulting them in the following order:

### 0. Introduction
This post presents the context, the objectives of the study, and the structure of the project.

### 1. Data Preparation
This section explains how the dataset was imported, cleaned and transformed before the study.

### 2.1 Alpha Diversity
This post explores within-sample diversity using diversity indices and visualizations.

### 2.2 Beta Diversity
This section investigates differences between samples and groups using distance metrics and ordination methods.

### 3. Metagenomic Network Analysis
The network analysis focuses on relationships between microbial taxa through co-occurrence network approaches.

### Additional sections include:
Several resources are available through the navigation bar at the top of the website:

#### Navigation Guide: How to navigate the blog.
Presentation: slides summarizing the project and main results.

#### Sources: references and resources used throughout the analyses.

#### About: information about the contributors and the context of the project.

#### Social Links: direct access to the contributors’ e-mail addresses, GitHub and LinkedIn profiles.

#Repository Structure

The repository is structured like that :

miso_m1_metagenomics/
│
├── Blog_template/     # Quarto website structure and templates
├── data/              # Datasets
├── docs/              # Rendered website files for GitHub Pages
├── posts/             # Blog posts
├── img/               # Images and figures used in the blog
├── _quarto.yml        # Quarto website configuration
├── resources/         # The different sources
└── README.md          # README.md

# Contributors

This project was carried out collaboratively, drawing on skills in biological data analysis, statistics, and bioinformatics.

The contributors are Ines Bakli, Paul Lemonnier, and Alden Sneath. We would also like to thank our referent, Clément Poupelin, for his guidance throughout this project.

