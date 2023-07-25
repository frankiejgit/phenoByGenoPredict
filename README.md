# Phenotype-By-Genotype Prediction Tool

__Outline__
- About
- Installation
- Requirements
- Usage Instructions

## About

This tool was created by Francisco Gonzalez, graduate student at the University of Florida and member of the Jarquin Lab within the Agronomy Department. The purpose of this tools is to help plant breeders and researchers estimate the expression of phenotypic traits by providing environmental and genomic line data. 

Using cross validation techniques and BGLR, the tool can handle sparse datasets to provide predictions with reliable accuracy despite being fed limited data points. 

## Installation

First, you must verify that you have Docker and git installed in your computer

```
docker -v
```

```
git -v
```

If `Docker` is not installed, please follow the guidelines to install for your particular OS [here](docs.docker.com/engine/install). Similarly, if `git` is not installed in your device, follow this [guide](github.com/git-guides/install-git) to install it.


Once you have `Docker` and `git` successfully installed, follow these steps to run the platform.

1. In your terminal, clone this repository to your local computer

```
git clone https://github.com/frankiejgit/phenoByGenoPredict.git
```

2. Navigate to the local repository
```
cd phenoByGenoPredict/
```

3. Use `Docker` to build and configure the tool, wait for it to finish building completely
```
docker build --no-cache -t pheno_by_geno .
```

4. Once the build is done, you cna now access the tool via a UI developed using RShiny
```
docker run -p 3838:3838 pheno_by_geno
```

## Requirements

In your local computer, you should have access to a terminal to download the tool from the GitHub repository. Additionally, `docker` and `git` should be installed in your computer. 

## Usage Instructions

The tool's model accepts two data files in either CSV or RDA format. The first file should be the training data that contains the phenotypic traits of the crop varieties while the second file should have information of the molecular markers. There are some data formatting requirements for the first file:

- An environment ID/name column must be specified: This can be either in character or integer format
- The column with the data for the desired phenotypic trait should be specified
- A column containing the genotypic line or ID must be specified as well

For model tuning, you have the option of configuring the model types you want to use, the type of data processing you want to do prior to training the model with your data, and the model parameters.

The models are as follows:
- `E+L` = environment
- `E+L+G` = environment + genotypic line
- `E+L+G+GE` = environment + genotypic line + GxE interaction
