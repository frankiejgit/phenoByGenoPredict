# Phenotype-By-Genotype Prediction Tool


_Still in progress_

__Outline__
- About
- Installation
- Requirements
- Usage Instructions

## About

This tool was created by Francisco Gonzalez from the Jarquin Lab at the University of Florida, Department of Agronomy. Its purpose is to help estimate the expression of phenotypic traits by using sparse data collected by plant breeders. 

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

