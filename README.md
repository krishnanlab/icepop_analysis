# icepop_analysis

This repository reproduces the results from the ICePop paper ([link](https://doi.org/10.64898/2026.04.01.715877)).  
For installation instructions and a tutorial for the ICePop package, see: https://github.com/krishnanlab/icepop

We assume the following directory structure when running the analysis:
```
data/
results/ (directory to save results generated from codes in src/)
src/
paper/
notebook/
run/ (working directory for running scripts)
```

All data required to reproduce the analyses are available on [Zenodo](https://doi.org/10.5281/zenodo.20146708).

We assume all data have been downloaded and extracted under `data/`.

## Repository Structure
- notebook/: Jupyter notebooks used to generate figures and tables in the paper  
- src/: Scripts for reproducing the results  
- paper/: Output figures and plots included in the manuscript  

## Documentation
- Detailed descriptions of scripts are available in [src/README.md](src/README.md)  
- Notebook-specific explanations are available in [notebook/README.md](notebook/README.md)