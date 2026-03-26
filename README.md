# icepop_analysis

This repository reproduces the results from the ICePop paper.  
For installation instructions and a tutorial for the ICePop package, see: https://github.com/krishnanlab/icepop_public/tree/main  

We assume the following directory structure when running the analysis:
```
data/  
results/  
src/  
notebook/  
run/ (working directory for running scripts)
```

All data required to reproduce the analyses are available on [Zenodo](https://doi.org/10.5281/zenodo.19238928).

We assume all data have been downloaded and extracted under `data/`.

The `magmaz` directory should be placed under `data/TM_FACS/`.

## Repository Structure
- notebook/: Jupyter notebooks used to generate figures and tables in the paper  
- src/: Scripts for reproducing the results  
- paper/: Output figures and plots included in the manuscript  

## Documentation
- Detailed descriptions of scripts are available in [src/README.md](src/README.md)  
- Notebook-specific explanations are available in [notebook/README.md](notebook/README.md)