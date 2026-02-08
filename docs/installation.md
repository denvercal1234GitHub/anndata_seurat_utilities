## Installation

Current version: 0.1.1


### Install official release from PyPI

```bash

pip install anndata-seurat-utilities
## Quick test; if printing version then succeeded
python -c "import anndata_seurat_utilities; print(anndata_seurat_utilities.__version__)"

```


### Development install

```bash
git clone https://github.com/denvercal1234GitHub/anndata_seurat_utilities
cd anndata_seurat_utilities
pip install -e .
pytest -q
```


## Requirements
anndata >= 0.9
pandas
numpy
scipy
