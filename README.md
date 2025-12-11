# anndata_Seurat_Utilities

**Quick install**

```bash
pip install anndata-seurat-utilities
```


Utilities to prepare **AnnData (.h5ad)** files for **Seurat / SeuratDisk** conversion and safe exchange between **Python** and **R** single-cell analysis workflows.

<!-- Uncomment these once CI & PyPI are active -->
<!-- 
[![PyPI](https://img.shields.io/pypi/v/anndata-seurat-utilities.svg)](https://pypi.org/project/anndata-seurat-utilities/)
[![CI](https://github.com/denvercal1234GitHub/anndata_Seurat_Utilities/actions/workflows/ci.yml/badge.svg)](https://github.com/denvercal1234GitHub/anndata_Seurat_Utilities/actions)
-->

---

## What this package does?

This library provides the robust function:

`prepare_adata_for_seurat_drop_empty_v3()`

which:

- **Drops empty `obs` / `var` columns** (with provenance stored in `.uns`)
- **Sanitizes column names** (avoids `/`, spaces, commas, semicolons, etc.)  
  → prevents accidental HDF5 path splitting (common SeuratDisk failure mode)
- **Converts pandas extension dtypes**  
  (`boolean`, `Int64`, `string`, `Categorical`, masked arrays)  
  → into safe atomic dtypes so HDF5 does *not* create `{mask, values}` groups
- **Optionally drops heavy layers** like `counts_mouse`
- **Ensures all metadata is stored in a SeuratDisk-compatible form**
- **Writes `.h5ad` cleanly** using `convert_strings_to_categoricals=False` when possible

The output h5ad can be used as input to Convert and LoadH5Seurat. See below for an example. 

## Installation

Current version: 0.1.1


### Install official release from PyPI

```bash

pip install anndata-seurat-utilities
## Quick test; if printing version then succeeded
python -c "import anndata_seurat_utils; print(anndata_seurat_utils.__version__)"

```


### Development install

```bash
git clone https://github.com/denvercal1234GitHub/anndata_Seurat_Utilities
cd anndata_Seurat_Utilities
pip install -e .
pytest -q
```


## Requirements
anndata >= 0.9
pandas
numpy
scipy




## Minimal usage example 

```python
## In Python 
from anndata_seurat_utils.prepare_for_seurat import (
    prepare_adata_for_seurat_drop_empty_v3
)
import anndata as ad

# Load an AnnData object (raw or processed)
adata = ad.read_h5ad("input.h5ad")

# Clean, sanitize, and write to a Seurat-compatible h5ad
## Note that if keep_layers is not None, it treats keep_layers as a whitelist and removes every layer whose name is not in keep_layers. drop_layers is ignored in this case.
# Note that if keep_layers is None, it instead looks at drop_layers and deletes every layer listed in drop_layers that exists in the AnnData. If drop_layers is None or empty, nothing is removed.
prepare_adata_for_seurat_drop_empty_v3(
    adata=adata,
    out_h5ad_path="/Users/......./cleaned_for_seurat.h5ad",
    drop_layers=['nonSCVI_normalize','counts'],   # or None or empty will not remove anything when keep_layers=None
    keep_layers=None, # if None will delete those defined in drop_layers 
    convert_strings_to_categoricals_on_write=False, # recommended for SeuratDisk
    sanitize_var_names=True,
    sanitize_obs_names=True,
    drop_empty_obs_and_var=True,
    verbose=True,
    keep_obsm=["X_pca", "X_umap"], # if None will remove all obsm; if 'all' will keep all 
    keep_obsp=["connectivities"], # if None will remove all obsp; if 'all' will keep all 
    keep_obs_cols=["sample_id", "nCount_RNA", "nFeature_RNA"], # if None will remove all obsp; if 'all' will keep all 
    keep_uns=["_scvi_uuid"], # if None will remove all obsp; if 'all' will keep all 
)

print("Wrote cleaned Seurat-ready file: cleaned_for_seurat.h5ad")

## Now proceed in R with SeuratDisk::Convert and LoadH5Seurat

```

The resulting `.h5ad` file loads reliably in:

```r

SeuratDisk::Convert("file.h5ad", dest = "h5seurat", overwrite = TRUE)
LoadH5Seurat("file.h5seurat", verbose = TRUE,
                    assays = "RNA",
                    reductions = NULL,
                    graphs = NULL,
                    neighbors = NULL,
                    images = FALSE) 

### If you run into memory issues, try running DietSeurat next before further downstream processing, e.g., Azimuth etc.
```

_Note from https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html#converting-from-anndata-to-seurat-via-h5seurat-1 that "The final main parameter is the images parameter; this parameter controls which spatial image data is loaded. All spatial image data are marked global by default, so they are loaded whether or not their associated assays are loaded as well. The images parameter has three special values: NULL for all spatial image data (the default), NA for global spatial image data (typically the same as NULL), or FALSE for no spatial image data."_


**Full tutorial with more utilities will be added at https://denvercal1234GitHub.github.io/anndata_Seurat_Utilities/**


## Why using this package?

Direct conversion of AnnData → Seurat using `LoadH5Seurat` often fails due to:

### Pandas extension dtypes
boolean, Int64, or masked arrays serialize as:
```bash
/meta.data/<column>/mask  
/meta.data/<column>/values  
```

SeuratDisk expects {levels, values}, not {mask, values}.

### Column names containing /

AnnData writes / as HDF5 paths, causing:

```
/meta.data/Kir3dl1
    2__nonSCVI_normlog
```

instead of one metadata column → SeuratDisk throws missing-column errors.


### Empty columns
Fully-NA metadata columns become incomplete HDF5 groups →
SeuratDisk error: "Missing required datasets 'levels' and 'values'".

This utility fixes all of these issues.



## Contributing
Contributions are welcome!
Please open an issue or pull request if you’d like to improve the utilities.

## License
MIT License — see the LICENSE file.

## Citation

If this package contributes to your work, please cite the repository:

```bibtex
@software{Nguyen2025anndata_Seurat_Utilities,
  author = {Quang Nguyen},
  title = {anndata_Seurat_Utilities},
  year = {2025},
  url = {https://github.com/denvercal1234GitHub/anndata_Seurat_Utilities},
  version = {0.1.0},
}
