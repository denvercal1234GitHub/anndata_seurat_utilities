# anndata_seurat_utilities

Welcome to the documentation for **anndata_seurat_utilities**,  
a Python package designed to ensure seamless interoperability between  
**AnnData (.h5ad)** objects and downstream **Seurat / SeuratDisk** workflows in R.

## Features

- Sanitizes and normalizes `obs` / `var` metadata
- Converts problematic pandas dtypes (boolean, Int64, categorical)
- Fixes HDF5 structures that cause SeuratDisk conversion errors
- Removes empty metadata columns
- Sanitizes column names to avoid HDF5 path issues
- Optional removal of heavy layers (e.g., `counts_mouse`)
- Writes `.h5ad` files compatible with `SeuratDisk::Convert()`

## Quick Start

```bash
pip install anndata_seurat_utilities
```

Please see the Usage Guide for an example at https://github.com/denvercal1234GitHub/anndata_Seurat_Utilities

More utilities (including advanced visualisation) will be provided soon. Cheers!