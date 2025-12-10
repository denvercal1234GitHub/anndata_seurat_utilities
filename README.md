# anndata_Seurat_Utilities

Utilities to prepare AnnData files for Seurat/SeuratDisk conversion and safe exchange between Python and R.

Includes the robust `prepare_adata_for_seurat_drop_empty_v3()` function which:

- drops empty obs/var columns (records them in `.uns`),
- sanitizes obs/var names (avoids `/` and other HDF5 path-breaking chars),
- converts pandas extension dtypes (nullable booleans/Int64/categorical) into atomic types,
- optionally drops heavy layers like `counts_mouse`,
- writes `.h5ad` safe for `SeuratDisk::Convert()` and `LoadH5Seurat()`.

## Quick install

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
