import tempfile
import anndata as ad
import numpy as np
import pandas as pd
from sage_broccoliqn.prepare_for_seurat import prepare_adata_for_seurat_drop_empty_v3

def test_write_simple(tmp_path):
    adata = ad.AnnData(X=np.zeros((5,3)))
    adata.obs['col1'] = pd.Series([True, False, None, True, False], dtype="boolean")
    out = tmp_path / "test.h5ad"
    prepare_adata_for_seurat_drop_empty_v3(adata, out, drop_layers=(), verbose=False)
    assert out.exists()
