"""
Utilities to prepare AnnData for SeuratDisk conversion.
Contains prepare_adata_for_seurat_drop_empty_v3().
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Iterable, Optional
import anndata as ad

def _is_all_missing_series(ser: pd.Series) -> bool:
    if pd.api.types.is_string_dtype(ser) or ser.dtype == object:
        vals = ser.astype(object).values
        mask_nonmissing = [ (v is not None and not (isinstance(v, float) and np.isnan(v)) and str(v).strip() != "" ) for v in vals ]
        return not any(mask_nonmissing)
    else:
        return ser.isna().all()

def _is_risky_extension_dtype(ser: pd.Series) -> bool:
    dt = str(ser.dtype)
    if dt.startswith("boolean") or dt.startswith("Int64") or "Nullable" in dt:
        return True
    if pd.api.types.is_categorical_dtype(ser) or pd.api.types.is_object_dtype(ser):
        return True
    if dt == "string":
        return True
    return False

def _serializable_series(ser: pd.Series, missing_to_empty_str: bool = True) -> pd.Series:
    if pd.api.types.is_integer_dtype(ser) or pd.api.types.is_float_dtype(ser) or pd.api.types.is_bool_dtype(ser):
        return ser.copy()
    if _is_risky_extension_dtype(ser):
        out_vals = []
        for v in ser.values:
            if v is None or (isinstance(v, float) and np.isnan(v)):
                out_vals.append("" if missing_to_empty_str else None)
            elif isinstance(v, (list, tuple)):
                try:
                    out_vals.append(";".join(map(str, v)))
                except Exception:
                    out_vals.append(json.dumps(v))
            elif isinstance(v, dict):
                out_vals.append(json.dumps(v))
            else:
                try:
                    out_vals.append(str(v))
                except Exception:
                    out_vals.append(json.dumps(v))
        return pd.Series(out_vals, index=ser.index)
    return ser.astype(str).fillna("" if missing_to_empty_str else "nan")

def _clean_dataframe_for_h5ad(df: pd.DataFrame, missing_to_empty_str: bool = True, drop_empty: bool = True):
    df_clean = df.copy(deep=False)
    dropped = []
    converted = []
    for col in list(df_clean.columns):
        ser = df_clean[col]
        if drop_empty and _is_all_missing_series(ser):
            dropped.append(col)
            del df_clean[col]
            continue
        if pd.api.types.is_categorical_dtype(ser):
            df_clean[col] = ser.astype(str).fillna("" if missing_to_empty_str else "nan").values
            converted.append(col)
            continue
        if _is_risky_extension_dtype(ser):
            df_clean[col] = _serializable_series(ser, missing_to_empty_str=missing_to_empty_str).values
            converted.append(col)
            continue
        if pd.api.types.is_string_dtype(ser) and missing_to_empty_str:
            df_clean[col] = ser.fillna("").values
        else:
            df_clean[col] = ser.values
    return df_clean, dropped, converted

def _sanitize_names(names):
    mapping = {}
    seen = {}
    out = []
    for name in names:
        s = str(name)
        s2 = s.replace("/", "_").replace(" ", "_").replace(",", "_").replace(";", "_").replace(":", "_")
        s2 = "_".join([p for p in s2.split("_") if p != ""])
        if s2 == "":
            s2 = "NA_col"
        base = s2
        i = 1
        while s2 in seen:
            i += 1
            s2 = f"{base}__dup{i}"
        seen[s2] = True
        out.append(s2)
        mapping[s2] = s
    return out, mapping

def _final_force_convert_mask_values(ad: ad.AnnData, missing_to_empty_str: bool = True, verbose: bool = True):
    forced_obs = []
    forced_var = []
    for col in list(ad.obs.columns):
        ser = ad.obs[col]
        if _is_risky_extension_dtype(ser):
            ad.obs[col] = ser.astype(str).fillna("" if missing_to_empty_str else "nan")
            forced_obs.append(col)
    for col in list(ad.var.columns):
        ser = ad.var[col]
        if _is_risky_extension_dtype(ser):
            ad.var[col] = ser.astype(str).fillna("" if missing_to_empty_str else "nan")
            forced_var.append(col)
    if verbose:
        print("Final pass forced conversion of risky dtypes -> strings")
        print("  obs forced-converted cols (sample 50):", forced_obs[:50])
        print("  var forced-converted cols (sample 50):", forced_var[:50])
    return forced_obs, forced_var

def prepare_adata_for_seurat_drop_empty_v3(
    adata: ad.AnnData,
    out_h5ad_path: str | Path,
    *,
    drop_layers: Optional[Iterable[str]] = ('counts_mouse',),
    keep_layers: Optional[Iterable[str]] = None,
    compression: Optional[str] = None,
    missing_to_empty_str: bool = True,
    convert_index_to_str: bool = True,
    convert_strings_to_categoricals_on_write: bool = False,
    sanitize_var_names: bool = True,
    sanitize_obs_names: bool = True,
    drop_empty_obs_and_var: bool = True,
    remove_uns_refs: bool = False,
    verbose: bool = True
) -> str:
    out_path = Path(out_h5ad_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if verbose:
        print("prepare_adata_for_seurat_drop_empty_v3: preparing copy and cleaning...")
    ad_cp = adata.copy()
    if keep_layers is not None:
        keep_set = set(keep_layers)
        for L in list(ad_cp.layers.keys()):
            if L not in keep_set:
                if verbose: print(f"  removing layer '{L}'")
                del ad_cp.layers[L]
    else:
        if drop_layers:
            for L in list(drop_layers):
                if L in ad_cp.layers:
                    if verbose: print(f"  dropping layer '{L}'")
                    del ad_cp.layers[L]
    if verbose: print("  cleaning adata.obs ...")
    obs_clean, obs_dropped, obs_converted = _clean_dataframe_for_h5ad(ad_cp.obs, missing_to_empty_str=missing_to_empty_str, drop_empty=drop_empty_obs_and_var)
    ad_cp.obs = obs_clean
    if verbose:
        print(f"    dropped {len(obs_dropped)} obs cols; converted {len(obs_converted)} obs cols")
    if verbose: print("  cleaning adata.var ...")
    var_clean, var_dropped, var_converted = _clean_dataframe_for_h5ad(ad_cp.var, missing_to_empty_str=missing_to_empty_str, drop_empty=drop_empty_obs_and_var)
    ad_cp.var = var_clean
    if verbose:
        print(f"    dropped {len(var_dropped)} var cols; converted {len(var_converted)} var cols")
    if convert_index_to_str:
        try:
            ad_cp.obs.index = ad_cp.obs.index.astype(str)
        except Exception:
            ad_cp.obs.index = ad_cp.obs.index.map(str)
    if sanitize_var_names:
        orig_var_names = list(ad_cp.var_names)
        ad_cp.var['orig_var_name'] = orig_var_names
        new_var_names, var_map = _sanitize_names(orig_var_names)
        ad_cp.var_names = new_var_names
        ad_cp.uns['orig_var_names_map'] = var_map
        if verbose: print("  sanitized var names; mapping saved in uns['orig_var_names_map']")
    if sanitize_obs_names:
        orig_obs_cols = list(ad_cp.obs.columns)
        new_obs_cols, obs_map = _sanitize_names(orig_obs_cols)
        ad_cp.obs.columns = new_obs_cols
        ad_cp.uns['orig_obs_names_map'] = obs_map
        if verbose: print("  sanitized obs column names; mapping saved in ad_cp.uns['orig_obs_names_map']")
    ad_cp.uns['prepare_adata_for_seurat_dropped'] = {
        'obs_dropped': obs_dropped,
        'var_dropped': var_dropped,
        'obs_converted_initial': obs_converted,
        'var_converted_initial': var_converted
    }
    forced_obs, forced_var = _final_force_convert_mask_values(ad_cp, missing_to_empty_str=missing_to_empty_str, verbose=verbose)
    ad_cp.uns['prepare_adata_for_seurat_dropped']['obs_forced_converted_final'] = forced_obs
    ad_cp.uns['prepare_adata_for_seurat_dropped']['var_forced_converted_final'] = forced_var
    if verbose:
        print("  final sanity: n_obs:", ad_cp.n_obs, "n_vars:", ad_cp.n_vars)
        print("  sample obs columns:", list(ad_cp.obs.columns)[:20])
        print("  sample var names:", list(ad_cp.var_names)[:20])
    if verbose:
        print(f"  writing .h5ad to: {out_path} (convert_strings_to_categoricals={convert_strings_to_categoricals_on_write})")
    try:
        ad_cp.write_h5ad(str(out_path), compression=compression, convert_strings_to_categoricals=convert_strings_to_categoricals_on_write)
    except TypeError:
        ad_cp.write_h5ad(str(out_path), compression=compression)
    if verbose: print("  write complete.")
    return str(out_path)
