# prepare_adata_for_seurat_drop_empty_v3.py
"""
Robust writer to produce an .h5ad safe for Seurat/SeuratDisk:
  - drops empty columns
  - converts risky pandas extension dtypes to safe strings/numeric
  - sanitizes var/obs names and records mapping
  - optionally subsets/drops adata.obsm / adata.obsp / adata.obs columns / adata.uns
  - final forced conversion to avoid mask+values groups in h5ad serialization
  - logs removed/converted keys into adata.uns for auditing
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Iterable, Optional, Union
import anndata as ad


def _is_all_missing_series(ser: pd.Series) -> bool:
    if pd.api.types.is_string_dtype(ser) or ser.dtype == object:
        vals = ser.astype(object).values
        mask_nonmissing = [
            (v is not None and not (isinstance(v, float) and np.isnan(v)) and str(v).strip() != "")
            for v in vals
        ]
        return not any(mask_nonmissing)
    else:
        return ser.isna().all()


def _is_risky_extension_dtype(ser: pd.Series) -> bool:
    """
    True if ser has a pandas extension dtype or other dtype that anndata may serialize
    as mask+values (e.g., 'boolean' nullable bool, 'Int64' nullable int, pandas string/categorical, object).
    """
    dt = str(ser.dtype)
    if dt.startswith("boolean") or dt.startswith("Int64") or "Nullable" in dt:
        return True
    if pd.api.types.is_categorical_dtype(ser) or pd.api.types.is_object_dtype(ser):
        return True
    # also catch pandas 'string' extension dtype
    if dt == "string":
        return True
    return False


def _serializable_series(ser: pd.Series, missing_to_empty_str: bool = True) -> pd.Series:
    """
    Convert risky extension dtype series into plain Python strings or safe numeric arrays.
    Keeps numeric/bool numpy dtypes unchanged.
    """
    # safe numeric/bool: keep as-is
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

    # fallback: cast to str (should not happen for numeric/bool)
    return ser.astype(str).fillna("" if missing_to_empty_str else "nan")


def _clean_dataframe_for_h5ad(df: pd.DataFrame, missing_to_empty_str: bool = True, drop_empty: bool = True):
    """
    Clean a DataFrame (used for adata.obs and adata.var) so that
    serialization to .h5ad will not create mask+values groups.

    Returns (clean_df, dropped_columns_list, converted_columns_list)
    """
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
        # safe numeric/bool/string
        if pd.api.types.is_string_dtype(ser) and missing_to_empty_str:
            df_clean[col] = ser.fillna("").values
        else:
            df_clean[col] = ser.values
    return df_clean, dropped, converted


def _sanitize_names(names):
    """
    Turn a list of names into sanitized unique names (no slashes/spaces/punctuation)
    and return (sanitized_list, mapping_sanitized_to_original).
    """
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
    """
    Final safety pass: detect any ad.obs/ad.var columns that still have risky extension dtypes
    and force-convert to plain strings ("" for missing).
    Returns lists of columns forced-converted in obs and var.
    """
    forced_obs = []
    forced_var = []

    # OBS
    for col in list(ad.obs.columns):
        ser = ad.obs[col]
        if _is_risky_extension_dtype(ser):
            ad.obs[col] = ser.astype(str).fillna("" if missing_to_empty_str else "nan")
            forced_obs.append(col)

    # VAR
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
    # obsm / obsp / obs control:
    # - iterable => keep only these keys
    # - 'all' => keep everything (default)
    # - None => drop all
    keep_obsm: Union[Iterable[str], str, None] = 'all',
    keep_obsp: Union[Iterable[str], str, None] = 'all',
    keep_obs_cols: Union[Iterable[str], str, None] = 'all',
    # NEW: control which uns keys to KEEP:
    # - iterable => keep only listed keys
    # - 'all' => keep everything (default)
    # - None => drop all entries in adata.uns
    keep_uns: Union[Iterable[str], str, None] = 'all',
    compression: Optional[str] = None,
    missing_to_empty_str: bool = True,
    convert_index_to_str: bool = True,
    convert_strings_to_categoricals_on_write: bool = False,
    sanitize_var_names: bool = True,
    sanitize_obs_names: bool = True,
    drop_empty_obs_and_var: bool = True,
    verbose: bool = True
) -> str:
    """
    Prepare an AnnData copy and write to .h5ad in a robust, Seurat-friendly way.

    Parameters (high level):
      - keep_obsm / keep_obsp / keep_obs_cols / keep_uns control subsetting behavior.
        semantics: iterable -> keep only those items; 'all' -> keep everything (default); None -> drop all.
      - drop_layers / keep_layers similar to previous behavior for adata.layers.
      - Final metadata about changes recorded in ad_cp.uns['prepare_adata_for_seurat_dropped'] and
        removed uns keys logged in ad_cp.uns['prepare_adata_for_seurat_removed_uns_keys'] when applicable.

    Returns: path string of written file.
    """
    out_path = Path(out_h5ad_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if verbose:
        print("prepare_adata_for_seurat_drop_empty_v3: preparing copy and cleaning...")

    ad_cp = adata.copy()

    # -------------------------
    # Manage .uns according to keep_uns policy
    # -------------------------
    if keep_uns != 'all':
        if verbose:
            print("  applying keep_uns policy:", keep_uns)
        existing_uns_keys = list(ad_cp.uns.keys())

        if keep_uns is None:
            keep_set = set()
        else:
            keep_set = set(map(str, keep_uns))

        keys_to_remove = [k for k in existing_uns_keys if k not in keep_set]

        removed_keys_logged = []
        if keys_to_remove:
            for k in keys_to_remove:
                try:
                    del ad_cp.uns[k]
                    removed_keys_logged.append(k)
                    if verbose:
                        print(f"    removed uns['{k}']")
                except KeyError:
                    pass
        else:
            if verbose:
                print("    no uns keys to remove under keep_uns policy.")

        # always write the audit log (even if empty)
        ad_cp.uns['prepare_adata_for_seurat_removed_uns_keys'] = removed_keys_logged

    # -------------------------
    # Manage layers
    # -------------------------
    if keep_layers is not None:
        keep_set = set(keep_layers)
        for L in list(ad_cp.layers.keys()):
            if L not in keep_set:
                if verbose:
                    print(f"  removing layer '{L}'")
                del ad_cp.layers[L]
    else:
        if drop_layers:
            for L in list(drop_layers):
                if L in ad_cp.layers:
                    if verbose:
                        print(f"  dropping layer '{L}'")
                    del ad_cp.layers[L]

    # -------------------------
    # Manage obsm
    # -------------------------
    if keep_obsm != 'all':
        if verbose:
            print("  processing obsm keep policy:", keep_obsm)
        existing_obsm = list(ad_cp.obsm.keys())
        if keep_obsm is None:
            for k in existing_obsm:
                if verbose:
                    print(f"    dropping obsm['{k}']")
                del ad_cp.obsm[k]
        else:
            keep_set = set(map(str, keep_obsm))
            for k in existing_obsm:
                if k not in keep_set:
                    if verbose:
                        print(f"    dropping obsm['{k}']")
                    del ad_cp.obsm[k]

    # -------------------------
    # Manage obsp
    # -------------------------
    if keep_obsp != 'all':
        if verbose:
            print("  processing obsp keep policy:", keep_obsp)
        existing_obsp = list(ad_cp.obsp.keys())
        if keep_obsp is None:
            for k in existing_obsp:
                if verbose:
                    print(f"    dropping obsp['{k}']")
                del ad_cp.obsp[k]
        else:
            keep_set = set(map(str, keep_obsp))
            for k in existing_obsp:
                if k not in keep_set:
                    if verbose:
                        print(f"    dropping obsp['{k}']")
                    del ad_cp.obsp[k]

    # -------------------------
    # Manage obs columns subset
    # -------------------------
    if keep_obs_cols != 'all':
        if verbose:
            print("  processing obs columns keep policy:", keep_obs_cols)
        if keep_obs_cols is None:
            if verbose:
                print("    dropping all obs columns (keep_obs_cols=None)")
            ad_cp.obs = pd.DataFrame(index=ad_cp.obs.index)
        else:
            keep_set = set(map(str, keep_obs_cols))
            present_keep = [c for c in ad_cp.obs.columns if c in keep_set]
            if verbose:
                missing = list(keep_set - set(present_keep))
                if missing:
                    print(f"    warning: requested obs cols not present and will be skipped: {missing}")
                print(f"    keeping obs columns: {present_keep}")
            ad_cp.obs = ad_cp.obs.loc[:, present_keep]

    # -------------------------
    # Clean obs
    # -------------------------
    if verbose:
        print("  cleaning adata.obs ...")
    obs_clean, obs_dropped, obs_converted = _clean_dataframe_for_h5ad(
        ad_cp.obs,
        missing_to_empty_str=missing_to_empty_str,
        drop_empty=drop_empty_obs_and_var
    )
    ad_cp.obs = obs_clean
    if verbose:
        print(f"    dropped {len(obs_dropped)} obs cols; converted {len(obs_converted)} obs cols")

    # -------------------------
    # Clean var
    # -------------------------
    if verbose:
        print("  cleaning adata.var ...")
    var_clean, var_dropped, var_converted = _clean_dataframe_for_h5ad(
        ad_cp.var,
        missing_to_empty_str=missing_to_empty_str,
        drop_empty=drop_empty_obs_and_var
    )
    ad_cp.var = var_clean
    if verbose:
        print(f"    dropped {len(var_dropped)} var cols; converted {len(var_converted)} var cols")

    # -------------------------
    # Ensure obs index is string
    # -------------------------
    if convert_index_to_str:
        try:
            ad_cp.obs.index = ad_cp.obs.index.astype(str)
        except Exception:
            ad_cp.obs.index = ad_cp.obs.index.map(str)

    # -------------------------
    # Sanitize var names
    # -------------------------
    if sanitize_var_names:
        orig_var_names = list(ad_cp.var_names)
        ad_cp.var['orig_var_name'] = orig_var_names
        new_var_names, var_map = _sanitize_names(orig_var_names)
        ad_cp.var_names = new_var_names
        ad_cp.uns['orig_var_names_map'] = var_map
        if verbose:
            print("  sanitized var names; mapping saved in uns['orig_var_names_map']")

    # -------------------------
    # Sanitize obs column NAMES
    # -------------------------
    if sanitize_obs_names:
        orig_obs_cols = list(ad_cp.obs.columns)
        new_obs_cols, obs_map = _sanitize_names(orig_obs_cols)
        ad_cp.obs.columns = new_obs_cols
        ad_cp.uns['orig_obs_names_map'] = obs_map
        if verbose:
            print("  sanitized obs column names; mapping saved in uns['orig_obs_names_map']")

    # -------------------------
    # Record dropped/converted summary
    # -------------------------
    ad_cp.uns['prepare_adata_for_seurat_dropped'] = {
        'obs_dropped': obs_dropped,
        'var_dropped': var_dropped,
        'obs_converted_initial': obs_converted,
        'var_converted_initial': var_converted
    }

    # -------------------------
    # Final force convert remaining risky dtypes -> strings
    # -------------------------
    forced_obs, forced_var = _final_force_convert_mask_values(
        ad_cp,
        missing_to_empty_str=missing_to_empty_str,
        verbose=verbose
    )
    ad_cp.uns['prepare_adata_for_seurat_dropped']['obs_forced_converted_final'] = forced_obs
    ad_cp.uns['prepare_adata_for_seurat_dropped']['var_forced_converted_final'] = forced_var

    # -------------------------
    # Diagnostics / Logging
    # -------------------------
    if verbose:
        print("  final sanity: n_obs:", ad_cp.n_obs, "n_vars:", ad_cp.n_vars)
        print("  sample obs columns:", list(ad_cp.obs.columns)[:20])
        print("  sample var names:", list(ad_cp.var_names)[:20])
        if hasattr(ad_cp, "obsm"):
            print("  obsm keys:", list(ad_cp.obsm.keys()))
        if hasattr(ad_cp, "obsp"):
            print("  obsp keys:", list(ad_cp.obsp.keys()))
        print("  remaining uns keys (sample 50):", list(ad_cp.uns.keys())[:50])

    # -------------------------
    # Write out to .h5ad with requested flags
    # -------------------------
    if verbose:
        print(f"  writing .h5ad to: {out_path} (convert_strings_to_categoricals={convert_strings_to_categoricals_on_write})")
    try:
        ad_cp.write_h5ad(str(out_path), compression=compression, convert_strings_to_categoricals=convert_strings_to_categoricals_on_write)
    except TypeError:
        # older anndata versions may not accept convert_strings_to_categoricals kw
        ad_cp.write_h5ad(str(out_path), compression=compression)
    if verbose:
        print("  write complete.")
    return str(out_path)