import os
import argparse

import mudata

import numpy as np
from scipy import sparse
from sklearn.metrics import explained_variance_score, r2_score

from joblib import Parallel, delayed
from tqdm.auto import tqdm

import logging
logging.basicConfig(level = logging.INFO)


def to_dense_block(X, j0, j1):
    '''Dense (cells x g) slice for columns [j0:j1).'''
    Xb = X[:, j0:j1]
    return Xb.toarray() if sparse.issparse(Xb) else np.asarray(Xb)

def sum_var_cols_dense(A, ddof=1):
    '''Sum of per-gene variances in a dense block.'''
    return float(np.sum(np.var(A, axis=0, ddof=ddof)))

def var_X_sum_blockwise(X, *, block=1024, ddof=1):
    '''Sum_j Var(X[:, j]) computed block-wise.'''
    _, n_genes = X.shape
    total = 0.0
    for j0 in range(0, n_genes, block):
        j1 = min(j0 + block, n_genes)
        total += sum_var_cols_dense(to_dense_block(X, j0, j1), ddof=ddof)
    return total

def l1_normalize_columns(H):
    '''L1-normalize each program column of H (genes x K).'''
    H = np.asarray(H, dtype=float)
    colsum = H.sum(axis=0, keepdims=True)
    colsum[colsum == 0.0] = 1.0
    return H / colsum

def matmul_X_h_blockwise(X, h, *, block=1024):
    '''Compute B = X @ h (cells,) without densifying X.'''
    _, n_genes = X.shape
    B = None
    for j0 in range(0, n_genes, block):
        j1 = min(j0 + block, n_genes)
        Xb = to_dense_block(X, j0, j1)
        part = Xb @ h[j0:j1]
        B = part if B is None else (B + part)
    return B

def residual_var_sum_global(X, W, Hn, *, block=1024, ddof=1, match_totals=False):
    '''
    Sum_j Var((X - X̂)[:, j]) for global reconstruction X̂ = W @ Hn^T, block-wise.
    Hn should already be normalized as intended (e.g., L1).
    '''
    _, n_genes = X.shape
    total = 0.0
    for j0 in range(0, n_genes, block):
        j1 = min(j0 + block, n_genes)
        Hb   = Hn[j0:j1, :]             # (g x K)
        Xhat = W @ Hb.T                 # (cells x g), dense
        Xb   = to_dense_block(X, j0, j1)
        if match_totals:
            s_true = Xb.sum(axis=1, keepdims=True)
            s_hat  = np.clip(Xhat.sum(axis=1, keepdims=True), 1e-12, None)
            Xhat   = (Xhat.T * (s_true / s_hat).squeeze()).T
        total += sum_var_cols_dense(Xb - Xhat, ddof=ddof)
    return total

def residual_var_sum_component(X, H, k, B, *, block=1024, ddof=1):
    '''
    Sum_j Var((X - B ⊗ h_k)[:, j]) for a single component k, block-wise.
    H is genes x K, B is (cells,).
    '''
    _, n_genes = X.shape
    total = 0.0
    for j0 in range(0, n_genes, block):
        j1 = min(j0 + block, n_genes)
        hk  = H[j0:j1, k]               # (g,)
        Xb  = to_dense_block(X, j0, j1) # (cells x g)
        Xkb = np.outer(B, hk)           # (cells x g)
        total += sum_var_cols_dense(Xb - Xkb, ddof=ddof)
    return total

def _compute_global_variance_explained(X, W, H, *, l1_norm=True, match_totals=False, block=1024, ddof=1):
    '''
    Global explained variance for full reconstruction X̂ = W @ H^T (variance_weighted),
    computed block-wise over genes.
    '''
    W  = np.asarray(W, dtype=float)
    Hn = l1_normalize_columns(H) if l1_norm else np.asarray(H, dtype=float)

    denom = var_X_sum_blockwise(X, block=block, ddof=ddof)
    num   = residual_var_sum_global(X, W, Hn, block=block, ddof=ddof, match_totals=match_totals)
    return 1.0 - (num / denom)

def _compute_variance_explained(X, H, *, mode='ls', block=1024, ddof=1):
    '''
    Per-component EV using ONLY a single program vector h_k each time.
      mode='kang' -> B = (X h)/||h||
      mode='ls'     -> B = (X h)/(h·h)
    Returns: (K,) array of EVs.
    '''
    H   = np.asarray(H, dtype=float)
    K   = H.shape[1]
    denom = var_X_sum_blockwise(X, block=block, ddof=ddof)

    h_sq   = np.sum(H * H, axis=0)
    h_norm = np.sqrt(np.maximum(h_sq, 0.0))
    h_norm[h_norm == 0.0] = 1.0

    ev = np.empty(K, dtype=float)
    it = tqdm(range(K), desc=f'Computing variance explained - ({mode})', unit='programs')
    for k in it:
        h = H[:, k]
        B = matmul_X_h_blockwise(X, h, block=block)
        if mode == 'kang':
            B = B / h_norm[k]
        elif mode == 'ls':
            if h_sq[k] == 0.0:
                ev[k] = 0.0
                continue
            B = B / h_sq[k]
        else:
            raise ValueError('mode must be kang or ls')

        num = residual_var_sum_component(X, H, k, B, block=block, ddof=ddof)
        ev[k] = 1.0 - (num / denom)
    return ev

def compute_explained_variance_ratio(mdata, prog_key='prog', data_key='rna',
                                     data_layer='X', prog_layer='X', n_jobs=1, inplace=True, *,
                                     l1_normalize_H=True, match_totals=False, block=1024,
                                     compute_per_component='kang', ddof=1):
    '''
    Computes the proportion of variance in the data explained by program usage x spectra
    reconstructions (global) and optionally by each individual program (per-component).

    ARGS
        mdata : MuData
            MuData object containing program usages/loadings and expression data,
            or a path to a .h5mu file.
        prog_key: str (default: 'prog')
            Key for the program-level AnnData object inside the MuData.
        data_key: str (default: 'rna')
            Key for the expression data AnnData object inside the MuData.
        data_layer: str (default: 'X')
            Layer of the data_key AnnData to use as the expression matrix (cells x genes).
        prog_layer: str (default: 'X')
            Layer of the prog_key AnnData to use as the program usage matrix (cells x K).
        n_jobs: int (default: 1)
            Number of threads to run processes on. (Currently unused, kept for compatibility.)
        inplace: bool (default: True)
            Update the MuData object in place. If False, return a copy.

    KWARGS
        l1_normalize_H: bool (default: True)
            If True, L1-normalize columns of the spectra/loadings matrix before computing EV,
            matching the Engreitz cNMF variance_explained_v2.py script.
        match_totals: bool (default: False)
            If True, rescale reconstructed X̂ per cell to match total counts of X
            over the gene set, neutralizing differences in usage scaling.
        block: int (default: 1024)
            Block size for iterating over genes to keep memory usage manageable.
        compute_per_component: tuple (default: ('kang','ls'))
            Which per-component EVs to compute:
              - 'kang' : heuristic rank-1 EV as in variance_explained_v2.py
              - 'ls'     : least-squares optimal rank-1 EV (recommended)
        ddof: int (default: 1)
            Degrees of freedom for variance calculations. Use 1 to match the variance_explained_v2.py script,
            or 0 to match sklearn defaults.

    RETURNS
        If inplace=True:
            Results are written into mdata in place:
              - Global EV: mdata[prog_key].uns['metrics']['global_explained_variance_<data_layer>_<prog_layer>']
              - Per-component EVs: mdata[prog_key].var['explained_variance_ratio_<data_layer>_<prog_layer>']
        If inplace=False:
            Returns a DataFrame view of mdata[prog_key].var containing the column
            'explained_variance_ratio_<data_layer>_<prog_layer>'.
    '''
    # Read in mudata if it is provided as a path
    frompath=False
    if isinstance(mdata, str):
        if os.path.exists(mdata):
            mdata = mudata.read(mdata)
            if inplace:
                logging.warning('Changed to inplace=False since path was provided')
                inplace=False
            frompath=True
        else: raise ValueError('Incorrect mudata specification.')
    
    if not inplace and not frompath:
        mdata = mudata.MuData({prog_key: mdata[prog_key].copy(),
                               data_key: mdata[data_key].copy()})

    # Check if same number of vars are present
    try: assert mdata[prog_key].varm['loadings'].shape[1]==mdata[data_key].var.shape[0]
    except: raise ValueError('Different number of features present in data and program loadings')

    # Check input matrices
    X = mdata[data_key].X if data_layer == 'X' else mdata[data_key].layers[data_layer]
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
        if data_layer == 'X': mdata[data_key].X = X
        else:                 mdata[data_key].layers[data_layer] = X

    W = mdata[prog_key].X if prog_layer == 'X' else mdata[prog_key].layers[prog_layer]
    if sparse.issparse(W): W = W.toarray()

    H = mdata[prog_key].varm['loadings'].T  # genes x K
    if H.shape[0] != mdata[data_key].var.shape[0]:
        raise ValueError('Different number of features present in data and program loadings')

    # Global explained variance
    ev_total = _compute_global_variance_explained(
        X, W, H, l1_norm=l1_normalize_H, match_totals=match_totals, block=block, ddof=ddof
    )
    if 'metrics' not in mdata[prog_key].uns:
        mdata[prog_key].uns['metrics'] = {}
    key = f'global_explained_variance_{data_layer}_{prog_layer}'
    mdata[prog_key].uns['metrics'][key] = float(ev_total)

    # per component explained variance
    if compute_per_component:
        Hn = l1_normalize_columns(H) if l1_normalize_H else np.asarray(H, dtype=float)
        if 'kang' in compute_per_component:
            ev_kang = _compute_variance_explained(
                X, Hn, mode='kang', block=block, ddof=ddof,
            )
            mdata[prog_key].var[f'explained_variance_ratio_{data_layer}_{prog_layer}'] = ev_kang
        elif 'ls' in compute_per_component:
            ev_ls = _compute_variance_explained(
                X, Hn, mode='ls', block=block, ddof=ddof,
            )
            mdata[prog_key].var[f'explained_variance_ratio_{data_layer}_{prog_layer}'] = ev_ls
        else:
            raise ValueError('Provide a valid component wise EV calculation option')

    if not inplace: return (mdata[prog_key].var.loc[:, [f'explained_variance_ratio_{data_layer}_{prog_layer}']])

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Compute global and per-component explained variance for programs."
    )
    parser.add_argument('mudataObj_path', help='Path to .h5mu (MuData) object')

    parser.add_argument('--data_layer', default='X', type=str)
    parser.add_argument('--prog_layer', default='X', type=str)
    parser.add_argument('-pk', '--prog_key', default='prog', type=str)
    parser.add_argument('-dk', '--data_key', default='rna', type=str)
    parser.add_argument('-n', '--n_jobs', default=1, type=int)  # kept for API compat
    parser.add_argument('--output', action='store_false',
                        help='If set, run NOT in-place (kept for backwards compatibility)')

    # Optional
    parser.add_argument('--no_l1_norm', action='store_true',
                        help='Disable L1 column normalization of H (default: enabled)')
    parser.add_argument('--match_totals', action='store_true',
                        help='Per-cell scale X̂ to match X totals over HVGs (default: off)')
    parser.add_argument('--block', default=1024, type=int,
                        help='Gene block size for block-wise computation (default: 1024)')
    parser.add_argument('--per_component', choices=['kang', 'ls'],
                        default='both',
                        help="Per-component EVs to compute (default: both)")
    parser.add_argument('--ddof', default=1, type=int,
                        help='Degrees of freedom for variance (default: 1 to mirror script)')

    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    # Load and run
    md = mudata.read(args.mudataObj_path)
    compute_explained_variance_ratio(
        md,
        prog_key=args.prog_key,
        data_key=args.data_key,
        data_layer=args.data_layer,
        prog_layer=args.prog_layer,
        n_jobs=args.n_jobs,
        inplace=args.output,          
        l1_normalize_H=not args.no_l1_norm,
        match_totals=args.match_totals,
        block=args.block,
        compute_per_component=args.per_component,
        ddof=args.ddof
    )
