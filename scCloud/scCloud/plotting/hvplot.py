import holoviews as hv

hv.extension('bokeh')
import hvplot.pandas
import numpy as np
import pandas as pd
import scipy.sparse


def hv_violin(adata, keys, by=None, log=False, width=200, cols=3, cmap='Category20', **kwds):
    plots = []
    keywords = {'padding': 0.02}
    keywords.update(kwds)
    if not isinstance(keys, (list, tuple)):
        keys = [keys]
    for key in keys:
        if key in adata.var.index:
            X = adata[:, key].X
            if scipy.sparse.issparse(X):
                X = X.toarray()
            df = pd.DataFrame(X, columns=[key])
            if by is not None:
                df[by] = adata.obs[by].values
        else:
            df = adata.obs
        if by is not None and str(df[by].dtype) == 'category':
            df[by] = df[by].astype(str)
        p = df.hvplot.violin(key, width=width, logy=log, by=by, violin_color=by, cmap=cmap, **keywords)
        plots.append(p)

    return hv.Layout(plots).cols(cols)


def hv_heatmap(adata, keys, by, reduce_function=np.mean, **kwds):
    if not isinstance(keys, (list, tuple)):
        keys = [keys]
    df = None
    keywords = {'colorbar': True, 'xlabel': str(by), 'ylabel': ''}
    keywords.update(kwds)
    for key in keys:
        X = adata[:, key].X
        if scipy.sparse.issparse(X):
            X = X.toarray()
        _df = pd.DataFrame(X, columns=['value'])
        _df['feature'] = key
        _df[by] = adata.obs[by].values
        df = _df if df is None else pd.concat((df, _df))

    return df.hvplot.heatmap(x=by, y='feature', C='value', reduce_function=reduce_function,
                             **keywords)


def hv_scatter_matrix(adata, keys, c=None, **kwds):
    if not isinstance(keys, (list, tuple)):
        keys = [keys]

    df = pd.DataFrame(index=adata.obs.index)
    if c is not None:
        df[c] = adata.obs[c].values
    for key in keys:
        if key in adata.var.index:
            X = adata[:, key].X
            if scipy.sparse.issparse(X):
                X = X.toarray()
            df[key] = X
        else:
            df[key] = adata.obs[key].values

    return hvplot.scatter_matrix(df, c=c, **kwds)


def hv_scatter(adata, basis, keys, alpha=1, cmap='viridis', sort=True, size=12, width=300, height=300, cols=3, **kwds):
    # color can be obs (by) or var_name (c)
    if not isinstance(keys, (list, tuple)):
        keys = [keys]
    keywords = {'fontsize': {'title': 9}, 'padding': 0.02, 'xaxis': False, 'yaxis': False}
    keywords.update(kwds)
    df = pd.DataFrame(adata.obsm['X_' + basis][:, 0:2], columns=[basis + c for c in ['1', '2']])
    plots = []
    for key in keys:
        if key in adata.var.index:
            X = adata[:, key].X
            if scipy.sparse.issparse(X):
                X = X.toarray()
            df[key] = X
        else:
            df[key] = adata.obs[key].values
        is_numeric = pd.api.types.is_numeric_dtype(df[key])
        if sort and is_numeric:
            df.sort_values(by=key, inplace=True)

        p = df.hvplot.scatter(x=basis + '1', y=basis + '2',
                              title=str(key), by=key if not is_numeric else None, c=key if is_numeric else None,
                              alpha=alpha, cmap=cmap, size=size, colorbar=is_numeric,
                              width=width, height=height, **keywords)
        plots.append(p)
    return hv.Layout(plots).cols(cols)
