from anndata import AnnData
from typing import Union, Optional
import pandas as pd
import numpy as np
from .distance import Distance
from .align import DistanceAlignment
from . import logger

def selfmatch(X: Union[AnnData, pd.DataFrame],
              columns: Union[list, tuple, np.ndarray, pd.Series, pd.Index],
              calculate_distance: bool = False,
              use_rep: Optional[str] = None, metric: Optional[str] = None,
              normalize: bool = True, Gaussian_kernel: bool = False,
              minimum_unique_percents: Union[list, tuple, np.ndarray, pd.Series, pd.Index, float] = (0.4, 0.5, 0.6, 0.7, 0.8),
              minimum_divide_percents: Union[list, tuple, np.ndarray, pd.Series, pd.Index, float] = (0.1, 0.15, 0.2),
              reannotate: bool = True, prefix: str = '') -> DistanceAlignment:
    """
    Match different versions of cell type annotations (e.g., different resolutions of clustering) for cells from a single dataset.

    Parameters
    ----------
    X
        An :class:`~anndata.AnnData` or :class:`~pandas.DataFrame` object containing information of different cell type annotations as multiple columns of cell metadata.
    columns
        Column names (keys) of cell metadata representing cell type annotations or clusterings.
    calculate_distance
        Whether to calculate the cell-by-cell-type distance matrix. This is usually not necessary as all annotations are in place for a single dataset.
        (Default: `False`)
    use_rep
        Representation used to calculate distances. This can be `'X'` or any representations stored in `.obsm`.
        This argument will be ignored when `calculate_distance = False` (the default).
        Default to the PCA coordinates if present (if not, use the expression matrix `X`).
    metric
        Metric to calculate the distance between each cell and each cell type. Can be `'euclidean'`, `'cosine'`, `'manhattan'` or any metrics applicable to :func:`sklearn.metrics.pairwise_distances`.
        This argument will be ignored when `calculate_distance = False` (the default).
        Default to `'euclidean'` if latent representations are used for calculating distances, and to `'correlation'` if the expression matrix is used.
    normalize
        Whether to normalize the distance matrix.
        This argument will be ignored when `calculate_distance = False` (the default).
        (Default: `True`)
    Gaussian_kernel
        Whether to apply the Gaussian kernel to the distance matrix.
        This argument will be ignored when `calculate_distance = False` (the default).
        (Default: `False`)
    minimum_unique_percents
        The minimum cell assignment fraction(s) to claim two cell types as uniquely matched.
        By default, five values will be tried (0.4, 0.5, 0.6, 0.7, 0.8) to find the one that produces least alignments in each harmonization iteration.
    minimum_divide_percents
        The minimum cell assignment fraction(s) to claim a cell type as divisible into two or more cell types.
        By default, three values will be tried (0.1, 0.15, 0.2) to find the one that produces least alignments in each harmonization iteration.
    reannotate
        Whether to reannotate cells into harmonized cell types.
        (Default: `True`)
    prefix
        Column prefix for the reannotation data frame.

    Returns
    ----------
    DistanceAlignment
        A :class:`~cellhint.align.DistanceAlignment` object. Four important attributes within this class are:
        1) :attr:`~cellhint.align.DistanceAlignment.base_distance`, within-dataset distances between all cells and all cell types.
        2) :attr:`~cellhint.align.DistanceAlignment.relation`, the harmonization table.
        3) :attr:`~cellhint.align.DistanceAlignment.groups`, high-hierarchy cell types categorizing rows of the harmonization table.
        4) :attr:`~cellhint.align.DistanceAlignment.reannotation`, reannotated cell types and cell type groups.
    """
    #input
    if not isinstance(X, (AnnData, pd.DataFrame)):
        raise TypeError(
                f"🛑 Please provide a correct input - an `anndata.AnnData` or `pandas.DataFrame`")
    df = X.obs if isinstance(X, AnnData) else X
    #columns
    if isinstance(columns, str) or len(columns) == 1:
        raise TypeError(
                f"🛑 Please provide at least two columns")
    columns = np.array(columns)
    non_columns = set(columns).difference(df.columns)
    if len(non_columns) >= 1:
        raise ValueError(
                f"🛑 The following column(s) are not found: {non_columns}")
    #meta
    cell = pd.concat([pd.DataFrame(dict(dataset = column, ID = df.index, cell_type = df[column].astype(str))) for column in np.unique(columns)], axis = 0, ignore_index = True)
    cell['ID'] = cell.dataset + '__' + cell.ID
    cell_type = pd.concat([pd.DataFrame(dict(dataset = column, cell_type = np.unique(df[column]))) for column in np.unique(columns)], axis = 0, ignore_index = True)
    #dist mat
    if calculate_distance and not isinstance(X, AnnData):
        raise TypeError(
                f"🛑 To calculate the distance, please provide an AnnData as input")
    if calculate_distance:
        distances = []
        for column in np.unique(columns):
            X.obs['__dataset__'] = '__constant__'
            distances.append(Distance.from_adata(X, dataset = '__dataset__', cell_type = column, use_rep = use_rep, metric = metric, n_jobs = -1, check_params = True))
        distance = distances[0].concatenate(distances[1:], by = 'cell_type', check = False)
        if normalize:
            distance.normalize(Gaussian_kernel = Gaussian_kernel, rank = True, normalize = True)
        combined_distance = distance.concatenate([distance] * (len(columns) - 1), by = 'cell', check = False)
        dist_mat = combined_distance.dist_mat
        X.obs.drop(columns = ['__dataset__'], inplace = True)
    else:
        dist_mat = np.zeros((cell.shape[0], cell_type.shape[0]))
    #distance
    combined_distance = Distance(dist_mat, cell, cell_type)
    #assignment
    combined_distance.assignment = pd.concat([pd.concat([df[column].astype(str) for column in np.unique(columns)], axis = 1)] * len(columns), axis = 0)
    combined_distance.assignment.index = combined_distance.cell.index
    combined_distance.assignment.columns = np.unique(columns)
    #alignment
    alignment = DistanceAlignment(combined_distance, check = False, dataset_order = columns, row_normalize = True, maximum_novel_percent = -1)
    alignment.best_align(dataset_order = None, minimum_unique_percents = minimum_unique_percents, minimum_divide_percents = minimum_divide_percents)
    if reannotate:
        logger.info(f"🖋️ Reannotating cells")
        alignment.reannotate(prefix = prefix)
    logger.info(f"✅ Harmonization done!")
    #return
    return alignment
