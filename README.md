<p align="left"><img src="https://github.com/Teichlab/cellhint/blob/main/docs/source/_static/img/logo_cellhint.png" width="250" height="85"></p>

[![Python Versions](https://img.shields.io/badge/python-3.6+-brightgreen.svg)](https://pypi.org/project/cellhint) [![Documentation Status](https://readthedocs.org/projects/cellhint/badge/?version=latest)](https://cellhint.readthedocs.io/en/latest/?badge=latest)

CellHint is an automated tool for cell type harmonisation and integration.
- _harmonisation_: match and harmonise cell types defined by independent datasets
- _integration_: integrate cell and cell types with supervision from harmonisation

# Interactive tutorials
### Harmonisation
[Using CellHint for cell type harmonisation ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/cellhint/blob/main/docs/notebook/cellhint_tutorial_harmonisation.ipynb)
### Integration
[Using CellHint for annotation-aware data integration ![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Teichlab/cellhint/blob/main/docs/notebook/cellhint_tutorial_integration.ipynb)

# Install CellHint
### Using pip [![PyPI](https://img.shields.io/pypi/v/cellhint.svg?color=brightgreen&style=flat)](https://pypi.org/project/cellhint)
```console
pip install cellhint
```

### Using conda [![install with bioconda](https://img.shields.io/conda/vn/bioconda/cellhint.svg?color=brightgreen&style=flat)](https://anaconda.org/bioconda/cellhint)
```console
conda install -c bioconda -c conda-forge cellhint
```

# Usage (harmonisation)

<details>
<summary><strong>1. Cross-dataset cell type harmonisation</strong></summary>

+ <details>
  <summary><strong>1.1. Cell type harmonisation</strong></summary>

  The input [AnnData](https://anndata.readthedocs.io/en/latest/) needs two columns in `.obs` representing dataset origin and cell original annotation respectively. The aim is to harmonise cell types across datasets using [cellhint.harmonize](https://cellhint.readthedocs.io/en/latest/cellhint.harmonize.html).  
    
  Internally, transcriptional distances between cells and cell types (denoted here as the cell centroid) will first be calculated. Since cell type is usually defined at the cluster level and no cluster is 100% pure, you can set `filter_cells = True` (default to `False`) to filter out cells whose gene expression profiles do not correlate most with the cell type they belong to. This will speed up the run as only a subset of cells are used, but will render these filtered cells unannotated (see `2.2.`). Distances are calculated at either gene or low-dimensional space. The latter is preferred to denoise the data by providing a latent representation via the argument `use_rep` (default to PCA coordinates).
  ```python
  #`use_rep` can be omitted here as it defaults to 'X_pca'.
  alignment = cellhint.harmonize(adata, dataset = 'dataset_column', cell_type = 'celltype_column', use_rep = 'X_pca')
  ```
  If `X_pca` is not detected in `.obsm` and no other latent representations are provided via `use_rep`, gene expression matrix in `.X` will be used to calculate the distances. In such case, subsetting the AnnData to informative genes (e.g. highly variable genes) is suggested and `.X` should be log-normalised (to a constant total count per cell).  
    
  The resulting `alignment` is an instance of the class [DistanceAlignment](https://cellhint.readthedocs.io/en/latest/cellhint.align.DistanceAlignment.html) as defined by CellHint, and can be written out as follows.
  ```python
  #Save the harmonisation output.
  alignment.write('/path/to/local/folder/some_name.pkl')
  ```
  </details>

+ <details>
  <summary><strong>1.2. Cell type harmonisation with PCT</strong></summary>

  Inferring cell type relationships based on directly calculated distances will suffice in most cases due to a normalisation procedure applied to the derived distances. If a very strong batch effect exists across datasets, you can turn on `use_pct = True` (default to `False`) to predict instead of calculate these distances. Through this parameter, a predictive clustering tree (PCT) is built for each dataset, and distances between cells in query datasets and cell types in the reference dataset are predicted, often resulting in unbiased distance measures.
  ```python
  #Use PCT to predict transcriptional cell-cell distances across datasets.
  alignment = cellhint.harmonize(adata, dataset = 'dataset_column', cell_type = 'celltype_column', use_rep = 'X_pca', use_pct = True)
  ```
  Due to the nonparametric nature of PCT, the format of the expression `.X` in the AnnData is flexible (normalised, log-normalised, z-scaled, etc.), but subsetting the AnnData to highly variable genes is always suggested. To avoid overfitting, each PCT is pruned at nodes where no further splits are needed based on F-test, which is turned on by default (`F_test_prune = True`). You can increase the p-value cutoff (default to 0.05, `p_thres = 0.05`) to prune fewer nodes for improved accuracy at the cost of reduced generalisability.
  </details>

+ <details>
  <summary><strong>1.3. Specify the dataset order</strong></summary>

  In CellHint, datasets are iteratively incorporated and harmonised. The order of datasets can be specified by providing a list of dataset names to the argument `dataset_order`. Otherwise, the order will be determined by CellHint through iteratively adding a dataset that is most similar (i.e., more shared cell types) to the datasets already incorporated. This behaviour can be disabled by setting `reorder_dataset = False` (default to `True`) and an alphabetical order of datasets will be used.
  ```python
  #Specify the order of datasets to be harmonised.
  alignment = cellhint.harmonize(adata, dataset = 'dataset_column', cell_type = 'celltype_column', use_rep = 'X_pca', dataset_order = a_list_of_datasets)
  ```
  </details>

+ <details>
  <summary><strong>1.4. Categories of harmonised cell types</strong></summary>

  Four kinds of harmonisations are anchored with [cellhint.harmonize](https://cellhint.readthedocs.io/en/latest/cellhint.harmonize.html):
     1) Novel cell types as determined by `maximum_novel_percent` (default to `0.05`). In each harmonisation iteration, a cell type (or meta-cell-type) whose maximal alignment fraction is < `maximum_novel_percent` with any cell types in any other datasets is designated as a novel cell type (`NONE`).
     2) One-to-one aligned cell types as determined by `minimum_unique_percents` and `minimum_divide_percents`. If the alignments (in both directions) between two cell types from two respective datasets are greater than `minimum_unique_percents`, plus that these alignments are not one-to-many (see the third point below), this will be an 1:1 (`=`) match. Dynamic thresholds of `minimum_unique_percents` (default to 0.4, 0.5, 0.6, 0.7, 0.8) and `minimum_divide_percents` (default to 0.1, 0.15, 0.2) are exhaustively tested until the least number of alignments is found between datasets.
     3) One-to-many (or many-to-one) aligned cell types as determined by `minimum_unique_percents` and `minimum_divide_percents`. If one cell type has more than two cell types aligned in the other dataset with a match proportion greater than `minimum_divide_percents`, and these matched cell types have a back-match proportion greater than `minimum_unique_percents`, this will be an 1:N (`∋`) or N:1 (`∈`) match. Dynamic thresholds of `minimum_unique_percents` (default to 0.4, 0.5, 0.6, 0.7, 0.8) and `minimum_divide_percents` (default to 0.1, 0.15, 0.2) are exhaustively tested until the least number of alignments is found between datasets.
     4) Unharmonised cell types. If after the above categorisation, a cell type remains unharmonised, then this cell type will be an unharmonised cell type (`UNRESOLVED`).  
    
  |If there are many datasets to harmonise and each dataset has many cell types, harmonisation may take longer time. You can restrict the test scope of `minimum_unique_percents` and `minimum_divide_percents` to reduce runtime. The default is a 15 (5X3) combo test; setting the two parameters to, for example a 3X2 combo, can decrease 60% of the runtime.|
  |:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
  ```python
  #`minimum_unique_percents` is set to three values (default is 0.4, 0.5, 0.6, 0.7, 0.8).
  #`minimum_divide_percents` is set to two values (default is 0.1, 0.15, 0.2).
  alignment = cellhint.harmonize(adata, dataset = 'dataset_column', cell_type = 'celltype_column', use_rep = 'X_pca', minimum_unique_percents = [0.5, 0.6, 0.7], minimum_divide_percents = [0.1, 0.15])
  ```
  </details>
</details>

<details>
<summary><strong>2. Inspection of the harmonisation result</strong></summary>

+ <details>
  <summary><strong>2.1. Harmonisation table</strong></summary>

  The previously saved harmonisation object can be loaded using `cellhint.DistanceAlignment.load`.
  ```python
  alignment = cellhint.DistanceAlignment.load('/path/to/local/folder/some_name.pkl')
  ```
  In `alignment`, the harmonisation table, which summarises cell types across datasets into semantically connected ones, is stored as the attribute `.relation` (`alignment.relation`). One illustrative example is:
  <div align="center">

  |D1   |relation|D2   |relation|D3        |
  |:---:|:---:   |:---:|:---:   |:---:     |
  |A    |=       |B    |=       |C         |
  |D    |=       |NONE |=       |UNRESOLVED|
  |E    |∈       |G    |=       |H         |
  |F    |∈       |G    |=       |I         |
  |J    |=       |K    |∋       |L         |
  |J    |=       |K    |∋       |M         |
  </div>

  The table columns are the dataset1 name, relation, dataset2 name, ..., all the way to the name of the last dataset. Accordingly, each row of the table is a list of cell types connected by predefined symbols of `=`, `∈`, and `∋`. In addition to cell type names, there are two extra definitions of `NONE` and `UNRESOLVED` in the table, representing two levels of novelties (see `1.4.`).  
    
  The table should be interpreted from left to right. For example, for the first row `A = B = C`, although it may look like an 1:1 match between A and B plus an 1:1 match between B and C, a correct interpretation should be an 1:1 match between A and B, resulting in a meta cell type of `A = B`. This meta cell type, as a whole, has an 1:1 match with C, further leading to `A = B = C`. Similarly, for the second row `D = NONE = UNRESOLVED`, instead of a novel cell type D in dataset1, this cell type should be read as a dataset1-specific cell type not existing in dataset2 (`D = NONE`), which as a whole is unharmonised when aligning with dataset3 (`D = NONE = UNRESOLVED`).  
    
  Extending this interpretation to the third and fourth rows, they denote two cell types (E and F) in dataset1 collectively constituting the cell type G in dataset2. The resulting subtypes (`E ∈ G` and `F ∈ G`) are 1:1 matched with H and I in dataset3, respectively. For the last two rows, they describe the subdivision of a meta cell type (`J = K`) into L and M in dataset3, being more than a subdivision of K.  
    
  In the table, each row corresponds to a harmonised low-hierarchy cell type, in other words, the most fine-grained level of annotation that can be achieved by automatic alignment. At a high hierarchy, some cell types such as `E ∈ G = H` and `F ∈ G = I` belong to the same group. CellHint defines a high-hierarchy cell type as fully connected rows in the harmonisation table. As a result, each high-hierarchy cell type is a cell type group independent of each other. This information can be accessed in the attribute `.groups` which is an array/vector with an length of the number of rows in the harmonisation table.
  ```python
  #Access the high-hierarchy cell types (cell type groups).
  alignment.groups
  ```
  </details>

+ <details>
  <summary><strong>2.2. Cell reannotation</strong></summary>

  After cell type harmonisation, each cell can be assigned a cell type label corresponding to a given row of the harmonisation table, denoted as the process of cell reannotation. By default, reannotation is enabled (`reannotate = True`) when using [cellhint.harmonize](https://celltypist.readthedocs.io/en/latest/celltypist.harmonize.html) and information of reannotated cell types is already in place as the attribute `.reannotation`.
  ```python
  #Access the cell reannotation information.
  alignment.reannotation
  ```
  This is a data frame with an example shown below. Unless `filter_cells = True` is set (see `1.1.`), all cells in the AnnData will be present in this data frame.
  <div align="center">

  |     |dataset|cell_type|reannotation         |group |
  |:---:|:---:  |:---:    |:---:                |:---: |
  |cell1|D1     |A        |A = B = C            |Group1|
  |cell2|D1     |D        |D = NONE = UNRESOLVED|Group2|
  |cell3|D2     |G        |E ∈ G = H            |Group3|
  |cell4|D2     |G        |F ∈ G = I            |Group3|
  |cell5|D3     |L        |J = K ∋ L            |Group4|
  |cell6|D3     |M        |J = K ∋ M            |Group4|
  </div>

  The four columns represent information of dataset origin, original author annotation, reannotated low- and high-hierarchy annotation, respectively. For the last column, it contains grouping (high-hierarchy) information, and each group corresponds to a subset of the harmonisation table. You can check this correspondence by coupling the table (`alignment.relation`) with the grouping (`alignment.groups`) (see `2.1.`).
  </details>

+ <details>
  <summary><strong>2.3. Meta-analysis</strong></summary>

  A distance matrix-like instance, which is from the class [Distance](https://celltypist.readthedocs.io/en/latest/celltypist.contro.distance.Distance.html) as defined by CellHint, is also stashed in `alignment` as the attribute `.base_distance`.
  ```python
  #Access the distance object.
  alignment.base_distance
  ```
  The main content of this object is the distance matrix (`alignment.base_distance.dist_mat`) between all cells (rows) and all cell types (columns). Values in this matrix are either calculated (the default) or inferred (if `use_pct` is `True`) by `cellhint.harmonize`, and after a normalisation procedure, lie between 0 and 1. If there are strong cross-dataset batches, an inferred distance matrix obtained from the PCT algorithm is usually more accurate. Metadata of cells and cell types for this matrix can be found in `alignment.base_distance.cell` and `alignment.base_distance.cell_type`, which record raw information such as the dataset origin and original author annotation.  
    
  During the internal harmonisation process, each cell is assigned the most similar cell type from each dataset. This result is stored in the assignment matrix (`alignment.base_distance.assignment`), with rows being cells (cell metadata can be found in `alignment.base_distance.cell` as mentioned above), columns being datasets, and elements being the assigned cell types in different datasets. This matrix can be interpreted as a summary of multi-data label transfers.
  ```python
  #Access the cell type assignment result.
  alignment.base_distance.assignment
  ```
  Each column (corresponding to one dataset) of the assignment matrix can be thought as a unified naming schema when all cells are named by this given dataset.  
    
  CellHint provides a quick way to summarise the above information including cells' distances and assignments into meta-analysis at the cell type level. Specifically, a distance matrix among all cell types can be obtained by:
  ```python
  #Get the cell-type-to-cell-type distance matrix.
  alignment.base_distance.to_meta()
  ```
  An optional `turn_binary = True` (default to `False`) can be added to turn the distance matrix into a cell membership matrix before meta-analysis, showing how cell types are assigned across datasets.
  ```python
  #Get the cell-type-to-cell-type membership matrix.
  alignment.base_distance.to_meta(turn_binary = True)
  ```
  </details>
</details>

<details>
<summary><strong>3. Reharmonisation</strong></summary>

+ <details>
  <summary><strong>3.1. Change the dataset order</strong></summary>

  The order of datasets used by `cellhint.harmonize` can be found in the attribute `.dataset_order` (`alignment.dataset_order`), which is either auto-determined by CellHint or specified by the user (via the `dataset_order` parameter in `cellhint.harmonize`). This order is also reflected by the column order of the harmonisation table.  
    
  Along the order of datasets, optimal choices of `minimum_unique_percents` and `minimum_divide_percents` (see `1.4.`) in each iteration can be found in `alignment.minimum_unique_percents` and `alignment.minimum_divide_percents`. For instance, harmonising five datasets requires four iterations, and thus both `.minimum_unique_percents` and `.minimum_divide_percents` have a length of four.  
    
  CellHint provides a method [best_align](https://celltypist.readthedocs.io/en/latest/celltypist.contro.align.DistanceAlignment.html#celltypist.contro.align.DistanceAlignment.best_align) to change the order of datasets post-harmonisation. Through this, datasets will be reharmonised in a different order (this post-harmonisation adjustment is more efficient than re-running `cellhint.harmonize` with a new order).
  ```python
  #Reharmonise cell types across datasets with a different dataset order.
  alignment.best_align(dataset_order = a_list_of_new_dataset_order)
  ```
  As in `cellhint.harmonize`, the combos of `minimum_unique_percents` and `minimum_divide_percents` will be tested to find the best alignment in each iteration. Importantly, as well as a full dataset list, you can provide a subset of datasets for reharmonisation. This is useful in terms of focusing on part of the data for inspection or visualisation (see `4.`).
  ```python
  #Reharmonise cell types across datasets with part of datasets.
  alignment.best_align(dataset_order = a_subset_of_dataset_names)
  ```
  A new harmonisation table will be generated in `alignment.relation`, which only includes datasets specified in `.best_align`. `.minimum_unique_percents` and `.minimum_divide_percents` are also overridden by new values used during reharmonisation.
  </details>

+ <details>
  <summary><strong>3.2. Reannotation</strong></summary>

  After changing the dataset order and reharmonising cell types, cells need to be reannotated based on the newly generated harmonisation table using the method [reannotate](https://celltypist.readthedocs.io/en/latest/celltypist.contro.align.DistanceAlignment.html#celltypist.contro.align.DistanceAlignment.reannotate).
  ```python
  #Reannotate cells based on the new harmonisation table.
  alignment.reannotate()
  ```
  Similarly, information of reannotated cells is stored in `alignment.reannotation`.
  </details>
</details>

<details>
<summary><strong>4. Visualisation</strong></summary>

+ <details>
  <summary><strong>4.1. Tree plot</strong></summary>

  The most intuitive way to visualise the harmonised cell types is the tree plot using the function [cellhint.treeplot](https://cellhint.readthedocs.io/en/latest/cellhint.treeplot.html).
  ```python
  #Visualise the harmonisation result with a tree plot.
  cellhint.treeplot(alignment)
  ```
  Alternatively, since only the harmonisation table (`alignment.relation`) is used when plotting this tree, `cellhint.treeplot` also accepts the input directly from the table. This is more convenient as a table is easier to manipulate, such as writing it out as a csv file and loading it later for tree plot.
  ```python
  #Write out the harmonisation table as a csv file.
  #Note - if cell type names contain commas, set a different `sep` here.
  alignment.relation.to_csv('/path/to/local/folder/HT.csv', sep = ',', index = False)
  ```
  ```python
  #Read the harmonisation table.
  HT = pd.read_csv('/path/to/local/folder/HT.csv', sep = ',')
  #Visualise the harmonisation result with a tree plot.
  cellhint.treeplot(HT)
  #Visualise the harmonisation result only for cell types (rows) of interest.
  cellhint.treeplot(HT[row_flag])
  ```
  In a tree plot, each column is a dataset and cell types are connected across datasets. By default, cell types belonging to one low hierarchy (one row in the harmonisation table) are in the same color. You can change the color scheme by providing a data frame to the `node_color` parameter, with three consecutive columns representing dataset, cell type, and color (in hex code), respectively. `node_color` can also be a data frame with columns of dataset, cell type, and numeric value (for mapping color gradient in combination with `cmap`). Other parameters controlling the appearance of the tree plot (node shape, line width, label size, figure size, etc.) are detailed in [cellhint.treeplot](https://celltypist.readthedocs.io/en/latest/celltypist.treeplot.html).  
  |The tree plot considers all pairs of reference-to-query assignments. Therefore, a restricted representation in two dimensionalities may overlay some cell types when they have complex 1:1 and 1:N intersections. These cross-connections are usually not solvable at 2D space; you may need to revisit the harmonisation table in some cases.|
  |:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
    
  By changing the dataset (column) order in each high-hierarchy cell type, broader (more divisible) cell types can be positioned to the left, followed by fine-grained cell types to the right. The resulting plot shows how different authors group these cell types, thereby being more characteristic of the potential underlying biological hierarchy. This hierarchy can be generated and visualised by adding `order_dataset = True`.
  ```python
  #Visualise the cell type hierarchy.
  #Again, the input can also be a harmonisation table.
  cellhint.treeplot(alignment, order_dataset = True)
  ```
  Because each high-hierarchy cell type is independent of each other, the new orders of datasets will be different across groups. To recognise the dataset origin of each cell type within the hierarchy, you can assign the same color or shape to cell types from the same dataset using the parameter `node_color` or `node_shape`. An example is:
  ```python
  #Cell types from the same dataset are in the same shape.
  #`node_shape` should be the same length as no. datasets in the harmonisation table.
  cellhint.treeplot(alignment, order_dataset = True, node_shape = list_of_shapes)
  ```
  Export the plot if needed.
  ```python
  cellhint.treeplot(alignment, show = False, save = '/path/to/local/folder/some_name.pdf')
  ```
  </details>

+ <details>
  <summary><strong>4.2. Sankey plot</strong></summary>

  The other way to visualise harmonised cell types is the Sankey plot by [cellhint.sankeyplot](https://celltypist.readthedocs.io/en/latest/celltypist.sankeyplot.html). CellHint builds this plot on the [plotly](https://pypi.org/project/plotly) package. `plotly` is not mandatory when installing CellHint, so you need to install it first if you want a visualisation form of Sankey diagram (and engines for exporting images such as [kaleido](https://pypi.org/project/kaleido)).
  ```python
  #Visualise the harmonisation result with a Sankey plot.
  #As with the tree plot, the input can also be a harmonisation table.
  cellhint.sankeyplot(alignment)
  ```
  Similar to the tree plot, this diagram shows how cell types are connected across datasets. Parameters controlling the appearance of the Sankey plot (node color, link color, figure size, etc.) are detailed in [cellhint.sankeyplot](https://celltypist.readthedocs.io/en/latest/celltypist.sankeyplot.html).  
    
  Different from the tree plot where novel (`NONE`) and unharmonised (`UNRESOLVED`) cell types are blank, in the Sankey plot they are colored in white and light grey, respectively. You can adjust these by changing the values of `novel_node_color` and `remain_node_color`.  
    
  Export the plot if needed.
  ```python
  #Export the image into html.
  cellhint.sankeyplot(alignment, show = False, save = '/path/to/local/folder/some_name.html')
  #Export the image into pdf.
  cellhint.sankeyplot(alignment, show = False, save = '/path/to/local/folder/some_name.pdf')
  ```
  </details>
</details>

# Usage (integration)

<details>
<summary><strong>1. Supervised data integration</strong></summary>

+ <details>
  <summary><strong>1.1. Specify batch and biological covariates</strong></summary>

  The input [AnnData](https://anndata.readthedocs.io/en/latest/) needs two columns in `.obs` representing the batch confounder and unified cell annotation respectively. The aim is to integrate cells by correcting batches and preserving biology (cell annotation) using [cellhint.integrate](https://celltypist.readthedocs.io/en/latest/celltypist.integrate.html).
  ```python
  #Integrate cells with `cellhint.integrate`.
  cellhint.integrate(adata, batch = 'a_batch_key', cell_type = 'a_celltype_key')
  ```
  With this function, CellHint will build the neighborhood graph by searching neighbors across matched cell groups in different batches, on the basis of a low-dimensional representation provided via the argument `use_rep` (default to PCA coordinates).
  ```python
  #`use_rep` can be omitted here as it defaults to 'X_pca'.
  cellhint.integrate(adata, batch = 'a_batch_key', cell_type = 'a_celltype_key', use_rep = 'X_pca')
  ```
  The batch confounder can be the dataset origin, donor ID, or any relevant covariate. For the biological factor, it is the consistent annotation across cells, such as manual annotations of all cells, transferred cell type labels from a single reference model, and as an example here, the harmonised cell types from the CellHint harmonisation pipeline (see the harmonisation section). Specifically, you can add two extra columns in the `.obs` of the input AnnData using the reannotation information from `alignment.reannotation`.
  ```python
  #Insert low- and high-hierarchy annotations into the AnnData.
  adata.obs[['harmonized_low', 'harmonized_high']] = alignment.reannotation.loc[adata.obs_names, ['reannotation', 'group']]
  ```
  Perform data integration using either of the two annotation columns.
  ```python
  #Integrate cells using the reannotated high-hierarchy cell annotation.
  cellhint.integrate(adata, batch = 'a_batch_key', cell_type = 'harmonized_high')
  #Not run; integrate cells using the reannotated low-hierarchy cell annotation.
  #cellhint.integrate(adata, batch = 'a_batch_key', cell_type = 'harmonized_low')
  ```
  Finally, generate a UMAP based on the reconstructed neighborhood graph.
  ```python
  sc.tl.umap(adata)
  ```
  </details>

+ <details>
  <summary><strong>1.2. Adjust the influence of annotation on integration</strong></summary>

  Influence of cell annotation on the data structure can range from forcibly merging the same cell types to a more lenient cell grouping. This is achieved by adjusting the parameter `n_meta_neighbors`.
  ```python
  #Actually the default value of `n_meta_neighbors` is 3.
  cellhint.integrate(adata, batch = 'a_batch_key', cell_type = 'a_celltype_key', n_meta_neighbors = 3)
  ```
  With `n_meta_neighbors` of 1, each cell type only has one neighboring cell type, that is, itself. This will result in strongly separated cell types in the final UMAP. Increasing `n_meta_neighbors` will loosen this restriction. For example, a `n_meta_neighbors` of 2 allows each cell type to have, in addition to itself, one nearest neighboring cell type based on the transcriptomic distances calculated by CellHint. This parameter defaults to 3, meaning that a linear spectrum of transcriptomic structure can possibly exist for each cell type.
  </details>
</details>

<details>
<summary><strong>2. Tips for data integration</strong></summary>

+ <details>
  <summary><strong>2.1. Partial annotation</strong></summary>

  Partial annotation (an `.obs` column combining annotated and unannotated cells) is allowed as the `cell_type` parameter of `cellhint.integrate`. You need to explicitly name unannotated cells as `'UNASSIGNED'` for use in CellHint (definition of symbols can be found [here](https://github.com/Teichlab/celltypist/blob/main/celltypist/contro/symbols.py)).
  </details>

+ <details>
  <summary><strong>2.2. Rare cell types</strong></summary>

  When an abundant cell type is annotated/distributed across multiple batches (e.g., datasets), sometimes not all batches can harbour adequate numbers. This leads to a rare cell type defined within the context of a specific batch. During neighborhood construction, if this batch cannot provide enough neighboring cells for this cell type, search space will be expanded to all cells in this batch.  
    
  Although this represents a safe solution in CellHint to anchor nearest neighbors for rare cell types, runtime of the algorithm will be increased and cells from this cell type may not be robustly clustered. Keeping them is fine for CellHint, but you can also remove such rare cell types in associated batches before running `celltypist.integrate` (a cell type with only a small number in a given batch naturally means that this batch may not be qualified for hosting this cell type). Example code is:
  ```python
  #Remove cells from cell types that have <=5 cells in a batch.
  combined = adata.obs['a_batch_key'].astype(str) + adata.obs['a_celltype_key'].astype(str)
  combined_counts = combined.value_counts()
  remove_combn = combined_counts.index[combined_counts <= 5]
  adata = adata[~combined.isin(remove_combn)].copy()
  ```
  </details>

+ <details>
  <summary><strong>2.3. Use CellTypist models for annotation and integration</strong></summary>

  `cellhint.integrate` requires cell annotation to be stored in the AnnData. This information can be obtained by different means. One quick way is to use available CellTypist models to annotate the data of interest (see the CellTypist model list [here](https://www.celltypist.org/models)).
  ```python
  #Annotate the data with a relevant model (immune model as an example here).
  adata = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True).to_adata()
  ```
  Then integrate cells on the basis of the predicted cell types.
  ```python
  #`cell_type` can also be 'majority_voting'.
  cellhint.integrate(adata, batch = 'a_batch_key', cell_type = 'predicted_labels')
  ```
  Even the model does not exactly match the data (e.g., using an immune model to annotate a lung data), this approach can be still useful as cells from the same cell type will probably be assigned the same identity by the model, therefore containing information with respect to which cells should be placed together in the neighborhood graph.
  </details>
</details>

# Citation
Xu et al., Automatic cell type harmonization and integration across Human Cell Atlas datasets. bioRxiv (2023). [Preprint](https://doi.org/10.1101/2023.05.01.538994)  
