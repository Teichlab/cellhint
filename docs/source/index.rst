.. cellhint documentation master file, created by
   sphinx-quickstart on Fri Aug  4 16:25:59 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to CellHint's documentation!
======================================

.. include:: ../../README.md
   :parser: myst_parser.sphinx_
   :start-line: 2
   :end-line: 7

.. include:: ../../README.md
   :parser: myst_parser.sphinx_
   :start-line: 14
   :end-line: 24

.. include:: ../../CHANGELOG.md
   :parser: myst_parser.sphinx_

.. toctree::
   :maxdepth: 2
   :caption: Tutorials:
   :hidden:

   notebook/cellhint_tutorial_harmonisation
   notebook/cellhint_tutorial_integration

.. toctree::
   :maxdepth: 2
   :caption: API (harmonisation):
   :hidden:

   cellhint.harmonize
   cellhint.treeplot
   cellhint.sankeyplot

.. toctree::
   :maxdepth: 2
   :caption: API (integration):
   :hidden:

   cellhint.integrate

.. toctree::
   :maxdepth: 2
   :caption: Package organization:
   :hidden:

   cellhint.align.DistanceAlignment
   cellhint.distance.Distance
   cellhint.pct.PredictiveClusteringTree
