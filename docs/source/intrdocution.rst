Intrdocution
==============

.. _installation:

Installation
------------

The latest stable release is available on PyPi (and on `GitHub`_), and you can install it by saying

.. code-block:: console

   (.venv) $ pip install bacteria

To build Bacteria from source, say ``python setup.py build``.
Then, to install Bacteria, say ``python setup.py install``.
If all went well, you should be able to execute the demo scripts (read docs)
(OS X users should follow the installation guide given below).

Alternatively, you can download or clone the repository and use ``pip`` to handle dependencies:

.. code-block:: console

   (.venv) $ unzip bacteria.zip
   (.venv) $ pip install -e bacteria

or

.. code-block:: console

   (.venv) $ git clone https://github.com/tuliofalmeida/bacteria
   (.venv) $ pip install -e bacteria

By calling ``pip list`` you should see ``bacteria`` now as an installed package:
``bacteria (1.x.x, /path/to/bacteria)``

.. _IO:

IO
---

The library's input is the SuperSegger data, the clist, and the organization used in 
pandas DataFrame. Its possible to load one clist and extract the 2D and 3D data and 
also concatenate clists from different fields of view into a single dataframe.

Single data
-----------

To load data from a single field of view

    >>> import bacteria as bac
    >>> path_2d = "data/2d_clist.mat"
    >>> path_3d = "data/3d_clist.mat"    
    >>> df2d = bac.df_data2d(path_2d,fps = 3)
    >>> df3d = bac.df_data3d(path_3d,fps = 3)

.. note::
   The fps parameter is to create the time column in minutes.

Concatenate data
----------------

To load and concatenate multiple field of view data from a single experiment,
if the data is organized inside the SuperSegger folder, use:

    >>> path = "/SuperSegger Data/Exp001/"
    >>> df2d,df3d,_ = bac.concatenate_clist(path)

.. note::
   By concatenating the data, it is possible to identify the cells of each FOV and the 
   Cell ID is modified to remain unique, the IDs are a continuous count at each FOV.

Or you can pass it the path to a folder with all the clists of the experiment

    >>> path_clists = "/SuperSegger Data/clists folder/"
    >>> df2d,df3d,_ = bac.concatenate_clist(path_clists, direct = True)

.. note::
   Is it possible to save this arrangement in a .mat file (this arrangement is different
   from the concatenation done by SuperSegger)

Filters
--------

The library has 2 filters implemented along with the SuperSegger filter. The data is 
filtered using the "stat0" column of SupperSegger (applied in the functions of reading 
data directly from the .mat file), size filter (excludes very long cells) and the mother
and daughter filter, which will exclude cells that were born smaller than 40% the mother's
size or greater than 60% of the mother's size.

    >>> df2d_f,df3d_f,_ = bac.combined_filters(df2d,df3d)

.. note::
   The function returns a histogram with the effect of each filter on the data and it is
   possible to adjust some parameters and select whether all filters will be applied.

Hint!
-----

After loading, concatenating and filtering the data, it is interesting to save them in a 
".csv" file to optimize the analysis, leaving the original data intact and making data 
loading faster. For that use:

    >>> import pandas as pd
    >>> df2d_f.to_csv('/data/"2D_filtered.csv",index=False)
    >>> df3d_f.to_csv('/data/"3D_filtered.csv",index=False)

Colab Tutorial
---------------


Check out the tutorial on Google Colab!

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/tuliofalmeida/bacteria/blob/main/notebooks/Tutorial_Concatenate_Filters.ipynb

.. _GitHub: https://github.com/tuliofalmeida/bacteria

