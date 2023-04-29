API
===

This section has the code examples to do a complete analysis using the library.

Imports and Data
-------------------

For execution we will use the following imports:

    >>> import bacteria as bac
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt

We will use the example data available on GitHub:

    >>> !git clone https://github.com/tuliofalmeida/bacteria    # we use '!' to execute linux commands in notebooks
    >>> df2d = pd.read_csv('/content/bacteria/data_tutorial/df2d_filtered.csv')
    >>> df3d = pd.read_csv('/content/bacteria/data_tutorial/df3d_filtered.csv')

.. note::
   This data is not representative about any experiment, it has a shift in the growth medium.

Experiment in time
-------------------

Calculating the instantaneous volume change (Growth Rate) and fluorescence.

    >>> gr_df,_ = bac.derivative(df3d,column='Volume')
    >>> fluo_df,_ = bac.derivative(df3d,column='Fluor1 sum')

Calculating the average of all cells for different columns from derivative and 3D data

    >>> time_gr,mean_gr,ci_gr = bac.column_mean(gr_df,'Derivative')
    >>> time_fluo,mean_fluo,ci_fluo = bac.column_mean(fluo_df,'Derivative')
    >>> time_div,mean_div,ci_div = bac.column_mean(df3d,'Long axis (L)')
    >>> time_ratio,mean_ratio,ci_ratio = bac.column_mean(df3d,'F/V')

Rearranging 2D data based on division time to plot the data over time

    >>> df_2d_temp = pd.melt(df2d, id_vars=['Cell ID','Cell age','Area death','Time Division'], value_vars=['Vd-Vb']).sort_values(by=['Time Division'])

Entire time plot

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/experiment_time_1.png

Cell Cycle Analysis
-------------------


Lineages Analysis
-----------------


Minima Analysis
---------------


Colab Tutorial
---------------

Check out the tutorial on Google Colab!

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/tuliofalmeida/bacteria/blob/main/notebooks/Tutorial_Concatenate_Filters.ipynb

.. _GitHub: https://github.com/tuliofalmeida/bacteria