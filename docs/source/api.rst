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
    
    >>> # we use '!' to execute linux commands in notebooks
    >>> !git clone https://github.com/tuliofalmeida/bacteria
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

    >>> time_gr,mean_gr,ci_gr = bac.column_mean(gr_df,column='Derivative')
    >>> time_fluo,mean_fluo,ci_fluo = bac.column_mean(fluo_df,column='Derivative')
    >>> time_div,mean_div,ci_div = bac.column_mean(df3d,column='Long axis (L)')
    >>> time_ratio,mean_ratio,ci_ratio = bac.column_mean(df3d,column='F/V')

Rearranging 2D data based on division time to plot the data over time

    >>> df_2d_temp = pd.melt(df2d, id_vars=['Cell ID','Cell age','Area death','Time Division'], value_vars=['Vd-Vb']).sort_values(by=['Time Division'])

Entire time plot

.. note::
   The plot codes are available on Google Colab.

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/experiment_time_1.png?raw=true

Pre and pos shift plots

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/experiment_time_2.png?raw=true

Transition plots

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/experiment_time_3.png?raw=true

Cell Cycle Analysis
-------------------

First we will extract the Cell ID's pre and pos shift

    >>> cells_pre = bac.cells_pre_shift(df3d,pre=600)
    >>> cells_pos = bac.cells_pos_shift(df3d,pos=800)

Bin the data for all the entire time

    >>> # Volume derivative 
    >>> bins_fluo_vol_norm,ci_fluo_vol_norm = bac.derivative_binning(fluo_df,derivative_column='Derivative/V',sort_by='Volume',print_bins=True)
    >>> bins_fluo_vol,ci_fluo_vol = bac.derivative_binning(fluo_df,derivative_column='Derivative',sort_by='Volume',print_bins=False)                                                            
    >>> # Fluorescence derivative
    >>> bins_fluo_cc_norm,ci_fluo_cc_norm = bac.derivative_binning(fluo_df,derivative_column='Derivative/V',sort_by='Cell Cycle',print_bins=False)
    >>> bins_fluo_cc,ci_fluo_cc = bac.derivative_binning(fluo_df,derivative_column='Derivative',sort_by='Cell Cycle',print_bins=False)
    >>> # Fluorescence/Volume Ratio
    >>> bins_ratio_vol, ci_ratio_vol = bac.bin_column(df3d,column = 'F/V',sort_by = 'Volume',print_bins=False)
    >>> bins_ratio_cc, ci_ratio_cc = bac.bin_column(df3d,column = 'F/V',sort_by = 'Cell Cycle',print_bins=False)

Bin plots entire time

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/cell_cycle_1.png?raw=true

Pre plots

    >>> # Volume derivative pre
    >>> bins_fluo_vol_norm,ci_fluo_vol_norm = bac.derivative_binning(fluo_df[fluo_df['Cell ID'].isin(cells_pre)],derivative_column='Derivative/V',sort_by='Volume',print_bins=False)
    >>> bins_fluo_vol,ci_fluo_vol = bac.derivative_binning(fluo_df[fluo_df['Cell ID'].isin(cells_pre)],derivative_column='Derivative',sort_by='Volume',print_bins=False)                                                         
    >>> # Fluorescence derivative pre
    >>> bins_fluo_cc_norm,ci_fluo_cc_norm = bac.derivative_binning(fluo_df[fluo_df['Cell ID'].isin(cells_pre)],derivative_column='Derivative/V',sort_by='Cell Cycle',print_bins=False)
    >>> bins_fluo_cc,ci_fluo_cc = bac.derivative_binning(fluo_df[fluo_df['Cell ID'].isin(cells_pre)],derivative_column='Derivative',sort_by='Cell Cycle',print_bins=False)
    >>> # Fluorescence/Volume Ratio pre
    >>> bins_ratio_vol, ci_ratio_vol = bac.bin_column(df3d[df3d['Cell ID'].isin(cells_pre)],column = 'F/V',sort_by = 'Volume',print_bins=False)
    >>> bins_ratio_cc, ci_ratio_cc = bac.bin_column(df3d[df3d['Cell ID'].isin(cells_pre)],column = 'F/V',sort_by = 'Cell Cycle',print_bins=False)

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/cell_cycle_2.png?raw=true

Pos plots

    >>> # Volume derivative pos
    >>> bins_fluo_vol_norm,ci_fluo_vol_norm = bac.derivative_binning(fluo_df[fluo_df['Cell ID'].isin(cells_pos)],derivative_column='Derivative/V',sort_by='Volume',print_bins=False)
    >>> bins_fluo_vol,ci_fluo_vol = bac.derivative_binning(fluo_df[fluo_df['Cell ID'].isin(cells_pos)],derivative_column='Derivative',sort_by='Volume',print_bins=False)                                                         
    >>> # Fluorescence derivative pos
    >>> bins_fluo_cc_norm,ci_fluo_cc_norm = bac.derivative_binning(fluo_df[fluo_df['Cell ID'].isin(cells_pos)],derivative_column='Derivative/V',sort_by='Cell Cycle',print_bins=False)
    >>> bins_fluo_cc,ci_fluo_cc = bac.derivative_binning(fluo_df[fluo_df['Cell ID'].isin(cells_pos)],derivative_column='Derivative',sort_by='Cell Cycle',print_bins=False)
    >>> # Fluorescence/Volume Ratio pos
    >>> bins_ratio_vol, ci_ratio_vol = bac.bin_column(df3d[df3d['Cell ID'].isin(cells_pos)],column = 'F/V',sort_by = 'Volume',print_bins=False)
    >>> bins_ratio_cc, ci_ratio_cc = bac.bin_column(df3d[df3d['Cell ID'].isin(cells_pos)],column = 'F/V',sort_by = 'Cell Cycle',print_bins=False)

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/cell_cycle_3.png?raw=true


Lineages Analysis
-----------------

Calculate all lineages in the dataset

    >>> reverse_lineages, large_lineages, large_nodes , all_nodes = bac.lineages(df2d)

Extract the longer lineages in dict and array

    >>> filtered_lineage = bac.lineage_corr_filter(df3d,reverse_lineages)

Extract lineages start pre and pos the shift

    >>> lineage_pre = bac.lineages_pre_shift(df3d,reverse_lineages, shift = 700)
    >>> lineage_pos = bac.lineages_pos_shift(df3d,reverse_lineages, shift = 600)

Plot one cell mother lineage and the correlation between the nodes

    >>> bac.plot_corr_lineage(df3d,reverse_lineages,filtered_lineage[5][0])

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/lineage_1.png?raw=true

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/lineage_1_1.png?raw=true

The order is a paramater that will influence the amount of noise accepet to estimate the smoothed derivative
Fluorescence concentration pre shift (using order = 5)

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/lineage_2.png?raw=true

Fluorescence concentration pos shift (using order = 5)

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/lineage_3.png?raw=true

Fluorescence rate pre shift (using order = 5)

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/lineage_4.png?raw=true

Fluorescence rate pos shift (using order = 5)

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/lineage_5.png?raw=true


Minima Analysis
---------------

    >>> diff_dict = bac.diff_minima(df3d,df2d,pre = 600,pos = 800,order = 3, adjust = False)

Plot Minima Volume by Volume at birth

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/minima_1.png?raw=true

Plot Daughther Minima Volume - Mother Minima Volume

    >>> diff_dict = bac.plot_diff_minima(diff_dict,color='purple')

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/minima_2.png?raw=true

Plot minima lineage by volume (order = 3)

.. image:: https://github.com/tuliofalmeida/bacteria/blob/main/plots/minima_3.png?raw=true


Colab Tutorial
---------------

All codes and plots are available on a notebook at Google Colab!

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/tuliofalmeida/bacteria/blob/main/notebooks/Tutorial_Analysis.ipynb

.. _GitHub: https://github.com/tuliofalmeida/bacteria