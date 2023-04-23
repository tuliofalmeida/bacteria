import os
import sys
import time
import shutil
import graphviz
import scipy.io
# import matlab.engine
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# from IPython import display
from scipy.stats import sem
from natsort import natsorted # pip install natsort
from datetime import date,timedelta
from matplotlib.lines import Line2D
from sklearn import preprocessing as pre
from sklearn.linear_model import LinearRegression
from scipy.signal import savgol_filter, argrelextrema
from anytree import Node, RenderTree, PostOrderIter # pip install anytree


def interactive_table(df_list, cell = None, labels = None):
    """
    Creates an interactive table with statistics using a DataFrame as input. 
    
    Parameters
    --------------
    df_list: list
        list of DataFrame
    cell: int (optional)
        The number of the specific cell (just for 3D data)
    labels: list
        list with the label for each DataFrame
        
    Returns
    --------------
    None
    
    Requirements
    --------------
    facets-overview==1.0.0
    https://stackoverflow.com/questions/71759248/importerror-cannot-import-name-builder-from-google-protobuf-internal
    
    """
    import pandas as pd
    from facets_overview.feature_statistics_generator import FeatureStatisticsGenerator
    from IPython.core.display import display, HTML
    import base64

    fsg = FeatureStatisticsGenerator()
    dataframes = []

    if isinstance(df_list, list) == False: raise TypeError("input data must be a list")
    if labels != None and len(labels) != len(df_list): raise TypeError("labels length must be equal to df_list")
    if cell is None:
        for idx,df in enumerate(df_list):
            if isinstance(df, pd.DataFrame) == False: raise TypeError("input data must be a pandas DataFrame from df_data2d() or df_data3d()")
            else:
                if labels is None:
                    dataframes.append({'table': df, 'name': 'Data_'+str(idx)})
                else:
                    dataframes.append({'table': df, 'name': labels[idx]})
    else:
        for idx,df in enumerate(df_list):    
            if df.shape[1] != 15: raise TypeError("input data must be a pandas DataFrame from df_data3d()")
            else:
                if labels is None:
                    dataframes.append({'table': df, 'name': 'Data_'+str(idx)})
                else:
                    dataframes.append({'table': df, 'name': labels[idx]})
        
    censusProto = fsg.ProtoFromDataFrames(dataframes)
    protostr = base64.b64encode(censusProto.SerializeToString()).decode("utf-8")

    HTML_TEMPLATE = """<script src="https://cdnjs.cloudflare.com/ajax/libs/webcomponentsjs/1.3.3/webcomponents-lite.js"></script>
            <link rel="import" href="https://raw.githubusercontent.com/PAIR-code/facets/1.0.0/facets-dist/facets-jupyter.html">
            <facets-overview id="elem"></facets-overview>
            <script>
            document.querySelector("#elem").protoInput = "{protostr}";
            </script>"""
    html = HTML_TEMPLATE.format(protostr=protostr)
    display(HTML(html))

def _off_set(df):
    """
    Create a pandas DataFrame of the 2D data using the headers pre-determined. 
    
    Parameters
    --------------
    df: pd.DataFrame
        2D or 3D data
        
    Returns
    --------------
    off_sets : list
        list valued to add
    intervals : list
        list with tuple of wich interval ([0] init , [1] end)
    """
    # getting the offsets for each fov in the data 
    count = 0
    add_value = [0]
    idx_off_set = []
    for i in range(len(df)-1):
        if df['Cell ID'][i+1] < df['Cell ID'][i]:
            count += df['Cell ID'][i]
            add_value.append(count) 
            idx_off_set.append(i+1)
    idx_off_set.append(len(df))
    off_sets = dict(zip(idx_off_set, add_value))

    # create the intervals for each fov
    intervals = []
    count = 0
    for idx,_ in off_sets.items():
        intervals.append((count, idx))
        count = idx

    return off_sets,intervals

def _volume(df,const = 0.1067):
    """
    Create a the volume column in 3d DataFrame. 
    v = np.pi*(l-w) * ((w/2)**2) + (4/3) * np.pi * ((w/2)**3)
    
    Parameters
    --------------
    df : DataFrame
        df from 3D data 
        
    Returns
    --------------
    df : pandas DataFrame
        df with the Volume column
    """
    vol_temp = []
    
    for cell in natsorted(set(df['Cell ID'].values)):
        l = df[df['Cell ID']==cell]['Long axis length'].values*const
        w = df[df['Cell ID']==cell]['Short axis'].values*const
        v = np.pi*(l-w) * ((w/2)**2) + (4/3) * np.pi * ((w/2)**3)
        vol_temp.extend(v) # *0.1067

    df['Volume'] = vol_temp

    return df

def _volume2d(df2d,const = 0.1067):
    """
    Create a the volume at birth and division column in 2d DataFrame. 
    v = np.pi*(l-w) * ((w/2)**2) + (4/3) * np.pi * ((w/2)**3)
    
    Parameters
    --------------
    df : DataFrame
        df from 2D data 
        
    Returns
    --------------
    df : pandas DataFrame
        df with the Volume columns
    """
    vol_b,vol_d = [],[]
    
    for cell in natsorted(set(df2d['Cell ID'].values)):
        l = df2d[df2d['Cell ID']==cell]['Long axis (L) birth'].values*const
        w = df2d[df2d['Cell ID']==cell]['Short axis birth'].values*const
        v = np.pi*(l-w) * ((w/2)**2) + (4/3) * np.pi * ((w/2)**3)
        vol_b.extend(v)
        
        l = df2d[df2d['Cell ID']==cell]['Long axis (L) death'].values*const
        w = df2d[df2d['Cell ID']==cell]['Short axis death'].values*const
        v = np.pi*(l-w) * ((w/2)**2) + (4/3) * np.pi * ((w/2)**3)
        vol_d.extend(v)
        
    df2d['Volume birth'] = vol_b
    df2d['Volume division'] = vol_d
    df2d['Vd-Vb'] = df2d['Volume division'].values - df2d['Volume birth'].values

    return df2d

def fluor_volume_ratio(df, good_cells = None):
    """
    Calculate Fluor/Volume ratio.
    
    Parameters
    --------------
    df : DataFrame
        DataFrame of 3D data
    good_cells : List (optional)
        List with cells to use 
        if None good_cells = df['Cell ID'].values
        
    Returns
    --------------
    df : pandas DataFrame
        df with the F/V column
    """
    fv_ratio = []

    if good_cells is None: good_cells = natsorted(set(df['Cell ID'].values)) 

    for cell in good_cells:
        fv_ratio.extend(df[df['Cell ID']==cell]['Fluor1 sum'].values/df[df['Cell ID']==cell]['Volume'].values)

    df['F/V'] = fv_ratio

    return df

def cell_cycle(df, good_cells = None):
    """
    Calculate Fluor/Volume ratio.

    Parameters
    --------------
    df : DataFrame
        DataFrame of 3D data
    good_cells : List (optional)
        List with cells to use 
        if None good_cells = df['Cell ID'].values
        
    Returns
    --------------
    df : pandas DataFrame
        df with the F/V column
    """
    cell_cycles = []

    if good_cells is None: good_cells = natsorted(set(df['Cell ID'].values)) 

    for cell in good_cells:
        time_frames = df[(df['Cell ID'] == cell) & (df['Time (fps)'])]['Time (fps)'].values
        cell_cycles.extend(pre.MinMaxScaler((0,1)).fit_transform(np.arange(0,len(time_frames),1).reshape(-1, 1)).flatten())

    df['Cell Cycle'] = cell_cycles

    return df

def df_data3d(data, fps = 3, filters = True, comp = False):
    """
    Create a pandas DataFrame of the 2D data using the headers pre-determined. 
    
    Parameters
    --------------
    data : clist.mat
        Supersegger clist already loaded
    fps : int (optional)
        Frames per second used in the experiment, usual 3
    filters : bool (optional)
        Filter data based on data2d['stat0'] = 2 and non NaN's
    comp : bool
        Input data are the compiled clists (clist for all xy/fovs)
        
    Returns
    --------------
    df : pandas DataFrame
        3D data in df
    gc : nd.array
        array with good cells ID
    """
    # vector of good cells
    gc = _good_cells(data, comp = comp)

    # creating the header of the dataframe removing the Focus and Fluor2 columns
    header_3d = {idx:column[0] for idx,column in enumerate(data['def3D'][0]) if column[0][:5] != 'Focus' and column[0][:6] != 'Fluor2'} 
    df_temp = []
    # for each cell take the time vector data (for each frames)
    for cell in range(data['data3D'].shape[0]):
        if filters:
            if cell in gc:
                temp = {}
                for column in header_3d:
                    temp[header_3d[column]] = data['data3D'][cell,column]
                    temp['Cell ID'] = cell # using the same numbers than matlab - done in gc
                # add the time vector in fps
                temp['Time (fps)'] = temp['Time (Frames)'] * fps
                temp['Age (fps)'] = temp['Age (Frames)'] * fps
                df_temp.append(pd.DataFrame(temp))     
        else:
            temp = {}
            for column in header_3d:
                temp[header_3d[column]] = data['data3D'][cell,column]
                temp['Cell ID'] = cell+1 # using the same numbers than matlab
            # add the time vector in fps
            temp['Time (fps)'] = temp['Time (Frames)'] * fps
            temp['Age (fps)'] = temp['Age (Frames)'] * fps
            df_temp.append(pd.DataFrame(temp))

    df = pd.concat(df_temp)
    if filters:
        # remove the NaN values
        df = df[df['Age (Frames)'].isna()== False]
    # remove trailing space
    df.columns = df.columns.str.rstrip()

    df = _volume(df)
    
    return df,gc

def df_data2d(data, fps = 3, comp = False):
    """
    Create a pandas DataFrame of the 2D data using the headers pre-determined. 
    
    Parameters
    --------------
    data: clist.mat
        Supersegger clist already loaded
    fps: int (optional)
        Frames per second used in the experiment, usual 3
    comp : bool
        Input data are the compiled clists (clist for all xy/fovs)
        
    Returns
    --------------
    df: pandas DataFrame
    
    """
    # header pre-determined by Bianca
    header = ['Cell ID','Cell birth time','Cell age','stat0','Long axis (L) birth','Long axis (L) death','Short axis birth','Short axis death',
              'Area birth','Area death','Fluor1 sum','Fluor1 mean','Mother ID','Daughter1 ID','Daughter2 ID','L death / L birth',
              'Fluor1 sum death','Fluor1 mean death','Long axis/Short axis birth','Long axis/Short axis death','Growth Rate']

    # creating a full header from clist with the according idx exluding the Focus data
    header_full = {idx:column[0] for idx,column in enumerate(data['def'][0]) if column[0][:5] != 'Focus'}

    # taking the idx from the supersegger for the pre-determined header
    header_idx = {idx:column for idx in header_full for column in header if header_full[idx] == column}

    # creating a dict with the 2D data and the predetermined columns 
    df_dict = {header_idx[column]:data['data'][:,column] for idx, column in enumerate(header_idx)}

    # create a pandas dataframe
    df = pd.DataFrame(df_dict)

    # multiply the frame by the fps
    df['Cell birth time'] = df['Cell birth time'] * fps
    df['Cell age'] = df['Cell age'] * fps

    # add the volume columns
    df = _volume2d(df)

    if comp:
        off_sets , intervals = _off_set(df)
        # **this lasts step are unnecessary because the Cell ID follow the number of points in the data, it is just to assure to find this off sets**

        # update the intervals for each fov with the offsets
        for idx in range(len(intervals)):
            df.loc[intervals[idx][0]:intervals[idx][1]-1,'Cell ID'] = [i+off_sets[intervals[idx][1]] for i in df.iloc[intervals[idx][0]:intervals[idx][1]]['Cell ID']]
            df.loc[intervals[idx][0]:intervals[idx][1]-1,'Daughter1 ID'] = [i+off_sets[intervals[idx][1]] for i in df.iloc[intervals[idx][0]:intervals[idx][1]]['Daughter1 ID']]
            df.loc[intervals[idx][0]:intervals[idx][1]-1,'Daughter2 ID'] = [i+off_sets[intervals[idx][1]] for i in df.iloc[intervals[idx][0]:intervals[idx][1]]['Daughter2 ID']]
            df.loc[intervals[idx][0]:intervals[idx][1]-1,'Mother ID'] = [i+off_sets[intervals[idx][1]] if i != 0 else i for i in df.iloc[intervals[idx][0]:intervals[idx][1]]['Mother ID']]

    return df

def _good_cells(data, fps = 3, comp = False):
    """
    Return a array with ID of good cells. 
    
    Parameters
    --------------
    data: clist.mat
        Supersegger clist already loaded
    fps: int (optional)
        Frames per second used in the experiment, usual 3
    comp : bool
        Input data are the compiled clists (clist for all xy/fovs)
        
    Returns
    --------------
    good_cells: array
    
    """
    # use the 2d data to get the ID of good cells
    data2d = df_data2d(data,fps,comp)
    
    # vector of good cells
    good_cells = data2d[data2d['stat0']==2]['Cell ID'].values

    return good_cells

def concatenate_clist(path, fps = 3, filters = True, save_mat = False, path2save = None , direct = False):
    """
    Concatenate all c_list.mat from one experiment in DataFrames.
    
    Parameters
    --------------
    path : str
        Experiment path
    fps : int (optional)
        Frames per second used in the experiment, usual 3
    filters : bool (optional)
        Filter data based on data2d['stat0'] = 2 and is.nan()
    save_mat : bool (optional)
        save the dataframes in .mat format.
    path2save : str (optional)
        If false, the data will be saved into the current directory
        else, should pass a path to save.
    direct : bool
        If you have all the clists in one folder.
        
    Returns
    --------------
    data2D : pandas DataFrame
        with 2D data
    data3D : pandas DataFrame
        with 3D data
    gc : nd.array
        array with good cells ID
    """
    if direct:
        paths_xy = natsorted([os.path.join(path, files) for files in os.listdir(path) if files.endswith(".mat") == True])
    else:
        paths_xy = natsorted([os.path.join(path, files+os.sep+'clist.mat') for files in os.listdir(path) if files.startswith("xy") == True])

    df_temp_2 = []
    df_temp_3 = []

    for idx,mat in enumerate(paths_xy,1):
        data_temp = scipy.io.loadmat(mat)
        data2D_temp = df_data2d(data_temp,fps)
        data3D_temp,_ = df_data3d(data_temp,fps,False)

        if direct:
            data2D_temp['fov'] = 'xy' + str(idx)
            data3D_temp['fov'] = 'xy' + str(idx)
        else:
            data2D_temp['fov'] = mat.split('/')[-2]
            data3D_temp['fov'] = mat.split('/')[-2]

        df_temp_2.append(data2D_temp)
        df_temp_3.append(data3D_temp)

    data2D = pd.concat(df_temp_2,ignore_index=True)
    data3D = pd.concat(df_temp_3,ignore_index=True)

    off2d,int2d = _off_set(data2D)
    off3d,int3d = _off_set(data3D)

    # update the intervals for each fov with the offsets
    for idx in range(len(int2d)):
        data2D.loc[int2d[idx][0]:int2d[idx][1]-1,'Cell ID'] = [i+off2d[int2d[idx][1]] for i in data2D.iloc[int2d[idx][0]:int2d[idx][1]]['Cell ID']]
        data2D.loc[int2d[idx][0]:int2d[idx][1]-1,'Daughter1 ID'] = [i+off2d[int2d[idx][1]] for i in data2D.iloc[int2d[idx][0]:int2d[idx][1]]['Daughter1 ID']]
        data2D.loc[int2d[idx][0]:int2d[idx][1]-1,'Daughter2 ID'] = [i+off2d[int2d[idx][1]] for i in data2D.iloc[int2d[idx][0]:int2d[idx][1]]['Daughter2 ID']]
        data2D.loc[int2d[idx][0]:int2d[idx][1]-1,'Mother ID'] = [i+off2d[int2d[idx][1]] if i != 0 else i for i in data2D.iloc[int2d[idx][0]:int2d[idx][1]]['Mother ID']]

    for idx in range(len(int3d)):
        data3D.loc[int3d[idx][0]:int3d[idx][1]-1,'Cell ID'] = [i+off3d[int3d[idx][1]] for i in data3D.iloc[int3d[idx][0]:int3d[idx][1]]['Cell ID']]

    gc = data2D[data2D['stat0']==2]['Cell ID'].values
    if filters:
        data2D = data2D[data2D['stat0']==2]
        data3D = data3D[data3D['Cell ID'].isin(gc)]
        data3D = data3D[data3D['Age (Frames)'].isna() == False]
        # test
        data2D['Time Division'] = data2D['Cell birth time'].values + data2D['Cell age'].values
         
        data3D = fluor_volume_ratio(data3D)
        data3D = cell_cycle(data3D)

    if save_mat:
        # data dictionary
        mat_data2D = {}
        mat_data3D = {} 

        # convert DF to dictionary
        mat_data2D['clist2D'] = data2D.to_dict('list')
        mat_data3D['clist3D'] = data3D.to_dict('list')

        if path2save is None:
            path_2d = os.getcwd() + os.sep + '2d_clist.mat'
            path_3d = os.getcwd() + os.sep + '3d_clist.mat'

            # save matlab file
            scipy.io.savemat(path_2d, mat_data2D)
            scipy.io.savemat(path_3d, mat_data3D)
        else:
            # save matlab file
            scipy.io.savemat(path2save + os.sep +'2d_clist', mat_data2D)
            scipy.io.savemat(path2save + os.sep +'3d_clist', mat_data3D)

    return data2D,data3D,gc

def combined_filters(df_2d, df_3d, std = 1, daughter = True, lentgh = True, plot=False):
    """
    Plot the mean fluo for the lineage for the entire experiment.
    Check also fluor_lineage().

    Parameters
    --------------
    df_2d : DataFrame
        DataFrame of 2D data
    df_3d : DataFrame
        DataFrame of 3D data
    std : int (optional)
        Number of standard deviations desired to use as 
        reference for the cell size filter
    daughter : bool (optional)
        If True the daughter filter will be used. If True
        alongised with 'lentgh' the two filters will be used
    lentgh : bool (optional)
        If True the lentgh filter will be used. If True
        alongised with 'daughter' the two filters will be used
    plot : bool (optional)
        Plot the histogram pre and pos filter the data
        
    Returns
    --------------
    df_2d : DataFrame
        Filtered DataFrame of 2D data
    df_3d : DataFrame
        Filtered DataFrame of 3D data
    gc_filter : nd.array
        Array with good cells ID after filter
    """
    gc = df_2d[df_2d['stat0'] == 2]['Cell ID'].values
    wm_cell = [cell for cell in gc if df_2d[df_2d['Cell ID']==cell]['Mother ID'].values in gc]
    good_daughters = []
    good_cells_d = []

    gc_mean_d = np.mean(df_2d['Area death'].values)
    gc_std_d = np.std(df_2d['Area death'].values)
    size_ref_d = gc_mean_d + gc_std_d*std

    # daughter filter
    for cell in wm_cell:
        mae = df_2d[df_2d['Cell ID']==cell]['Mother ID'].values[0]
        mae_d = df_2d[df_2d['Cell ID']==mae]['Area death'].values[0]
        filha_b = df_2d[df_2d['Cell ID']==cell]['Area birth'].values[0]
        if filha_b > .4*mae_d and filha_b < .6*mae_d:
            good_daughters.append(cell)

    # lentgh filter
    for cell in df_2d['Cell ID'].values:
        if (df_2d[df_2d['Cell ID']==cell]['Area death'].values) < size_ref_d:
            good_cells_d.append(cell)

    # combined filters
    gc_filter = natsorted(list(set(good_cells_d) & set(good_daughters)))

    if plot:
        fig,ax = plt.subplots(2,2, figsize = (12,8),sharex=True, sharey= True)
        yx = ax[0][0].hist(df_2d['Area death'], bins=20, range = [0,1000],label = 'Pre-Filter',histtype='stepfilled', facecolor='blue', edgecolor='k')
        ax[0][0].set_title('Stat0 == 2 - '+ str(len(df_2d['Cell ID'])) + " Cells", fontsize = 17)
        ax[0][0].axvline(size_ref_d, linestyle='-', color='red',label = 'Size ref')
        ax[0][0].axvline(gc_mean_d, linestyle='-', color='k',label = 'Pre-Filter Mean')
        ax[0][0].fill_between(np.arange(gc_mean_d-gc_std_d,gc_mean_d+gc_std_d,20),max(yx[0]),
                            linestyle='--', color='k',label = 'Pre-Filter Std', alpha = .1)
        ax[0][0].legend(loc='upper right')

        ax[0][1].hist(df_2d[df_2d['Cell ID'].isin(gc_filter)]['Area death'], bins=20, range = [0,1000],label = "Combined Filters",histtype='stepfilled', facecolor='yellow', edgecolor='k')
        ax[0][1].set_title("Combined Filters - "+ str(len(gc_filter)) + " Cells", fontsize = 17)
        ax[0][1].axvline(size_ref_d, linestyle='-', color='red',label = 'Size ref')
        ax[0][1].axvline(gc_mean_d, linestyle='-', color='k',label = 'Pre-Filter Mean')
        ax[0][1].fill_between(np.arange(gc_mean_d-gc_std_d,gc_mean_d+gc_std_d,20),max(yx[0]),
                            linestyle='--', color='k',label = 'Pre-Filter Std', alpha = .1)
        ax[0][1].legend(loc='upper right')

        ax[1][0].hist(df_2d[df_2d['Cell ID'].isin(good_daughters)]['Area death'], bins=20, range = [0,1000],label = "Daughters Filter",histtype='stepfilled', facecolor='green', edgecolor='k')
        ax[1][0].set_title("Daughters Filter - " + str(len(good_daughters)) + " Cells", fontsize = 17)
        ax[1][0].axvline(size_ref_d, linestyle='-', color='red',label = 'Size ref')
        ax[1][0].axvline(gc_mean_d, linestyle='-', color='k',label = 'Pre-Filter Mean')
        ax[1][0].fill_between(np.arange(gc_mean_d-gc_std_d,gc_mean_d+gc_std_d,20),max(yx[0]),
                            linestyle='--', color='k',label = 'Pre-Filter Std', alpha = .1)
        ax[1][0].legend(loc='upper right')

        ax[1][1].hist(df_2d[df_2d['Cell ID'].isin(good_cells_d)]['Area death'], bins=20, range = [0,1000],label = "Length Filter",histtype='stepfilled', facecolor='orange', edgecolor='k')
        ax[1][1].set_title('Length Filter - '+ str(len(good_cells_d)) + " Cells", fontsize = 17)
        ax[1][1].axvline(size_ref_d, linestyle='-', color='red',label = 'Size ref')
        ax[1][1].axvline(gc_mean_d, linestyle='-', color='k',label = 'Pre-Filter Mean')
        ax[1][1].fill_between(np.arange(gc_mean_d-gc_std_d,gc_mean_d+gc_std_d,20),max(yx[0]),
                            linestyle='--', color='k',label = 'Pre-Filter Std', alpha = .1)
        ax[1][1].legend(loc='upper right')

        fig.supxlabel('Division size (px)', fontsize = 15)
        fig.supylabel('Number of cells', fontsize = 15)
        plt.tight_layout()

    # Filters
    if daughter == True and lentgh == True:
        gc_filter = natsorted(list(set(good_cells_d) & set(good_daughters)))
    elif daughter == True and lentgh == False:
        gc_filter = natsorted(set(good_daughters))
    elif lentgh == False and lentgh == True:
        gc_filter = natsorted(set(good_cells_d))

    df_3d = df_3d[df_3d['Cell ID'].isin(gc_filter)]
    df_2d = df_2d[df_2d['Cell ID'].isin(gc_filter)]

    return df_2d,df_3d,np.asarray(gc_filter)

def plot_2d_data(df_2d):
    """
    Plot The Growth Rate, Cell Age and
    Long axis (L) death by 'Birth Time +
    Age'.

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    reverse_lineage : dict
        Dict with the lineage Cell ID's as key,
        output from 'lineages()'
    shift : int
        The time value of the shift
            
    Returns
    --------------
    None
    """	
    df_2d['Cell birth time + age'] = df_2d['Cell birth time'].values + df_2d['Cell age'].values
    _,ax = plt.subplots(1,3, figsize =(18,4))
    df = pd.melt(df_2d, id_vars=['Cell ID','Cell birth time + age'], value_vars=['Growth Rate']).sort_values(by=['Cell birth time + age'])
    ax[0].plot(df['Cell birth time + age'],df.value.rolling(150).mean())
    ax[0].set_ylabel('Growth Rate', fontsize = 15)
    ax[0].set_xlabel('Cell birth time + age', fontsize = 15)
    ax[0].set_title('Growth Rate by Birth Time + Age', fontsize = 17)

    df = pd.melt(df_2d, id_vars=['Cell ID','Cell birth time + age'], value_vars=['Cell age']).sort_values(by=['Cell birth time + age'])
    ax[1].plot(df['Cell birth time + age'],df.value.rolling(150).mean())
    ax[1].set_ylabel('Cell age', fontsize = 15)
    ax[1].set_xlabel('Cell birth time + age', fontsize = 15)
    ax[1].set_title('Cell Age by Birth Time + Age', fontsize = 17)
    
    df = pd.melt(df_2d, id_vars=['Cell ID','Cell birth time + age'], value_vars=['Long axis (L) death']).sort_values(by=['Cell birth time + age'])
    ax[2].plot(df['Cell birth time + age'],df.value.rolling(150).mean())
    ax[2].set_ylabel('Long axis (L) death', fontsize = 15)
    ax[2].set_xlabel('Cell birth time + age', fontsize = 15)
    ax[2].set_title('Long axis (L) death by Birth Time + Age', fontsize = 17)
    plt.show()

def roll_mean_plot2D(ax, df, y, x = 'Cell birth time', window = 50, fsize = 15):
    """
    A helper function to make a roll mean plot for 2D data

    Parameters
    ----------
    ax : Axes
       The axes to draw to
    df : pandas DataFrame
       df from 2D data
    y : str
       The y data
    x : str (optional)
       The x data, using 'Cell birth time'
    window : int (optional)
       Window value to roll mean, using 50
    fsize = int (optional)
       Plot font size

    Returns
    -------
    out : list
        list of artists added
    """
    ci = 1.96 * np.std(df[y])/np.sqrt(len(df[x]))
    out = ax.plot(df[x],df[y].rolling(window).mean()),\
          ax.fill_between(df[x], (df[y].rolling(window).mean()-ci), (df[y].rolling(window).mean()+ci), color='b', alpha=.1),\
          ax.set_title(y + ' roll mean', fontsize= fsize),\
          ax.set_xlabel(x, fontsize = fsize-1),\
          ax.set_ylabel(y + ' average', fontsize = fsize-1)
    
    return out

def roll_mean_plot3D(ax, df, y, id, x = 'Time (fps)', window = 1, fsize = 15):
    """
    A helper function to make a roll mean plot for 3D data

    Parameters
    ----------
    ax : Axes
       The axes to draw to
    df : pandas DataFrame
       df from 2D data
    y : str
       The y data
    id : int
       Cell ID for plot (check good cells using 'good_cells()')
    x : str (optional)
       The x data, using 'Cell birth time'
    window : int (optional)
       Window value to roll mean, using 1
    fsize = int (optional)
       Plot font size

    Returns
    -------
    out : list
        list of artists added
    """
    ci = 1.96 * np.std(df[df['Cell ID']==id][y])/np.sqrt(len(df[df['Cell ID']==id][x]))
    out = ax.plot(df[df['Cell ID']==id][x],df[df['Cell ID']==id][y].rolling(window).mean()),\
          ax.fill_between(df[df['Cell ID']==id][x], (df[df['Cell ID']==id][y].rolling(window).mean()-ci), 
                         (df[df['Cell ID']==id][y].rolling(window).mean()+ci), color='b', alpha=.1),\
          ax.set_title(y + ' roll mean Cell ' + str(id), fontsize= fsize),\
          ax.set_xlabel(x, fontsize = fsize-1),\
          ax.set_ylabel(y + ' average', fontsize = fsize-1)
    
    return out

def cells_length(df,limit_inf,limit_sup, ax = None):
    """
    Plot the cell's lenght over time, filtering by size. The data
    will be filted based on the size of the cells (< mean + 1 std)
    
    Parameters
    --------------
    df : DataFrame
        DataFrame of 3D data
    limit_inf : int
        Start moment in time (fps) of the stationary phase to filter
        the data (use growth_rate() to see it)
    limit_sup : int
        End moment in time (fps) of the stationary phase to filter
        the data (use growth_rate() to see it)
    ax : nd.array (optional)
        the ax position to plot
        
    Returns
    --------------
    None
    """
    cells_stats = list(set(df[(df['Time (fps)']>limit_inf) & (df['Time (fps)']<limit_sup)]['Cell ID']))
    size_ref = np.mean(df['Long axis (L)'].values) + np.std(df['Long axis (L)'].values)
    if ax is None:
        for cell in cells_stats:
            if np.mean(df[df['Cell ID']==cell]['Long axis (L)'].values) < size_ref:
                plt.plot(df[df['Cell ID']==cell]['Time (fps)'],df[df['Cell ID']==cell]['Long axis (L)'])
        plt.title('Cell’s length over time',fontsize = 17)
        plt.xlabel('Time (min)',fontsize = 15)
        plt.ylabel('Length (px)',fontsize = 15)
        plt.show()
    
    else:
        out = []
        for cell in cells_stats:
            if np.mean(df[df['Cell ID']==cell]['Long axis (L)'].values) < size_ref:
                out.append(ax.plot(df[df['Cell ID']==cell]['Time (fps)'],df[df['Cell ID']==cell]['Long axis (L)']))
        out.append(ax.set_xlabel('Time (min)',fontsize = 15))
        out.append(ax.set_ylabel('Length (px)',fontsize = 15))

        return out

def scatter_linear_regression(df,x = 'Area birth',y = 'Area death', ax = None):
    """
    Scatter plot with linear regression.
    
    Parameters
    --------------
    df : DataFrame
        DataFrame of 2D data
    x : str (optional)
        X column for the Linear Regression model, usual 'Area birth'
    y: str (optional)
        Y column for the Linear Regression model, usual 'Area death'
    ax : nd.array (optional)
        the ax position to plot

    Returns
    --------------
    None
    """
    diff = df[y]-df[x]
    model = LinearRegression()
    model.fit(df[x].values.reshape(len(df), 1), diff.values.reshape(len(df)))

    if ax is None:
        plt.scatter(df[x],diff,s=3)
        plt.plot(df[x], model.predict(df[x].values.reshape(len(df), 1)), color='k', linewidth=1, label = 'Coef: ' + str(model.coef_[0])[:5])
        plt.ylabel(r'$%s_{%s} - %s_{%s}$'%(y.split()[0],y.split()[-1],x.split()[0],x.split()[-1]), fontsize = 15)
        plt.xlabel(r'$%s_{%s}$'%(x.split()[0],x.split()[-1]), fontsize = 15)
        plt.title('Adder', fontsize = 17)
        plt.legend()
    else:
        out = ax.scatter(df[x],diff,s=3),\
              ax.plot(df[x], model.predict(df[x].values.reshape(len(df), 1)), color='k', linewidth=1, label = 'Coef: ' + str(model.coef_[0])[:5]),\
              ax.set_xlabel(r'$%s_{%s}$'%(x.split()[0],x.split()[-1]), fontsize = 15),\
              ax.set_ylabel(r'$%s_{%s} - %s_{%s}$'%(y.split()[0],y.split()[-1],x.split()[0],x.split()[-1]), fontsize = 15),\
              ax.set_title('Adder', fontsize = 17),\
              ax.legend()
    
        return out

def plot_fluor_lineage_single_cell(df_2d, df_3d, cell_id, ax = None):
    """
    Plot the cell's fluo in the cell cycle from the lineage,
    mother and 1 or 2 daughter's. If just plot the mother data
    it's because this cell has no daughter's, one can check this
    on the dataframe. 
    
    Parameters
    --------------
    df_2d : DataFrame
        DataFrame of 2D data
    df_3d : DataFrame
        DataFrame of 3D data
    cell_id: int
        Cell to plot the fluo
    ax : nd.array (optional)
        the ax position to plot
        
    Returns
    --------------
    None
    """
    from sklearn import preprocessing as pre

    t_mother = pre.MinMaxScaler((0,1)).fit_transform(df_3d[df_3d['Cell ID'] == cell_id]['Time (Frames)'].values.reshape(-1, 1))
    # Plot mother fluorescence:
    plt.plot(t_mother,df_3d[df_3d['Cell ID'] == cell_id]['Fluor1 mean'],label = 'Mother ' + str(int(cell_id)))

    # Plot daughters fluorescence:
    color_add = 0
    if ax is None:
        for cell in df_2d[df_2d['Mother ID'] == cell_id]['Cell ID']:
            color_add += 1
            color_array = np.array([.2,.2,.2]) * color_add
            # df_2d[df_2d['Mother ID'] == cell_id]
            t_daughter = pre.MinMaxScaler((1,2)).fit_transform(df_3d[df_3d['Cell ID'] == cell]['Time (Frames)'].values.reshape(-1, 1))
            plt.plot(t_daughter,df_3d[df_3d['Cell ID'] == cell]['Fluor1 mean'],color=color_array, label ='Daughter ' + str(int(cell)))

        plt.xlabel('Cell cycle',fontsize = 15)
        plt.ylabel('[GFP]',fontsize = 15)
        plt.title('Fluor mean - Cell ' + str(int(cell_id)),fontsize = 17)
        plt.legend()
        plt.show()
    
    else:
        out = []
        for cell in df_2d[df_2d['Mother ID'] == cell_id]['Cell ID']:
            color_add += 1
            color_array = np.array([.2,.2,.2]) * color_add
            # df_2d[df_2d['Mother ID'] == cell_id]
            t_daughter = pre.MinMaxScaler((1,2)).fit_transform(df_3d[df_3d['Cell ID'] == cell]['Time (Frames)'].values.reshape(-1, 1))
            out.append(ax.plot(t_daughter,df_3d[df_3d['Cell ID'] == cell]['Fluor1 mean'],color=color_array, label ='Daughter ' + str(int(cell))))

            return out
        
def _bin_vector(df_3d,t_cell,d_cell):
    """
    Binarize the data and make the mean for each bin.
    Suport function for fluor_lineage(). Bins are
    calculated based on the time data (previous reshaped
    and sclaed) using np.histogram(), default values for bins
    is 10.
    
    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    t_cell : nd.array
        time frames reshaped
    d_cell : int
        cell id
        
    Returns
    --------------
    mean_vector : list
        list with the mean value for each bin
    bins : nd.array
        array with bins used
    """
    count= 0
    intervals = []
    mean_vector = []
    hist_intervals,bins = np.histogram(t_cell)
    for idx in hist_intervals:
        intervals.append((count, idx+count))
        count = idx+count
    for i in range(len(intervals)):
        mean_vector.append(np.mean(df_3d[df_3d['Cell ID'] == d_cell]['Fluor1 mean'][intervals[i][0]:intervals[i][1]]))
    
    bins = (bins[1:]+bins[:-1])/2 # normalize the bins
    
    return mean_vector,bins

def fluor_lineage(df_2d,df_3d,limit_inf,limit_sup):
    """
    Calculate the fluor mean for every mother cell and daugther.
    The bins here were not made using pandas.

    Parameters
    --------------
    df_2d : DataFrame
        DataFrame of 2D data
    df_3d : DataFrame
        DataFrame of 3D data
    limit_inf : int
        Start moment in time (fps) of the stationary phase to filter
        the data (use growth_rate() to see it)
    limit_sup: int
        End moment in time (fps) of the stationary phase to filter
        the data (use growth_rate() to see it)
        
    Returns
    --------------
    mean_m : nd.array
        vector with the mean value for every mother cell
    bins_m : nd.array
        vector with bins used for mother cell's (cell cycle)
    mean_d : nd.array
        vector with the mean value for every daugther cell
    bins_d : nd.array
        vector with bins used for daugthers cell's (cell cycle)
    cells_removed : list
        list with cell id for each cell without daugther
    """
    cells_removed = []
    cells_stats = list(set(df_3d[(df_3d['Time (fps)']>limit_inf) & (df_3d['Time (fps)']<limit_sup)]['Cell ID']))
    mother_cells = natsorted(list(set(df_2d['Cell ID'].values) & set(cells_stats)))
    mean_d = []
    mean_m = []

    for mother in range(len(mother_cells)):
        if len(df_2d[df_2d['Mother ID'] == mother_cells[mother]]) > 0:
            # mother cells
            t_mother = pre.MinMaxScaler((-0.05,1.05)).fit_transform(df_3d[df_3d['Cell ID'] == mother_cells[mother]]['Time (Frames)'].values.reshape(-1, 1))
            temp_vector_m,bins_m = _bin_vector(df_3d,t_mother,mother_cells[mother])
            mean_m.append(temp_vector_m)

            # daughter cells
            daughter_cells = df_2d[df_2d['Mother ID'] == mother_cells[mother]]['Cell ID']
            daughter_cells = natsorted(list(set(daughter_cells) & set(cells_stats)))
            for daughter in daughter_cells:   
                t_daughter = pre.MinMaxScaler((0.95,2.05)).fit_transform(df_3d[df_3d['Cell ID'] == daughter]['Time (Frames)'].values.reshape(-1, 1))
                temp_vector_d,bins_d = _bin_vector(df_3d,t_daughter,daughter)
                mean_d.append(temp_vector_d)

        else:
            cells_removed.append(mother_cells[mother])

    mean_m = np.array([np.array(vector) for vector in mean_m]) # transform list to vector
    mean_d = np.array([np.array(vector) for vector in mean_d]) # transform list to vector

    mean_m = mean_m[~np.isnan(mean_m).any(axis=1)] # remove nan values
    mean_d = mean_d[~np.isnan(mean_d).any(axis=1)] # remove nan values

    return mean_m,bins_m,mean_d,bins_d,cells_removed
    
def plot_fluor_lineage(df_2d,df_3d,limit_inf,limit_sup,conf = 0.95, method= np.mean, ax = None):
    """
    Plot the mean fluo for the lineage for the entire experiment.
    Check also fluor_lineage().
    
    Parameters
    --------------
    df_2d : DataFrame
        DataFrame of 2D data
    df_3d : DataFrame
        DataFrame of 3D data
    limit_inf : int
        Start moment in time (fps) of the stationary phase to filter
        the data (use growth_rate() to see it)
    limit_sup: int
        End moment in time (fps) of the stationary phase to filter
        the data (use growth_rate() to see it)
    conf : float (optional)
        confidence interval desired for bootstrap, default - 0.95
    method : function (optional)
        Method to use as basis for the boostrap,
        default its np.mean. Other options are np.std
        and stats.sem.
    ax : nd.array (optional)
        the ax position to plot
        
    Returns
    --------------
    None
    """
    mean_m,bins_m,mean_d,bins_d,_ = fluor_lineage(df_2d,df_3d,limit_inf,limit_sup)

    mother_ci = np.zeros((len(bins_m),2))
    daughther_ci = np.zeros((len(bins_m),2))

    for i in range(len(bins_m)):
        mother_ci[i,:] = ci_bootstrap(mean_m[:,i],conf = conf, method = method, plot = False)
        daughther_ci[i,:] = ci_bootstrap(mean_d[:,i],conf = conf, method = method, plot = False)

    mean_m = np.mean(mean_m,0)
    mean_d = np.mean(mean_d,0)

    if ax is None:
        for i in range(len(bins_m)):
            left_m,left_d = bins_m[i] - 0.30 / 15, bins_d[i] - 0.30 / 15
            top_m,top_d = mother_ci[i,1], daughther_ci[i,1]
            right_m,right_d = bins_m[i] + 0.30 / 15, bins_d[i] + 0.30 / 15
            bottom_m,bottom_d = mother_ci[i,0], daughther_ci[i,0] 

            plt.plot([bins_m[i],bins_m[i]],[mother_ci[i,0],mother_ci[i,1]],color = 'blue', alpha = .3)
            plt.plot([bins_d[i],bins_d[i]],[daughther_ci[i,0],daughther_ci[i,1]],color = 'orange', alpha = .3)

            # mother error bar
            plt.plot([bins_m[i], bins_m[i]], [top_m, bottom_m], color='k' )
            plt.plot([left_m, right_m], [top_m, top_m], color='k')
            plt.plot([left_m, right_m], [bottom_m, bottom_m], color='k')

            # daughther error bar
            plt.plot([bins_d[i], bins_d[i]], [top_d, bottom_d], color='k' )
            plt.plot([left_d, right_d], [top_d, top_d], color='k')
            plt.plot([left_d, right_d], [bottom_d, bottom_d], color='k')

        plt.plot(bins_m,mean_m, 'o-',label = "Mother cell's",color='blue')
        plt.plot(bins_d,mean_d, 'o-',label = "Daughter cell's", color= 'orange')

        plt.title('Fluor Mean Lineage', fontsize = 17)
        plt.xlabel('Cell cycle', fontsize = 15)
        plt.ylabel('[GFP]', fontsize = 15)
        plt.legend()
        plt.show()

    else:
        out = []
        for i in range(len(bins_m)):

            left_m,left_d = bins_m[i] - 0.30 / 15, bins_d[i] - 0.30 / 15
            top_m,top_d = mother_ci[i,1], daughther_ci[i,1]
            right_m,right_d = bins_m[i] + 0.30 / 15, bins_d[i] + 0.30 / 15
            bottom_m,bottom_d = mother_ci[i,0], daughther_ci[i,0] 

            out.append(ax.plot([bins_m[i],bins_m[i]],[mother_ci[i,0],mother_ci[i,1]],color = 'blue', alpha = .3))
            out.append(ax.plot([bins_d[i],bins_d[i]],[daughther_ci[i,0],daughther_ci[i,1]],color = 'orange', alpha = .3))
    
            # mother error bar
            out.append(ax.plot([bins_m[i], bins_m[i]], [top_m, bottom_m], color='k' ))
            out.append(ax.plot([left_m, right_m], [top_m, top_m], color='k'))
            out.append(ax.plot([left_m, right_m], [bottom_m, bottom_m], color='k'))

            # daughther error bar
            out.append(ax.plot([bins_d[i], bins_d[i]], [top_d, bottom_d], color='k' ))
            out.append(ax.plot([left_d, right_d], [top_d, top_d], color='k'))
            out.append(ax.plot([left_d, right_d], [bottom_d, bottom_d], color='k'))

        out.append(ax.plot(bins_m,mean_m, 'o-',label = "Mother cell's",color='blue'))
        out.append(ax.plot(bins_d,mean_d, 'o-',label = "Daughter cell's", color= 'orange'))    
        out.append(ax.set_title('Fluor Mean Lineage', fontsize = 17))
        out.append(ax.set_xlabel('Cell cycle', fontsize = 15))
        out.append(ax.set_ylabel('[GFP]', fontsize = 15))
        out.append(ax.set_title('Fluor Mean Lineage', fontsize = 17))
        out.append(ax.legend())

        return out

def ci_bootstrap(arr, conf = 0.95, method = np.mean, plot= False):
    """
    Calculate the confidence interval (CI) based on 
    different methods.By default the bootstrap is done 
    for 10000 times. 
    
    Parameters
    --------------
    arr: nd.array
        data to calculate the CI
    conf : float (optional)
        confidence interval desired, default - 0.95
    method : function (optional)
        Method to use as basis for the boostrap,
        default its np.mean. Other options are np.std
        and stats.sem.
    plot : bool (optional)
        If true plot the distribution based of the boostrap
        and the CI.
        
    Returns
    --------------
    ci_low : float
        Value for the lower CI (default 2.5%)
    ci_high : float
        Value for the higher CI (default 97.5%)
    """
    from scipy.stats import bootstrap, norm, sem
    
    rng = np.random.default_rng()
    data = (arr,) 
    res = bootstrap(data, method, confidence_level= conf,random_state=rng,n_resamples=10000)
    
    if method is np.mean:
        sample = np.mean(data)
        ci_low  = res.confidence_interval[0] 
        ci_high = res.confidence_interval[1]
    elif method is np.std:
        sample = np.std(data)
        ci_low  = np.mean(data) - res.confidence_interval[0] 
        ci_high = np.mean(data) + res.confidence_interval[1]
    elif method is sem:
        sample = sem(arr)
        ci_low  = np.mean(data) - res.confidence_interval[0] 
        ci_high = np.mean(data) + res.confidence_interval[0] 

    x = np.linspace(res.bootstrap_distribution.min(), res.bootstrap_distribution.max())
    pdf = norm.pdf(x, loc=sample, scale=res.standard_error)

    if plot:
        fig, ax = plt.subplots(1,2, figsize = (12,5))
        ax[0].hist(res.bootstrap_distribution, bins=25)
        ax[0].set_title('Bootstrap Distribution', fontsize = 17)
        ax[0].set_xlabel('Statistic value', fontsize = 15)
        ax[0].set_ylabel('Frequency', fontsize = 15)

        ax[1].hist(res.bootstrap_distribution, bins=25, density=True)
        ax[1].axvline(res.confidence_interval[0],linestyle='--', color='red',label = 'CI Low')
        ax[1].axvline(res.confidence_interval[1],linestyle='--', color='k',label = 'CI High')
        ax[1].plot(x, pdf)
        ax[1].set_title('Normal Approximation of the Bootstrap Distribution', fontsize = 17)
        ax[1].set_xlabel('Statistic value', fontsize = 15)
        ax[1].set_ylabel('PDF', fontsize = 15)
        ax[1].legend()

        plt.tight_layout()
        plt.show()

    return ci_low,ci_high

def mean_vector_np(df, column_y ,column_x = 'Time (fps)', fps = 3, good_cells = None):
    """
    Calculate the mean vector for a column in time or volume
    using numpy logic.
    
    Parameters
    --------------
    df : DataFrame
        DataFrame of 3D data
    column_y : str
        Column to calculate the mean (y axis data)
    column_x: str (optional)
        x axis data, 'Volume' or 'Time (fps)'
    fps : int (optional)
        frames per second used in the experiment, default = 3
    good_cells : List (optional)
        List with cells to use 
        if None good_cells = df['Cell ID'].values
        
    Returns
    --------------
    time :  nd.array
        time array using column_x
    vector : nd.array
        mean of the vector from column_y
    """
    if good_cells is None: 
        good_cells = list(set(df['Cell ID'].values))    
    times = np.arange(df[column_x].min(),df[column_x].max()+fps,fps)
    vector = np.zeros((len(good_cells),len(times)))

    for idx,cell in enumerate(good_cells):
        for _,data in enumerate(zip(df[df['Cell ID']==cell][column_x],df[df['Cell ID']==cell][column_y])):
            vector[idx,np.where(times==data[0])[0][0]] = data[1]

    return times, np.mean(vector,0)

def mean_vector_df(df,id_vars = 'Time (fps)', value_vars = 'F/V', conf = 0.95, method = np.mean,fps = 3):
    """
    Calculate the mean vector for a column in time or volume
    using pd.melt.
    
    Parameters
    --------------
    df : DataFrame
        DataFrame of 3D data
    id_vars : str (optional)
        Column to do the melt, default 'Time (fps)'
    value_vars : str (optional)
        Column to calculate the mean, default 'F/V'
    conf : float (optional)
        confidence interval desired, default - 0.95
    method : function (optional)
        Method to use as basis for the boostrap,
        default its np.mean. Other options are np.std
        and stats.sem.
    fps : int (optional)
        frames per second used in the experiment, default = 3

    Returns
    --------------
    time : nd.array
        time array for the experiment (using 'Time (fps)')
    mean_nd.array
        mean values array
    ci : nd.array
        matrix with CI low and high (ci[:,0] = low, ci[:,1] = high)
    """
    
    df_ratio = pd.melt(df, id_vars=[id_vars], value_vars=[value_vars]).sort_values(by=[id_vars])
    time_frames = list(set(df_ratio['Time (fps)'].values))
    mean_list = []
    times = []
    ci = np.zeros((len(range(1,len(time_frames)-2)),2))

    for i in range(1,len(time_frames)-2):
        mean_list.append(np.mean(df_ratio[df_ratio['Time (fps)'].isin(time_frames[i-1:i+2])]['value'].values))
        times.append(time_frames[i])
        ci[i-1,:] = ci_bootstrap(df_ratio[df_ratio['Time (fps)'].isin(time_frames[i-1:i+2])]['value'].values,
                                 conf = conf, method = method, plot = False)

    mean_list = np.asarray(mean_list)
    times = np.asarray(times)

    return times,mean_list,ci

def plot_mean_vector_df(df, id_vars = 'Time (fps)', value_vars = 'F/V', conf = 0.95, method = np.mean, fps = 3, ax = None):
    """
    Plot the mean vector for the entire experiment.
    Check also mean_vector_df().
    
    Parameters
    --------------
    df : DataFrame
        DataFrame of 3D data
    id_vars : str (optional)
        Column to do the melt, default 'Time (fps)'
    value_vars: str (optional)
        Column to calculate the mean, default 'F/V'
    conf : float (optional)
        confidence interval desired, default - 0.95
    method : function (optional)
        Method to use as basis for the boostrap,
        default its np.mean. Other options are np.std
        and stats.sem.
    fps : int (optional)
        frames per second used in the experiment, default = 3
    ax : nd.array (optional)
        the ax position to plot
        
    Returns
    --------------
    None
    """

    times,mean_vector,ci = mean_vector_df(df,id_vars=id_vars,value_vars=value_vars,conf=conf,method=method,fps=fps)

    if ax is None:
        plt.plot(times,mean_vector)
        plt.fill_between(times,ci[:,0],ci[:,1], alpha=.3)
        plt.title(value_vars + ' mean', fontsize = 17)
        plt.xlabel(id_vars, fontsize = 15)
        plt.ylabel(value_vars + ' (a.u.)', fontsize = 15)
        plt.show()

    else:
        out = ax.plot(times,mean_vector),\
              ax.fill_between(times,ci[:,0],ci[:,1], alpha=.3),\
              ax.set_xlabel(id_vars, fontsize = 15),\
              ax.set_ylabel(value_vars + ' (a.u.)', fontsize = 15),\
              ax.set_title(value_vars + ' mean', fontsize = 17)
        
        return out

def instantaneous_measuraments(df,good_cells, column, fps = 3, minus = 1, plus = 2):
    """
    Calculate the instantaneous mesasurament for each cell 
    in the dataset using a window of 3 data points. 
    
    Parameters
    --------------
    df : DataFrame
        DataFrame of 3D data
    good_cells : nd.array
        array with good cells (output of df_data3d() )
    column = str
        df column that you want to calculate
    fps : int (optional)
        frames per second used in the experiment
    minus : int
        value to change the window to the convolution
        to adapt the range for the loop 
        (e.g, default = i in range(1,len(x)-2))
    plus : int
        value to change the window to the convolution
        to adapt the range for the loop 
        (e.g, default = i in range(1,len(x)-2))
        
    Returns
    --------------
    im : nd.array
        growth rate for each cell
    time : nd.array
        time array in fps
    """
    gr_inst = {}
    time_dict = {}
    for cell in good_cells:
        y = np.log(df[df['Cell ID']==cell][column].values)
        x = df[df['Cell ID']==cell]['Time (fps)'].values
        gr_temp = []
        t_temp = []
        for i in range(minus,len(x)-plus):
            model = LinearRegression()
            if plus == 0 and minus == 0:
                model.fit(x.reshape(len(x), 1), y.reshape(len(y)))
            else:
                x_temp = x[i-minus:i+plus]
                y_temp = y[i-minus:i+plus]
                model.fit(x_temp.reshape(len(x_temp), 1), y_temp.reshape(len(y_temp)))
            gr_temp.append(model.coef_[0])
            t_temp.append(x[i])
        gr_inst[cell] = np.array(gr_temp)
        time_dict[cell] = np.array(t_temp)

    times = np.arange(df['Time (fps)'].min(),df['Time (fps)'].max()+fps,fps)
    im = np.zeros((len(gr_inst),len(times)))

    for idx,cell in enumerate(good_cells):
        for _,data in enumerate(zip(time_dict[cell],gr_inst[cell])):
            im[idx,np.where(times==data[0])[0][0]] = data[1]

    return im,times

def simple_bootstrap(arr, conf = 0.95,times = 10000):
    """
    Calculate the confidence interval (CI).By default 
    the bootstrap is done for 10000 times. 
    
    Parameters
    --------------
    arr: nd.array
        data to calculate the CI
    conf : float (optional)
        confidence interval desired, default - 0.95
    times : int (optional)
        Number of times to run the bootstrap

    Returns
    --------------
    ci_low : float
        Value for the lower CI (default 2.5%)
    ci_high : float
        Value for the higher CI (default 97.5%)
    """
    values = [np.random.choice(arr,size=len(arr),replace=True).mean() for i in range(times)] 
    ci_low,ci_high = np.percentile(values,[100*(1-conf)/2,100*(1-(1-conf)/2)])
    
    return ci_low,ci_high

def progressbar(it, prefix="", size=60, out=sys.stdout):
    """
    Binarize the data and make the mean for each bin.
    Suport function for fluor_lineage(). Bins are
    calculated based on the time data (previous reshaped
    and sclaed) using np.histogram(), default values for bins
    is 10.
    
    Parameters
    --------------
    it : int
        range of iterations (for size)
    prefix : str
        Term before the computation
    size : int (optional)
        progress bar size
    out : function
        show the bar in the console
        
    Returns
    --------------
    None

    Author
    --------------
    https://stackoverflow.com/questions/3160699/python-progress-bar
    """
    import sys
    count = len(it)
    def show(j):
        x = int(size*j/count)
        print("{}[{}{}] {}/{}".format(prefix, "#"*x, "."*(size-x), j, count), 
                end='\r', file=out, flush=True)
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)

def derivative(df_3d, column = 'Fluor1 sum',minus = 1, plus = 2, bar = True):
    """
    Calculate the derivatine of the fluor sum using
    a 3 slide window. Return a dataframe with the
    derivative, derivative normalized by volume,
    logrith

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    column : str
        column to calculate the derivative
    minus : int (optional)
        value to change the window to the convolution
        to adapt the range for the loop 
        (e.g, default = i in range(1,len(x)-2))
    plus : int (optional)
        value to change the window to the convolution
        to adapt the range for the loop 
        (e.g, default = i in range(1,len(x)-2))
        
    Returns
    --------------
    df_derivative : DataFrame
        DataFrame with the deravtives data alongside with 
        the time and volume 
    dict_derivative : dict
        A dict with other dicts that contains the 
        [volume,fluor,time,F/V and all derivatives] 
        of each cell at the slide window (cell number = key).
    """
    import sys
    if 'tqdm' in sys.modules:
        from tqdm.notebook import tqdm
    else:
        bar = False
    if column == 'Volume':
        df_dev = pd.melt(df_3d, id_vars=['Cell ID','Time (fps)','F/V'], value_vars=[column]).sort_values(by=['Cell ID','Time (fps)'])
    else:
        df_dev = pd.melt(df_3d, id_vars=['Cell ID','Time (fps)','F/V','Volume'], value_vars=[column]).sort_values(by=['Cell ID','Time (fps)'])

    dict_derivative,vol_dict,time_dict,bin_dict,ratio_dict = {},{},{},{},{}
    dev_dict,dev_dict_norm,dev_dict_log,dev_dict_log_norm = {},{},{},{}

    cells = natsorted(set(df_dev['Cell ID']))
    count = 0

    if bar:
        for idx in tqdm(range(len(cells))):
            cell = cells[count]
            count += 1
            time_frames = df_dev[(df_dev['Cell ID'] == cell) & (df_dev['Time (fps)'])]['Time (fps)'].values
            vol,deriv,deriv_log,deriv_norm,deriv_log_norm,time,ratio = [],[],[],[],[],[],[]
            model_deriv = LinearRegression()
            model_deriv_log = LinearRegression()

            for i in range(minus,len(time_frames)-plus):
                x_temp = time_frames[i-minus:i+plus]
                y_temp = df_dev[(df_dev['Cell ID'] == cell) & (df_dev['Time (fps)'].isin(x_temp))].value.values
                y_temp_log = np.log(y_temp)
                model_deriv.fit(x_temp.reshape(len(x_temp), 1), y_temp.reshape(len(y_temp)))
                model_deriv_log.fit(x_temp.reshape(len(x_temp), 1), y_temp_log.reshape(len(y_temp_log)))
                if column == 'Volume':
                    deriv_norm.append(model_deriv.coef_[0]/df_dev[df_dev['Cell ID'] == cell].value.values[i])
                    deriv_log_norm.append(model_deriv_log.coef_[0]/df_dev[df_dev['Cell ID'] == cell].value.values[i])
                    vol.append(df_dev[df_dev['Cell ID'] == cell].value.values[i])
                else:
                    deriv_norm.append(model_deriv.coef_[0]/df_dev[df_dev['Cell ID'] == cell]['Volume'].values[i])
                    deriv_log_norm.append(model_deriv_log.coef_[0]/df_dev[df_dev['Cell ID'] == cell]['Volume'].values[i])
                    vol.append(df_dev[df_dev['Cell ID'] == cell]['Volume'].values[i])

                deriv.append(model_deriv.coef_[0])
                deriv_log.append(model_deriv_log.coef_[0])
                time.append(time_frames[i])
                ratio.append(df_dev[df_dev['Cell ID'] == cell]['F/V'].values[i])
            
            vol_dict[cell] = np.array(vol)
            dev_dict[cell] = np.array(deriv)
            dev_dict_norm[cell] = np.array(deriv_norm) 
            dev_dict_log[cell] = np.array(deriv_log)
            dev_dict_log_norm[cell] = np.array(deriv_log_norm)
            time_dict[cell] = np.array(time)
            ratio_dict[cell] = np.asarray(ratio)
            bin_dict[cell] = pre.MinMaxScaler((0,1)).fit_transform(np.arange(0,len(time),1).reshape(-1, 1)).flatten()
    else:
        for idx in range(len(cells)):
            cell = cells[count]
            count += 1
            time_frames = df_dev[(df_dev['Cell ID'] == cell) & (df_dev['Time (fps)'])]['Time (fps)'].values
            vol,deriv,deriv_log,deriv_norm,deriv_log_norm,time,ratio = [],[],[],[],[],[],[]
            model_deriv = LinearRegression()
            model_deriv_log = LinearRegression()

            for i in range(minus,len(time_frames)-plus):
                x_temp = time_frames[i-minus:i+plus]
                y_temp = df_dev[(df_dev['Cell ID'] == cell) & (df_dev['Time (fps)'].isin(x_temp))].value.values
                y_temp_log = np.log(y_temp)
                model_deriv.fit(x_temp.reshape(len(x_temp), 1), y_temp.reshape(len(y_temp)))
                model_deriv_log.fit(x_temp.reshape(len(x_temp), 1), y_temp_log.reshape(len(y_temp_log)))

                if column == 'Volume':
                    deriv_norm.append(model_deriv.coef_[0]/df_dev[df_dev['Cell ID'] == cell].value.values[i])
                    deriv_log_norm.append(model_deriv_log.coef_[0]/df_dev[df_dev['Cell ID'] == cell].value.values[i])
                    vol.append(df_dev[df_dev['Cell ID'] == cell].value.values[i])
                else:
                    deriv_norm.append(model_deriv.coef_[0]/df_dev[df_dev['Cell ID'] == cell]['Volume'].values[i])
                    deriv_log_norm.append(model_deriv_log.coef_[0]/df_dev[df_dev['Cell ID'] == cell]['Volume'].values[i])
                    vol.append(df_dev[df_dev['Cell ID'] == cell]['Volume'].values[i])

                deriv.append(model_deriv.coef_[0])
                deriv_log.append(model_deriv_log.coef_[0])
                time.append(time_frames[i])
                ratio.append(df_dev[df_dev['Cell ID'] == cell]['F/V'].values[i])

            vol_dict[cell] = np.array(vol)
            dev_dict[cell] = np.array(deriv)
            dev_dict_norm[cell] = np.array(deriv_norm) 
            dev_dict_log[cell] = np.array(deriv_log)
            dev_dict_log_norm[cell] = np.array(deriv_log_norm)
            time_dict[cell] = np.array(time)
            ratio_dict[cell] = np.asarray(ratio)
            bin_dict[cell] = pre.MinMaxScaler((0,1)).fit_transform(np.arange(0,len(time),1).reshape(-1, 1)).flatten()

    df_temp = []
    for cell in vol_dict.keys():
        ref_list = [cell]*len(vol_dict[cell])
        df_temp.append(pd.DataFrame(list(zip(ref_list,time_dict[cell],ratio_dict[cell],vol_dict[cell],dev_dict[cell],dev_dict_norm[cell],dev_dict_log[cell], dev_dict_log_norm[cell],bin_dict[cell])),
                                    columns=['Cell ID','Time (fps)','F/V','Volume','Derivative','Derivative/V','Derivative Log','Derivative Log/V','Cell Cycle']))
             
    df_derivative = pd.concat(df_temp, ignore_index=True)
    dict_derivative['Volume'],dict_derivative['Time (fps)'],dict_derivative['F/V'] = vol_dict,time_dict,ratio_dict
    dict_derivative['Cell Cycle'],dict_derivative['Derivative'],dict_derivative['Derivative/V'] = bin_dict,dev_dict,dev_dict_norm
    dict_derivative['Derivative Log'], dict_derivative['Derivative Log/V'] = dev_dict_log,dev_dict_log_norm

    return df_derivative,dict_derivative

def column_mean(df,column,conf = .95, method = np.mean, plot_hist = False):
    """
    Estimate the column mean for all cells using
    the 'pd.melt()' method, sorting by time.

    Parameters
    --------------
    df : DataFrame
        DataFrame with the deravtive data alongside 
        with the time and volume or the 3D df
    column : str
        column to calculate the mean   
    conf : float (optional)
        confidence interval desired for bootstrap,
        default - 0.95
    method : function (optional)
        Method to use as basis for the boostrap,
        default its np.mean. Other options are np.std
        and stats.sem.
    plot_hist: bool (optional)
        Plot a histogram with the data distribution
        
    Returns
    --------------
    time : nd.array
        1D time vector with the uniques time points (fps)
    mean : nd.array
        1D time vector with the mean for 
        each time point
    ci : nd.array
        array(2,bins) with inferior Confidence
        Interval in the first column and the 
    """
    len_data = []
    time = natsorted(set(df['Time (fps)'].values))
    mean = np.zeros((len(time),))
    ci = np.zeros((len(time),2))

    for idx,ts in enumerate(time):
        temp = df[df['Time (fps)']==ts][column].values
        mean[idx] = np.mean(temp)
        len_data.append(len(df[df['Time (fps)']==ts][column].values))
        if len(temp) > 1:
            ci[idx,:] = ci_bootstrap(temp,conf = conf, method = method, plot = False)
        else:
            ci[idx,:] = [-temp[0],+temp[0]]

    if plot_hist:
        plt.hist(len_data)
        plt.ylabel('Amount of time points')
        plt.xlabel('len(derivative)')
        plt.show()

    return time,mean,ci

def plot_column_mean(time,mean,ci,column,color = 'Orange', ax = None):
    """
    Plot the colunm mean for the entire experiment.
    Check also 'column_mean()'.

    Parameters
    --------------
    time : nd.array
            1D time vector with the uniques time points (fps)
            output of 'column_mean()'
    mean : nd.array
            1D time vector with the mean growth rate for 
            each time point output of 'column_mean()'
    ci : nd.array
            array(2,bins) with inferior Confidence
            Interval in the first column and the 
            upper one in the second output of 
            'column_mean()'
    column : str
            DataFrame column to estimate the column_mean,
            default is 'Volume'
    color: str (optional)
            Color to plot the growth rate.
    ax : nd.array (optional)
            the ax position to plot
        
    Returns
    --------------
    None
    """
    if ax is None:
        plt.plot(time,mean, color = color)
        plt.fill_between(time,ci[:,0],ci[:,1],color=color, alpha = .3)
        plt.title(column+ ' Mean', fontsize = 17)
        plt.xlabel('Time (min)', fontsize = 15)
        plt.ylabel(column, fontsize = 17)
        plt.show()
    else:
        out = []
        out.append(ax.plot(time,mean, color = color))
        out.append(ax.fill_between(time,ci[:,0],ci[:,1],color=color, alpha = .3))
        out.append(ax.set_title(column+ ' Mean', fontsize = 17))
        out.append(ax.set_xlabel('Time (min)', fontsize = 15))
        out.append(ax.set_ylabel(column, fontsize = 17))
        
        return out

def	derivative_binning(df_deriv, derivative_column = 'Derivative', bins = 10 , plot_params = None, sort_by = 'Cell Cycle', print_bins = True, conf = 0.95, method = np.mean):
    """
    Perform the binning of the derivative of one column.
    Please check 'derivative()'.
    
    Parameters
    --------------
    df_deriv : DataFrame
        DataFrame with the derivative data, output
        from 'derivative()'
    derivative_column : str (optional)
        Derivative column to bin. Columns: 'Derivative',
        'Derivative/V' is the derivative normalized by 
        volume, 'Log Derivative' logarithm of the derivative,
        'Log Derivative/V' logarithm of the derivative
        normalized by volume.
    bins : int or list (optional)
        Amount of bins to use, default = 10. This
        parameter can also receive a list/nd.array
        with the desired bins intervals. 
        e.g., [0,500,1000,1500,2000,2500,3000,
        3500,4000,4500].
    plot_params : dict
        Is the output of 'growth_rate()', this func
        will add the bin parameter to the dict, if None
        the output will a dict with just one item 'bin'.
    sort_by : str (optional)
        The parameter to bin the data (x axis), the 
        default is 'Cell Cycle', but is also possible 
        to bin by 'Volume' or 'Time (fps)' columns.
    print_bins : bool
        If true will plot a report of bins
    conf : float (optional)
        confidence interval desired for bootstrap,
		default - 0.95
    method : function (optional)
        Method to use as basis for the boostrap,
        default its np.mean. Other options are np.std
        and stats.sem.
        
    Returns
    --------------
    bins_mean : nd.array
        array(2,bins) with bin value in the first
		column and the mean value in the second
    ci : dict
        array(2,bins) with inferior Confidence
		Interval in the first column and the 
		upper one in the second
    """
    if sort_by == 'Time (fps)':
        df_bin = pd.melt(df_deriv, id_vars=['Cell ID','Time (fps)','Volume','Cell Cycle'], value_vars=[derivative_column]).sort_values(by=[sort_by])
        x = df_bin['Time (fps)'].value_counts(bins = bins, sort = False)
    elif sort_by == 'Cell Cycle':
        df_bin = pd.melt(df_deriv, id_vars=['Cell ID','Time (fps)','Volume','Cell Cycle'], value_vars=[derivative_column]).sort_values(by=[sort_by])
        x = df_bin['Cell Cycle'].value_counts(bins = 10, sort = False)
    else:
        df_bin = pd.melt(df_deriv, id_vars=['Cell ID','Time (fps)','Volume','Cell Cycle'], value_vars=[derivative_column]).sort_values(by=[sort_by])
        x = df_bin['Volume'].value_counts(bins = bins, sort = False)
        
    if isinstance(bins, int):
        bins_array = np.zeros((bins,2))
    else:
        bins_array = np.zeros((len(bins)-1,2))

    for idx,interval in enumerate(x.keys()):
        bins_array[idx][0] = interval.left
        bins_array[idx][1] = interval.right
        if print_bins:
            if sort_by == 'Cell Cycle':
                print('Bin: {:>2}|Bin Min Value: {:>6}|Bin Max Value: {:>4}|Unique Points: {:>5}'.format((idx+1),float(interval.left),float(interval.right),x.values[idx]))
            else:
                print('Bin: {:>2}|Bin Min Value: {:>4}|Bin Max Value: {:>4}|Unique Points: {:>5}'.format((idx+1),int(interval.left),int(interval.right),x.values[idx]))
    if print_bins:
        print('\n')
    bins_mean = np.zeros((len(bins_array),2))
    ci = np.zeros((len(bins_array),2))
    for i in range(len(bins_array)):
        temp = df_bin[df_bin[sort_by].between(bins_array[i][0], bins_array[i][1])].value
        bins_mean[i,0] = bins_array[i][1]
        bins_mean[i,1] = np.mean(temp)
        if len(temp) > 1:
            ci[i,:] = ci_bootstrap(temp,conf = conf, method = method, plot = False)

    return bins_mean,ci

def bin_column(df_3d,column = 'F/V',bins = 10,sort_by = 'Cell Cycle',print_bins = True,conf = .95,method = np.mean):
    """
    Perform the binning of one column.

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame with the 3D data
    column : str
        DataFrame column to bin
    bins : int or list (optional)
        Amount of bins to use, default = 10. This
        parameter can also receive a list/nd.array
        with the desired bins intervals. 
        e.g., [0,500,1000,1500,2000,2500,3000,
        3500,4000,4500].
    sort_by : str (optional)
        The parameter to bin the data (x axis), the 
        default is 'Cell Cycle', but is also possible 
        to bin by 'Volume' or 'Time (fps)' columns.
    print_bins : bool
        If true will plot a report of bins
    conf : float (optional)
        confidence interval desired for bootstrap,
        default - 0.95
    method : function (optional)
        Method to use as basis for the boostrap,
        default its np.mean. Other options are np.std
        and stats.sem.
        
    Returns
    --------------
    bins_mean : nd.array
        array(2,bins) with bin value in the first
        column and the mean value in the second
    ci : dict
        array(2,bins) with inferior Confidence
        Interval in the first column and the 
        upper one in the second
    """
    if sort_by == 'Time (fps)':
        df_bin = pd.melt(df_3d, id_vars=['Cell ID','Time (fps)','Volume','Cell Cycle'], value_vars=[column]).sort_values(by=[sort_by])
        x = df_bin['Time (fps)'].value_counts(bins = bins, sort = False)
    elif sort_by == 'Cell Cycle':
        df_bin = pd.melt(df_3d, id_vars=['Cell ID','Time (fps)','Volume','Cell Cycle'], value_vars=[column]).sort_values(by=[sort_by])
        x = df_bin['Cell Cycle'].value_counts(bins = 10, sort = False)
    else:
        df_bin = pd.melt(df_3d, id_vars=['Cell ID','Time (fps)','Volume','Cell Cycle'], value_vars=[column]).sort_values(by=[sort_by])
        x = df_bin['Volume'].value_counts(bins = bins, sort = False)

    if isinstance(bins, int):
        bins_array = np.zeros((bins,2))
    else:
        bins_array = np.zeros((len(bins)-1,2))

    for idx,interval in enumerate(x.keys()):
        bins_array[idx][0] = interval.left
        bins_array[idx][1] = interval.right
        if print_bins:
            if sort_by == 'Cell Cycle':
                print('Bin: {:>2}|Bin Min Value: {:>6}|Bin Max Value: {:>4}|Unique Points: {:>5}'.format((idx+1),float(interval.left),float(interval.right),x.values[idx]))
            else:
                print('Bin: {:>2}|Bin Min Value: {:>4}|Bin Max Value: {:>4}|Unique Points: {:>5}'.format((idx+1),int(interval.left),int(interval.right),x.values[idx]))
    if print_bins:
        print('\n')
    bins_mean = np.zeros((len(bins_array),2))
    ci = np.zeros((len(bins_array),2))
    for i in range(len(bins_array)):
        temp = df_bin[df_bin[sort_by].between(bins_array[i][0], bins_array[i][1])].value
        bins_mean[i,0] = bins_array[i][1]
        bins_mean[i,1] = np.mean(temp)
        if len(temp) > 1:
            ci[i,:] = ci_bootstrap(temp,conf = conf, method = method, plot = False)
            
    return bins_mean,ci

def plot_derivative_mean(time,mean,ci,column = 'Fluor1 sum',derivative_column = 'Derivative',color = 'Orange', ax = None):
    """
    Plot the derivative mean for the entire experiment.
    Check also derivative() and column_mean().

    Parameters
    --------------
    time : nd.array
        1D time vector with the uniques time points (fps)
        output of 'derivative()'
    mean : nd.array
        1D time vector with the mean derivative for 
        each time point output of 'derivative()'
    ci : nd.array
        array(2,bins) with inferior Confidence
        Interval in the first column and the 
        upper one in the second output of 
        'derivative()'
    column : str (optional)
        The feature used to calculate the derivative,
        yaxis label
    derivative_column : str (optional)
        Derivative column to bin. Columns: 'Derivative',
        'Derivative/V' is the derivative normalized by 
        volume, 'Log Derivative' logarithm of the derivative,
        'Log Derivative/V' logarithm of the derivative
        normalized by volume.
    color: str (optional)
        Color to plot the growth rate.
    ax : nd.array (optional)
        the ax position to plot
        
    Returns
    --------------
    None
    """
    if ax is None:
        plt.plot(time,mean, color = color)
        plt.fill_between(time,ci[:,0],ci[:,1],color=color, alpha = .3)
        if column == 'Volume':
            plt.title('{} Derivative'.format(column), fontsize = 17)
        else:
            plt.title('Growth Rate'.format(column), fontsize = 17)
        plt.xlabel('Time (min)', fontsize = 15)
        if derivative_column == 'Log Derivative/V':
            plt.ylabel(r'$\frac{\frac{\partial (log(%s))}{\partial t}}{Volume}$'%(column), fontsize = 15)
        elif derivative_column == 'Derivative/V':
            plt.ylabel(r'$\frac{\frac{\partial %s}{\partial t}}{Volume}$'%(column), fontsize = 15)
        elif derivative_column == 'Log Derivative':
            plt.ylabel(r'$\frac{\partial (log(%s))}{\partial t}$'%(column), fontsize = 15)
        elif derivative_column == 'Derivative':
            plt.ylabel(r'$\frac{\partial %s}{\partial t}$'%(column), fontsize = 15)
        plt.show()
    else:
        out = []
        out.append(ax.plot(time,mean, color = color))
        out.append(ax.fill_between(time,ci[:,0],ci[:,1],color=color, alpha = .3))
        if column == 'Volume':
            out.append(ax.set_title('Growth Rate', fontsize = 17))
        else:
            out.append(ax.set_title('{} Derivative'.format(column), fontsize = 17))
        out.append(ax.set_xlabel('Time (min)', fontsize = 15))
        if derivative_column == 'Log Derivative/V':
            out.append(ax.set_ylabel(r'$\frac{\frac{\partial (log(%s))}{\partial t}}{Volume}$'%(column), fontsize = 15))
        elif derivative_column == 'Derivative/V':
            out.append(ax.set_ylabel(r'$\frac{\frac{\partial %s}{\partial t}}{Volume}$'%(column), fontsize = 15))
        elif derivative_column == 'Log Derivative':
            out.append(ax.set_ylabel(r'$\frac{\partial (log(%s))}{\partial t}$'%(column), fontsize = 15))
        elif derivative_column == 'Derivative':
            out.append(ax.set_ylabel(r'$\frac{\partial %s}{\partial t}$'%(column), fontsize = 15))

        return out

def plot_bins(arr,ci,column = 'Fluor1 sum',derivative_column = 'Derivative',sort_by = 'Cell Cycle',cmap = 'rainbow',line = 'k',ax =None):
    """
    Plot the bins that are organized into arrays.

    Parameters
    --------------
    arr : nd.array
        2D array containing bins and the data to plot
    ci : nd.array
        2D array containing the inferior and superior
        confidence interval values
    column : str (optional)
        The feature used to calculate the derivative,
        yaxis label
    derivative_column : str (optional)
        Derivative column to bin. Columns: 'Derivative',
        'Derivative/V' is the derivative normalized by 
        volume, 'Log Derivative' logarithm of the derivative,
        'Log Derivative/V' logarithm of the derivative
        normalized by volume.
    sort_by : str (optional)
        The parameter used to bin the data. Default
        is 'Cell Cycle', xaxis legend.
    cmap : str (optional)
        The desired colormap to plot the interval
        between bins. Default is 'rainbow', to
        remove this pass 'None'. Check matplotlib
        colormaps for all options.
    line : str (optional)
        The color of the line, default is black. 
        Can be adjusted accordingly to the color
        map.
    ax : nd.array (optional)
        the ax position to plot
        
    Returns
    --------------
    None
    """
    interval = np.diff(arr[:,0])[0]
    n = len(arr)
    data_line = arr[:,1]
    max_v = 10e10 # max(ci[:,1])
    min_v = -10e10 # min(ci[:,0])
    bins_xthick = arr[:,0]
    bins_thick = arr[:,0] - interval/2
    color = iter(plt.get_cmap(cmap)(np.linspace(0, 1, n)))
    if ax is None:
        plt.plot(bins_thick,data_line,'-',color = line)
        for i in range(n):
            c = next(color)

            bins = arr[:,0][i] - interval/2
            up = ci[:,1][i]
            down = ci[:,0][i]
            left  = bins - (0.30*(interval/2)) / 2
            right = bins + (0.30*(interval/2)) / 2
            data = arr[:,1][i]
            right_fill = bins + interval/2
            left_fill = bins - interval/2

            plt.plot([bins,bins], [up,down],'-',color = line)
            plt.plot([left, right], [up, up], color= line)
            plt.plot([left, right], [down, down], color= line)
            # plt.plot([left_fill,right_fill],[data,data],'-',color = line, alpha = .3,linestyle='--');
            plt.plot(bins,data,'o',color = line)
            if cmap is not None:
                plt.fill_between([left_fill,right_fill],min_v,max_v, color=c,alpha=.3)
        
        if sort_by == 'Cell Cycle':
            plt.xlabel('Cell Cycle', fontsize = 15)
        elif sort_by == 'Time (fps)':
            plt.xlabel('Time (min)', fontsize = 15)
        elif sort_by == 'Volume':
            plt.xlabel(r'Volume [$\mu m^3$]', fontsize = 15)
        if derivative_column == 'Log Derivative/V':
            plt.ylabel(r'$\frac{\frac{\partial (log(%s))}{\partial t}}{Volume}$'%(column), fontsize = 15)
        elif derivative_column == 'Derivative/V':
            plt.ylabel(r'$\frac{\frac{\partial %s}{\partial t}}{Volume}$'%(column), fontsize = 15)
        elif derivative_column == 'Log Derivative':
            plt.ylabel(r'$\frac{\partial (log(%s))}{\partial t}$'%(column), fontsize = 15)
        elif derivative_column == 'Derivative':
            plt.ylabel(r'$\frac{\partial %s}{\partial t}$'%(column), fontsize = 15)
        plt.xticks(list(bins_xthick).insert(0, bins_xthick[0]-bins_xthick[0]))
            
    else:
        out = []
        out.append(ax.plot(bins_thick,data_line,'-',color = line))
        for i in range(n):
            c = next(color)

            bins = arr[:,0][i] - interval/2
            up = ci[:,1][i]
            down = ci[:,0][i]
            left  = bins - (0.30*(interval/2)) / 2
            right = bins + (0.30*(interval/2)) / 2
            data = arr[:,1][i]
            right_fill = bins + interval/2
            left_fill = bins - interval/2
            out.append(ax.plot([bins,bins], [up,down],'-',color = line))
            out.append(ax.plot([left, right], [up, up], color= line))
            out.append(ax.plot([left, right], [down, down], color= line))
            out.append(ax.plot(bins,data,'o',color = line))
            if cmap is not None:
                out.append(ax.fill_between([left_fill,right_fill],min_v,max_v, color=c,alpha=.3))

        if sort_by == 'Cell Cycle':
            out.append(ax.set_xlabel('Cell Cycle', fontsize = 15))
        elif sort_by == 'Time (fps)':
            out.append(ax.set_xlabel('Time (min)', fontsize = 15))
        elif sort_by == 'Volume':
            out.append(ax.set_xlabel(r'Volume [$\mu m^3$]',fontsize = 15))
        if derivative_column == 'Log Derivative/V':
            out.append(ax.set_ylabel(r'$\frac{\frac{\partial (log(%s))}{\partial t}}{Volume}$'%(column), fontsize = 15))
        elif derivative_column == 'Derivative/V':
            out.append(ax.set_ylabel(r'$\frac{\frac{\partial %s}{\partial t}}{Volume}$'%(column), fontsize = 15))
        elif derivative_column == 'Log Derivative':
            out.append(ax.set_ylabel(r'$\frac{\partial (log(%s))}{\partial t}$'%(column), fontsize = 15))
        elif derivative_column == 'Derivative':
            out.append(ax.set_ylabel(r'$\frac{\partial %s}{\partial t}$'%(column), fontsize = 15))
        out.append(ax.set_xticks(bins_xthick))
        
        return out

def lineages(df_2d, nodes_min=5):
    """
    Extract the lineage for the entire dataset using
    AnyTree lib, data is organized using the node 
    structure. The output is in dicts and anytree.node
    types. The last value of reverse lineage is the ID
    of a cell that we don't have data, cells without 
    daughters are removed by Stat0 == 2 filter.

    Parameters
    --------------
    df_2d : pd.DataFrame
        2D data
    nodes_min : int
        The minimum value of nodes in the lineage 
        
    Returns
    --------------
    reverse_lineages : dict
        Key values are the cell id of cells with more than
        X nodes (nodes_min), and the data values is a 
        list with the nodes to reach the key cell.
    large_lineages : dict
        Key values are the cell id of the root of the lineage,
        the data values are the all cells ids that are part
        of this lineage (mother and daughters).
    large_nodes : dict
        Key values are the cell id of the root of the lineage,
        and the values are the Node with all the information
        about this lineage, structure in anytree.Node dataype
    all_nodes : dict
        Key values are all the cell id's in the dataset,
        and the values are the Node with all the information
        about this lineage, structure in anytree.Node dataype 
    """
    reverse_lineages = {}
    large_lineages = {}
    large_nodes = {}
    all_nodes = {}

    for cell in df_2d['Cell ID'].values:
        if cell not in df_2d['Daughter1 ID'].values and cell not in df_2d['Daughter2 ID'].values:
            d1 = int(df_2d[df_2d['Cell ID'] == cell]['Daughter1 ID'].values[0])
            d2 = int(df_2d[df_2d['Cell ID'] == cell]['Daughter2 ID'].values[0])

            vars()["cell"+str(int(cell))] = Node(int(cell))
            all_nodes[int(cell)] = vars()["cell"+str(int(cell))]

            vars()["cell"+str(d1)] = Node(d1,parent = all_nodes[int(cell)])
            all_nodes[d1] = vars()["cell"+str(d1)]

            vars()["cell"+str(d2)] = Node(d2,parent = all_nodes[int(cell)])
            all_nodes[d2] = vars()["cell"+str(d2)]

        elif cell in df_2d['Daughter1 ID'].values:
            d1 = int(df_2d[df_2d['Cell ID'] == cell]['Daughter1 ID'].values[0])
            d2 = int(df_2d[df_2d['Cell ID'] == cell]['Daughter2 ID'].values[0])
            mom = int(df_2d[df_2d['Daughter1 ID'] == cell]['Cell ID'].values[0])

            vars()["cell"+str(int(cell))] = Node(int(cell), parent = all_nodes[mom])

            vars()["cell"+str(d1)] = Node(d1, parent = all_nodes[int(cell)])
            all_nodes[d1] = vars()["cell"+str(d1)]
            
            vars()["cell"+str(d2)] = Node(d2, parent = all_nodes[int(cell)])
            all_nodes[d2] = vars()["cell"+str(d2)]

        elif cell in df_2d['Daughter2 ID'].values:
            d1 = int(df_2d[df_2d['Cell ID'] == cell]['Daughter1 ID'].values[0])
            d2 = int(df_2d[df_2d['Cell ID'] == cell]['Daughter2 ID'].values[0])
            mom = int(df_2d[df_2d['Daughter2 ID'] == cell]['Cell ID'].values[0])

            vars()["cell"+str(int(cell))] = Node(int(cell), parent = all_nodes[mom])

            vars()["cell"+str(d1)] = Node(d1, parent = all_nodes[int(cell)])
            all_nodes[d1] = vars()["cell"+str(d1)]

            vars()["cell"+str(d2)] = Node(d2, parent = all_nodes[int(cell)])
            all_nodes[d2] = vars()["cell"+str(d2)]

    for cell in df_2d['Cell ID'].values:
        if all_nodes[int(cell)].is_root:
            if all_nodes[int(cell)].height >= nodes_min:
                large_lineages[int(cell)] = sorted(list(set([node.name for node in PostOrderIter(all_nodes[int(cell)])])))
                large_nodes[int(cell)] = all_nodes[int(cell)]

    for cell in df_2d['Cell ID'].values:
        for node in PostOrderIter(all_nodes[int(cell)]):
            temp = []
            if len(node.anchestors) > nodes_min:
                for name in node.anchestors:
                    temp.append(name.name)
                reverse_lineages[node.name] = temp
                reverse_lineages[node.name].extend([node.name])

    return reverse_lineages, large_lineages, large_nodes, all_nodes

def plot_reverse_lineages(df_2d, lineage_list, column = 'Long axis (L) death', label = 'Cell Size ',for_ref = 2,return_plot = False):
    """
    Plot the reverse lineage data with one column value.
    This plot remove the last value, as by the filter 
    we dont have the data from cells without daughters.
    This plot is made using graphviz, the variable can
    be used to save the plot in pdf.

    Parameters
    --------------
    df_2d : DataFrame
        2D data
    lineage_list : list
        list with the Cell ID's to plot
    column : str (optional)
        column value to plot alongside with the ID's
    label : str (optional)
        Label to use alongside with the column value,
        as usual the column name is to big.
    for_ref : int (optional)
        The limit of the for range (in range - for_ref)
    return_plot : bool (optional)
        If true return the graphviz variable.
        
    Returns
    --------------
    None
    """
    from IPython import display
    graph = graphviz.Digraph('unix', filename='unix.gv',
                        node_attr={'color': 'white', 'style': 'filled'})
    graph.attr(size='1920,1080',rankdir='LR')

    for i in range(len(lineage_list)-for_ref):
        a = 'ID '+str(lineage_list[i]) +'\n'+ label + str(int(df_2d[df_2d['Cell ID'] == lineage_list[i]][column].values[0])) 
        b = 'ID '+str(lineage_list[i+1]) +'\n'+ label + str(int(df_2d[df_2d['Cell ID'] == lineage_list[i+1]][column].values[0])) 
        graph.edge(a,b)

    if return_plot:
        return graph
    else:
        display.display_svg(graph)
  
def plot_large_lineages(df_2d, lineage_list, return_plot = False):
    """
    Plot the large lineage data without prune the data.
    This plot is made using graphviz, the variable can
    be used to save the plot in pdf.

    Parameters
    --------------
    df_2d : DataFrame
        2D data
    lineage_list : list
        list with the Cell ID's to plot
    return_plot : bool (optional)
        If true return the graphviz variable.
        
    Returns
    --------------
    None
    """
    from IPython import display
    graph = graphviz.Digraph('unix', filename='unix.gv',
                        node_attr={'color': 'white', 'style': 'filled'})
    graph.attr(size='1920,1080',rankdir='LR')

    for i in df_2d[df_2d['Cell ID'].isin(lineage_list)]['Cell ID'].values:
        id_d1 = int(df_2d[df_2d['Cell ID']==i]['Daughter1 ID'].values[0])
        id_d2 = int(df_2d[df_2d['Cell ID']==i]['Daughter2 ID'].values[0])
        parent = 'ID '+ str(int(i)) 
        d1 = 'ID '+str(id_d1)
        d2 = 'ID '+str(id_d2)
        graph.edge(parent, d1)
        graph.edge(parent, d2)

    if return_plot:
        return graph
    else:
        display.display_svg(graph)

def plot_cell_size_lineage(df_2d,reverse_lineage, ax = None):
    """
    Plot the cell size lineage. X-axis is the Cell ID and the
    Y-axis is the Cell Size ('Long axis (L) death') in px, 
    the size of the circle is proportionald to the cell size.

    Parameters
    --------------
    df_2d : DataFrame
        2D data
    reverse_lineage : list
        list with the Cell ID's to plot, check 'lineages()'
    ax : nd.array (optional)
        Pass the ax for subplot.
        
    Returns
    --------------
    None
    """	
    dict_values = {value:int(df_2d[df_2d['Cell ID']== value]['Long axis (L) death'].values[0]) for _,value in enumerate(reverse_lineage[:-1])}
    sizes = [i**1.99 for i in dict_values.values()]

    if ax is None:
        sc = plt.scatter(list(dict_values.keys()),list(dict_values.values()),sizes,zorder=2,c = list(dict_values.values()),cmap=plt.get_cmap('rainbow'),edgecolor="k")
        plt.plot(list(dict_values.keys()),list(dict_values.values()),color = 'k',zorder=1)
        for i, txt in dict_values.items(): plt.annotate(txt, (i, txt),horizontalalignment='center',verticalalignment='center')
        plt.xlim(min(dict_values.keys())-50,max(dict_values.keys())+30)
        plt.ylim(0,max(dict_values.values())+15)
        plt.ylabel('Cell Size (px)', fontsize = 15)
        plt.xlabel('Cell ID', fontsize = 15)
        plt.title('Cell Size Lineage {}-{}'.format(list(dict_values.keys())[0],list(dict_values.keys())[-1]), fontsize = 17)
        plt.colorbar(sc)
    else:
        out = []
        out.append(ax.scatter(list(dict_values.keys()),list(dict_values.values()),sizes,zorder=2,c = list(dict_values.values()),cmap=plt.get_cmap('rainbow'),edgecolor="k"))
        out.append(ax.plot(list(dict_values.keys()),list(dict_values.values()),color = 'k',zorder=1))
        for i, txt in dict_values.items(): out.append(ax.annotate(txt, (i, txt),horizontalalignment='center',verticalalignment='center'))
        out.append(ax.set_xlim(min(dict_values.keys())-50,max(dict_values.keys())+30))
        out.append(ax.set_ylim(0,max(dict_values.values())+15))
        out.append(ax.set_xlabel('Cell ID', fontsize = 15))
        out.append(ax.set_ylabel('Cell Size (px)', fontsize = 15))
        out.append(ax.set_title('Cell Size Lineage {}-{}'.format(list(dict_values.keys())[0],list(dict_values.keys())[-1]), fontsize = 17))
        out.append(plt.colorbar(out[0]))
        
        return out

def plot_cell_lineage(df_2d,reverse_lineage,x_axis,y_axis,z_axis,c_axis ,colormap = 'rainbow',ax = None):
    """
    Plot the cell size lineage. It's possible to pick different
    varibales for each axis. The function is already removing
    the last value of the reverse_lineage, because the last
    cell has no daughter and is removed by our filters. 

    Parameters
    --------------
    df_2d : DataFrame
        2D data
    reverse_lineage : list
        list with the Cell ID's to plot, check 'lineages()'
    x_axis : str
        Data to plot in the X axis, use the df_2d terms.
    y_axis : str
        Data to plot in the Y axis, use the df_2d terms.
    z_axis : str
        Data to use as Z axis, use the df_2d terms.
        In this case the Z is the size of each scatter point.
        It's intrensting to use cell size as parameter, at
        division or at birth.
    c_axis : str
        Data to use in the color of each saccater point and 
        in the color bar, use the df_2d terms. It's intresting
        to use fluorescence data.
    colormap : str (optional)
        Colormap to use. Check the matplotlib documentation.
        Default = rainbow.
    ax : nd.array (optional)
        Pass the ax for subplot.
        
    Returns
    --------------
    None
    """	
    x = [np.round(df_2d[df_2d['Cell ID']== cell][x_axis].values[0],2) for cell in reverse_lineage[:-1]]
    y = [np.round(df_2d[df_2d['Cell ID']== cell][y_axis].values[0],2) for cell in reverse_lineage[:-1]]
    z = [df_2d[df_2d['Cell ID']== cell][z_axis].values[0] for cell in reverse_lineage[:-1]]
    c = [df_2d[df_2d['Cell ID']== cell][c_axis].values[0] for cell in reverse_lineage[:-1]]
    annotations = {np.round(df_2d[df_2d['Cell ID']==c][x_axis].values[0],2):np.round(df_2d[df_2d['Cell ID']==c][y_axis].values[0],2) for c in reverse_lineage[:-1]}
    sizes = [i**1.99 for i in z]

    if ax is None:
        sc = plt.scatter(x,y,sizes,zorder=2,c = c,cmap=plt.get_cmap(colormap),edgecolor="k")
        plt.plot(x,y,color = 'k',zorder=1)
        for i, txt in annotations.items(): plt.annotate(txt, (i, txt),horizontalalignment='center',verticalalignment='center')
        plt.tight_layout()
        plt.xlim(min(x)-min(x)*.60,max(x)+max(x)*.1)
        plt.ylim(min(y)-min(y)*.6,max(y)+max(y)*.25)
        plt.ylabel(y_axis, fontsize = 15)
        plt.xlabel(x_axis, fontsize = 15)
        plt.title(y_axis+' Lineage {}-{}'.format(reverse_lineage[0],reverse_lineage[-1]), fontsize = 17)
        plt.colorbar(sc).set_label(label=c_axis,size=15)

    else:
        out = []
        out.append(ax.scatter(x,y,sizes,zorder=2,c = c,cmap=plt.get_cmap(colormap),edgecolor="k"))
        out.append(ax.plot(x,y,color = 'k',zorder=1))
        for i, txt in annotations.items(): out.append(ax.annotate(txt, (i, txt),horizontalalignment='center',verticalalignment='center'))
        out.append(plt.tight_layout())
        out.append(ax.set_xlim(min(x)-min(x)*.60,max(x)+max(x)*.1))
        out.append(ax.set_ylim(min(y)-min(y)*.6,max(y)+max(y)*.25))
        out.append(ax.set_xlabel(x_axis, fontsize = 15))
        out.append(ax.set_ylabel(y_axis, fontsize = 15))
        out.append(ax.set_title(y_axis+' Lineage {}-{}'.format(reverse_lineage[0],reverse_lineage[-1]), fontsize = 17))
        out.append(plt.colorbar(out[0]).set_label(label=c_axis,size=15))
        
        return out

def plot_lineages_check(df_2d,reverse_lineages,y_axis,colormap = 'rainbow'):
    """
    Plot all lineages using the length as x-axis. It's possible to chose
    different columns to inspect. The function is already removing
    the last value of the reverse_lineage, because the last
    cell has no daughter and is removed by our filters. 

    Parameters
    --------------
    df_2d : DataFrame
        2D data
    reverse_lineages : dict
        Dict with all lineages, output from 'lineages()'
    y_axis : str
        Data to plot in the Y axis, use the df_2d terms.
    colormap : str (optional)
        Colormap to use. Check the matplotlib documentation.
        Default = rainbow.
        
    Returns
    --------------
    None
    """	
    import matplotlib
    parameters = np.linspace(min(reverse_lineages.keys()),max(reverse_lineages.keys()),len(reverse_lineages))
    colormap = plt.get_cmap(colormap)
    norm = matplotlib.colors.Normalize(vmin=min(reverse_lineages.keys()),vmax=max(reverse_lineages.keys()))
    s_m = matplotlib.cm.ScalarMappable(cmap=colormap, norm=norm)
    s_m.set_array([])

    plt.figure(figsize = (14,9))
    for idx,lineage in enumerate(reverse_lineages):
        y = [np.round(df_2d[df_2d['Cell ID']== cell][y_axis].values[0],2) for cell in reverse_lineages[lineage][:-1]]
        x = np.arange(1,len(y)+1)
        plt.scatter(x,y,zorder=2,edgecolor="k", color =  s_m.to_rgba(parameters[idx]))
        plt.plot(x,y,color = s_m.to_rgba(parameters[idx]),zorder=1,label = lineage) 

    plt.ylabel(y_axis, fontsize=15)
    plt.xlabel('Lineage Length', fontsize=15)
    plt.title('Lineage Check', fontsize=17)
    plt.legend(bbox_to_anchor=(0, -.15),loc="upper left", borderaxespad=0.,ncol=11)
    plt.colorbar(s_m)
    plt.tight_layout()
    plt.show()

def plot_total_fluo_volume(df_3d, reverse_lineage, ax = None):
    """
    Plot the total fluorescence by the volume. 

    Parameters
    --------------
    df_3d : DataFrame
        3D data
    reverse_lineage : list
        List with the lineage Cell ID's, output 
        from 'lineages()'
    ax : nd.array
        axis to plot
            
    Returns
    --------------
    None
    """	
    df_oder = pd.melt(df_3d, id_vars=['Cell ID','Time (fps)','Volume'], value_vars=['Fluor1 sum']).sort_values(by=['Cell ID','Volume'])
    if ax is None:
        for i in reverse_lineage[:-1]:
            plt.plot(df_oder[df_oder['Cell ID']==i]['Volume'].values,df_oder[df_oder['Cell ID']==i].value.values,'.',label='Cell {}'.format(i))
        plt.title('Total Fluorescene Lineage {}-{}'.format(reverse_lineage[:-1][0],reverse_lineage[:-1][-1]),fontsize = 17)
        plt.xlabel(r'Volume [$\mu m^3$]',fontsize = 15)
        plt.ylabel('Fluor1 sum',fontsize = 15)
        plt.legend()
        plt.show()
    else:
        out = []
        for i in reverse_lineage[:-1]:
            out.append(ax.plot(df_oder[df_oder['Cell ID']==i]['Volume'].values,df_oder[df_oder['Cell ID']==i].value.values,'.',label='Cell {}'.format(i)))
        out.append(ax.set_title('Total Fluorescene Lineage {}-{}'.format(reverse_lineage[:-1][0],reverse_lineage[:-1][-1]),fontsize = 17))
        out.append(ax.set_xlabel(r'Volume [$\mu m^3$]',fontsize = 15))
        out.append(ax.set_ylabel('Fluor1 sum',fontsize = 15))
        out.append(ax.legend())

        return out

def _adjust_order(order,data,perc):
    """
    Function to adapt the order of the smooth
    filter for small data. 
    
    Parameters
    --------------
    order: int
        Initial order desired
    data : nd.array
        Data that will be filter
    perc : float
        Percentage of the length of the
        data to use as window.
        
    Returns
    --------------
    None
    """
    if order < int(len(data)*perc):
        return order
    else:
        order = order-1
        return _adjust_order(order,data,perc)

def lineage_derivative(df_3d, reverse_lineage, column = 'Fluor1 sum',derivative_column = 'Derivative',std = 2,window = None, order = 5, perc = .7, mode = 'nearest'):
    """
    Plot a single lineages derivative. It's possible to chose
    different columns to inspect. The function is already removing
    the last value of the reverse_lineage, because the last
    cell has no daughter and is removed by our filters. In this code
    the volume on the X axis are filtred using a standard window 
    (lenght*perc) and summed with the last value of the previous
    cell and with the mean diff inside the cell. To create this 
    volume array the data is filtered using the parameters 'perc' and
    'order', in case of the order is biggter than the window the order
    will be reduced automatically and saved in 'plot_params' (outout of
    the function). This manipulation isn't used to normalize the 
    derivative data. Also, the derivative data is filtered to remove
    the outliers based on the 'std' parameter.

    Parameters
    --------------
    df_3d : DataFrame
        3D data
    reverse_lineage : list
        List with the lineage Cell ID's, output 
        from 'lineages()'
    column : str (optional)
        The df column to compute the derivative, 
        default is 'Fluor1 sum'.
    derivative_column : str (optional)
        Derivative column to bin. Columns: 'Derivative',
        'Derivative/V' is the derivative normalized by 
        volume, 'Log Derivative' logarithm of the derivative,
        'Log Derivative/V' logarithm of the derivative
        normalized by volume. Default is 'Derivative'.
    std : int (optional)
        Number of standard deviations to filter the data 
    window : int (optional)
        The window value to smooth the graph using 
        savgol_filter from scipy. Default 19
        based on previous tests.
    order : int (optional)
        The order of the polynomial used to fit the 
        savgol_filter from scipy. Default 5 based
        on previous tests. If the order is too large
        for one cell the order will be automatically 
        decreased using '_adjust_order()'.
    perc : float (optional)
        Percentage of the length of each cell to estimate
        the window to smooth the data for the volume adjust.
        Default is .7, should be between .1 to 1.
    mode : str (optional)
        This determines the type of extension to use 
        for the padded signal to which the filter is applied
        check 'scipy.signal.savgol_filter()'
            
    Returns
    --------------
    df_lineage : DataFrame
        Data frame with the columns 'Cell ID', 'Time (fps)', 
        'Volume', 'Derivative', 'Smoothed Volume',
        'Smoothed Derivative'. The columns 'Cell ID', 
        'Time (fps)', 'Volume' are related and should be 
        used to access the data relative for each cell. 
        The columns 'Smoothed Volume' and 'Smoothed Volume'
        are the data filtered and concatenated of the 
        entire lineage (these columns are not related with
        the 'Cell ID). The column 'Time (fps)' is used to 
        plot the time for each cell and for entire lineage.
    peaks : dict
        Dictionary data with the indexes of local minima and
        local maxima related to 'Smoothed Volume','Smoothed 
        Derivative' and 'Time (fps)'.
    order_dict : dict
        Dictionary with the order used in each cell
    """	
    if len(df_3d[df_3d['Cell ID']==reverse_lineage[-1]]) == 0:
        keys = natsorted(reverse_lineage[:-1])
    else:
        keys = natsorted(reverse_lineage)

    df_deriv,_ = derivative(df_3d[df_3d['Cell ID'].isin(keys)],column=column,bar = False)

    x_vol_sep,x_vol_order,cell_cycles = {},{},{}

    count = 0
    order_dict = {}

    for idx,cell in enumerate(keys):
        if idx == 0:
            temp = df_deriv[df_deriv['Cell ID']==cell]['Volume'].values
            e_order = _adjust_order(order=order,data=temp,perc=perc)
            temp_vol = savgol_filter(temp, int(len(temp)*perc), e_order)
            x_vol_sep[cell] = temp_vol + count
            cell_cycles[cell] = [x_vol_sep[cell][0],x_vol_sep[cell][-1]]
            x_vol_order[cell] = e_order
        else:
            temp1 = df_deriv[df_deriv['Cell ID']==cell]['Volume'].values
            temp2 = df_deriv[df_deriv['Cell ID']==keys[idx-1]]['Volume'].values

            e_order = _adjust_order(order=order,data=temp1,perc=perc)
            temp_vol = savgol_filter(temp1, int(len(temp1)*perc), e_order)
            x_vol_order[cell] = e_order
            e_order = _adjust_order(order=order,data=temp2,perc=perc)
            temp_fut = savgol_filter(temp2, int(len(temp2)*perc), e_order)

            count = temp_fut[-1] - temp_vol[0] + count
            diff = np.mean(np.diff(temp_vol)) + count
            x_vol_sep[cell] = temp_vol+diff
            cell_cycles[cell] = [x_vol_sep[cell][0],x_vol_sep[cell][-1]]
            
    order_dict['volume order'] = x_vol_order
    vol_dict,dev_dict,time_dict,ratio_dict = {},{},{},{}
    vol_concat,dev_concat,time_concat,ratio_concat,df_temp = [],[],[],[],[]

    for cell in keys:
        temp_vol = df_deriv[df_deriv['Cell ID']==cell][derivative_column].values
        low = np.mean(temp_vol) - (np.std(temp_vol) * std)
        up = np.mean(temp_vol) + (np.std(temp_vol) * std)
        dev_filt = temp_vol[(temp_vol > low) & (temp_vol < up)]
        vol_idx = np.where(np.in1d(temp_vol, temp_vol[(temp_vol > low) & (temp_vol < up)]))[0]

        vol_dict[cell] = x_vol_sep[cell][vol_idx]
        time_dict[cell] = df_deriv[df_deriv['Cell ID']==cell]['Time (fps)'].values[vol_idx]
        ratio_dict[cell] = df_deriv[df_deriv['Cell ID']==cell]['F/V'].values[vol_idx]
        dev_dict[cell] = dev_filt

        time_concat.extend(df_deriv[df_deriv['Cell ID']==cell]['Time (fps)'].values[vol_idx])
        vol_concat.extend(x_vol_sep[cell][vol_idx])
        dev_concat.extend(dev_filt)
        ratio_concat.extend(df_deriv[df_deriv['Cell ID']==cell]['F/V'].values[vol_idx])

        ref_list = [cell]*len(vol_dict[cell])
        df_temp.append(pd.DataFrame(list(zip(ref_list,time_dict[cell],ratio_dict[cell],vol_dict[cell],dev_dict[cell])),
                                    columns=['Cell ID','Time (fps)','F/V','Volume','Derivative']))

    if window is None:
        window = int(len(vol_concat)*.2)
            
    smooth_vol, smooth_dev = savgol_filter((vol_concat,dev_concat), window , order ,mode=mode)
        
    df_lineage = pd.concat(df_temp, ignore_index=True)
    df_lineage['Smoothed Volume'] = smooth_vol
    df_lineage['Smoothed Derivative'] = smooth_dev

    local_min = argrelextrema(smooth_dev, np.less) 
    local_max = argrelextrema(smooth_dev, np.greater)
    local_min_dict,local_max_dict = {},{}

    for _, key in enumerate(cell_cycles):
        local_min_dict[key] = []
        local_max_dict[key] = []
        for values in local_min[0]:
            if cell_cycles[key][0] < smooth_vol[values] < cell_cycles[key][1]:
                local_min_dict[key] += [values]
        for values in local_max[0]:
            if cell_cycles[key][0] < smooth_vol[values] < cell_cycles[key][1]:
                local_max_dict[key] += [values]

    local_min_filt = {key:np.where(smooth_dev == min(smooth_dev[local_min_dict[key]]))[0][0] if len(smooth_dev[local_min_dict[key]])>0 else 0 for key in local_min_dict}
    local_max_filt = {key:np.where(smooth_dev == max(smooth_dev[local_max_dict[key]]))[0][0] if len(smooth_dev[local_max_dict[key]])>0 else 0 for key in local_max_dict}
    peaks = {'minima':local_min_filt,'maxima':local_max_filt}

    return df_lineage,peaks,order_dict

def lineages_shift(df_3d,reverse_lineages, shift):
    """
    Filter the lineages based on time, used to extract
    lineages that pass through a time point (up shift
    or down shift). Save the Cell ID of every lineage
    that the first cell was born before the shift and
    the last was born after the shift.

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    reverse_lineage : dict
        Dict with the lineage Cell ID's as key,
        output from 'lineages()'
    shift : int
        The time value of the shift
            
    Returns
    --------------
    lineage_filter : list
        List with the Cell ID used as a key in
        'reverse_lineages' dict.
    """	
    lineage_filter = []
    for cell in list(reverse_lineages.keys()):
        first = reverse_lineages[cell][0]
        last = reverse_lineages[cell][-2]
        if df_3d[df_3d['Cell ID']==first]['Time (fps)'].values[0] < shift and df_3d[df_3d['Cell ID']==last]['Time (fps)'].values[-1] > shift:
            lineage_filter.append(cell)

    return natsorted(lineage_filter)

def lineages_pre_shift(df_3d,reverse_lineages, shift):
    """
    Filter the lineages based on time, used to extract
    lineages pre shift. Save the Cell ID of every 
    lineage that the first cell was born before the shift.

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    reverse_lineage : dict
        Dict with the lineage Cell ID's as key,
        output from 'lineages()'
    shift : int
        The time value of the shift
            
    Returns
    --------------
    lineage_filter : list
        List with the Cell ID used as a key in
        'reverse_lineages' dict.
    """	
    lineage_filter = []
    for cell in list(reverse_lineages.keys()):
        last = reverse_lineages[cell][-2]
        if df_3d[df_3d['Cell ID']==last]['Time (fps)'].values[-1] < shift:
            lineage_filter.append(reverse_lineages[cell])

    return natsorted(lineage_filter)

def lineages_pos_shift(df_3d,reverse_lineages, shift):
    """
    Filter the lineages based on time, used to extract
    lineages pre shift. Save the Cell ID of every lineage
    that the first cell was born after the shift.

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    reverse_lineage : dict
        Dict with the lineage Cell ID's as key,
        output from 'lineages()'
    shift : int
        The time value of the shift
            
    Returns
    --------------
    lineage_filter : list
        List with the Cell ID used as a key in
        'reverse_lineages' dict.
    """	
    lineage_filter = []
    for cell in list(reverse_lineages.keys()):
        first = reverse_lineages[cell][0]
        if df_3d[df_3d['Cell ID']==first]['Time (fps)'].values[0] > shift:
            lineage_filter.append(reverse_lineages[cell])

    return natsorted(lineage_filter)

def longer_lineages(df_3d, reverse_lineages, column = 'Fluor1 mean', method = np.var):
    """
    Filter the lineages based on length and return 
    the variance (default or other method - np.mean, 
    np.std, stats.sem and others) for one feature (as 
    default 'Fluor1 mean', but can be used any column 
    of df_3d). Return a dictionary with the root cell 
    as a key, each root cell have dictionaries with 
    the longest lineages and use the last cell ID as key. 
    The values in the dict is the method applied into 
    the column choosed. This function is used as support 
    for other functions.

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    reverse_lineage : dict
        Dict with the lineage Cell ID's as key,
        output from 'lineages()'
    column : str (optional)
        The column to estimate the
    method : function (optional)
        Method to use to extract feature of
        each cell, default its np.mean. Other 
        options are np.std and stats.sem.
            
    Returns
    --------------
    root_dict : dict
        Dict with the mother Cell ID used as a key and
        inside each mother cell another dict with the
        longest daughters, the key of this dict is the
        last cell of the lineage.
    lineage_array : nd.array
        Array with all lineages in arrays. Used in plot
        functions.
    len_dict : dict
        Dict with the mother Cell ID used as a key and
        inside the lenght of the longest lineage.
    """	
    lineage_array = np.array([np.array(lineage) for lineage in list(reverse_lineages.values())],dtype=object)
    fist_cells = natsorted(set([cell[0] for cell in lineage_array]))

    len_dict = {}
    for root_cell in fist_cells:
        temp_len = 0
        for cell in lineage_array:
            if cell[0] == root_cell:
                if len(cell) > temp_len:
                    temp_len = len(cell)
        len_dict[root_cell] = temp_len

    root_dict = {}
    for root_cell in fist_cells:
        count = 0
        temp_dict = {}
        for lineage in lineage_array: # other length?
            if lineage[0] == root_cell and len_dict[root_cell] == len(lineage):
                temp_list = []
                for cell in lineage:
                    temp_list.append(method(df_3d[df_3d['Cell ID']==cell][column].values))    
                temp_dict[cell] = temp_list
                count += 1
        root_dict[root_cell] = temp_dict

    return root_dict, lineage_array, len_dict

def lineage_corr_filter(df_3d,reverse_lineages,shift = None,threshold = .96,column = 'Fluor1 mean', method = np.var):
    """
    Filter the lineages based on the correlation
    between the longest lineage of a mother cell.
    Return the lineages that match the parameters.

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    reverse_lineage : dict
        Dict with the lineage Cell ID's as key,
        output from 'lineages()'
    shift : int (optional)
        Value to use to extract the lineages that
        pass through this time point, usually used
        to get the shift moment. If shift is passed
        automatically the output will be filtered 
        for these cells. Default is None.
    threshold : float (optional)
        Value to be used as threshold to determine
        if the lineage must be maintained. If the 
        correlation between the daughter lineages 
        is less than the threshold value, the 
        lineage is maintained. Default '.96'.
    column : str (optional)
        The column to extract the feature of
        the dataset. Default 'Fluor1 mean'.
    method : function (optional)
        Method to use to extract feature of
        each cell, default its np.mean. Other 
        options are np.std and stats.sem.	
            
    Returns
    --------------
    filtered_lineage : list
        List with the Cell ID for each lineage that
        match the paramters (threshold and shift).
        Longest lineages for each mother cell and
        different from each other in case of more 
        than one.
    """
    root_dict,_,_ = longer_lineages(df_3d = df_3d,
                                    reverse_lineages = reverse_lineages,
                                    column = column, 
                                    method = method)
    if shift is not None:
        lineage_filter = lineages_shift(df_3d = df_3d,
                                        reverse_lineages = reverse_lineages,
                                        shift = shift)

    lineages_corr = []
    for mother_cell in root_dict.keys():
        keys_lineage = list(root_dict[mother_cell].keys())
        corr_arr = np.array([np.array(root_dict[mother_cell][i][:-1]) for i in root_dict[mother_cell]])
        corr = np.corrcoef(corr_arr)
        if len(corr[::2,::2]) == 1:
            lineages_corr.append(keys_lineage[0])
        elif len(corr[::2,::2][0]) > 1:
            for i in range(1,len(corr[::2,::2][0])):
                if corr[::2,::2][0][i] < threshold:
                    lineages_corr.append(keys_lineage[i])

    filtered_lineage = []
    for cell in reverse_lineages:
        if shift is not None:
            if cell in natsorted(list(set(lineages_corr) & set(lineage_filter))):
                filtered_lineage.append(reverse_lineages[cell])
        else:
            if cell in natsorted(set(lineages_corr)):
                filtered_lineage.append(reverse_lineages[cell])

    return filtered_lineage

def plot_lineage_derivative(df_lineage,peaks,column = 'Fluor1 sum',derivative_column = 'Derivative',x_axis = 'Smoothed Volume',ax = None):
    """
    Plot a single lineages derivative with the
    local minima and local maxmia. The input parameters
    are the output of 'lineage_derivative()'. This 
    function can receive both volume ('Smoothed Volume')
    and time ('Time (fps)') as the X axis.

    Parameters
    --------------
    df_lineage : DataFrame
        DataFrame output of 'lineage_derivative()' with all
        data necesseary to plot the lineages.
    peaks : dict
        Dictionary data with the indexes of local minima and
        local maxima related to 'Smoothed Volume','Smoothed 
        Derivative' and 'Time (fps)'. Output of 
        'lineage_derivative()'
    x_axis : str (optional)
        Column to use as X axis in the plot. Should be 'Smoothed 
        Volume' (default) or 'Time (fps)'.
    ax : nd.array (optional)
        Pass the ax for subplot.
            
    Returns
    --------------
    None
    """	
    cells = natsorted(list(set(df_lineage['Cell ID'].values)))
    for_label = 'Cell Division'
    minima = np.asarray(list(peaks['minima'].values()))
    maxima = np.asarray(list(peaks['maxima'].values()))
    if ax is None:
        plt.figure(figsize=(17,6))
        for cell in cells:
            plt.plot(df_lineage[df_lineage['Cell ID']==cell][x_axis],df_lineage[df_lineage['Cell ID']==cell]['Derivative'],'.',label = 'Cell {}'.format(cell))
        plt.plot(df_lineage[x_axis].values,df_lineage['Smoothed Derivative'].values,'k',label = 'Smoothed Derivative')
        plt.scatter(df_lineage[x_axis].values[minima],df_lineage['Smoothed Derivative'].values[minima],
                    75,marker='^',edgecolor='black',facecolor='white',zorder = 3, label = 'Local Minima')
        plt.scatter(df_lineage[x_axis].values[maxima],df_lineage['Smoothed Derivative'].values[maxima],
                    75,marker='v',edgecolor='black',facecolor='white',zorder = 3, label = 'Local Maxima')
        for cell in cells:
            plt.axvline(df_lineage[df_lineage['Cell ID']==cell][x_axis].values[-1],linestyle='--',color="darkgray", zorder=1, label=for_label)
            for_label = "_nolegend_"
        plt.title('Lineage {}-{}'.format(cells[0],cells[-1]), fontsize = 17)
        plt.axhline(np.mean(df_lineage['Smoothed Derivative'].values),linestyle=':',color="gray",label='Mean')
        if x_axis == 'Time (fps)':
            plt.xlabel('Time (min)', fontsize = 15)
        elif x_axis == 'Smoothed Volume' or x_axis == 'Volume':
            plt.xlabel(r'Volume [$\mu m^3$]', fontsize = 15)
        else:
            plt.xlabel(x_axis, fontsize = 15)
        if derivative_column == 'Log Derivative/V':
            plt.ylabel(r'$\frac{\frac{\partial (log(%s))}{\partial t}}{Volume}$'%(column), fontsize = 15)
        elif derivative_column == 'Derivative/V':
            plt.ylabel(r'$\frac{\frac{\partial %s}{\partial t}}{Volume}$'%(column), fontsize = 15)
        elif derivative_column == 'Log Derivative':
            plt.ylabel(r'$\frac{\partial (log(%s))}{\partial t}$'%(column), fontsize = 15)
        elif derivative_column == 'Derivative':
            plt.ylabel(r'$\frac{\partial %s}{\partial t}$'%(column), fontsize = 15)
        plt.legend(bbox_to_anchor=(1.05, 1),loc="upper left", borderaxespad=0.,ncol=1)
        plt.show()
    else:
        out = []
        for cell in cells:
            out.append(ax.plot(df_lineage[df_lineage['Cell ID']==cell][x_axis],df_lineage[df_lineage['Cell ID']==cell]['Derivative'],'.',label = 'Cell {}'.format(cell)))
        out.append(ax.plot(df_lineage[x_axis].values,df_lineage['Smoothed Derivative'].values,'k',label = 'Smoothed Derivative'))
        out.append(ax.scatter(df_lineage[x_axis].values[minima],df_lineage['Smoothed Derivative'].values[minima],
                              75,marker='^',edgecolor='black',facecolor='white',zorder = 3, label = 'Local Minima'))
        out.append(ax.scatter(df_lineage[x_axis].values[maxima],df_lineage['Smoothed Derivative'].values[maxima],
                              75,marker='v',edgecolor='black',facecolor='white',zorder = 3, label = 'Local Maxima'))
        for cell in cells:
            out.append(ax.axvline(df_lineage[df_lineage['Cell ID']==cell][x_axis].values[-1],linestyle='--',color="darkgray", zorder=1, label=for_label))
            for_label = "_nolegend_"
        out.append(ax.set_title('Lineage {}-{}'.format(cells[0],cells[-1]), fontsize = 17))
        out.append(ax.axhline(np.mean(df_lineage['Smoothed Derivative'].values),linestyle=':',color="gray",label='Mean'))
        if x_axis == 'Time (fps)':
            out.append(ax.set_xlabel('Time (min)', fontsize = 15))
        elif x_axis == 'Smoothed Volume':
            out.append(ax.set_xlabel(r'Volume [$\mu m^3$]', fontsize = 15))
        else:
            out.append(ax.set_xlabel(x_axis, fontsize = 15))
        if derivative_column == 'Log Derivative/V':
            out.append(ax.set_ylabel(r'$\frac{\frac{\partial (log(%s))}{\partial t}}{Volume}$'%(column), fontsize = 15))
        elif derivative_column == 'Derivative/V':
            out.append(ax.set_ylabel(r'$\frac{\frac{\partial %s}{\partial t}}{Volume}$'%(column), fontsize = 15))
        elif derivative_column == 'Log Derivative':
            out.append(ax.set_ylabel(r'$\frac{\partial (log(%s))}{\partial t}$'%(column), fontsize = 15))
        elif derivative_column == 'Derivative':
            out.append(ax.set_ylabel(r'$\frac{\partial %s}{\partial t}$'%(column), fontsize = 15))
        out.append(ax.legend(bbox_to_anchor=(1.05, 1),loc="upper left", borderaxespad=0.,ncol=1))

        return out

def plot_lineage_ratio(df_lineage,x_axis = 'Smoothed Volume',ax = None):
    """
    Plot a single lineages F/V ratio. The input parameters
    are the output of 'lineage_derivative()'. This 
    function can receive both volume ('Smoothed Volume')
    and time ('Time (fps)') as the X axis.

    Parameters
    --------------
    df_lineage : DataFrame
        DataFrame output of 'lineage_derivative()' with all
        data necesseary to plot the lineages.
    x_axis : str (optional)
        Column to use as X axis in the plot. Should be 'Smoothed 
        Volume' (default) or 'Time (fps)'.
    ax : nd.array (optional)
        Pass the ax for subplot.
            
    Returns
    --------------
    None
    """	
    cells = natsorted(list(set(df_lineage['Cell ID'].values)))
    for_label = 'Cell Division'
    temp_x,temp_y = [],[]
    if ax is None:
        plt.figure(figsize=(17,6))
        for cell in cells:
            plt.plot(df_lineage[df_lineage['Cell ID']==cell][x_axis],df_lineage[df_lineage['Cell ID']==cell]['F/V'],'.',label = 'Cell {}'.format(cell))
            temp_x.extend(df_lineage[df_lineage['Cell ID']==cell][x_axis].values)
            temp_y.extend(df_lineage[df_lineage['Cell ID']==cell]['F/V'].values)
        x,y = savgol_filter((temp_x,temp_y), 15 , 5)
        plt.plot(x,y,color='k',label = 'Smoothed Ratio')
        for cell in cells:
            plt.axvline(df_lineage[df_lineage['Cell ID']==cell][x_axis].values[-1],linestyle='--',color="darkgray", zorder=1, label=for_label)
            for_label = "_nolegend_"
        plt.title('Lineage {}-{}'.format(cells[0],cells[-1]), fontsize = 17)
        if x_axis == 'Time (fps)':
            plt.xlabel('Time (min)', fontsize = 15)
        elif x_axis == 'Smoothed Volume':
            plt.xlabel(r'Volume [$\mu m^3$]', fontsize = 15)
        else:
            plt.xlabel(x_axis, fontsize = 15)
        plt.ylabel('F/V', fontsize = 17)
        plt.legend(bbox_to_anchor=(1.05, 1),loc="upper left", borderaxespad=0.,ncol=1)
        plt.show()
    else:
        out = []
        for cell in cells:
            out.append(ax.plot(df_lineage[df_lineage['Cell ID']==cell][x_axis],df_lineage[df_lineage['Cell ID']==cell]['F/V'],'.',label = 'Cell {}'.format(cell)))
            temp_x.extend(df_lineage[df_lineage['Cell ID']==cell][x_axis].values)
            temp_y.extend(df_lineage[df_lineage['Cell ID']==cell]['F/V'].values)
        x,y = savgol_filter((temp_x,temp_y), 15 , 5)
        out.append(ax.plot(x,y,color='k',label = 'Smoothed Ratio'))
        for cell in cells:
            out.append(ax.axvline(df_lineage[df_lineage['Cell ID']==cell][x_axis].values[-1],linestyle='--',color="darkgray", zorder=1, label=for_label))
            for_label = "_nolegend_"
        out.append(ax.set_title('Lineage {}-{}'.format(cells[0],cells[-1]), fontsize = 17))
        if x_axis == 'Time (fps)':
            out.append(ax.set_xlabel('Time (min)', fontsize = 15))
        elif x_axis == 'Smoothed Volume':
            out.append(ax.set_xlabel(r'Volume [$\mu m^3$]', fontsize = 15))
        else:
            out.append(ax.set_xlabel(x_axis, fontsize = 15))
        out.append(ax.set_ylabel('F/V', fontsize = 17))        
        out.append(ax.legend(bbox_to_anchor=(1.05, 1),loc="upper left", borderaxespad=0.,ncol=1))

        return out

def plot_corr_lineage(df_3d, reverse_lineages, mother_cell, column = 'Fluor1 mean', method = np.var, derivative_plot = False):
    """
    Function to produce three plots. The first
    one is the correlation between the lineages
    from the same mother cell (and same size),
    the second one present the variance of the
    column for each cell inside the lineages
    easy to see if the lineages are differents,
    and the third graph present the dV/dt, dF/dt,
    Volume and Fluor1 mean of the lineage, and
    if the 'derivative_plot = True', will plot the
    lineage oscilations using 'plot_lineage_derivative()'

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    reverse_lineage : dict
        Dict with the lineage Cell ID's as key,
        output from 'lineages()'
    mother_cell : int 
        Mother cell (first cell of the lineage),
        to plot the lineage
    column : str (optional)
        The column to extract the feature of
        the dataset. Default 'Fluor1 mean'.
    method : function (optional)
        Method to use to extract feature of
        each cell, default its np.mean. Other 
        options are np.std and stats.sem.
    derivative_plot : bool (optional)
        If true beside plot the dV/dt, dF/dt,
        Volume and Fluor1 mean, will plot the
        lineage oscilations using 'plot_lineage_derivative()'
            
    Returns
    --------------
    None
    """
    root_dict, lineage_array, len_dict = longer_lineages(df_3d = df_3d,
                                                            reverse_lineages = reverse_lineages,
                                                            column = column, 
                                                            method = method)

    # pearson correlation
    corr_arr = np.array([np.array(root_dict[mother_cell][i][:-1]) for i in root_dict[mother_cell]]) 
    corr = np.corrcoef(corr_arr)

    # plot params
    colors = ['red', 'k', 'purple', 'gray']
    lines = ['--','--','-.','-.']
    lines = [Line2D([0], [0], color=colors[i], linewidth=1, linestyle=lines[i]) for i in range(len(colors))]
    labels = [r'$\frac{dV}{dt}$', r'$\frac{dF}{dt}$', 'Volume', 'Fluor1 mean']

    temp_lineage = []
    for lineage in lineage_array:
        if lineage[0] == mother_cell and len_dict[mother_cell] == len(lineage):
            temp_lineage.append(lineage)
    plot_labels = [lineage[-1] for lineage in temp_lineage[::2]]
    lineage_plot_size = int(len(root_dict[mother_cell].keys())/2)
    corr_plot_size = int(len(root_dict[mother_cell].keys()))
    len_lineage = len_dict[mother_cell]
    lineage_plot_labels = ['Cell {}'.format(cell+1) for cell in range(len_lineage)]
    lineage_plot_keys = list(root_dict[mother_cell].keys())

    # corr heatmap plot
    fig,ax = plt.subplots(1,2, figsize = (16,5))
    sns.heatmap(corr[::2,::2],ax=ax[0],xticklabels = plot_labels,yticklabels = plot_labels)
    ax[0].set_title('Correlation Heatmap',fontsize = 17)
    ax[0].set_xlabel('Lineages',fontsize = 15)
    ax[0].set_ylabel('Lineages',fontsize = 15)

    # lineage corr plot
    ax[1].plot(root_dict[mother_cell][lineage_plot_keys[0]],'x:',color = colors[0],label='Lineage {} \nCorr: {}'.format(plot_labels[0],np.round(corr[:,0][0],2)))
    if corr_plot_size > 2:
        ax[1].plot(root_dict[mother_cell][lineage_plot_keys[3]],'o-.',color = colors[1],label='Lineage {} \nCorr: {}'.format(plot_labels[1],np.round(corr[:,0][3],2)))
    if corr_plot_size > 4:
        ax[1].plot(root_dict[mother_cell][lineage_plot_keys[5]],'*--',color = colors[2],label='Lineage {} \nCorr: {}'.format(plot_labels[2],np.round(corr[:,0][5],2)))
    if corr_plot_size > 6:
        ax[1].plot(root_dict[mother_cell][lineage_plot_keys[7]],'v--',color = colors[3],label='Lineage {} \nCorr: {}'.format(plot_labels[3],np.round(corr[:,0][7],2)))
    ax[1].set_title(column,fontsize = 17)
    ax[1].set_xticklabels(lineage_plot_labels,rotation=45)
    ax[1].set_xticks(range(len(lineage_plot_labels)-1))
    ax[1].set_ylabel('Fluor1 mean variance',fontsize = 15)
    ax[1].set_xlabel('Cells',fontsize = 15)
    ax[1].legend()
    plt.show()

    # lineage plot
    if derivative_plot:
        limit_plot = len(temp_lineage[::2])-1
        fig,ax = plt.subplots(lineage_plot_size,1, figsize = (12.5,6+lineage_plot_size),constrained_layout = True)
        fig.suptitle('Lineage of Cell {}'.format(mother_cell),fontsize = 17)
        for idx,lineage in enumerate(temp_lineage[::2]):
            if lineage_plot_size != 1:
                df_lineage,extremes,plot_params = lineage_derivative(df_3d,temp_lineage[::2][idx],column='Volume',order=4)
                plot_lineage_derivative(df_lineage,extremes,plot_params,x_axis='Time (fps)',ax=ax[idx])
                ax[idx].set_title('')
                if limit_plot != idx:
                    ax[idx].set_xlabel('')
                    ax[idx].set_xticks([])
            else:
                df_lineage,extremes,plot_params = lineage_derivative(df_3d,temp_lineage[::2][idx],column='Volume',order=4)
                plot_lineage_derivative(df_lineage,extremes,plot_params,x_axis='Time (fps)',ax=ax)
                ax.set_title('')
    else:
        fig,ax = plt.subplots(lineage_plot_size,1, figsize = (12.5,3+lineage_plot_size),constrained_layout = True)
        count,temp_min,temp_max = 0,100000,0
        for idx,lineage in enumerate(temp_lineage[::2]):
            if temp_max < df_3d[df_3d['Cell ID']==lineage[-2]]['Time (fps)'].values[-1]:
                temp_max = df_3d[df_3d['Cell ID']==lineage[-2]]['Time (fps)'].values[-1]
            if df_3d[df_3d['Cell ID']==lineage[0]]['Time (fps)'].values[0] < temp_min:
                temp_min = df_3d[df_3d['Cell ID']==lineage[0]]['Time (fps)'].values[0]

        for idx,lineage in enumerate(temp_lineage[::2]):
            if lineage_plot_size != 1:
                ax2 = ax[count].twinx()
                df_lineage,extremes,plot_params = lineage_derivative(df_3d,temp_lineage[::2][idx],column='Volume',order=4)
                
                y_axis = pre.MinMaxScaler((0,50)).fit_transform(df_lineage['Smoothed Derivative'].values.reshape(-1, 1))
                ax[count].plot(df_lineage['Time (fps)'].values, y_axis,'--', color='red')
                df_lineage,extremes,plot_params = lineage_derivative(df_3d,temp_lineage[::2][idx],column='Fluor1 sum',order=4)
                y_axis = pre.MinMaxScaler((0,50)).fit_transform(df_lineage['Smoothed Derivative'].values.reshape(-1, 1))

                ax[count].plot(df_lineage['Time (fps)'].values, y_axis,'--', color='k')
                ax[count].legend(lines, labels,loc='upper left',ncols = 4)
            else:
                ax2 = ax.twinx()
                df_lineage,extremes,plot_params = lineage_derivative(df_3d,temp_lineage[::2][idx],column='Volume',order=4)
                
                y_axis = pre.MinMaxScaler((0,50)).fit_transform(df_lineage['Smoothed Derivative'].values.reshape(-1, 1))
                ax.plot(df_lineage['Time (fps)'].values, y_axis,'--', color='red')
                
                df_lineage,extremes,plot_params = lineage_derivative(df_3d,temp_lineage[::2][idx],column='Fluor1 sum',order=4)
                y_axis = pre.MinMaxScaler((0,50)).fit_transform(df_lineage['Smoothed Derivative'].values.reshape(-1, 1))

                ax.plot(df_lineage['Time (fps)'].values, y_axis,'--', color='k')
                ax.legend(lines, labels,loc='upper left',ncols = 4)

            temp_x,temp_y,scale_y = [],[],[]
            for cell in lineage:
                temp_x.extend(df_3d[df_3d['Cell ID']==cell]['Time (fps)'].values)
                temp_y.extend(df_3d[df_3d['Cell ID']==cell]['Volume'].values)
                scale_y.extend(df_3d[df_3d['Cell ID']==cell]['Fluor1 mean'].values)
                # ax2.plot(df_3d[df_3d['Cell ID']==cell]['Time (fps)'].values,df_3d[df_3d['Cell ID']==cell]['Volume'].values,
                # 				':',color = 'brown') 
                ax2.plot(df_3d[df_3d['Cell ID']==cell]['Time (fps)'].values,df_3d[df_3d['Cell ID']==cell]['Fluor1 mean'].values,
                                ':',color = 'gray') 
                if lineage_plot_size != 1 and count != lineage_plot_size-1:
                    ax[count].set_xticks([])
                ax2.set_xlim(temp_min-15,temp_max+15)
            ax2.plot(temp_x,pre.MinMaxScaler((min(scale_y),max(scale_y))).fit_transform(np.asarray(temp_y).reshape(-1, 1)),':',color = 'purple')
            if count != lineage_plot_size-1:
                count +=1
                
        fig.text(-0.035, 0.5, r'$\frac{dV}{dt}\;and\;\frac{dF}{dt}$', va='center', rotation='vertical',fontsize = 15,zorder=2)
        fig.text(1.02, 0.5, 'Volume and Fluor1 mean', va='center', rotation=-90,fontsize = 15,zorder=2)
        if lineage_plot_size != 1:
            ax[count].set_xlabel('Time (min)',fontsize = 15)
        else:
            ax.set_xlabel('Time (min)', fontsize = 15)
        fig.suptitle('Lineage of Cell {}'.format(mother_cell),fontsize = 17)
        plt.show()

def plot_distance_minima_lineage(df_3d,df_2d,lineages_list,order=9,ax = None):
    """
    Function to plot the distance between the minimas for
    an specifics lineages.

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    df_2d : DataFrame
        DataFrame of 2D data
    lineages_list : list
        List with the lineages to plot, its better to use
        the output of 'lineage_corr_filter()' 
    order : int (optional)
        The order of the polynomial used to fit the 
        savgol_filter from scipy in 'lineage_derivative()'
        higher order is better to find more minimas.
    ax : nd.array (optional)
        Pass the ax for subplot.
            
    Returns
    --------------
    None
    """
    size_birth,minima_volume = [],[] 
    for lineage in lineages_list:
        df,peaks,_ = lineage_derivative(df_3d,lineage,order=order,column='Fluor1 sum')
        for cell,idx in peaks['minima'].items():
            if idx !=0:
                size_birth.extend(df_2d[df_2d['Cell ID']==cell]['Volume birth'].values)
                minima_volume.append(df['Smoothed Volume'].values[idx])
    
    model = LinearRegression()
    model.fit(np.asarray(size_birth).reshape(len(size_birth), 1),
              np.asarray(minima_volume).reshape(len(minima_volume)))

    if ax is None:
        plt.plot(size_birth,minima_volume,'.')
        plt.plot(size_birth,model.predict(np.asarray(size_birth).reshape(len(size_birth), 1)),'gray',label = 'Coef: ' + str(model.coef_[0])[:5])
        plt.title('Volume Minima by Birh Size',fontsize = 17)
        plt.xlabel(r'Volume at birth [$\mu m^3$]',fontsize = 15)
        plt.ylabel(r'Volume minima [$\mu m^3$]',fontsize = 15)
        plt.legend()
    else:
        out = []
        out.append(ax.plot(size_birth,minima_volume,'.'))
        out.append(ax.plot(size_birth,model.predict(np.asarray(size_birth).reshape(len(size_birth), 1)),'gray',label = 'Coef: ' + str(model.coef_[0])[:5]))
        out.append(ax.set_title('Minima Distance',fontsize = 17))
        out.append(ax.set_xlabel(r'Volume at birth [$\mu m^3$]',fontsize = 15))
        out.append(ax.set_ylabel(r'Volume minima [$\mu m^3$]',fontsize = 15))
        out.append(ax.legend())

        return out

def diff_minima(df_3d,df_2d,order = 3,pre = None,pos = None, adjust = False):
    """
    Function to plot the distance between the minimas of
    the daughter and mother by the mothe volume. The 'pre' 
    and 'pos' to be used, both must be different than 'None'.
    It is possible to use a higher order together with adjust 
    order for each cell, the higher the order, the higher 
    the noise acceptance. Higher orders take longer and can 
    reach 20 minutes to run. The value of the minimum used 
    comes from the column 'Smoothed Volume', output of
    'lineage_derivative()'. Also, plot the histogram of 
    the diff.

    Parameters
    --------------
    df_3d : DataFrame
        DataFrame of 3D data
    df_2d : DataFrame
        DataFrame of 2D data
    order : int (optional)
        The order of the polynomial used to fit the 
        savgol_filter from scipy in 'lineage_derivative()'
        higher order is better to find more minimas, accept
        more noise. Default is 3.
    pre : int (optional)
        Save the data pre this moment in minutes
    pos : int (optional)
        Save the data pos this moment in minutes
    adjust : bool (optional)
        If true it's possible to pass a bigger order value
        and the algorithm will adjust for each cell if 
        needed and the order used for each cell is stored
        in the dictionary vol_dict['order']['cell id'].
        Default false
    plot : bool (optional)
        If true will plot the result of the analysis, if
        'pre' or 'pos' is None will plot just all the data.

    Returns
    --------------
    vol_dict : dict
        Dictionary with the following keys 'all', 'pre',
        'pos' with subkeys 'diff' and 'minima', 'diff' is 
        the difference between daughter and mother volume at
        minima and 'minima' the minima volume of the mother.
        The key 'order' store the order used for each cell and
        the value of 'pre' and 'pos' for plot purpose.
    """
    vol_dict = {}
    vol_dict['all'],vol_dict['pre'],vol_dict['pos'],vol_dict['order'] = {},{},{},{}
    volume_minima,volume_minima_pos,volume_diff_pos = [],[],[]
    volume_diff,volume_minima_pre,volume_diff_pre = [],[],[]
    size_birth,size_birth_pre,size_birth_pos = [],[],[]
    vol_dict['order']['pre'],vol_dict['order']['pos'] = pre,pos
    if adjust:
        for_range = 9
    else:
        for_range = 2

    for cell in natsorted(df_2d['Cell ID'].values):
        temp_order = order
        daughter_cell1 = df_2d[df_2d['Cell ID']==cell]['Daughter1 ID'].values[0]
        daughter_cell2 = df_2d[df_2d['Cell ID']==cell]['Daughter2 ID'].values[0]
        
        if len(df_3d[df_3d['Cell ID']==daughter_cell1]) > 0:
            for attempt in range(for_range):
                try:            
                    df,peaks,_= lineage_derivative(df_3d,[cell,daughter_cell1],order = temp_order)
                    vol_dict['order'][cell] = temp_order
                    if 0 not in list(peaks['minima'].values()):
                        temp = df['Smoothed Volume'].values[list(peaks['minima'].values())]
                        time_pre = df_3d[df_3d['Cell ID']==daughter_cell1]['Time (fps)'].values[-1]
                        time_pos = df_3d[df_3d['Cell ID']==cell]['Time (fps)'].values[0]
                        size_birth.extend(df_2d[df_2d['Cell ID']==cell]['Volume birth'].values)
                        diff_size = temp[1]-temp[0]
                        volume_minima.append(temp[0])
                        volume_diff.append(diff_size)
                        if pre is not None and pos is not None:
                            if time_pre < pre:
                                volume_minima_pre.append(temp[0])
                                size_birth_pre.extend(df_2d[df_2d['Cell ID']==cell]['Volume birth'].values)
                                volume_diff_pre.append(diff_size)
                            elif time_pos > pos:
                                volume_minima_pos.append(temp[0])
                                size_birth_pos.extend(df_2d[df_2d['Cell ID']==cell]['Volume birth'].values)
                                volume_diff_pos.append(diff_size) 
                except ValueError:
                    if adjust:
                        temp_order = temp_order-1
                    else:
                        continue
            # else:
            # 	pass

        elif len(df_3d[df_3d['Cell ID']==daughter_cell2]) > 0:
            for attempt in range(for_range):
                try:
                    df,peaks,_= lineage_derivative(df_3d,[cell,daughter_cell2],order = temp_order)
                    vol_dict['order'][cell] = temp_order
                    if 0 not in list(peaks['minima'].values()):
                        temp = df['Smoothed Volume'].values[list(peaks['minima'].values())]
                        time_pre = df_3d[df_3d['Cell ID']==daughter_cell2]['Time (fps)'].values[-1]
                        time_pos = df_3d[df_3d['Cell ID']==cell]['Time (fps)'].values[0]
                        size_birth.extend(df_2d[df_2d['Cell ID']==cell]['Volume birth'].values)
                        diff_size = temp[1]-temp[0]
                        volume_minima.append(temp[0])
                        volume_diff.append(diff_size)
                        if pre is not None and pos is not None:
                            if time_pre < pre:
                                volume_minima_pre.append(temp[0])
                                size_birth_pre.extend(df_2d[df_2d['Cell ID']==cell]['Volume birth'].values)
                                volume_diff_pre.append(diff_size)
                            elif time_pos > pos:
                                volume_minima_pos.append(temp[0])
                                size_birth_pos.extend(df_2d[df_2d['Cell ID']==cell]['Volume birth'].values)
                                volume_diff_pos.append(diff_size) 
                except ValueError:
                    if adjust:
                        temp_order = temp_order-1
                    else:
                        continue
            # else:
            # 	pass

    vol_dict['all']['diff'],vol_dict['all']['min'],vol_dict['all']['birth'] = volume_diff,volume_minima,size_birth
    vol_dict['pre']['diff'],vol_dict['pre']['min'],vol_dict['pre']['birth'] = volume_diff_pre,volume_minima_pre,size_birth_pre
    vol_dict['pos']['diff'],vol_dict['pos']['min'],vol_dict['pos']['birth'] = volume_diff_pos,volume_minima_pos,size_birth_pos

    return vol_dict

def plot_diff_minima(vol_dict,key = None, ax = None):
    """
    Function to plot the difference between the minimas.
    The function will plot the 'pre' and 'pos' moments
    automatically if present in the dictionary. Also,
    plot the histogram of the diff. It's possible to 
    choose the key to plot ('all','pre' or 'pos'), if
    you pass just the key will plot the cloud and the
    histogram, if you pass the 'key' and the 'ax' will
    return just the cloud plot. If the 'pre' os 'pos'
    are less than 2, the function will plot just the 
    all data vector, to plot those moments you should
    pass the 'key'.

    Parameters
    --------------
    vol_dict : dict
        Dict with the values of the diff minima, output
        of 'diff_minima()'.
    key : str (optional)
        key used to plot just one moment, 'all', 'pre'
        or 'por'
    ax : nd.array (optional)
        the ax position to plot

    Returns
    --------------
    None
    """
    if key is not None:
        model = LinearRegression()
        model.fit(np.asarray(vol_dict[key]['min']).reshape(len(vol_dict[key]['min']), 1),
                    np.asarray(vol_dict[key]['diff']).reshape(len(vol_dict[key]['diff'])))
        if ax is None:
            fig,ax = plt.subplots(1,2, figsize = (13,5))


            ax[0].plot(vol_dict[key]['min'],vol_dict[key]['diff'],'.')
            ax[0].plot(vol_dict[key]['min'],model.predict(np.asarray(vol_dict[key]['min']).reshape(len(vol_dict[key]['min']), 1)),
                    'gray',label = 'Coef: ' + str(model.coef_[0])[:5])
            ax[0].set_title('{} the time'.format(key.capitalize()),fontsize = 17)
            ax[0].set_xlabel(r'Volume minima [$\mu m^3$]',fontsize = 15)
            ax[0].set_ylabel(r'Daughter Vm - Mother Vm [$\mu m^3$]',fontsize = 15)
            ax[0].set_ylim(0,5)
            ax[0].set_xlim(.5,5.5)
            ax[0].legend()
            
            ax[1].hist(vol_dict[key]['diff'])
            if key == 'all':
                ax[1].set_title('{} the time'.format(key.capitalize()),fontsize = 17)
            else:
                ax[1].set_title('{} {} min'.format(key.capitalize(),vol_dict['order'][key]),fontsize = 17)
            ax[1].set_xlabel(r'Daughter Vm - Mother Vm [$\mu m^3$]',fontsize = 15)
        else:
            out = []
            out.append(ax.plot(vol_dict[key]['min'],vol_dict[key]['diff'],'.'))
            out.append(ax.plot(vol_dict[key]['min'],model.predict(np.asarray(vol_dict[key]['min']).reshape(len(vol_dict[key]['min']), 1)),
                        'gray',label = 'Coef: ' + str(model.coef_[0])[:5]))
            out.append(ax.set_title('{} the time'.format(key.capitalize()),fontsize = 17))
            out.append(ax.set_xlabel(r'Volume minima [$\mu m^3$]',fontsize = 15))
            out.append(ax.set_ylabel(r'Daughter Vm - Mother Vm [$\mu m^3$]',fontsize = 15))
            out.append(ax.set_ylim(0,5))
            out.append(ax.set_xlim(.5,5.5))
            out.append(ax.legend())	

            return out
    else:	
        if len(vol_dict['pre']['diff']) < 2 or len(vol_dict['pos']['diff']) < 2:
            fig,ax = plt.subplots(1,2, figsize = (13,5))
            model = LinearRegression()
            model.fit(np.asarray(vol_dict['all']['min']).reshape(len(vol_dict['all']['min']), 1),
                    np.asarray(vol_dict['all']['diff']).reshape(len(vol_dict['all']['diff'])))

            ax[0].plot(vol_dict['all']['min'],vol_dict['all']['diff'],'.')
            ax[0].plot(vol_dict['all']['min'],model.predict(np.asarray(vol_dict['all']['min']).reshape(len(vol_dict['all']['min']), 1)),
                    'gray',label = 'Coef: ' + str(model.coef_[0])[:5])
            ax[0].set_title('All the time',fontsize = 17)
            ax[0].set_xlabel(r'Volume minima [$\mu m^3$]',fontsize = 15)
            ax[0].set_ylabel(r'Daughter Vm - Mother Vm [$\mu m^3$]',fontsize = 15)
            ax[0].set_ylim(0,5)
            ax[0].set_xlim(.5,5.5)
            ax[0].legend()
            
            ax[1].hist(vol_dict['all']['diff'])
            ax[1].set_title('All the time Histogram',fontsize = 17)
            ax[1].set_xlabel(r'Daughter Vm - Mother Vm [$\mu m^3$]',fontsize = 15)	 
        else:
            y_min,y_max = [],[]
            fig,ax = plt.subplots(2,3, figsize = (17,10))

            model = LinearRegression()
            model.fit(np.asarray(vol_dict['pre']['min']).reshape(len(vol_dict['pre']['min']), 1),
                    np.asarray(vol_dict['pre']['diff']).reshape(len(vol_dict['pre']['diff'])))

            ax[0][0].plot(vol_dict['pre']['min'],vol_dict['pre']['diff'],'.',color = 'k')
            ax[0][0].plot(vol_dict['pre']['min'],model.predict(np.asarray(vol_dict['pre']['min']).reshape(len(vol_dict['pre']['min']), 1)),
                    'gray',label = 'Pre Coef: ' + str(model.coef_[0])[:5])
            ax[0][0].set_title('Pre {} min'.format(vol_dict['order']['pre']),fontsize = 17)
            ax[0][0].set_ylabel(r'Daughter Vm - Mother Vm [$\mu m^3$]',fontsize = 15)
            ax[0][0].set_xlabel(r'Volume minima [$\mu m^3$]',fontsize = 15)
            ax[0][0].set_ylim(0,5)
            ax[0][0].set_xlim(.5,5.5)
            ax[0][0].legend()

            ax[1][0].hist(vol_dict['pre']['diff'],color = 'k')
            ax[1][0].set_title('Pre {} Histogram'.format(vol_dict['order']['pre']),fontsize = 17)
            ax[1][0].set_xlabel(r'Daughter Vm - Mother Vm [$\mu m^3$]',fontsize = 15)
            y_min.append(ax[1][0].get_ylim()[0])
            y_max.append(ax[1][0].get_ylim()[1])

            model = LinearRegression()
            model.fit(np.asarray(vol_dict['pos']['min']).reshape(len(vol_dict['pos']['min']), 1),
                    np.asarray(vol_dict['pos']['diff']).reshape(len(vol_dict['pos']['diff'])))

            ax[0][1].plot(vol_dict['pos']['min'],vol_dict['pos']['diff'],'.',color = 'red')
            ax[0][1].plot(vol_dict['pos']['min'],model.predict(np.asarray(vol_dict['pos']['min']).reshape(len(vol_dict['pos']['min']), 1)),
                    'gray',label = 'Pos Coef: ' + str(model.coef_[0])[:5])
            ax[0][1].set_title('Pos {} min'.format(vol_dict['order']['pos']),fontsize = 17)
            ax[0][1].set_xlabel(r'Volume minima [$\mu m^3$]',fontsize = 15)
            ax[0][1].set_ylim(0,5)
            ax[0][1].set_xlim(.5,5.5)
            ax[0][1].legend()

            ax[1][1].hist(vol_dict['pos']['diff'],color = 'red')
            ax[1][1].set_title('Pos {} Histogram'.format(vol_dict['order']['pos']),fontsize = 17)
            ax[1][1].set_xlabel(r'Daughter Vm - Mother Vm [$\mu m^3$]',fontsize = 15)
            y_min.append(ax[1][1].get_ylim()[0])
            y_max.append(ax[1][1].get_ylim()[1])

            model = LinearRegression()
            model.fit(np.asarray(vol_dict['all']['min']).reshape(len(vol_dict['all']['min']), 1),
                    np.asarray(vol_dict['all']['diff']).reshape(len(vol_dict['all']['diff'])))

            ax[0][2].plot(vol_dict['all']['min'],vol_dict['all']['diff'],'.')
            ax[0][2].plot(vol_dict['all']['min'],model.predict(np.asarray(vol_dict['all']['min']).reshape(len(vol_dict['all']['min']), 1)),
                    'gray',label = 'Coef: ' + str(model.coef_[0])[:5])
            ax[0][2].set_title('All the time',fontsize = 17)
            ax[0][2].set_xlabel(r'Volume minima [$\mu m^3$]',fontsize = 15)
            ax[0][2].set_ylim(0,5)
            ax[0][2].set_xlim(.5,5.5)
            ax[0][2].legend()

            ax[1][2].hist(vol_dict['all']['diff'])
            ax[1][2].set_title('All the time Histogram',fontsize = 17)
            ax[1][2].set_xlabel(r'Daughter Vm - Mother Vm [$\mu m^3$]',fontsize = 15)
            y_min.append(ax[1][2].get_ylim()[0])
            y_max.append(ax[1][2].get_ylim()[1])

            ax[1][0].set_ylim(min(y_min),max(y_max)+100)
            ax[1][1].set_ylim(min(y_min),max(y_max)+100)
            ax[1][2].set_ylim(min(y_min),max(y_max)+100)
            ax[1][0].set_xlim(0,5.5)
            ax[1][1].set_xlim(0,5.5)
            ax[1][2].set_xlim(0,5.5)

            plt.tight_layout()

def plot_distance_minima(vol_dict,key,color = None,ax=None):
    """
    Function to plot the distance between the minimas.
    The function will plot the volume minima by the 
    volume at birth. It's necessary to pass the 'key' 
    perform the cloud plot ('pre','pos','all').

    Parameters
    --------------
    vol_dict : dict
        Dict with the values of the diff minima, output
        of 'diff_minima()'.
    key : str (optional)
        key used to plot just one moment, 'all', 'pre'
        or 'por'
    color : str (optional)
        Cloud color
    ax : nd.array (optional)
        the ax position to plot

    Returns
    --------------
    None
    """
    x = np.asarray(vol_dict[key]['birth'])
    y = np.asarray(vol_dict[key]['min'])
    model = LinearRegression()
    model.fit(x.reshape(len(x), 1),y.reshape(len(y)))
    if ax is None:
        plt.plot(x,y,'.',color = color)
        plt.plot(x,model.predict(x.reshape(len(x), 1)),'gray',label = 'Coef: ' + str(model.coef_[0])[:5])
        plt.xlabel(r'Volume at birth [$\mu m^3$]',fontsize = 15)
        plt.ylabel(r'Volume minima [$\mu m^3$]',fontsize = 15)
        if key == 'all':
            plt.title('{} cells'.format(key.capitalize()),fontsize = 17)
        else:
            plt.title('{} {} min cells'.format(key.capitalize(),vol_dict['order'][key]),fontsize = 17)
        plt.legend()
    else:
        out = []
        out.append(ax.plot(x,y,'.',color = color))
        out.append(ax.plot(x,model.predict(x.reshape(len(x), 1)),'gray',label = 'Coef: ' + str(model.coef_[0])[:5]))
        out.append(ax.set_xlabel(r'Volume at birth [$\mu m^3$]',fontsize = 15))
        out.append(ax.set_ylabel(r'Volume minima [$\mu m^3$]',fontsize = 15))
        if key == 'all':
            out.append(ax.set_title('{} cells'.format(key.capitalize()),fontsize = 17))
        else:
            out.append(ax.set_title('{} {} min cells'.format(key.capitalize(),vol_dict['order'][key]),fontsize = 17))
        out.append(ax.legend())
        return out

# !pip install session-info
# !pip install pipreqs
# !pip install nbconvert
# %cd code
# !jupyter nbconvert --output-dir="./reqs" --to script 200707.ipynb
# %cd reqs
# !pipreqs
# import session_info
# session_info.show()