# Pipeline Documentation

This pipeline was developed to use SuperSegger + Omnipose more fluidly in Python. The pipeline transitions SuperSegger(MATLAB) -> Omnipose (Python) -> SuperSegger (MATLAB) automatically. The pipeline has 3 routines, which allow you to analyze only one Field of View (FOV), analyze multiple FOV (of the same experiment or not, it will analyze all FOV in the folder) and analyze multiple experiments.

## Requirements

1. [Python](https://www.python.org/downloads/) (extensively tested in version 3.9 inside conda env)
   * [Anaconda](https://www.anaconda.com/download/)
   * [Omnipose](https://github.com/kevinjohncutler/omnipose)
   * Basic libraries as: os,sys,time,shutil (native in python)
2. [MATLAB](https://fr.mathworks.com/products/matlab.html) (extensively tested on version R2022b, but the version must be with compatible python 3.9 our earlier)
   * [SuperSegger-Omnipose](https://github.com/tlo-bot/supersegger-omnipose)
3. [matlabengine](https://pypi.org/project/matlabengine/) (version that match with your Python and MATLAB)
4. Download the Python and MATLAB scripts here LINK (you can do this cloning this repo, also check the Installation topic)
5. OS, the pipeline should work on Windows and Linux (extensively tested in Ubuntu 20.04.5 LTS) 


## Installation

1. Install [Omnipose](https://github.com/kevinjohncutler/omnipose) (check their tutorial)
2. Inside the Omnipose environment you should install the matlabengine
3. Clone [SuperSegger](https://github.com/tlo-bot/supersegger-omnipose) repo and add it into your MATLAB [path - Add folder and subfolders](https://fr.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html)
4. Configure [Python inside MATLAB](https://fr.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
5. Add the scripts _pipeline.m,pipeline2.m_ and _pipelineSeg.m_ in the _batch_ folder
6. Add the script _loadConstants_pipeline.m_ in the _settings_ folder
7. Create a folder called _pipeline_ (or whatever name you want, you can also use the pipeline folder when downloading the repositories)
8. Put the _pipeline.py_ and _functions.py_ in the same the folder that you choosed
