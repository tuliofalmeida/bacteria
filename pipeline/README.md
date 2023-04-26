# Pipeline Documentation

This pipeline was developed to use SuperSegger + Omnipose more fluidly in Python. The pipeline transitions SuperSegger(MATLAB) -> Omnipose (Python) -> SuperSegger (MATLAB) automatically. The pipeline has 3 routines, which allow you to analyze only one Field of View (FOV), analyze multiple FOV (of the same experiment or not, it will analyze all FOV in the folder) and analyze multiple experiments.

## Requirements

1. [Python](https://www.python.org/downloads/) (extensively tested in version 3.9 inside conda env)
   * [Anaconda](https://www.anaconda.com/download/)
   * [Omnipose](https://github.com/kevinjohncutler/omnipose)
   * Basic libraries as: os,sys,time,shutil (native in python)
2. [MATLAB](https://fr.mathworks.com/products/matlab.html) (extensively tested on version R2022b, but the version must be compatible with python 3.9 our earlier)
   * [SuperSegger-Omnipose](https://github.com/tlo-bot/supersegger-omnipose)
3. [matlabengine](https://pypi.org/project/matlabengine/) (version that match with your Python and MATLAB)
4. Download the Python and MATLAB scripts here **LINK** (you can do this cloning this repo, also check the Installation topic)
5. OS, the pipeline should work on Windows and Linux (extensively tested in Ubuntu 20.04.5 LTS) 

## Installation

1. Install [Omnipose](https://github.com/kevinjohncutler/omnipose) (check their tutorial)
2. Inside the Omnipose environment you should install the matlabengine `pip install matlabengine`
3. Clone [SuperSegger](https://github.com/tlo-bot/supersegger-omnipose) repo and add it into your MATLAB [path - Add folder and subfolders](https://fr.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html)
4. Configure [Python inside MATLAB](https://fr.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
5. Add the scripts _pipeline.m,pipeline2.m_ and _pipelineSeg.m_ in the _batch_ folder
6. Add the script _loadConstants_pipeline.m_ in the _settings_ folder
7. Create a folder called _pipeline_ (or whatever name you want, you can also use the pipeline folder when downloading the repositories)
8. Put the _pipeline.py_ and _functions.py_ in the same the folder that you choosed

## Running the pipeline

This is a step by step guide to using the pipeline after installation, but before you run your data read the section "Before you start".

1. Open python terminal and navigate to the folder _pipeline_ (_cd_ command [tutorial](https://fernando-mc.github.io/python3-workshop/navigating-with-a-terminal.html))
2. Run the following code
```python
conda run -n omnipose --live-stream python pipeline.py
```
3. You will see the following message in the terminal and choose the routine:
```python
"""
Pipeline - SuperSegger-Omnipose

Chose the code to run: 
0 - One FOV Analysis
1 - Multiple FOV's Analysis 
2 - Multiple Experiments Analysis 
3 - To see the documentation
"""
```
* Input 0: Run the Pipeline for one Field of View (FOV) - One folder with images (FOV).
* Input 1: Run the Pipeline for multiple Field of Views (FOV's) - One folder (Experiment) with other folders (FOV's) with images.
* Input 2: Run the Pipeline for multiple experiments - One folder (All data) with other folders (Experiments) with folders (FOV's) with images.
* Input 3: Will print the documentation in the terminal.
4. After choosing the routine, you will need to pass the folder path:
```python
"""
Set the folder to run the analysis:
"""
```
You should pass as example: /home/documents/experiment_test/ (don't need to add quotation marks '')

5. The pipeline you ask if you want to format the images:
```python
"""
Format images (0 = No, 1 = Yes): 
"""
```
* Input 0: It will not format the images, they should already be organized with the two separate channels and with the names following the SuperSegger pattern.
* Input 1: The pipeline will modify the name of the file according to SuperSegger and will duplicate the images in channel 1 and channel 2 (it is doing this automatically because we do not have the data in contract phase, only fluorescence. This can be changed from the parameter _'double'_ to _False_, in the function _format_images()_ in the file _functions.py_).
<details>
  <summary>If you have more than one fluorescence channel!</summary>
You will have to modify this function, it is possible to use the same logical structure, continuing to fold the image with GFP and adding a line with a condition and just change the name of the mCherry images to channel 3 (variable "ref_cherry" in function _format_image()_).</details>

6. Finally, the pipeline will ask if you want to align the images:
```python
"""
Align images to the first (0 = No, 1 = Yes):
"""
```
This is an internal SuperSegger function to align all frames with the first frame (or you can change the frame to align in SuperSegger's settings in the _loadConstants_pipeline.m_ file). We didn't perform extensive tests for this parameter, but we got slightly better results aligning the images.
* Input 0: SuperSegger will not align images 
* Input 1: SuperSegger will align the images 

7. After that the pipeline will start processing the images (Check the section Before you start)

## Before you start!




