# Pipeline Documentation

This pipeline was developed to use SuperSegger + Omnipose more fluidly in Python. The pipeline transitions SuperSegger(MATLAB) -> Omnipose (Python) -> SuperSegger (MATLAB) automatically. The pipeline has 3 routines, which allow you to analyze only one Field of View (FOV), analyze multiple FOV (of the same experiment or not, it will analyze all FOV in the folder) and analyze multiple experiments.

## Contents

- [Requirements](#requirements) 
- [Installation](#installation) 
- [Running the pipeline](#running-the-pipeline) 
- [Pipeline Flow](#pipeline-flow)
- [Before you start!](#before-you-start)

## Requirements

1. [Python](https://www.python.org/downloads/) (extensively tested in version 3.9 inside conda env)
   * [Anaconda](https://www.anaconda.com/download/)
   * [Omnipose](https://github.com/kevinjohncutler/omnipose)
   * Basic libraries as: os,sys,time,shutil,datetime (native in python) and natsort (to order/sort files)
2. [MATLAB](https://fr.mathworks.com/products/matlab.html) (extensively tested on version R2022b, but the version must be compatible with python 3.9 our earlier)
   * [SuperSegger-Omnipose](https://github.com/tlo-bot/supersegger-omnipose)
3. [matlabengine](https://pypi.org/project/matlabengine/) (version that match with your Python and MATLAB)
4. Download the Python and MATLAB scripts from this folder (you can do this cloning this repo, also check the Installation topic)
5. OS, the pipeline should work on Windows and Linux (extensively tested in Ubuntu 20.04.5 LTS) 

## Installation

1. Install [Omnipose](https://github.com/kevinjohncutler/omnipose) (check their tutorial)
2. Inside the Omnipose environment you should install the matlabengine `pip install matlabengine` and narsort `pip install natsort`
3. Clone [SuperSegger](https://github.com/tlo-bot/supersegger-omnipose) repo and add it into your MATLAB [path - Add folder and subfolders](https://fr.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html)
4. Configure [Python inside MATLAB](https://fr.mathworks.com/help/matlab/matlab_external/install-supported-python-implementation.html)
5. Add the scripts _pipeline.m,pipeline2.m_ and _pipelineSeg.m_ in the _batch_ folder
6. Add the script _loadConstants_pipeline.m_ in the _settings_ folder
7. Create a folder called _pipeline_ (or whatever name you want, you can also use the pipeline folder (this one) when downloading the repositories)
8. Put the _pipeline.py_ and _functions_pipeline.py_ in the same the folder that you choosed

## Running the pipeline

This is a step by step guide to using the pipeline after installation, but before you run your data read the section "Before you start".

1. Open python terminal and navigate to the folder _pipeline_ (_cd_ command [tutorial](https://fernando-mc.github.io/python3-workshop/navigating-with-a-terminal.html))
2. Run the following code
```python
conda run -n omnipose --live-stream python pipeline.py
```
3. You will see the following message in the terminal and choose the routine:

        Pipeline - SuperSegger-Omnipose

        Chose the code to run: 
        0 - One FOV Analysis
        1 - Multiple FOV's Analysis 
        2 - Multiple Experiments Analysis 
        3 - To see the documentation
  
 * Options:
    * Input 0: Run the Pipeline for one Field of View (FOV) - One folder with images (FOV).
    * Input 1: Run the Pipeline for multiple Field of Views (FOV's) - One folder (Experiment) with other folders (FOV's) with images.
    * Input 2: Run the Pipeline for multiple experiments - One folder (All data) with other folders (Experiments) with folders (FOV's) with images.
    * Input 3: Will print the documentation in the terminal (check functions documentation).

4. After choosing the routine, you will need to pass the folder path:

        Set the folder to run the analysis:

You should pass as example: /home/documents/experiment_test/ (don't need to add quotation marks '')

5. The pipeline you ask if you want to format the images:

        Format images (0 = No, 1 = Yes): 
 * Options:
      * Input 0: It will not format the images, they should already be organized with the two separate channels and with the names following the SuperSegger pattern.
      * Input 1: The pipeline will modify the name of the file according to SuperSegger and will duplicate the images in channel 1 and channel 2 (it is doing this automatically because we do not have the data in contract phase, only fluorescence. This can be changed from the parameter _'double'_ to _False_, in the function _format_images()_ in the file _functions_pipeline.py_).

<details>
  <summary>If you have more than one fluorescence channel!</summary>
You will have to modify this function, it is possible to use the same logical structure, continuing to fold the image with GFP and adding a line with a condition and just change the name of the mCherry images to channel 3 (variable "ref_cherry" in function _format_image()_).</details>

6. Also, the pipeline will ask if you want to align the images:

        Align images to the first (0 = No, 1 = Yes):

This is an internal SuperSegger function to align all frames with the first frame (or you can change the frame to align in SuperSegger's settings in the _loadConstants_pipeline.m_ file). We didn't perform extensive tests for this parameter, but we got slightly better results aligning the images.
 * Options:
      * Input 0: SuperSegger will not align images 
      * Input 1: SuperSegger will align the images 

7. Finally, you will have to choose which model you want to use:

        Which model do you want to use? (0 = bact_phase_omni, 1 = bact_fluor_omni):
 * Options:
      * Input 0: to use _bact_phase_omni_ 
      * Input 1: to use _bact_fluor_omni_ 

The bacht_phase_omni model is trained using phase contrast images and can work well in fluoresence data using inverted images, the code already perform the invertion. Also we need to use a higher value for mask parameter (using 2.4), this model can work well in experiments were you have a big change in fluorescence intensity and/or in focus. The bact_fluor_omni is trained using fluorescence images, and work well in experiments with little variation in focus and fluorescence.

8. After these steps, the pipeline will start processing the images (Check the section Before you start)

## Pipeline Flow

<p align="center">
  <img width="500" height="600" src="https://github.com/tuliofalmeida/bacteria/blob/main/pipeline/pipeline.png">
</p>

## Before you start!

This section is to give more information about the pipeline, how to organize the data, possible errors, how to improve the result...

<details>
  <summary>How to organize the data?</summary>
    
  Initially the FOV folders should be named in sequence as "01,02,03,04,05...10,11,12..."to facilitate the writing of the log by the pipeline and for the user to have a notion of which FOV is being analyzed. I believe this has already been solved and the names can be passed as "1,2,3,4,5...10,11,12..." this organization in numbers is important for SuperSegger and the FOV can be checked after the library processing in the 'FOV' column as "xy01,xy02,xy03...". Recently we had an error when analyzing an experiment with problems, each FOV had a different number of images. This led to errors. So before starting the process we probably need to adjust the number of images in each FOV. Unfortunately, we didn`t have time to debug this problem.
  
  To analyze only one FOV
  
    .
    ├── ...
    ├── experiment_folder             # Experiment folder with only one FOV
    │   ├── 01                        # FOV folder --> you must pass THIS FOLDER to the pipeline <--
    │       ├── img0001.tiff          # Inside the folder you should have the images in sequence
    |       ├── img0002.tiff  
    |       ├── img0003.tiff
    |       └── ...                                                                                                  
    └── ...
    
  To analyze multiple FOVs

    .
    ├── ...
    ├── experiment_folder             # Experiment folder with many FOVs --> you must pass THIS FOLDER to the pipeline <--
    │   ├── 01                        # FOV Folder
    │       ├── img0001.tiff          # Inside the folder you should have the images in sequence
    |       ├── img0002.tiff  
    |       └── ...  
    |   ├── 02                        # FOV Folder
    │       ├── img0001.tiff          # Inside the folder you should have the images in sequence
    |       ├── img0002.tiff  
    |       └── ...   
    |   ├── 03                        # FOV Folder
    │       ├── img0001.tiff          # Inside the folder you should have the images in sequence
    |       ├── img0002.tiff  
    |       └── ...  
    └── ...

   To analyze multiple Experiments          
 
    .                                                                                                                            
    ├── all_experiments_folder        # Folder with all experiments --> you must pass THIS FOLDER to the pipeline <--
    |   ├── experiment_1_folder       # Experiment folder with many FOVs
    │       ├── 01                    # FOV Folder
    │           ├── img0001.tiff      # Inside the folder you should have the images in sequence
    |           ├── img0002.tiff  
    |           └── ...  
    |       ├── 02                    # FOV Folder
    │           ├── img0001.tiff      # Inside the folder you should have the images in sequence
    |           ├── img0002.tiff  
    |           └── ...   
    |       ├── 03                    # FOV Folder
    │           ├── img0001.tiff      # Inside the folder you should have the images in sequence
    |           ├── img0002.tiff  
    |           └── ...   
    |   ├── experiment_2_folder       # Experiment folder with many FOVs 
    │       ├── 01                    # FOV Folder
    │           ├── img0001.tiff      # Inside the folder you should have the images in sequence
    |           ├── img0002.tiff  
    |           └── ...  
    |       ├── 02                    # FOV Folder
    │           ├── img0001.tiff      # Inside the folder you should have the images in sequence
    |           ├── img0002.tiff  
    |           └── ...   
    |       ├── 03                    # FOV Folder
    │           ├── img0001.tiff      # Inside the folder you should have the images in sequence
    |           ├── img0002.tiff  
    |           └── ...
    |   ├── experiment_3_folder       # Experiment folder with many FOVs 
    │       ├── 01                    # FOV Folder
    │           ├── img0001.tiff      # Inside the folder you should have the images in sequence
    |           ├── img0002.tiff  
    |           └── ...  
    |       ├── 02                    # FOV Folder
    │           ├── img0001.tiff      # Inside the folder you should have the images in sequence
    |           ├── img0002.tiff  
    |           └── ...   
    |       ├── 03                    # FOV Folder
    │           ├── img0001.tiff      # Inside the folder you should have the images in sequence
    |           ├── img0002.tiff
    |           └── ... 
    └── ...
</details>

<details>
  <summary>Possible erros</summary>
  
  Several things can produce an error in this code, because it has not been tested with much variability. If you get an error following this tutorial, check if the data is organized correctly, if all paths are correct, if the choices made in the pipeline are correct (the value of the inputs). Then, try to run it again ! Another thing that can be involved is different versions of packages needed by python/matlab and the OS used, the installation should be done carefully. If errors persist please open an [Issue here on GitHub](https://github.com/tuliofalmeida/bacteria/issues/new)! When creating an Issue try to be as specific as possible, put the complete error (copy and paste), tell how it happened, add a screenshot of the terminal.
  I believe SuperSegger is designed to analyze a single FOV at a time, so when trying to run multiple FOVs the folders must have the same amount of images.
  
</details>

<details>
  <summary>Improving Pipeline Results (SuperSegger-Omnipose)</summary>
  
  1. There are different ways to improve the pipeline result. The more controlled the experiment is (focus, stability, fluorescence level...) the easier it is for both software to perform well. A classic rule of machine learning and deep learning is GIGO (Garbage In Garbage Out), if the data isn't good, the model cannot perform well.
  2. For reasons of time, it was not possible to train a specific network for use. But, it is possible to adjust some parameters to improve the Omnipose segmentation process. If the two nets with pre-set parameters don't perform well on your data, you can test different networks (for E.coli basically _bact_phase_omni_ and _bact_fluor_omni_) and parameters (for us the most important was the _mask_) using the [Omnipose GUI](https://github.com/kevinjohncutler/omnipose). Using the interface you can add an image (by dragging it to the interface) and adjust the parameters for segmentation. After the segmentation if you like the results, in the bottom left corner of the screen you will see the code line to use these parameters. These parameters must be added manually in the pipeline code (I advise to test them in one FOV using the One FOV mode). To modify the network parameters used, you must change in the _functions_pipeline.py_ file the function you want to use (_one_fov()_,_multiple_fovs()_ or _multiple_experiments()_) and look for the line that contains a condition `if model == 0` (for the _bact_phase_omni_ model) and `elif model == 1` (for the _bact_fluor_omni_ model). After this condition you will have a command `os.system('python -m omnipose ...')` and inside this function you should make the changes. Type this to open the GUI:
  
          conda activate omnipose   # or the name that you choosed to the omnipose env
          python -m omnipose        # this will open the interface

  3. If you adjust the Omnipose parameters it does not improve the results. You may need to retrain the SuperSegger algorithm and an Omnipose network. For computational power reasons and convenience, I suggest you start by trying to retrain SuperSegger for the data ([how to train SuperSegger](https://github.com/wiggins-lab/SuperSegger/wiki/Creating-your-own-constants-file)), if not, you may need to train an Omnipose network. Training an Omnipose network is a little more complex, as you will need to create the reference masks for your images (labeling) and without a GPU this process can be time consuming. To do this, it will be necessary to interact with the Omnipose codes to train the data and the necessary organization of the images ([here](https://github.com/kevinjohncutler/omnipose#how-to-train-omnipose)). Also it might be interesting to read about the [CellPose](https://github.com/MouseLand/cellpose) process and there is also a project/repository that has tutorials on how to create masks and train networks that might be useful [ZeroCostDL4Mic](https://github.com/HenriquesLab/ZeroCostDL4Mic/wiki) (An alternative to the lack of GPU is to carry out the training using GoogleColab, they explain how to do this in this repository).
  4. Something that can be added is a pre-processing step to remove noise and improve the focus of images [paper](https://www.sciencedirect.com/science/article/pii/S2001037022001192). 

</details>

<details>
  <summary>GPU</summary>
  
  If you have a GPU, just configure it (check the [Omnipose Tutorial](https://github.com/kevinjohncutler/omnipose#gpu-support)). So I believe that the omnipose will automatically use it without the need to change the code, if necessary it should be changed the same way you select the network parameters (same line of code).
  
</details>

<details>
  <summary>Pipeline Log</summary>
  
  After running the analyses the pipeline saves in the project folder a .txt file with the time spent by the code in each FOV.
  
</details>

