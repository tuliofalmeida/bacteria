# Pipeline Documentation

This pipeline was developed to use SuperSegger + Omnipose more fluidly in Python. The pipeline transitions SuperSegger(MATLAB) -> Omnipose (Python) -> SuperSegger (MATLAB) automatically.

## Requirements

1. Python  (extensively tested in version 3.9)
   * Anaconda
   * [Omnipose](https://github.com/kevinjohncutler/omnipose)
2. MATLAB (extensively tested on version R2022b, but the version must be with compatible python 3.9 our earlier)
   * [SuperSegger-Omnipose](https://github.com/tlo-bot/supersegger-omnipose)
3. [matlabengine](https://pypi.org/project/matlabengine/) (version that mach with your Python and MATLAB)
