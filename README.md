# MATLAB code to replicate main results from "Functional ultrasound neuroimaging reveals mesoscopic organization of saccades in the lateral intraparietal area"


## Project Organization
------------
    
    ├── LICENSE                         <- license
    ├── LICENSE.md                      <- license
    ├── main_analysis.m                 <- Top-level MATLAB script that generates the data figures used in paper
    ├── preview_data_as_movie.m         <- MATLAB script that creates GUI showing the acquired fUSI doppler data
    ├── README.md                       <- The top-level README for developers using this project
    ├── setup.m                         <- MATLAB script to initialize git repository. Only needs to be run once when repo is first created.

------------
## Associated dataset
Available from CaltechData: [Click here to download](https://doi.org/10.22002/p5jan-02r60). 

Packaged as .zip file in data hierarchy used by MATLAB functions. Once downloaded, extract contents from the .zip file.
   
    ├── ProjectRecord_paper.json                            JSON file storing metadata for each fUSI experimental session
    ├── doppler                                             fUSI Power Doppler Data for each session
    |   ├──doppler_S{*A}_R{*B}+normcorre.mat                Task-aligned Doppler data where {*A} is the session number and {*B} is the run number.
    |   ├──doppler_S{*A}_R{*B}_allTrials+normcorre.mat      Task-aligned Doppler data from [real-time decoding project](https://doi.org/10.1038/s41593-023-01500-7) where {*A} is the session number and {*B} is the run number.
    |   └──dopplerContinuous_S{*A}_R{*B}+normcorre.mat      Continuous Doppler data where {*A} is the session number and {*B} is the run number.
    ├── output                                              Where processed data, figures, and movies will be saved
    |   ├──across session analyses                          Where all processed data for across-session analysis will be stored
    |   ├──decoding                                         Where all processed data for decoding analysis will be stored
    |   ├──glm                                              Where all processed data for GLM analysis will be stored
    |   └──glm                                              Where all processed data for Dice-Sørensen analysis will be stored
    ├── ROIs                                                Where all Regions Of Interest (ROIs) will be saved
    |   ├──Anatomical Polygons                              Where Anatomical Polygon files for specific sessions/runs will be stored
    |   └──SulcusMaps                                       Where sulcus maps for specific sessions/runs will be stored



------------
## Other information
Tested on MATLAB R2024b on Mac Sequoia 15.3.2
Please send feedback and suggestions to: [wsgriggs@gmail.com](mailto:wsgriggs@gmail.com)

Zenodo archive - 
[![DOI](https://zenodo.org/badge/755770495.svg)](https://zenodo.org/doi/10.5281/zenodo.15122175)

==============================

### In publications, please reference:
Griggs, W.S., Norman, S.L., Tanter, M., Liu, C., Christopoulos, V., Shapiro, M.G., and Andersen, R.A. Functional ultrasound neuroimaging reveals mesoscopic organization of saccades in the lateral intraparietal area. Nature Communications. _In review._ (See BioRxiv version [here](https://www.biorxiv.org/content/10.1101/2024.06.28.600796v1))

------------
## Getting set up

1. Clone the repository 

    ```bash
    $ git clone https://github.com/wsgriggs2/directional_tuning
    ```
    *note: If you are interested in contributing, please get in touch with the authors (Whitney or Sumner). In the meantime, feel free to submit issues and/or feature requests.*

2. Download paired dataset from [CaltechData](https://doi.org/10.22002/p5jan-02r60) and unzip in a location that is convenient.

3. Run `setup.m` within the root folder. This will set up any necessary pathing.

4. Try out any of the analyses included within the `main_analysis.m` file or look at Power Doppler movies using `preview_data_as_movie.m`!

## Dependencies

This repository makes use of multiple MATLAB toolboxes including
* Bioinformatics Toolbox
* Curve Fitting Toolbox
* Image Processing Toolbox
* Medical Imaging Toolbox
* Navigation Toolbox
* Optimization Toolbox
* Parallel Processing Toolbox
* Signal Processing Toolbox
* Statistics and Machine Learning Toolbox
* System Identification Toolbox
 

## Authors
* **Whitney Griggs** - [GitHub](https://github.com/wsgriggs2) | [wsgriggs@gmail.com](mailto:wsgriggs@caltech.edu)
* **Sumner Norman** - [GitHub](https://github.com/sumner15) | [sumnern@caltech.edu](mailto:sumnern@caltech.edu)
* **Vasileios Christopolous** - [vchristo@usc.edu](mailto:vchristo@usc.edu) 


