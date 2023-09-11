# SensorResponse-inSilico
SensinSilico is a simulation software published together with the respective peer-reviewed article *xXXx published in the peer-reviewed journal _Analyst_. It studies the error propagation when combining two (electrochemical) sensors with different response times. This document is meant as an installation and operation manual to use the software and explains the relevant application’s functions.

### Background
In environmental sensing, especially in total parameter monitoring, it is required to monitor not one but at least two analytes simultaneously and ideally at the same location. However, due to differences in sensor characteristics, e.g., acquisition rate, sensor response, etc., it can be cumbersome to do this for multiple sensors.
Although this bottleneck is well known, it is common to adapt the measurement procedure to the faster sensor, accepting a (small) inaccuracy for the total parameter. However, no one has yet investigated the extent of the inaccuracy. Is the inaccuracy acceptable or does the error propagation that is initially accepted falsify the overall result and make it unusable? 

This software application allows the user to examine a given dual-sensor system, study the influence of different sensor responses and thereby optimize the measurement procedure. 

### Download GUI
Download the *application* folder including the subfolder *(picture)*. To run the software application, open the file named *[SensinSilico.py](SensinSilico.py)* with a Python editor or IDE of your choice and run the script.

For the software application to run properly, leave all the files in the folder together without restructuring or removing any parts. For software development and testing, we tested the software application with PyCharm and VS Code 1.74.2 (Universal).

### Exploring different response curves
In addition to the GUI which only allows a sensor response defined by the Gompertz response curve, we included a Jupyter Notebook *Code_snippes.ipynb* in which we included a step-by-step explanation of the simulation and allow the user to add other response curves for explorational reason. 
This Notebook and as well as the related python script *dbs_CodeSnippes.py* can be found in the folder *[scripts](scripts/)*. To be able to run the Notebook, both scripts have to be stored in the same folder. 


### Requirements
When using the source code, python3 and the following Python packages are required for execution: 
```
numpy (version 1.24.2), pandas (version 1.5.3), scipy (version 1.10.1), and seaborn (version 0.12.2)
``` 
In addition, when using the GUI, the following packages are required:
```
PyQt5 (version 5.15.9), and pyqtgraph (version 0.13.3)
```

The package versions listed are the versions used when the so+ware is released. Older versions might cause errors and bugs, whereas newer versions should not make much of a difference. If you have difficulties running the software, please do not hesitate to reach out via email at info@silviazieger.com.

### Manual and Instructions on how to use the software
Next to the readme file, a manual *Instruction.pdf* is uploaded. Here you will find a step-by-step guide on how to use the simulation software including a summary of the context and theory.

### Disclaimer
The software has been tested for both MacOS (MacBook Pro – Apple M1, Version 13.4) and Windows (Version 10.0.19045). The original published version is linked to the following DOI: [![DOI](https://zenodo.org/badge/417093693.svg)](https://zenodo.org/badge/latestdoi/417093693)
