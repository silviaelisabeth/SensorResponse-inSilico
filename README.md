# SensorResponse-inSilico
This is an in-silico study to simulate bias due to delayed sensor response and the error propagation when working with two (electrochemical) sensors.

### Background
In environmental sensing, especially in total parameter monitoring, it is required to monitor not one but at least two analytes simultaneously and ideally at the same location. However, due to differences in sensor characteristics, e.g., acquisition rate, sensor response, etc., it can be cumbersome to do this for multiple sensors.
Although this bottleneck is well known, it is common to adapt the measurement procedure to the faster sensor, accepting a (small) inaccuracy for the total parameter. However, no one has yet investigated the extend of the inaccuracy. Is the inaccuracy acceptable or does the error propagation that is initially accepted falsify the overall result and make it unusable. 

This software application allows the user to examine a given dual sensor system, study the influence of different sensor responses and thereby optimize the measurement procedure. 

### Download 
Download the *application* folder including the subfolder (picture)*.* To run the software application, open the file named *[SensinSilico.py](http://SensinSilico.py)* with a Python editor or IDE of your choice and run the script.

For the software application to run properly, leave all the files in the folder together without restructuring or removing any parts. For software development and testing, we tested the software application with PyCharm and VS Code 1.74.2 (Universal).

### Requirements
When using the source code, python3 and the following Python packages are required for execution: 
Matplotlib (version 3.7.1), numpy (version 1.24.2), pandas (version 1.5.3), PyQt5 (version 5.15.9), scipy (version 1.10.1), and seaborn (version 0.12.2).

The package versions listed are the versions used when the software is released. Older versions might cause errors and bugs. If you have difficulties running the software, please do not hesitate to contact me via email under info@silviazieger.com.

### Manual and Instructions on how to use the software
Next to the readme file, a manual *Instruction.pdf* is uploaded. Here you will find a step by step guide on how to use the simulation software including a summary of the context and theory.

### Disclaimer
The application was tested on windows and mac (M1 Ventura 13.3).
