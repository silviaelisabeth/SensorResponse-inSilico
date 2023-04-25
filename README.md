# SensorResponse-inSilico
This is an in-silico study to simulate bias due to delayed sensor response and the error propagation when working with two (electrochemical) sensors.

### Background
In environmental sensing, especially in total parameter monitoring, it is required to monitor not one but at least two analytes simultaneously and ideally at the same location. However, due to differences in sensor characteristics, e.g., acquisition rate, sensor response, etc., it can be cumbersome to do this for multiple sensors.
Although this bottleneck is well known, it is common to adapt the measurement procedure to the faster sensor, accepting a (small) inaccuracy for the total parameter. However, no one has yet investigated the extend of the inaccuracy. Is the inaccuracy acceptable or does the error propagation that is initially accepted falsify the overall result and make it unusable. 

This software application allows the user to examine a given dual sensor system, study the influence of different sensor responses and thereby optimize the measurement procedure. 

### Disclaimer
The application was tested on windows and mac (M1 Ventura 13.3).
