# canvas_demeter

This is an unfinished repo containing code that was used during design and development of the CANVAS cubesat. 

There are three folders:
1. capacitance_req: code to calculate the crossover between capacitive and resistive coupling between the antenna and LEO plasma. 
We used this code to set our maximum allowable capacitance between the antenna housing and antenna. 
This code requires pyglow to get approximate density values of LEO (more details to come)

2. demeter_data_analysis: used to parse large quantities of DEMETER burst mode data to size the gain for the CANVAS instrument. 
I have several TB of DEMETER data, contact me if you want to use it. 

3. SVD: directory in progress, I'm working on using SVD methods to anticipate uncreatinites in the CANVAS data processing etc. 
