# A drone-delivered aid simulation
A simulation to discover how effective a drone delivered aid system could be in a conflict or disaster
## Requirements
This project now runs using python 3.9.1

In order to run the program you will need to install the correct python version (remember to change the environment path), and also the required packages.
The required packages are found in *requirements.txt* You will also need to install the two *.whl* files which are included in the files

## How to set up
Start by copying all the files contained within this project into a folder on your machine. You then should check you have the correct version of python installed and the required packages as mentioned earlier. To install the packages, you need to open up the command prompt and install in this order (go [here](https://geoffboeing.com/2014/09/using-geopandas-windows/) if you are having issues as we are using geopandas):-

First remember to change directory to the folder where you copied the files

```pip install GDAL-3.2.1-cp39-cp39-win_amd64.whl```

and

```pip install Fiona-1.8.18-cp39-cp39-win_amd64.whl```

to install the *.whl* files, then

```pip install -r requirements.txt```

to install the packages required


Once this is complete, you should be able to run the program by opening up command prompt and following these steps:-

1. Change directory to the folder where you copied the files (file path given as an example)
* ```cd "C:/Users/Hugh-PC/Documents/Drone_Simulation"```
2. Run this command to start the program
* ```bokeh serve Drone_Simulation_Program.py --show```

This should start a bokeh server and a browser window should open. The console will start loading shapefiles and setting up the map of Syria
(one only goes this far to avoid paying google for dev access to their maps), you will know when the program is done as the map will pop up on 
the browser and the last readout from the console will look like this:-

![alt text](https://i.imgur.com/sA2SF8r.png)

## How to use the simulation
The interactive map is a full map of Syrian population centres (grey), roads and admin boundaries (black) and conflict locations (red dots) gathered from the ACLED database. 
By double clicking on the map you can place a "drone station" with it's range of coverage shown by a light circle (in this simulation you can place 10 stations):-

![alt text](https://i.imgur.com/x3OidzP.png)

Your aim should be to deploy the drone stations in range of the active conflict areas in order to cover as many incidents as possible. Once the stations are deployed, you can start the simulation by clicking on the button.
Drones then deploy from the stations to the incidents as they happen (in this version of simulation there are 101 drones per station), and you can see the live readout of the simulation from the console. 

Once the simulation has finished running, it will output how well your placement of stations has performed.
The task you have is to find an ideal placement of drone stations in order to cover as many of the red dots as possible, delivering as much aid as possible.
