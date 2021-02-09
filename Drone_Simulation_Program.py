"""
 This model should take information from conflict using the ACLED data tables based on a country,
 allow the user to place stations in that country, and simulate a drone response if the event is within range of the drones. 
 The outcome should be that aid is delivered to the location of the incident and increase the probability of 
 immediate survival/provide aid for those involved.

 Counters will allow you to see the amount of aid deployed and other information gathered from the simulation.
"""




#----------------------Here's the import list, quite a lot of imports as I had to keep adding them to get the plots to function-------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#
import pandas as pd
from bokeh.io import output_file, show
import geopandas as gpd
from bokeh.plotting import figure, save, gmap, show, curdoc
from bokeh.models import ColumnDataSource, HoverTool, LogColorMapper, LinearColorMapper, GeoJSONDataSource, GMapOptions, Column
import numpy as np
from osgeo import osr
from bokeh.models import GeoJSONDataSource
from bokeh.palettes import Viridis256 as palette
import bokeh as bokeh
from bokeh.events import DoubleTap, ButtonClick
from bokeh.layouts import widgetbox
from bokeh.models.widgets import Button
import mpu
import time
import geopy.distance
from datetime import datetime

#
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------




#----------------------------------------This method takes a .prj file and outputs the EPSG code for CRS-------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#
def esriprj2standards(shapeprj_path):
   prj_file = open(shapeprj_path, 'r')
   prj_txt = prj_file.read()
   srs = osr.SpatialReference()
   srs.ImportFromESRI([prj_txt])
   print 'Shape prj is: %s' % prj_txt
   print 'WKT is: %s' % srs.ExportToWkt()
   print 'Proj4 is: %s' % srs.ExportToProj4()
   srs.AutoIdentifyEPSG()
   print 'EPSG is: %s' % srs.GetAuthorityCode(None)
#
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------




#All of the following methods are used to get coordinates and line coordinates from polygons, points, multipolygons, multilines and such from shapefile geometry
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#
def getXYCoords(geometry, coord_type):
    """ Returns either x or y coordinates from  geometry coordinate sequence. Used with LineString and Polygon geometries."""
    if coord_type == 'x':
        return geometry.coords.xy[0]
    elif coord_type == 'y':
        return geometry.coords.xy[1]

def getPolyCoords(row, geom, coord_type):
    """Returns the coordinates ('x' or 'y') of edges of a Polygon exterior"""

    # Parse the exterior of the coordinate
    exterior = row[geom].exterior

    if coord_type == 'x':
        # Get the x coordinates of the exterior
        return list( exterior.coords.xy[0] )
    elif coord_type == 'y':
        # Get the y coordinates of the exterior
        return list( exterior.coords.xy[1] )

def getLineCoords(geometry, coord_type):
    """ Returns Coordinates of Linestring object."""
    return getXYCoords(geometry, coord_type)

def getPointCoords(geometry, coord_type):
    """ Returns Coordinates of Point object."""
    if coord_type == 'x':
        return geometry.x
    elif coord_type == 'y':
        return geometry.y

def multiGeomHandler(multi_geometry, coord_type, geom_type):
    """
    Function for handling multi-geometries. Can be MultiPoint, MultiLineString or MultiPolygon.
    Returns a list of coordinates where all parts of Multi-geometries are merged into a single list.
    Individual geometries are separated with np.nan which is how Bokeh wants them.
    # Bokeh documentation regarding the Multi-geometry issues can be found here (it is an open issue)
    # https://github.com/bokeh/bokeh/issues/2321
    """

    for i, part in enumerate(multi_geometry):
        # On the first part of the Multi-geometry initialize the coord_array (np.array)
        if i == 0:
            if geom_type == "MultiPoint":
                coord_arrays = np.append(getPointCoords(part, coord_type), np.nan)
            elif geom_type == "MultiLineString":
                coord_arrays = np.append(getLineCoords(part, coord_type), np.nan)
            elif geom_type == "MultiPolygon":
                coord_arrays = np.append(getPolyCoords(part, coord_type), np.nan)
        else:
            if geom_type == "MultiPoint":
                coord_arrays = np.concatenate([coord_arrays, np.append(getPointCoords(part, coord_type), np.nan)])
            elif geom_type == "MultiLineString":
                coord_arrays = np.concatenate([coord_arrays, np.append(getLineCoords(part, coord_type))], np.nan)
            elif geom_type == "MultiPolygon":
                coord_arrays = np.concatenate([coord_arrays, np.append(getPolyCoords(part, coord_type), np.nan)])

    # Return the coordinates
    return coord_arrays

def getCoords(row, geom_col, coord_type):
    """
    Returns coordinates ('x' or 'y') of a geometry (Point, LineString or Polygon) as a list (if geometry is LineString or Polygon).
    Can handle also MultiGeometries.
    """
    # Get geometry
    geom = row[geom_col]

    # Check the geometry type
    gtype = geom.geom_type

    # "Normal" geometries
    # -------------------

    if gtype == "Point":
        return getPointCoords(geom, coord_type)
    elif gtype == "LineString":
        return list( getLineCoords(geom, coord_type) )
    elif gtype == "Polygon":
        return list( getPolyCoords(geom, coord_type) )

    # Multi geometries
    # ----------------

    else:
        return list( multiGeomHandler(geom, coord_type, gtype) )
#
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------




#-------------------------------------------------Reading file location for ACLED data--------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#
print('Reading Conflict data:')
conflict_file_loc = pd.read_csv('Syria_2016-2019_Jul13_1430.csv')
print('done')  
#
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------




#---------------------------------------------Reading the shapefile data into the program-----------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#
print("Reading roads shapefile:")
roads_fp = r"hotosm_syr_roads_lines.shp"
roads = gpd.read_file(roads_fp)
print("done")

print("Reading admin boundaries shapefile:")
admins_fp = r"syr_admin1.shp"
admins = gpd.read_file(admins_fp)
print("done")

print('Reading populated areas shapefile:')
populated_areas_fp = r'syr_pplp_adm4_unocha.shp'
populated_areas = gpd.read_file(populated_areas_fp)
print('done')
#
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------Finished reading shapefile data-----------------------------------------------------




#-------------This transforms the populated area shapefile into x and y coordinate lists, makes it easier to plot in bokeh--------------
#---------------------------------------------------------------------------------------------------------------------------------------
#
print('Transforming populated areas data:')
populated_areas['x'] = populated_areas.apply(getCoords, geom_col='geometry', coord_type='x', axis=1)
populated_areas['y'] = populated_areas.apply(getCoords, geom_col='geometry', coord_type='y', axis=1)
print('done')
#
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------




#---------Here we need to make sure the CRS of all the files is equal so that they share the same space in the plot---------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#
print('Equalising CRS:')
CRS = admins.crs
roads = roads.to_crs(crs=CRS)
populated_areas = populated_areas.to_crs(crs=CRS)
print('done')
#
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------




#-------------------We need to convert the shapefiles to GeoJSON format so that we can plot them as maps in bokeh-----------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#
print('Converting all data to GEOJSON:')
admins_source = GeoJSONDataSource(geojson=admins.to_json())
roads_source = GeoJSONDataSource(geojson=roads.to_json())
p_df = populated_areas.drop('geometry', axis=1).copy()
psource = ColumnDataSource(p_df)
print('done')
#
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------




#Tools that we will need for the interactive map
TOOLS = "pan, wheel_zoom, box_zoom, reset, hover, save, tap"
  
#This creates an empty list for the drone station coordinates to be saved into
coordList=[]

#Empty column data source created for the drone stations
source = ColumnDataSource(data=dict(x=[], y=[]))  
 
#These were for testing
#init_coords = (36.929,35.427)
#coordList.append(init_coords)




#----------------------------------------Here's where the plot of the interactive map is made-------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#
print("Plotting figure..")
p = figure(title="An interactive map for using drone-delivered aid in Syria", tools=TOOLS, active_scroll = "wheel_zoom")
p.grid.grid_line_color = None
print('plotting admin boundaries..')
geoplot = p.patches('xs', 'ys', source=admins_source, fill_color='#c2b280', fill_alpha=0.9, line_width=0.5, line_color='black')
print('plotting habitated areas..')
c = p.circle('x', 'y', source=psource, size=4, color="#808080")
print('plotting road network..')
r = p.multi_line('xs', 'ys', source=roads_source, color='black', line_alpha=0.5, line_width=1.5)
print('done')

#I'm kind of cheating so that the usefulness of the program can be demonstrated, but this plots all locations of conflict on top of the map
print('plotting conflict locations "CHEATING"..')
p.circle(conflict_file_loc["LONGITUDE"], conflict_file_loc["LATITUDE"], size = 3, color = "red",alpha=0.2)
print('done')

#---------------------------------This is the part that plots the coordinates of the drone stations----------------------------------------
p.circle(source=source,x='x',y='y', radius=0.2351, fill_color="#FFFFFF", fill_alpha=0.4, line_color="#FFFFFF") 
#------------------------------------------------------------------------------------------------------------------------------------------




#------This is where the hovertooltip comes into play, where we can actively see the coordinates on the plot as we move the mouse----------
#------------------------------------------------------------------------------------------------------------------------------------------
#
hover = p.select_one(HoverTool)
hover.point_policy = "follow_mouse"
hover.renderers =[geoplot]
hover.tooltips = [("(Long, Lat)", "($x, $y)")]
#
#------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------
button_widget = Button(label = "Run Drone Sim")




#---------------------------------Distance calcualtion between two points using latitude and longitude----------------------------------
#------------------------------------------taking into account the curvature of the earth-----------------------------------------------
#
def distanceCalc(lat1,lon1,lat2,lon2):
    dist = mpu.haversine_distance((lat1, lon1), (lat2, lon2))
    #Distance is returned in km
    return(dist)
#
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------




#--------------------------Flight time calculation for one round trip of one drone, rounds to nearest int-------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#
def droneFlightTime(dist):
    #Edit the drone speed to see what happens if we use better drones
    droneSpeed = 140.0 #km/h
    timetaken = dist/droneSpeed
    #We need the time in seconds to work with a UNIX clock
    secs = timetaken*60*60
    return(secs)
#
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------




#------------------------------This initialises a list of drone numbers as long as the list of the stations-----------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#
def droneStationInit(coordList):
    drones_per_station = []
    for i in range(len(coordList)):
        drones_per_station.append(101)
    #print("drones per station:",drones_per_station)
    return(drones_per_station)
#
#---------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------




#--------------------------------------------This returns a list of stations in range of the incident-----------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------
#
def stationsInRange(lat1,lon1,coordList):
    stations_in_range = []
    for x,y in coordList:
        #print(x,y)
        if distanceCalc(y,x,lat1,lon1) < 26:
            stations_in_range.append([x,y])
        else:
            continue
    if stations_in_range == []:
        return None
    else:
        return(stations_in_range)
#
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------




#-----------------------------------This chooses the station in range with the maximum number of drones available------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#
def maximumStationDrones(stations_available, coordList, drones_available):
    drones_list=[]
    stations_deploying = []
    for w,x in coordList:
        for y,z in stations_available:
            if w == y and x == z:
                #print("match found")
                i = coordList.index((w,x))
                drones_list.append(drones_available[i])
                stations_deploying.append([w,x])
    #print("drones list is:", drones_list)
    for  w,x in coordList:
        for y,z in stations_deploying:
            if w == y and x == z:
                maxdrone = max(drones_list)
                i = coordList.index((w,x))
                if maxdrone == drones_available[i]:
                    #print([w,x])
                    return([w,x])
                else:
                    continue
#
#----------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------




#-----------------------------------------Event triggering function for conflict data and drone response---------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------
#
def eventTrigger(coordList):
    #Input start date and end date required for the simulation in unix time
    start_date = 1552576403
    end_date =   1552577224
    
    #Initialise counters for the simulation
    drones_used = 0
    crisis_out_out_of_range = 0
    drones_unavailable = 0
    #Initialise the drones at each station
    drones_available = droneStationInit(coordList)
    #Initialise the lists to store drone mission data
    response_timers = []
    stations_used = []
    anoth_list = []
    station_list = ['ALPHA','BRAVO','CHARLIE','DELTA','ECHO','FOXTROT','GOLF','HOTEL','INDIA','JULIET','KILO','LIMA','MIKE','NOVEMBER','OSCAR','PAPA','QUEBEC','ROMEO','SIERRA','TANGO','UNIFORM','VICTOR','WHISKEY','XRAY','YANKEE']
    
    #This iterates through the clock until the end date is found
    while start_date != end_date:
        start_date = start_date
        print(datetime.utcfromtimestamp(start_date).strftime('%Y-%m-%d %H:%M:%S'))
        #This iterates through the conflict file timestamps
        for index, row in conflict_file_loc.iterrows():
            #If the current time is equal to the timestamp, an incident has occurred
            if start_date == row["TIMESTAMP"]:
                print('')
                print('Incident Reported')
                print('')
                #Find the stations in range of the incident
                stations_available = stationsInRange(row["LATITUDE"], row["LONGITUDE"], coordList)
                #If there are stations in range
                if  stations_available is not None:
                    print('')
                    print('Drone Deploying')
                    #Chooses a station in range with the highest number of drones available
                    station_chosen = maximumStationDrones(stations_available, coordList, drones_available)

                    #This calculates the distance between the event and the chosen station
                    dist = (distanceCalc(station_chosen[1],station_chosen[0],row["LATITUDE"],row["LONGITUDE"]))*2
                    print(dist)
                    #This calculates the drone response time to the event (how long it needs to be out of action for)
                    response_time = int(droneFlightTime(dist))
                    print(response_time)
                
                    #This will update the drones available list by removing a drone from the chosen station
                    if station_chosen is not None:
                        a_list = []
                        a_list.append(station_chosen)
                        for w,x in coordList:
                            for y,z in a_list:
                                if w == y and x == z:
                                    #Finds the index of the station chosen in the corrdList
                                    i = coordList.index((w,x))
                                    #If the station has more than one drone
                                    if drones_available[i] != 1:
                                        drones_available[i] = drones_available[i]-1
                                        #Increase the counter of drones deployed
                                        drones_used = drones_used + 1
                                        #Add which station was used to the a list
                                        stations_used.append(station_chosen)
                                        #Add how long the drone will be out of service for to a list
                                        response_timers.append([start_date + response_time])
                                        print('Drone at station ',station_list[i],' dispatched, drones available :', drones_available[i])
                                    #If the station is currently depleted of drones
                                    else:
                                        drones_available[i] = drones_available[i]
                                        #Add to the list of incidents missed because no drones left in station
                                        drones_unavailable = drones_unavailable + 1
                                        print('Drone at station ',station_list[i],' empty :', drones_available[i])
                                        continue
                    else:
                        continue
                #If the station is out of range of the incident
                else:
                    crisis_out_out_of_range = crisis_out_out_of_range + 1
                    print("Out of Range")
                    continue
            else:
                
                continue
        

        #This is the code for returning the drones
        #We make a list with one value so that it can be iterated over by another list
        anoth_list.append([start_date])
        #This loop checks the response timers list and sees if any match the current time (I.E. a drone has finished its mission)
        for i in response_timers:
            for j in anoth_list:
                if i == j:
                    print('')
                    print('Drone Returned')
                    #Finds the index of the time in response timers
                    h = response_timers.index(i)
                    an_list = []
                    #Create a new list with a single value with the station the response timer indexed to
                    #Because the timers are added to a list at the same time as stations used, their indexes will be the same.
                    an_list.append(stations_used[h])
                    #Again we iterate over the list with the single valie to find a match in another list
                    for v,b in coordList:
                        for n,m in an_list:
                            if v == n and b == m:
                                #print("match found")
                                #Finds the index of the station used in the list of all stations
                                l = coordList.index((v,b))
                                #Appends the number of drones that are at that station. If you deploy multiple drones you must return multiple drones here.
                                drones_available[l] = drones_available[l] + 1
                                print('Drone at station ', station_list[l], ' returned, drones available: ', drones_available[l])
                            else:
                                #print("Drone Deploying")
                                continue
                    #Destroy the index used after their drone has returned
                    response_timers.pop(h)
                    stations_used.pop(h)
                else:
                    #print("no match found")
                    continue

        #Increase the start date by 1
        start_date = start_date + 1
    #
    #--------------------------------------------------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------------------------------------------------

    #This calculates the amount of aid delivered, change the integer that multiplies drones_used to change the payload capacity of the drones
    kilos_aid_delivered = drones_used*3.6

    
    #The output after the simulation is run
    print('')
    print('')
    print('')
    print('')
    print('Conflicts out of range of stations: ',crisis_out_out_of_range)
    print('')
    print('Conflicts missed because station ran out of drones: ',drones_unavailable)
    print('')
    print('Total kilograms of aid delivered to conflict: ', kilos_aid_delivered, 'kg')
    print('')
    print('Total drones deployed to conflict: ',drones_used)
    print('')
    #print('Percentage of conflicts reached: ',"{:.0%}".format(drones_used))




#-----------------------This allows you to interact with the map, the callback function callbacks on mouse double click-------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#
def callback(event):
    #Here you can set the amount of drone stations available to deploy
    if len(coordList) < 10:
        Coords=(round(event.x, 3),round(event.y, 3))
        coordList.append(Coords)
        print(coordList)
        source.data = dict(x=[i[0] for i in coordList], y=[i[1] for i in coordList])
    else:
        #If the number of stations deployed is max, display message
        print("Max number of stations reached")
#
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------




#-------------------------------This allows you to interact with the drone sim button to start the drone simulation-----------------------
#-----------------------------------------------------------------------------------------------------------------------------------------
#
def callbackSimRun(event):
    eventTrigger(coordList)
#
#-----------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------




#These plot the new stations and the button to run the simulation
p.on_event(DoubleTap, callback)
button_widget.on_event(ButtonClick, callbackSimRun)

#This sets the layout of the plot
layout=Column(p,widgetbox(button_widget))

#This updates the plot to the server
curdoc().add_root(layout)