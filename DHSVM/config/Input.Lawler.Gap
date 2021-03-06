##########################################################################################################
# DHSVM INPUT FILE FORMAT
##########################################################################################################
# The file is organized in sections [...], which contain key = entry pairs.  
# The file is free format, in that correct reading of the file is not dependent 
# on spaces and/or the order of the key-entry pairs within a section.  
# The keys are not case-sensitive, but the entries are, because filenames on a 
# UNIX platform are case-sensitive.  
# Comments are preceded by a '#', and run from the occurrence of '#' till the 
# end of the line. You can comment out an entire line (like in this 
# header), or you can place a comment after an entry.  
# It is important to place the key-entry pair in the correct section, since it 
# will not be found if it is in another section.  
# The easiest way to make the input file is to fill out this default template.  
# Since DHSVM will only use the keys that it requires you do not have to worry 
# about empty entries for keys that are not needed. For example, if you are 
# running the model in point mode, you do not have to fill out the routing 
# section.  If you have already filled it out you can leave it, since DHSVM will 
# not use the information.  This allows easy switching between point and basin 
# mode.  
# For more information about the specific entries see the DHSVM web page
# the soil parameters are copied from Lan's lowland INPUT File (with extremely 
# low lateral conductivity $ changed the station from 1 to 2
##########################################################################################################
# OPTIONS SECTION 
##########################################################################################################

[OPTIONS]                                 # Model Options
Format               = BIN                # BIN or NETCDF
Extent               = BASIN              # POINT or BASIN
Gradient             = TOPOGRAPHY   	    # TOPOGRAPHY or WATERTABLE
Flow Routing         = NETWORK            # UNIT_HYDROGRAPH or NETWORK (StrEnv[flow_routing])
Sensible Heat Flux   = FALSE              # TRUE or FALSE
Overland Routing     = CONVENTIONAL	      # CONVENTIONAL or KINEMATIC (StrEnv[routing])
Infiltration	       = STATIC		          # Static or Dynamic
Interpolation        = INVDIST            # NEAREST or INVDIST or VARCRESS
MM5                  = FALSE		          # TRUE or FALSE
QPF                  = FALSE              # TRUE or FALSE
PRISM                = FALSE
PRISM data path      = 
PRISM data extension = 
Gridded Met data     = FALSE              # TRUE or FALSE
Canopy radiation attenuation mode = FIXED # FIXED or VARIABLE
Shading              = FALSE    	        # TRUE or FALSE
Shading data path    = 	
Shading data extension =  
Skyview data path    = 
Snotel               = FALSE              # TRUE or FALSE
Outside              = FALSE	            # TRUE or FALSE
Rhoverride           = FALSE              # TRUE or FALSE
Precipitation Source = STATION            # STATION or RADAR
Wind Source          = STATION            # STATION or MODEL
Temperature lapse rate   = CONSTANT       # CONSTANT or VARIABLE
Precipitation lapse rate = CONSTANT       # CONSTANT, MAP, or VARIABLE
Cressman radius      =                    # in model pixels
Cressman stations    =                    # the maximum stations for interpolation
Stream Temperature   = TRUE               # TRUE or FALSE
Riparian Shading       = TRUE             # used for stream tempeature only
Variable Light Transmittance = TRUE       # TRUE if light transmittance changes with solar altitude
Canopy Gapping = TRUE
Snow Sliding = FALSE                      # this function is not available - needs more tests
Precipitation Separation = FALSE          # TRUE if snow and rain are separately provided in meterological input data (e.g. WRF)

##########################################################################################################
# MODEL AREA SECTION
##########################################################################################################

[AREA]                                    # Model area
Coordinate System    = UTM                # UTM or USER_DEFINED
Extreme North        = 5189026.4          # Coordinate for northern edge of grid
Extreme West         = 521090.8           # Coordinate for western edge of grid
Center Latitude      = 46.854444          # Central parallel of basin 
Center Longitude     = -116.723333        # Central meridian of basin 
Time Zone Meridian   = -105.0             # Time zone meridian for area 
Number of Rows       = 3                  # Number of rows
Number of Columns    = 3                  # Number of columns
Grid spacing         = 90                 # Grid resolution in m 
Point North          =                    # North coordinate for point model if Extent = POINT 
Point East           =                    # East coordinate for point model if Extent = POINT 

##########################################################################################################
# TIME SECTION
##########################################################################################################

[TIME]                                    # Model period
Time Step            =  0.5               # Model time step (hours)
Model Start          =  10/1/2007-00:00   # Model start time (MM/DD/YYYY-HH)
Model End            =  5/15/2008-00:00   # Model end time (MM/DD/YYYY-HH) 
################################################################################
# CONSTANTS SECTION
################################################################################
[CONSTANTS]                               # Model constants
Ground Roughness     = 0.02               # Roughness of soil surface (m)
Snow Roughness       = 0.01               # Roughness of snow surface (m)
Snow Water Capacity  = 0.01               # Snow liquid water holding capacity(fraction)
Reference Height     = 70.0               # Reference height (m)
Rain LAI Multiplier  = 0.0001             # LAI Multiplier for rain interception
Snow LAI Multiplier  = 0.0002             # LAI Mulitplier for snow interception
Min Intercepted Snow = 0.001              # Intercepted snow that can only be melted (m)
Outside Basin Value  = 0                  # Value in mask that indicates outside the basin
Temperature Lapse Rate   = -0.0065        # Temperature lapse rate (C/m)
Precipitation Lapse Rate =  0.000    	    # Precipitation lapse rate (m/m). 
Max Surface Snow Layer Depth = 0.125      # maximum depth of the surface snow layer in water equivalent (m)
Gap Wind Adj Factor = 0.9                 # Adjustment for wind received by gap (0-1)
Snowslide Parameter1 = 700                # Empirical parameters used to calculate snow holding depth
Snowslide Parameter2 = 0.12               # Empirical parameters used to calculate snow holding depth

########################################################
# The following fields in this section can take either a single value (i.e. basin-wide constant parameter value)
# or the path to a spatial parameter binary file in the same format as DEM, soil and veg  

Precipitation Multiplier Map = ../input/multiplier.bin # no correction be made if leave blank, or a constant applied basin-wide, or a spatial input
Rain Threshold       = 0.0                             # Minimum temperature at which rain occurs (C)
Snow Threshold       = ../input/Ts.bin                 # Maximum temperature at which snow occurs (C)
Fresh Snow Albedo  = ../input/amax.bin                 # Fresh snow albedo (0.85 is the default value used in older versions)
Albedo Accumulation Lambda = 0.9 	                     # Lambda for accumulation period
Albedo Melting Lambda = 0.6 		                       # Lambda for melting period
Albedo Accumulation Min = 0.7 		                     # Min albedo value during accumulation period; set to -9999 if min value is not needed.
Albedo Melting Min = 0.55 		                         # Min albedo value during melting period; set to -9999 if min value is not needed.
                                                
##########################################################################################################
# TERRAIN INFORMATION SECTION
##########################################################################################################

[TERRAIN]                                 # Terrain information
DEM File    	  = ../input/dem.bin
Basin Mask File   = ../input/mask.bin
##########################################################################################################

# ROUTING SECTION
##########################################################################################################

[ROUTING]                                 # Routing information. This section is 
                                          # only relevant if the Extent = BASIN
############################### STREAM NETWORK ###########################################################

# The following three fields are only used if Flow Routing = NETWORK
Stream Map File      = ../input/stream.map.dat
Stream Network File  = ../input/stream.network.dat
Stream Class File    = ../input/adjust.classfile

# enter the path to riparian parameters if Canopy Gapping = TRUE
Riparian Veg File    = ../input/rveg.baseline.network

################################### ROAD NETWORK #########################################################

# The following three fields are only used if Flow Routing = NETWORK and there
# is a road network

#Road Map File        =                    # path for road map file
#Road Network File    =                    # path for road network file
#Road Class File      =                    # path for road network file

#################################### UNIT HYDROGRAPH #####################################################

# The following two fields are only used if Flow Routing = UNIT_HYDROGRAPH

Travel Time File     =                    # path for travel time file
Unit Hydrograph File =                    # path for unit hydrograph file

##########################################################################################################

# METEOROLOGY SECTION 
##########################################################################################################

[METEOROLOGY]     
# Meteorological station if [OPTION] Gridded Met data = FALSE
Number of Stations = 1                   # Number of meteorological stations
Station Name     1 = Lawler              # Name for station 1
North Coordinate 1 = 5188936             # North coordinate of station 1
East Coordinate  1 = 521180              # East coordinate of station 1
Elevation        1 = 1000.00000          # Elevation of station 1 in m
Station File     1 = ../met/subhourly_met

# Gridded Met if [OPTION] Gridded Met data = TRUE
Extreme North Lat =                  
Extreme South Lat = 
Extreme East Lon = 
Extreme West Lon = 
Number of Grids =         # estimate of total grids need for the basin for memory 
                          # allocation (must >= actual grids used)  
GRID_DECIMAL =            # Number of digits after decimal point in forcing file names
Met File Path = 
File Prefix = 

########################## MM5 ##########################

# The following block only needs to be filled out if MM5 = TRUE.  In that case
# This is the ONLY block that needs to be filled out

MM5 Start              =                  # Start of MM5 file (MM/DD/YYYY-HH),
MM5 Rows               =
MM5 Cols               =
MM5 Extreme North      =
MM5 Extreme West       = 
MM5 DY                 =

# MM5 met files
MM5 Temperature File   = 
MM5 Humidity File      = 
MM5 Wind Speed File    = 
MM5 Shortwave File     = 
MM5 Longwave File      = 
MM5 Pressure File      = 
MM5 Precipitation File = 
MM5 Terrain File       =
MM5 Temp Lapse File    =

# For each soil layer make a key-entry pair as below (n = 1, ..,
# Number of Soil Layers)

MM5 Soil Temperature File 0 = 
MM5 Soil Temperature File 1 = 
MM5 Soil Temperature File 2 = 

######################################### RADAR ##########################################################
# The following block only needs to be filled out if Precipitation Source = 
# RADAR.  
Radar Start            =
Radar File             =
Radar Extreme North    =
Radar Extreme West     =
Radar Number of Rows   =
Radar Number of Columns = 
Radar Grid Spacing     =

########################################## Wind ##########################################################
# The following block only needs to be filled out if Wind Source = MODEL
Number of Wind Maps    =
Wind File Basename     =
Wind Map Met Stations  =

########################################## Precipitation lapse rate ######################################
# The following block only needs to be filled out if Precipitation lapse rate 
# = MAP
Precipitation lapse rate =

##########################################################################################################

# SOILS INFORMATION SECTION
##########################################################################################################
[SOILS]                                   # Soil information
Soil Map File   = ../input/soil.bin
Soil Depth File = ../input/soild.bin
Number of Soil Types = 1           

######################################### SOIL 1 ###############################################
Soil Description       1 = Sandy Loam        # Description of soil type
Lateral Conductivity   1 = 5e-5              # Lateral saturated hydraulic conductivity (m/s)
Exponential Decrease   1 = 0.03              # Exponent for decrease in conductivity with depth
Depth Threshold        1 = 1.5
Maximum Infiltration   1 = 1e-5              # Maximum infiltration rate (m/s)
Capillary Drive        1 = 0.05
Surface Albedo         1 = 0.1               # Soil surface albedo
Number of Soil Layers  1 = 3                 # Number of soil layers
Porosity               1 = 0.42 0.40 0.35    # Soil porosity 
Pore Size Distribution 1 = 0.38 0.38 0.38    # Pore size distribution
Bubbling Pressure      1 = .146 .146 .146    # Bubbling Pressure (m)	
Field Capacity         1 = .207 .207 .207    # Soil moisture content at field capacity
Wilting Point          1 = .095 .095 .095    # Soil moisture content at wilting point
Bulk Density           1 = 1540 1590 1740    # Soil bulk density (kg/m3)
Vertical Conductivity  1 = 1e-5 1e-5 1e-5    # Vertical saturated hydraulic conductivity (m/s)
Thermal Conductivity   1 = 7.70 7.49 7.46    # Effective solids thermal conductivity (W/(m*K))
Thermal Capacity       1 = 1.4e6 1.4e6 1.4e6 # Dry soil thermal capacity (J/(m3*K))
Residual Water Content 1 = 0. 0. 0.

###########################################################################################################	
# VEGETATION INFORMATION SECTION
############################################################################################################
[VEGETATION]

Vegetation Map File = ../input/Veg.bin
Canopy Gap Map File = ../input/GapDiam.bin    # Gap diameter spatial input

Number of Vegetation Types = 1              # Number of different vegetation types

############################ Vegetation  1  ########################
Vegetation Description   1 = Broadleaf Evergreen Forest 
Impervious Fraction      1 = 0.0
Detention Fraction       1 = 0
Detention Decay          1 = 0
Overstory Present        1 = TRUE      
Understory Present       1 = FALSE             
Fractional Coverage      1 = 1.0         
Hemi Fract Coverage      1 = 
Clumping Factor          1 =
Leaf Angle A             1 =
Leaf Angle B             1 =
Scattering Parameter     1 =
Trunk Space              1 = 0.5         
Aerodynamic Attenuation  1 = 2.5     
Radiation Attenuation    1 = 0.2               
Diffuse Radiation Attenuation  1 = 0.215   
Max Snow Int Capacity    1 = 0.03 
Snow Interception Eff    1 = 0.6              
Mass Release Drip Ratio  1 = 0.4             
Height                   1 = 25
Overstory Monthly LAI    1 = 5.0 5.0 5.0 5.0 8.0 10.0 10.0 8.0 5.0 5.0 5.0 5.0 
Understory Monthly LAI   1 = 
Maximum Resistance       1 = 4000  
Minimum Resistance       1 = 460    
Moisture Threshold       1 = 0.33      
Vapor Pressure Deficit   1 = 4000      
Rpc                      1 = .108 
Overstory Monthly Alb    1 = 0.14 0.14 0.14 0.13 0.13 0.12 0.11 0.11 0.12 0.13 0.14 0.14
Understory Monthly Alb   1 = 
Number of Root Zones     1 = 3 		   
Root Zone Depths         1 = 0.10 0.15 0.2    
Overstory Root Fraction  1 = 0.20 0.40 0.40      
Understory Root Fraction 1 =          
Monthly Light Extinction 1 = 0.065 0.065 0.065 0.065 0.065 0.065 0.065 0.065 0.065 0.065 0.065 0.065
Canopy View Adj Factor   1 = 1.2

################################################################################
# MODEL OUTPUT SECTION

################################################################################
[OUTPUT]                                  # Information what to output when

Output Directory = ../output/
Initial State Directory   = ../input/modelstate/

################ PIXEL DUMPS ###################################################

Number of Output Pixels    = 0

# For each pixel make a key-entry pair as indicated below, varying the 
# number for the output pixel  (1, .. , Number of Output Pixel)



################ MODEL STATE ###################################################                
Number of Model States     =  0           # Number of model states to dump

# For each model state make a key-entry pair as indicated below, varying the 
# number for the model state dump (1, .. , Number of Model States)

State Date               1 = 1/1/1970-00


################ MODEL MAPS ####################################################

Number of Map Variables    = 0            # Number of different variables for
                                          # which you want to output maps

######################################### MODEL IMAGES #################################################################

Number of Image Variables  = 0            # Number of variables for which you 
                                          # would like to output images

# For each of the variables make a block like the one that follows, varying
# the number of the variable (n = 1, .. , Number of Image Variables)

Image Variable        1 = 801              # ID of the variable to output
Image Layer           1 = 1                # If the variable exists for a number
										   # of layers, specify the layers here with the top layer = 1


Image Start           1 = 12/10/1999-03    # First timestep for which to output an image										                

Image End             1 = 12/10/1999-06    # Last timestep for which to output an image

Image Interval        1 =                  # Time interval between images (hours)
Image Upper Limit     1 =                  # All values in the output equal to or
                                           # greater than this limit will be set 
                                           # to 255
Image Lower Limit     1 =                  # All values in the output equal to or
                                           # smaller than this limit will be set 
                                           # to 0
######################################### GRAPHIC IMAGES ################################################################


Number of Graphics      =  0               # Number of variables for which you 
                                           # would like to output images
Graphics ID           1 = 15               # ID of the variable to output
Graphics ID           1 = 22
Graphics ID           2 = 23
Graphics ID           3 = 24
Graphics ID           4 = 25
Graphics ID           5 = 43
Graphics ID           6 = 44
Graphics ID           7 = 8
Graphics ID           8 = 2
Graphics ID           9 = 50
Graphics ID           10 = 1

# 1	SWE (mm)
# 2	Water Table Depth (mm)
# 3	Digital Elevation Model (m)
# 4	Vegetation Class (index #)
# 5	Soil Class (index #)
# 6	Soil Depth (mm)
# 7	Precipitation at current time step (mm/time step)
# 8	Incoming Shortwave (Beam and Diffuse) (W/sqm)
# 9	Intercepted Snow (mm)
# 10	Snow Surface Temp (C)
# 11	Cold Content of snow entire snow pack (kJ)
# 12	Snow Melt (as Outflow minus Precip, can be negative) (mm/time step)
# 13	Snow Pack Outflow (mm/time step)
# 14	Saturated Subsurface Flow (mm/time step)
# 15	Overland Flow(mm)
# 16	Total Evapotranspiration (soil + all veg layers)
# 17	Ground Snow pack vapor flux (mm)
# 18	Intercepted snow pack vapor flux (mm)
# 19    Soil Moisture (Surface Layer) % of saturation (i.e. porosity)
# 20    Soil Moisture (2nd Layer) % of saturation (i.e. porosity)
# 21    Soil Moisture (3rd Layer) % of saturation (i.e. porosity)
# 22    Accumulated Precip (mm)
# 23    air temperature (C)
# 24    wind speed (m/s)
# 25    relative humidity
# 26    Prism Precip Field (mm)
# 31    Overstory Transpiration (mm)
# 32    Understory Transpirtation (mm)
# 33    Soil Evaporation (mm)
# 34    Overstory Evaporation (mm)
# 35    Understory Evaportation (mm)
# 41    Sky View Factor (%)
# 42    Shade Map (%)
# 43    Direct Beam Shortwave Rad (W/sqm)
# 44    Diffuse Beam Shortwave Rad (W/sqm)
# 45    Aspect (degrees)
# 46    Slope (percent)
# 50    Channel Subsurface Interception (mm)
# 51    Road Subsurface Interception (mm)
# WARNING Use soil mositure layers with caution, to minimize calculations during redraw
# DHSVM does not check to make sure that the assigned soil layer exists
#########################################################################################################


# END OF INPUT FILE
#########################################################################################################

[End]                                     # This is probably not needed, but 
                                          # just in case (to close the previous
                                          # section)
#########################################################################################################

