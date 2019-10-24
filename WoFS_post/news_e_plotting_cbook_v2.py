
#######################################################################################
#news_e_plotting_cbook.py - written 2016 by Pat Skinner and Jessie Choate
#
#This code is a collection of classes and subroutines used to create plots of NEWS-e 
#summary output for the real-time website:  
#http://www.nssl.noaa.gov/projects/wof/news-e/images.php? 
#
######################################################################################
#Classes available are: 
#
#cb_colors     -> Collection of RGB values for custom colormaps.  Colors taken from 
#                 the colorbrewer website:  http://colorbrewer2.org/
#web_plot      -> Object consisting of information needed (i.e. labels, contour levels)
#                 to create most of the NEWS-e website plots
#v_plot        -> Object consisting of information needed to create paintball plots
#ob_plot       -> Object consisting of information needed to create object matching plots
#
#Plotting Subroutines available are: 
#
#mymap         -> Creates basemap instance (Lambert Conformal projection) with 
#                 geographic boundaries for plotting
#create_fig    -> Creates a standardized figure template for NEWS-e website plots
#env_plot      -> Plots a NEWS-e environment or swath product, saves .png, then 
#                 removed plots and returns empty basemap (used for majority of plotting)
#paintqc_plot  -> Creates paintball plot of forecast and verification objects
#ob_ver_plot   -> Creates probability of object match vs. false alarm plot
#ob_cent_plot  -> Creates scatter plot of object centroid positions (no longer used for website)
#dot_plot      -> Same as env_plot, except includes scatter plot of OK Mesonet station error
#mrms_plot     -> Creates plots of MRMS DZ and Az. Shear over the NEWS-e domain
#                 (not used for website)
#plot_warn     -> Plots shapefiles of NWS SVR/TOR/FF warnings 
#plot_lsr      -> Adds scatter plots of hail/wind/tornado local storm reports
#remove_warn   -> Removes plots of NWS warnings
#rmove_lsr     -> Removes plots of LSRs
#pmm_rectangle -> Plots rectangle over domain where probability-matched mean is valid
#                 (no longer used) 
#rem_pmm_rect  -> Removes pmm_rectangle plot (no longer used) 
#
#######################################################################################
#Future needs: 
#
#
#######################################################################################

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib
import pylab as P
from scipy import signal
from scipy import *
from scipy import ndimage
import ctables 			#Separate file with Carbone42 (soloii) colortable amongst others

#######################################################################################

class cb_colors(object):

   #Defines RGB values for select colors from colorbrewer 2.0 (http://colorbrewer2.org/)

   #Input:

   #None

   #Returns:

   #None

#######################

#These colors are single color colortable RGB values from colorbrewer 2.0 (http://colorbrewer2.org/)

   orange1 = (255 / 255., 245 / 255., 235 / 255.)
   orange2 = (254 / 255., 230 / 255., 206 / 255.)
   orange3 = (253 / 255., 208 / 255., 162 / 255.)
   orange4 = (253 / 255., 174 / 255., 107 / 255.)
   orange5 = (253 / 255., 141 / 255., 60 / 255.)
   orange6 = (241 / 255., 105 / 255., 19 / 255.)
   orange7 = (217 / 255., 72 / 255., 1 / 255.)
   orange8 = (166 / 255., 54 / 255., 3 / 255.)
   orange9 = (127 / 255., 39 / 255., 4 / 255.)

   blue1 = (247/255., 251/255., 255/255.)
   blue2 = (222/255., 235/255., 247/255.)
   blue3 = (198/255., 219/255., 239/255.)
   blue4 = (158/255., 202/255., 225/255.)
   blue5 = (107/255., 174/255., 214/255.)
   blue6 = (66/255., 146/255., 198/255.)
   blue7 = (33/255., 113/255., 181/255.)
   blue8 = (8/255., 81/255., 156/255.)
   blue9 = (8/255., 48/255., 107/255.)

   purple1 = (252/255., 251/255., 253/255.)
   purple2 = (239/255., 237/255., 245/255.)
   purple3 = (218/255., 218/255., 235/255.)
   purple4 = (188/255., 189/255., 220/255.)
   purple5 = (158/255., 154/255., 200/255.)
   purple6 = (128/255., 125/255., 186/255.)
   purple7 = (106/255., 81/255., 163/255.)
   purple8 = (84/255., 39/255., 143/255.)
   purple9 = (63/255., 0/255., 125/255.)

   green1 = (247/255., 252/255., 245/255.)
   green2 = (229/255., 245/255., 224/255.)
   green3 = (199/255., 233/255., 192/255.)
   green4 = (161/255., 217/255., 155/255.)
   green5 = (116/255., 196/255., 118/255.)
   green6 = (65/255., 171/255., 93/255.)
   green7 = (35/255., 139/255., 69/255.)
   green8 = (0/255., 109/255., 44/255.)
   green9 = (0/255., 68/255., 27/255.)

   gray1 = (255/255., 255/255., 255/255.)
   gray2 = (240/255., 240/255., 240/255.)
   gray3 = (217/255., 217/255., 217/255.)
   gray4 = (189/255., 189/255., 189/255.)
   gray5 = (150/255., 150/255., 150/255.)
   gray6 = (115/255., 115/255., 115/255.)
   gray7 = (82/255., 82/255., 82/255.)
   gray8 = (37/255., 37/255., 37/255.)
   gray9 = (0/255., 0/255., 0/255.)

   red1 = (255/255., 245/255., 240/255.)
   red2 = (254/255., 224/255., 210/255.)
   red3 = (252/255., 187/255., 161/255.)
   red4 = (252/255., 146/255., 114/255.)
   red5 = (251/255., 106/255., 74/255.)
   red6 = (239/255., 59/255., 44/255.)
   red7 = (203/255., 24/255., 29/255.)
   red8 = (165/255., 15/255., 21/255.)
   red9 = (103/255., 0/255., 13/255.)

### Qualitative colors (pastels): 

   q1 = (141/255., 255/255., 199/255.)  #aqua
   q2 = (255/255., 255/255., 179/255.)  #pale yellow
   q3 = (190/255., 186/255., 218/255.)  #lavender
   q4 = (251/255., 128/255., 114/255.)  #pink/orange
   q5 = (128/255., 177/255., 211/255.)  #light blue
   q6 = (253/255., 180/255., 98/255.)   #light orange
   q7 = (179/255., 222/255., 105/255.)  #lime
   q8 = (252/255., 205/255., 229/255.)  #pink
   q9 = (217/255., 217/255., 217/255.)  #light gray
   q10 = (188/255., 128/255., 189/255.) #purple
   q11 = (204/255., 235/255., 197/255.) #pale green
   q12 = (255/255., 237/255., 111/255.) #yellow

### Qualitative colors (bright):

   b1 = (228/255., 26/255., 28/255.)   #red
   b2 = (55/255., 126/255., 184/255.)  #blue
   b3 = (77/255., 175/255., 74/255.)   #green
   b4 = (152/255., 78/255., 163/255.)  #purple
   b5 = (255/255., 127/255., 0/255.)   #orange
   b6 = (255/255., 255/255., 51/255.)  #yellow
   b7 = (166/255., 86/255., 40/255.)   #brown
   b8 = (247/255., 129/255., 191/255.) #pink

### NWS Reflectivity Colors (courtesy MetPy library):

   c5 =  (0.0,                 0.9254901960784314, 0.9254901960784314)
   c10 = (0.00392156862745098, 0.6274509803921569, 0.9647058823529412)
   c15 = (0.0,                 0.0,                0.9647058823529412)
   c20 = (0.0,                 1.0,                0.0)
   c25 = (0.0,                 0.7843137254901961, 0.0)
   c30 = (0.0,                 0.5647058823529412, 0.0)
   c35 = (1.0,                 1.0,                0.0)
   c40 = (0.9058823529411765,  0.7529411764705882, 0.0)
   c45 = (1.0,                 0.5647058823529412, 0.0)
   c50 = (1.0,                 0.0,                0.0)
   c55 = (0.8392156862745098,  0.0,                0.0)
   c60 = (0.7529411764705882,  0.0,                0.0)
   c65 = (1.0,                 0.0,                1.0)
   c70 = (0.6,                 0.3333333333333333, 0.788235294117647)
   c75 = (0.0,                 0.0,                0.0)

### Custom Colormaps: 

   cin_cmap = matplotlib.colors.ListedColormap([purple7, purple6, purple5, purple4, blue4, blue3, blue2, blue1])

   wz_cmap = matplotlib.colors.ListedColormap([blue2, blue3, blue4, red2, red3, red4, red5, red6, red7])

   dz_cmap_2 = matplotlib.colors.ListedColormap([blue5, blue3, green3, green5, green7, orange3, orange5, orange7, red7, red8, purple8, purple6])

   dz_cmap = matplotlib.colors.ListedColormap([green5, green4, green3, orange2, orange4, orange6, red6, red4, purple3, purple5])

   nws_dz_cmap = matplotlib.colors.ListedColormap([c20, c25, c30, c35, c40, c45, c50, c55, c60, c65, c70])
#   nws_dz_cmap = matplotlib.colors.ListedColormap([c5, c10, c15, c20, c25, c30, c35, c40, c45, c50, c55, c60, c65, c70])

#   wind_cmap = matplotlib.colors.ListedColormap([gray2, gray3, gray3, blue2, blue3, blue4, blue5, blue6, blue7, green3, green4, green5, orange3, orange4, orange5, red6, red7, red8, purple5, purple6, purple7 ])
#   wind_trop_cmap = matplotlib.colors.ListedColormap([gray1, gray2, gray3, blue1, blue2, blue3, green2, green3, green4, orange3, orange4, orange5, red4, red5, red6])

   wind_cmap = matplotlib.colors.ListedColormap([gray1, gray2, gray3, orange2, orange3, orange4, orange5, orange6, red7, red8])

   wz_cmap_extend = matplotlib.colors.ListedColormap([blue2, blue3, blue4, red2, red3, red4, red5, red6, red7, purple7, purple6, purple5])

   cape_cmap = matplotlib.colors.ListedColormap([blue2, blue3, blue4, orange2, orange3, orange4, orange5, red4, red5, red6, red7, purple7, purple6, purple5])

   td_cmap_ncar = matplotlib.colors.ListedColormap(['#ad598a', '#c589ac','#dcb8cd','#e7cfd1','#d0a0a4','#ad5960', '#8b131d', '#8b4513','#ad7c59', '#c5a289','#dcc7b8','#eeeeee', '#dddddd', '#bbbbbb', '#e1e1d3', '#e1d5b1','#ccb77a','#ffffe5','#f7fcb9', '#addd8e', '#41ab5d', '#006837', '#004529', '#195257', '#4c787c'])

   temp_cmap_ugly = matplotlib.colors.ListedColormap([blue6, blue4, blue2, green6, green4, green2, orange2, orange4, red5, red7, purple7])
#   temp_cmap_ugly = matplotlib.colors.ListedColormap([purple5, purple3, blue3, blue5, green5, green3, orange2, orange4, red5, red7, purple7])

   temp_cmap = matplotlib.colors.ListedColormap([purple4, purple5, purple6, purple7, blue8, blue7, blue6, blue5, blue4, blue3, green7, green6, green5, green4, green3, green2, orange2, orange3, orange4, orange5, red5, red6, red7, red8, purple6, purple5, purple4, purple3])

   blues_cmap = matplotlib.colors.ListedColormap([blue3, blue4, blue5, blue6, blue7])

   oranges_cmap = matplotlib.colors.ListedColormap([orange3, orange4, orange5, orange6, orange7])

   td_cmap_ugly = matplotlib.colors.ListedColormap([orange3, gray4, gray3, gray1, green3, green5, green7, blue3, blue5, blue7, purple3])
#   td_cmap_ugly = matplotlib.colors.ListedColormap([orange4, orange2, green3, green5, green7, blue3, blue5, blue7, purple3])

   td_cmap = matplotlib.colors.ListedColormap([gray6, gray5, gray4, gray3, gray2, gray1, green1, green2, green3, green4, green5, green6, blue3, blue4, blue5, purple3])

#   td_cmap = matplotlib.colors.ListedColormap([orange6, orange5, orange4, orange3, gray6, gray5, gray4, gray3, gray2, gray1, green1, green2, green3, green4, green5, green6, blue3, blue4, blue5, blue6, purple3, purple4, purple5 ])

   uv_cmap = matplotlib.colors.ListedColormap([purple5, purple4, purple3, purple2, purple1, orange1, orange2, orange3, orange4, orange5])

   diff_cmap = matplotlib.colors.ListedColormap([blue7, blue6, blue5, blue4, blue3, blue2, blue1, red1, red2, red3, red4, red5, red6, red7])

   mslp_cmap = matplotlib.colors.ListedColormap([purple7, purple6, purple5, red7, red6, red5,orange7, orange6, orange5, green8, green7, green6, blue6, blue5, blue4, gray3, gray2, gray1])

   paintball_colors = matplotlib.colors.ListedColormap([q1, b8, q3, q4, q5, q6, q7, q8, b6, q10, q11, b3, b2, purple5, red5, green5, blue5, orange5])
   paintball_colors_list = [q1, b8, q3, q4, q5, q6, q7, q8, b6, q10, q11, b3, b2, purple5, red5, green5, blue5, orange5]

   tl_paintball_colors = matplotlib.colors.ListedColormap([gray6, orange6, blue6, red6, green6, purple6, gray4, orange4, blue4, red4, green4, purple4, gray6, orange6, blue6, red6]) 
   tl_paintball_colors_list = [gray6, orange6, blue6, red6, green6, purple6, gray4, orange4, blue4, red4, green4, purple4, gray6, orange6, blue6, red6] 

#   tl_paintball_colors_list = [blue5, orange5, green5, purple5, red5, q1, b8, q3, q4, q5, q6, q7, q8, b6, b3, b2]

   mslp_paint_colors = matplotlib.colors.ListedColormap([purple8,purple6, purple4, red8, red6, red4, orange8, orange6, orange4, green8, green6, green4, blue8, blue6, blue4, gray3, gray2, gray1])

   all_blues_cmap =  matplotlib.colors.ListedColormap([blue1, blue2, blue3, blue4, blue5, blue6, blue7, blue8, blue9])
   all_greens_cmap =  matplotlib.colors.ListedColormap([green1, green2, green3, green4, green5, green6, green7, green8, green9])
   all_reds_cmap =  matplotlib.colors.ListedColormap([red1, red2, red3, red4, red5, red6, red7, red8, red9])

#   mfc_cmap = matplotlib.colors.ListedColormap([blue7, blue6, blue4, blue2, gray1, red2, red4, red6, red7])
   mfc_cmap = matplotlib.colors.ListedColormap([gray1, orange2, orange3, orange4, red5, red6, red7])

   corf_cmap = matplotlib.colors.ListedColormap([purple7, purple5, purple3, purple1, gray1, gray1, orange1, orange3, orange5, orange7])

   ul_dvg_cmap = matplotlib.colors.ListedColormap([green8, green6, green4, green2, gray1, gray1, purple2, purple4, purple6, purple8])

   rain_cmap = matplotlib.colors.ListedColormap([green2, green3, green5, blue4, blue5, blue6, purple6, purple5, purple4, red4, red5, red6])

   cp_cmap   = matplotlib.colors.ListedColormap([gray8,gray6,purple6,purple5,purple4, blue7, blue6, blue5, blue4,green6,green5,green4, green3, orange6, orange5,orange4, red5, red6,red8])

   cwp_cmap = matplotlib.colors.ListedColormap([gray2, gray3, gray4, gray5, gray6,  blue6, blue5, blue4, blue3,green6, green5, green4,green3, green2, orange2, orange3, orange4,orange5, red5, red6, red7, red8, purple6, purple7, purple8])
   mslp_cmap = matplotlib.colors.ListedColormap([purple8,purple6, purple4, red8, red6, red4, orange8, orange6, orange4, green8, green6, green4, blue8, blue6, blue4, gray3, gray2, gray1])

   pw_cmap = matplotlib.colors.ListedColormap([orange5, orange4, orange3, orange2,green1, green2, green3, green4,green5, green6, green7, green8, green9, blue4, blue5, blue6, blue7, blue8, purple4, purple5, purple6,purple7,red4,red5,red6,red7,red8])

   omega_cmap = matplotlib.colors.ListedColormap([purple7, purple5, purple3, purple2, gray1, green2, green3, green5, green7])
   
#######################################################################################

class web_plot:
   
   #Defines a plotting object used for the majority of the NEWS-e website plots: 
   #http://www.nssl.noaa.gov/projects/wof/news-e/images.php?

   #Input:

   #None

   #Returns:

   #None

#######################

   def __init__(self, name, var1_title, var2_title, var2_tcolor, var1_levels, var2_levels, var1_cont_levels, var1_cont_colors, var2_colors, over_color, under_color, cmap, extend, alpha, neighborhood):
      self.name			= name			#string name for plot (e.g. wz0to2rot90)
      self.var1_title		= var1_title		#string title for var1 plot
      self.var2_title		= var2_title		#string title for var2 plot
      self.var2_tcolor		= var2_tcolor		#title color for var2 plot
      self.var1_levels		= var1_levels		#contour levels for var1 plot
      self.var2_levels		= var2_levels		#contour levels for var2 plot
      self.var1_cont_levels	= var1_cont_levels	#levels for second contour plot with var1 (e.g. 60 deg isodrosotherm)
      self.var1_cont_colors	= var1_cont_colors	#colors for second contour plot with var1
      self.var2_colors		= var2_colors		#colors for var2 plot
      self.over_color		= over_color		#color for values > largest var1_level
      self.under_color		= under_color		#color for values < smallest var1_level
      self.cmap			= cmap			#colormap for var1 plot
      self.extend		= extend		#string for extend option in var1 plot (e.g. 'neither', 'max', ...)
      self.alpha		= alpha			#float value for transparency of var1 plot
      self.neighborhood		= neighborhood		#Neighborhood width for probability matched mean

#######################################################################################

class v_plot:

   #Defines a plotting object used for plotting paintball verification plots: 
   #http://www.nssl.noaa.gov/projects/wof/news-e/images.php?

   #Input:

   #None

   #Returns:

   #None

#######################

   def __init__(self, name, var1_title, var2_title, var1_tcolor, var2_tcolor, var1_levels, var2_levels, var1_colors, var2_colors, var1_units, var2_units, extend, alpha, neighborhood):
      self.name                 = name                  #string name for plot (e.g. wz0to2rot90)
      self.var1_title           = var1_title            #string title for var1 plot
      self.var2_title           = var2_title            #string title for var2 plot
      self.var1_tcolor          = var1_tcolor           #title color for var1 plot
      self.var2_tcolor          = var2_tcolor           #title color for var2 plot
      self.var1_levels          = var1_levels           #contour levels for var1 plot
      self.var2_levels          = var2_levels           #contour levels for var2 plot
      self.var1_colors          = var1_colors           #colors for var1 plot
      self.var2_colors          = var2_colors           #colors for var2 plot
      self.var1_units           = var1_units            #string of var1 unit label
      self.var2_units           = var2_units            #string of var2 unit label
      self.extend               = extend                #string for extend option in var1 plot (e.g. 'neither', 'max', ...)
      self.alpha                = alpha                 #float value for transparency of var1 plot
      self.neighborhood		= neighborhood		#Neighborhood width for probability matched mean

#######################################################################################

class ob_plot:

   #Defines a plotting object used for object matching verification plots:   
   #http://www.nssl.noaa.gov/projects/wof/news-e/images.php?

   #Input:

   #None

   #Returns:

   #None

#######################

   def __init__(self, name, var1_title, var2_title, var1_tcolor, var2_tcolor, var1_levels, var2_levels, var1_cmap, var2_cmap, var3_level, var3_color, extend, alpha, neighborhood): 
      self.name                 = name                  #string name for plot (e.g. wz0to2rot90)
      self.var1_title           = var1_title            #string title for var1 plot
      self.var2_title           = var2_title            #string title for var2 plot
      self.var1_tcolor          = var1_tcolor           #title color for var1 plot
      self.var2_tcolor          = var2_tcolor           #title color for var2 plot
      self.var1_levels          = var1_levels           #contour levels for var1 plot
      self.var2_levels          = var2_levels           #contour levels for var2 plot
      self.var1_cmap            = var1_cmap             #colormap for var1 plots
      self.var2_cmap            = var2_cmap             #colormap for var2 plots
      self.var3_level           = var3_level            #level for var3 plot
      self.var3_color           = var3_color            #color for var3 plot
      self.extend		= extend		#string for extend option in var1 plot (e.g. 'neither', 'max', ...)
      self.alpha		= alpha			#float value for transparency of var1 plot
      self.neighborhood		= neighborhood		#Neighborhood width for probability matched mean

#######################################################################################

def mymap(sw_lat, sw_lon, ne_lat, ne_lon, lat_1, lat_2, lat_0, lon_0, damage_files, resolution='i', area_thresh = 1000., counties='True'):

   #Creates basemap instance with Lambert Comformal projection for plotting

   #Dependencies:  

   #basemap, matplotlib, pyplot

   #Input:  

   #sw_lat - latitude (decimal deg.) of southwest corner of plotting domain
   #sw_lon - longitude (decimal deg.) of southwest corner of plotting domain
   #ne_lat - latitude (decimal deg.) of northeast corner of plotting domain
   #ne_lon - longitude (decimal deg.) of northeast corner of plotting domain
   #lat_1 - First true lat value for projection (typically 30.)
   #lat_2 - First true lat value for projection (typically 60.)
   #lat_0 - center lat value of projection
   #lon_0 - center lon value of projection
   #damage_files - array of paths to additional shapefiles to be plotted (assumed to be damage tracks)
   #resolution - resolution value for basemap instance (default = 'i' (intermediate))
   #area_thresh - area threshold for plotting water in basemap instance (default = 1000. km)

   #Returns: 
   #map - basemap instance for plotting

#######################

   map = Basemap(llcrnrlon=sw_lon, llcrnrlat=sw_lat, urcrnrlon=ne_lon, urcrnrlat=ne_lat, projection='lcc', lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0, resolution = resolution, area_thresh = area_thresh)

   if (counties == 'True'):
      map.drawcounties(linewidth=0.7, color=cb_colors.gray4)
   map.drawstates(linewidth=1., color=cb_colors.gray6)
   map.drawcoastlines(linewidth=1., color=cb_colors.gray5)
   map.drawcountries(linewidth=1., color=cb_colors.gray5)

   for f in range(0, len(damage_files)):
      map.readshapefile(damage_files[f], damage_files[f], linewidth=1.5) 

   return map

#######################################################################################

#######################################################################################

def mymap_boundaries(map, damage_files, counties='True'):

   #NOTE:  Used in 2017 processing for when the map instance has been pickled.  Speeds up plotting. 
   #Creates basemap instance with Lambert Comformal projection for plotting

   #Dependencies:  

   #basemap, matplotlib, pyplot

   #Input:  

   #map - basemap instance to draw boundaries on
   #damage_files - array of shapefiles (i.e. damage tracks) to be plotted
   #counties - boolean for whether or not to include county boundaries

   #Returns: 
   #map - updated basemap instance for plotting

#######################

   if (counties == 'True'):
      map.drawcounties(linewidth=0.7, color=cb_colors.gray4)
   map.drawstates(linewidth=1., color=cb_colors.gray6)
   map.drawcoastlines(linewidth=1., color=cb_colors.gray5)
   map.drawcountries(linewidth=1., color=cb_colors.gray5)

   for f in range(0, len(damage_files)):
      map.readshapefile(damage_files[f], damage_files[f], linewidth=3., color=cb_colors.gray9)

   return map

#######################################################################################

def create_fig(sw_lat, sw_lon, ne_lat, ne_lon, lat_1, lat_2, lat_0, lon_0, damage_files, resolution, area_thresh, verif='False', object='False'):

   #Return a standardized, empty figure for plotting.  Output figure should contain basemap with shapefile
   #plots, title, time initialized, time valid, and empty colortable axis
   
   #Dependencies:  

   #basemap, matplotlib, pyplot
   #mymap

   #Input:  

   #sw_lat - latitude (decimal deg.) of southwest corner of plotting domain
   #sw_lon - longitude (decimal deg.) of southwest corner of plotting domain
   #ne_lat - latitude (decimal deg.) of northeast corner of plotting domain
   #ne_lon - longitude (decimal deg.) of northeast corner of plotting domain
   #resolution - resolution value for basemap instance (default = 'c' (crude))
   #area_thresh - area threshold for plotting water in basemap instance (default = 10. km)
   #county_file - full path to a shapefile of counties to be plotted
   #damage_files - array of paths to additional shapefiles to be plotted (assumed to be damage tracks)

   #time_init - date/time string of initialization time
   #time_valid - date/time string of forecast time
   #title - string of plot title (1st value)
   #title_2 - string of plot title (2nd value)

   #Returns:
   #fig - empty, standardized figure for plotting

#######################

   fig = P.figure(figsize=(8.,9.))

   ax1 = fig.add_axes([0.01,0.09,0.98,0.85])
   ax1.spines["top"].set_alpha(0.9)
   ax1.spines["right"].set_alpha(0.9)
   ax1.spines["top"].set_linewidth(0.75)
   ax1.spines["right"].set_linewidth(0.75)
   ax1.spines["bottom"].set_alpha(0.9)
   ax1.spines["left"].set_alpha(0.9)
   ax1.spines["bottom"].set_linewidth(0.75)
   ax1.spines["left"].set_linewidth(0.75)

##### handle axes variations between paintball, object matching, and regular plots: #######

   if (verif == 'True'):   
      map = mymap(sw_lat, sw_lon, ne_lat, ne_lon, lat_1, lat_2, lat_0, lon_0, damage_files, resolution, area_thresh, counties='True') 
      ax2 = fig.add_axes([0.04,0.03,0.46,0.05])   
      ax3 = fig.add_axes([0.54,0.03,0.46,0.05])
   elif (object == 'True'):   
      map = mymap(sw_lat, sw_lon, ne_lat, ne_lon, lat_1, lat_2, lat_0, lon_0, damage_files, resolution, area_thresh, counties='True') 
      ax2 = fig.add_axes([0.02,0.06,0.46,0.02])   
      ax3 = fig.add_axes([0.52,0.06,0.46,0.02])
   else:
      map = mymap(sw_lat, sw_lon, ne_lat, ne_lon, lat_1, lat_2, lat_0, lon_0, damage_files, resolution, area_thresh, counties='True') 
      ax2 = fig.add_axes([0.04,0.06,0.6,0.02])   
      ax3 = fig.add_axes([0.66,0.03,0.3,0.05])

   return map, fig, ax1, ax2, ax3

#######################################################################################

def create_fig_nomap(verif='False', object='False'):

   #NOTE:  Used as of 2017 processing for when the map instance has been pickled.  Speeds up plotting. 

   #Return a standardized, empty figure for plotting.  Output figure should contain basemap with shapefile
   #plots, title, time initialized, time valid, and empty colortable axis

   #Dependencies:  

   #basemap, matplotlib, pyplot
   #mymap

   #Input:  

   #verif - boolean for creating a verification plot
   #object - boolean for creating an object-based verification plot 

   #Returns:
   #fig - empty, standardized figure for plotting
   #ax1, ax2, ax3 - three axes for plotting the data, colorbar, and contour intervals, respectively.

#######################

   fig = P.figure(figsize=(8.,9.))

   ax1 = fig.add_axes([0.01,0.09,0.98,0.85])
   ax1.spines["top"].set_alpha(0.9)
   ax1.spines["right"].set_alpha(0.9)
   ax1.spines["top"].set_linewidth(0.75)
   ax1.spines["right"].set_linewidth(0.75)
   ax1.spines["bottom"].set_alpha(0.9)
   ax1.spines["left"].set_alpha(0.9)
   ax1.spines["bottom"].set_linewidth(0.75)
   ax1.spines["left"].set_linewidth(0.75)

##### handle axes variations between paintball, object matching, and regular plots: #######

   if (verif == 'True'):
      ax2 = fig.add_axes([0.04,0.03,0.46,0.05])
      ax3 = fig.add_axes([0.54,0.03,0.46,0.05])
   elif (object == 'True'):
      ax2 = fig.add_axes([0.02,0.06,0.46,0.02])
      ax3 = fig.add_axes([0.52,0.06,0.46,0.02])
   else:
      ax2 = fig.add_axes([0.04,0.06,0.6,0.02])
      ax3 = fig.add_axes([0.66,0.03,0.3,0.05])

   return fig, ax1, ax2, ax3

#######################################################################################

#######################################################################################

def create_fig_tc(verif='False', object='False'):

   #NOTE:  Used as of 2017 processing for when the map instance has been pickled.  Speeds up plotting. 

   #Return a standardized, empty figure for plotting.  Output figure should contain basemap with shapefile
   #plots, title, time initialized, time valid, and empty colortable axis

   #Dependencies:  

   #basemap, matplotlib, pyplot
   #mymap

   #Input:  

   #verif - boolean for creating a verification plot
   #object - boolean for creating an object-based verification plot 

   #Returns:
   #fig - empty, standardized figure for plotting
   #ax1, ax2, ax3 - three axes for plotting the data, colorbar, and contour intervals, respectively.

#######################

   fig = P.figure(figsize=(8.,9.))

   ax1 = fig.add_axes([0.01,0.09,0.98,0.85])
   ax1.spines["top"].set_alpha(0.9)
   ax1.spines["right"].set_alpha(0.9)
   ax1.spines["top"].set_linewidth(0.75)
   ax1.spines["right"].set_linewidth(0.75)
   ax1.spines["bottom"].set_alpha(0.9)
   ax1.spines["left"].set_alpha(0.9)
   ax1.spines["bottom"].set_linewidth(0.75)
   ax1.spines["left"].set_linewidth(0.75)

##### handle axes variations between paintball, object matching, and regular plots: #######

   if (verif == 'True'):
      ax2 = fig.add_axes([0.04,0.03,0.46,0.05])
      ax3 = fig.add_axes([0.54,0.03,0.46,0.05])
   elif (object == 'True'):
      ax2 = fig.add_axes([0.02,0.06,0.46,0.02])
      ax3 = fig.add_axes([0.52,0.06,0.46,0.02])
   else:
      ax2 = fig.add_axes([0.04,0.0625,0.72,0.02])
      ax3 = fig.add_axes([0.04,0.0425,0.72,0.02])
      ax4 = fig.add_axes([0.04,0.0225,0.72,0.02])

   return fig, ax1, ax2, ax3, ax4

#######################################################################################

def env_plot(map, fig, ax1, ax2, ax3, x, y, plot, var1, var2, t, init_label, valid_label, domain, outdir, q_u, q_v, scale, t_star, t_plus, spec, quiv, showmax='False'):

   #This subroutine is used to create the majority of the NEWS-e plots used on the website
   #(all except the paintball, object-matching, and MRMS verification plots).  If passed an
   #empty basemap and figure object as well as the data and specifications for the plot, it
   #will create a plot over the model domain with appropriate labels, save the image us a 
   #.png using the appropriate timestamp, then remove the plots from the underlying basemap
   #template 

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted 
   #y - 2d numpy array of y axis, basemap-relative locations of data to be plotted 
   #plot - web_plot object to be plotted
   #var1 - 2d numpy array of the data to be plotted as a filled contour
   #var2 - 2d numpy array of the data to be plotted as a contour
   #t - time index associated with the plot (used for website indexing)
   #init_label - string of forecast initialization time
   #valid_label - string of time forecast valid
   #domain - string specifying subdomain over which data is plotted (currently always 'full') 
   #outdir - string of directory path to save .png image in
   #q_u - if quiv=True, u values to produce quiver plot
   #q_v - if quiv=True, v values to produce quiver plot
   #scale - if quiv=True, scaling factor for quiver vectors
   #t_star - multiplication factor for time index (used for website indexing)
   #t_plus - addition factor for time index (used for website indexing)
   #spec - boolean specifying if special contour level is to be plotted (e.g. 60 deg. dewpoint)
   #quiv - boolean specifying if quiver plot is to be overlain
   #showmax - boolean specifying if scatter of domain max value and associated label are to be plotted  

   #Returns:

   #Nothing, but all plots are removed from base template of the basemap after plot is saved

#######################

############# Calculate appropriate index for filename (e.g. 10-min forecast with 5-min output ###############
############# would be:  t_star = 5, t_plus = 0, t = 2, outtime = 'f010'                       ###############

   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   elif (len(temp_outtime) == 3): 
      outtime = 'f' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

   P.sca(ax1) #set plotting axis to 1

#   r1, r2, r3, r4 = pmm_rectangle(map, fig, x, y, plot) #plot pmm rectangle - deprecated

############# Create filled contour plot, with intervals highlighted for var1 data:    #################

   temp_cmap = plot.cmap
   temp_cmap.set_under(plot.under_color)
   temp_cmap.set_over(plot.over_color)

   p1 = map.contourf(x, y, var1, cmap=temp_cmap, levels = plot.var1_levels, alpha=plot.alpha, extend=plot.extend)
#   p1 = map.contourf(x, y, var1, cmap=plot.cmap, levels = plot.var1_levels, alpha=plot.alpha, extend=plot.extend)
   p2 = map.contour(x, y, var1, colors='k', levels = plot.var1_levels, linewidths=0.4, alpha=0.25, linestyles='solid')

#   p1.cmap.set_over(plot.over_color)
#   p1.cmap.set_under(plot.under_color)

############# Create contour plot for var2 data:  #################

   p3 = map.contour(x, y, var2, colors=plot.var2_colors, levels=plot.var2_levels, linewidths=[1., 1.5])

############# If special contour is to be plotted, plot it: ###############

   if (spec == 'True'):
      p4 = map.contour(x, y, var1, colors=plot.var1_cont_colors, levels=plot.var1_cont_levels, linewidths=1.25)

############ If quiver plot is to be included, include it: #################

   if (quiv == 'True'):
      quiv_x = x[0:-1:13,0:-1:13]
      quiv_y = y[0:-1:13,0:-1:13]
#      q1 = map.quiver(quiv_x, quiv_y, q_u, q_v, color='k', linewidth=0.3, minshaft=2.5, scale=scale, edgecolor='none', pivot='tail', alpha=0.35)
      barbs1, barbs2 = map.barbs(quiv_x, quiv_y, q_u, q_v, length=5, color='k', alpha=0.35, barb_increments=dict(half=5., full=10., flag=50.), fill_empty="True")

############ Add text labels of plot titles and initialization and valid forecast time: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold')
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)
   
############ If domain max value is to be plotted, plot it: #################

   if (showmax == 'True'): 
      mx, my = np.unravel_index(var1.argmax(), var1.shape)
      p4 = map.scatter(x[mx,my], y[mx,my], s=30, linewidth=0.85, marker='+', color='k', alpha=0.8)

########### Switch to axis 2 and set axis formatting: ################

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Get colormap information: #####################

   cmap = plot.cmap
   bounds = plot.var1_levels
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

########## Plot and format colormap for filled contour plot in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)
   cbar.set_label(plot.var1_title, fontsize=10)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=10)
   for label in cbar.ax.xaxis.get_ticklabels()[1::2]: #only label every other contour interval 
       label.set_visible(False)

########### Switch to axis 3 and set axis to be blank: ################

   P.sca(ax3)
   ax3.axis('off')   

########### Plot legend and description of contour plots (and max value if needed): ################

   for i in range(0, len(plot.var2_levels)):
      temp_label = str(plot.var2_levels[i])
      label = temp_label[0:2] + ' dBZ'

      P.axhline(y=0.75, xmin=(0.04 + (0.23 * i)), xmax=(0.17 + (0.23 * i)), linewidth=(1. + i), color=plot.var2_colors[i])
      P.text((0.025 + (0.235 * i)), 0.3, label, fontsize=10, color=plot.var2_colors[i])

   if (showmax == 'True'): 
      var1_max = str(np.max(var1))
      var1_label = 'Max Val.: ' + var1_max[0:6]
      P.text(0.55, 0.55, var1_label, fontsize=10) 

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight') 

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

   for coll in p1.collections:
      coll.remove()
   for coll in p2.collections:
      coll.remove()
   for coll in p3.collections:
      coll.remove()
   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()

#   rem_pmm_rect(r1, r2, r3, r4) #remove pmm rectangle plot - deprecated

   if (showmax == 'True'): 
      p4.remove()

   if (spec == 'True'):
      for coll in p4.collections: 
         coll.remove()
   if (quiv == 'True'):
      barbs1.remove()
      barbs2.remove()
#      q1.remove()

#######################################################################################

#######################################################################################

def multi_plot(map, fig, ax1, ax2, ax3, ax4, x, y, plot, var1, var2, var3, t, init_label, valid_label, domain, outdir, q_u, q_v, scale, t_star, t_plus, spec, quiv, showmax='False', showmin='False'):
   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   elif (len(temp_outtime) == 3):
      outtime = 'f' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

   P.sca(ax1) #set plotting axis to 1

############# Create filled contour plot, with intervals highlighted for var1 data:    #################

   #print plot.cmap

   temp_cmap = plot.cmap
   temp_cmap.set_under(plot.under_color)
   temp_cmap.set_over(plot.over_color)

#   p1 = map.contourf(x, y, var1, cmap=plot.cmap, levels = plot.var1_levels, alpha=plot.alpha, extend=plot.extend)

   #RAIN
   p3 = map.contourf(x, y, var3, cmap=cb_colors.all_greens_cmap, levels = plot.var1_levels, alpha=0.5, extend=plot.extend)
   # WIND
   p2 = map.contourf(x, y, var2, cmap=cb_colors.all_blues_cmap, levels = plot.var1_levels, alpha=0.5, extend=plot.extend)
   #TORN
   p1 = map.contourf(x, y, var1, cmap=cb_colors.all_reds_cmap, levels = plot.var1_levels, alpha=0.5, extend=plot.extend)

############ Add text labels of plot titles and initialization and valid forecast time: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold')
   t2 = fig.text(0.03, 0.945, 'Vert. Vort > 0.003 s$^{-1}$, Rain > 2 in., Wind > 64 kts', fontsize=10, fontweight='bold', color=cb_colors.gray7)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)

########### Switch to axis 2 and set axis formatting: ################

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Get colormap information: #####################

   cmap = cb_colors.all_reds_cmap
   #cmap = cb_colors.all_blues_cmap
   bounds = plot.var1_levels
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

########## Plot and format colormap for filled contour plot in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)
   #cbar.set_label(plot.var1_title, fontsize=10)
   #cl = P.getp(cbar.ax, 'xmajorticklabels')
   #P.setp(cl, fontsize=12)
   for label in cbar.ax.xaxis.get_ticklabels()[:]: #only label every other contour interval 
       label.set_visible(False)

########### Switch to axis 3 and set axis to be blank: ################

   P.sca(ax3)
   cmap = cb_colors.all_greens_cmap
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
   cbar = matplotlib.colorbar.ColorbarBase(ax3, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)
   for label in cbar.ax.xaxis.get_ticklabels()[:]: #only label every other contour interval 
       label.set_visible(False)

   P.sca(ax4)
   cmap = cb_colors.all_blues_cmap
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
   cbar = matplotlib.colorbar.ColorbarBase(ax4, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)

   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=12)
#   for label in cbar.ax.xaxis.get_ticklabels()[1::2]: #only label every other contour interval 
#       label.set_visible(False)

   #ax3.axis('off')
   t7 = fig.text(0.78, 0.070, '0-2 km Vorticity > 0.0015 s-1', fontsize=10, fontweight='bold', color=cb_colors.red9)
   t8 = fig.text(0.78, 0.048, 'Accumulated Rainfall > 2 in', fontsize=10, fontweight='bold', color=cb_colors.green9)
   t9 = fig.text(0.78, 0.026, 'Sustained wind > 64 kt', fontsize=10, fontweight='bold', color=cb_colors.blue9)


########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight')

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

   for coll in p1.collections:
      coll.remove()
   for coll in p2.collections:
      coll.remove()
   for coll in p3.collections:
      coll.remove()
   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()


#######################################################################################

def paintqc_plot(map, fig, ax1, ax2, ax3, x, y, vx, vy, plot, var1, var2, radmask, t, init_label, valid_label, domain, outdir, t_star, t_plus, blank="False"):

   #This subroutine is used to create paintball plots used on the website.  If passed an
   #empty basemap and figure object as well as the data and specifications for the plot, it
   #will create a plot over the model domain with appropriate labels, save the image us a
   #.png using the appropriate timestamp, then remove the plots from the underlying basemap
   #template

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, v_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted
   #y - 2d numpy array of y-axis, basemap-relative locations of data to be plotted
   #vx - 2d numpy array of x-axis, basemap-relative locations of verification data to be plotted
   #vy - 2d numpy array of y-axis, basemap-relative locations of verification data to be plotted
   #plot - web_plot object to be plotted
   #var1 - 2d numpy array of the data to be plotted as a filled contour
   #var2 - 2d numpy array of the data to be plotted as a contour
   #radmask - 2d, binary numpy array with regions near/distant to 88d sites set to 1 (set all to 0 to turn off)
   #t - time index associated with the plot (used for website indexing)
   #init_label - string of forecast initialization time
   #valid_label - string of time forecast valid
   #domain - string specifying subdomain over which data is plotted (currently always 'full')
   #outdir - string of directory path to save .png image in
   #t_star - multiplication factor for time index (used for website indexing)
   #t_plus - addition factor for time index (used for website indexing)
   #blank - boolean for when plot is during a period when rotation track objects are not available (i.e. < time_window from beginning or end)

   #Returns:

   #Nothing, but all plots are removed from base template of the basemap after plot is saved

#######################

############# Calculate appropriate index for filename (e.g. 10-min forecast with 5-min output ###############
############# would be:  t_star = 5, t_plus = 0, t = 2, outtime = 'f010'                       ###############

   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   elif (len(temp_outtime) == 3):
      outtime = 'f' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

   P.sca(ax1) #set current plotting axis to Axis 1

############ Add text labels of plot titles and initialization and valid forecast time: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold', color=plot.var1_tcolor)
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)

############ If forecast objects are present (i.e. not in blanking period at beginning/end of forecast), ###############
############ then shade masked locations (if present), plot each member's data in a different color, and  ##############
############ shade verification objects in dark gray.  If forecast objects are not present, add label.  ##############

   if (blank == 'False'): 
      p3 = map.contourf(vx, vy, radmask, colors=[cb_colors.gray2, cb_colors.gray2], levels=[0.1, 1000.], alpha=0.8)
      members = np.arange(0,var2.shape[0])
#      np.random.shuffle(members)
      for n, member in enumerate(members):
         p2 = map.contourf(x, y, var2[member,:,:], colors=[plot.var2_colors[member], plot.var2_colors[member]], levels=plot.var2_levels, alpha=0.15, extend=plot.extend)
         p2 = map.contour(x, y, var2[member,:,:], colors=[plot.var2_colors[member], plot.var2_colors[member]], levels=plot.var2_levels, alpha=0.9, linewidths=0.5, linestylyes='solid')
      p1 = map.contourf(vx, vy, var1, colors=plot.var1_colors, levels=plot.var1_levels, alpha=0.5, extend=plot.extend)
   else:
      t5 = fig.text(0.07,0.6, 'ROTATION TRACKS NOT AVAILABLE!', fontsize=24, fontweight='bold', color=plot.var1_tcolor)

############ Set Axes 2 and 3 to be blank: ###############

   P.sca(ax2)
   ax2.axis('off')

   P.sca(ax3)
   ax3.axis('off')

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight')

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

#   if (np.max(var1) > 0.): 
   if (blank == 'False'): 
      for coll in p1.collections:
         coll.remove()
      for coll in p2.collections:
         coll.remove()
      for coll in p3.collections:
         coll.remove()
   else:
      t5.remove()

   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()

#######################

#######################################################################################

def painttl_plot(map, fig, ax1, ax2, ax3, x, y, vx, vy, plot, var1, var2, radmask, t, init_label, valid_label, plot_cmap, fcst_times, domain, outdir, t_star, t_plus, blank="False"):

   #This subroutine is used to create paintball plots used on the website.  If passed an
   #empty basemap and figure object as well as the data and specifications for the plot, it
   #will create a plot over the model domain with appropriate labels, save the image us a
   #.png using the appropriate timestamp, then remove the plots from the underlying basemap
   #template

   #This differs from paintqc in that it is intended to plot paintball plots from multiple
   #forecasts valid at the same time on the same plot ('time-lagged paintball') 

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, v_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted
   #y - 2d numpy array of y-axis, basemap-relative locations of data to be plotted
   #vx - 2d numpy array of x-axis, basemap-relative locations of verification data to be plotted
   #vy - 2d numpy array of y-axis, basemap-relative locations of verification data to be plotted
   #plot - web_plot object to be plotted
   #var1 - 2d numpy array of the data to be plotted as a filled contour
   #var2 - 2d numpy array of the data to be plotted as a contour
   #radmask - 2d, binary numpy array with regions near/distant to 88d sites set to 1 (set all to 0 to turn off)
   #t - time index associated with the plot (used for website indexing)
   #init_label - string of forecast initialization time
   #valid_label - string of time forecast valid
   #domain - string specifying subdomain over which data is plotted (currently always 'full')
   #outdir - string of directory path to save .png image in
   #t_star - multiplication factor for time index (used for website indexing)
   #t_plus - addition factor for time index (used for website indexing)
   #blank - boolean for when plot is during a period when rotation track objects are not available (i.e. < time_window from beginning or end)

   #Returns:

   #Nothing, but all plots are removed from base template of the basemap after plot is saved

#######################

############# Calculate appropriate index for filename (e.g. 10-min forecast with 5-min output ###############
############# would be:  t_star = 5, t_plus = 0, t = 2, outtime = 'f010'                       ###############

   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   elif (len(temp_outtime) == 3):
      outtime = 'f' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

   P.sca(ax1) #set current plotting axis to Axis 1

############ Add text labels of plot titles and initialization and valid forecast time: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold', color=plot.var1_tcolor)
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)

############ If forecast objects are present (i.e. not in blanking period at beginning/end of forecast), ###############
############ then shade masked locations (if present), plot each member's data in a different color, and  ##############
############ shade verification objects in dark gray.  If forecast objects are not present, add label.  ##############

   if (blank == 'False'): 
      p3 = map.contourf(vx, vy, radmask, colors=[cb_colors.gray2, cb_colors.gray2], levels=[0.1, 1000.], alpha=0.8)
      members = np.arange(0,var2.shape[1])
      for f, fcst in enumerate(fcst_times): 
         for n, member in enumerate(members):
            p2 = map.contourf(x, y, var2[f,member,:,:], colors=[plot.var2_colors[f], plot.var2_colors[f]], levels=plot.var2_levels, alpha=0.25, extend=plot.extend)
#            p2 = map.contour(x, y, var2[f,member,:,:], colors=[plot.var2_colors[f], plot.var2_colors[f]], levels=plot.var2_levels, alpha=0.6, linewidths=0.5, linestylyes='solid')
         p1 = map.contourf(vx, vy, var1, colors=plot.var1_colors, levels=plot.var1_levels, alpha=0.5, extend=plot.extend)
   else:
      t5 = fig.text(0.07,0.6, 'ROTATION TRACKS NOT AVAILABLE!', fontsize=24, fontweight='bold', color=plot.var1_tcolor)

############ plot colormap in Axis 2: ###############

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Get colormap information: #####################

   cmap = plot_cmap
   bounds = np.arange(0,(len(fcst_times)+1))
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

########## Plot and format colormap for filled contour plot in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', alpha=plot.alpha)
   cbar.set_label('Forecast Initialization Time (UTC)', fontsize=10)
   cbar.ax.set_xticklabels(fcst_times)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=10)
   for label in cbar.ax.xaxis.get_ticklabels()[1::2]: #only label every other contour interval 
       label.set_visible(False)

############ Set Axis 3 to be blank: ###############

   P.sca(ax3)
   ax3.axis('off')

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight')

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

#   if (np.max(var1) > 0.): 
   if (blank == 'False'): 
      for coll in p1.collections:
         coll.remove()
      for coll in p2.collections:
         coll.remove()
      for coll in p3.collections:
         coll.remove()
   else:
      t5.remove()

   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()

#######################################################################################

def ob_ver_plot(map, fig, ax1, ax2, ax3, x, y, vx, vy, plot, var1, var2, var3, radmask, t, pod, far, bias, csi, init_label, valid_label, domain, outdir, t_star, t_plus, blank='False'):

   #This subroutine is used to create object-matching plots used on the website.  If passed an
   #empty basemap and figure object as well as the data and specifications for the plot, it
   #will create a plot over the model domain with appropriate labels, save the image us a
   #.png using the appropriate timestamp, then remove the plots from the underlying basemap
   #template

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, ob_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted
   #y - 2d numpy array of y-axis, basemap-relative locations of data to be plotted
   #vx - 2d numpy array of x-axis, basemap-relative locations of verification data to be plotted
   #vy - 2d numpy array of y-axis, basemap-relative locations of verification data to be plotted
   #plot - web_plot object to be plotted
   #var1 - 2d numpy array of the first variable (e.g. matches) to be plotted as a pcolormesh
   #var2 - 2d numpy array of the second variable (e.g. false alarms) to be plotted as a pcolormesh
   #var3 - 2d numpy array of the verification data to be shaded
   #radmask - 2d, binary numpy array with regions near/distant to 88d sites set to 1 (set all to 0 to turn off)
   #t - time index associated with the plot (used for website indexing)
   #ots - real value of the total interest-weighted, object-based threat score to display in the lower-right corner 
   #bin_ots - real value of the binary, object-based threat score to display in the lower-right corner 
   #init_label - string of forecast initialization time
   #valid_label - string of time forecast valid
   #domain - string specifying subdomain over which data is plotted (currently always 'full')
   #outdir - string of directory path to save .png image in
   #t_star - multiplication factor for time index (used for website indexing)
   #t_plus - addition factor for time index (used for website indexing)

   #Returns:

   #Nothing, but all plots are removed from base template of the basemap after plot is saved

#######################

############# Calculate appropriate index for filename (e.g. 10-min forecast with 5-min output ###############
############# would be:  t_star = 5, t_plus = 0, t = 2, outtime = 'f010'                       ###############

   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   elif (len(temp_outtime) == 3): 
      outtime = 'f' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

############ Convert OTS values to strings with labels: #################

   str_pod = str(pod)
   str_far = str(far)
   str_bias = str(bias)
   str_csi = str(csi)
#   str_ots = str(ots)
#   str_bin_ots = str(bin_ots)

   if (len(str_pod) > 4): 
      pod_text = 'POD: ' + str_pod[0:4]
   else: 
      pod_text = 'POD: ' + str_pod

   if (len(str_far) > 4): 
      far_text = 'FAR: ' + str_far[0:4]
   else: 
      far_text = 'FAR: ' + str_far

   if (len(str_bias) > 4): 
      bias_text = 'BIAS: ' + str_bias[0:4]
   else: 
      bias_text = 'BIAS: ' + str_bias

   if (len(str_csi) > 4): 
      csi_text = 'CSI: ' + str_csi[0:4]
   else: 
      csi_text = 'CSI: ' + str_csi

#   ots_text = 'OTS: ' + str_ots[0:4]
#   bin_ots_text = 'Binary OTS: ' + str_bin_ots[0:4]

   P.sca(ax1) #switch to Axis 1

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold', color=plot.var1_tcolor)
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)

########## Get colormap information for both colormaps: #####################

   bounds1 = plot.var1_levels
   cmap1 = plot.var1_cmap
   norm1 = matplotlib.colors.BoundaryNorm(bounds1, cmap1.N)

   bounds2 = plot.var2_levels
   cmap2 = plot.var2_cmap
   norm2 = matplotlib.colors.BoundaryNorm(bounds2, cmap2.N)

   if (blank == 'False'):

########## Plot radmask, var1 and var2 as a pcolormesh, var3 as a shaded contour: #####################

      p4 = map.contourf(vx, vy, radmask, colors=[cb_colors.gray2, cb_colors.gray2], levels=[0.1, 1000.], alpha=0.8)
      p1 = map.pcolormesh(x, y, var1, cmap=cmap1, norm=norm1, vmin=np.min(bounds1), vmax=np.max(bounds1), alpha=plot.alpha)  
      p2 = map.pcolormesh(x, y, var2, cmap=cmap2, norm=norm2, vmin=np.min(bounds2), vmax=np.max(bounds2), alpha=plot.alpha)  
      p3 = map.contourf(vx, vy, var3, colors=plot.var3_color, levels=plot.var3_level, alpha=0.6) 

############ Add text labels of plot titles, initialization and valid forecast time, and OTS info: ###############

      t5 = fig.text(0.09, 0.1125, pod_text, fontsize=16)
      t6 = fig.text(0.31, 0.1125, far_text, fontsize=16)
      t7 = fig.text(0.53, 0.1125, bias_text, fontsize=16)
      t8 = fig.text(0.75, 0.1125, csi_text, fontsize=16)

   else:
      t5 = fig.text(0.07,0.6, 'ROTATION TRACKS NOT AVAILABLE!', fontsize=24, fontweight='bold', color=plot.var1_tcolor)

########### Switch to axis 2 and set axis formatting: ################

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Plot and format colormap for pcolormesh plot for var1 in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap1, norm=norm1, ticks=bounds1, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)
   cbar.set_label(plot.var1_title, fontsize=10)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=10)
   for label in cbar.ax.xaxis.get_ticklabels()[1::2]:
       label.set_visible(False)

########### Switch to axis 3 and set axis formatting: ################

   P.sca(ax3)
   ax3.spines["top"].set_alpha(0.9)
   ax3.spines["right"].set_alpha(0.9)
   ax3.spines["top"].set_linewidth(0.5)
   ax3.spines["right"].set_linewidth(0.5)
   ax3.spines["bottom"].set_alpha(0.9)
   ax3.spines["left"].set_alpha(0.9)
   ax3.spines["bottom"].set_linewidth(0.5)
   ax3.spines["left"].set_linewidth(0.5)

########## Plot and format colormap for pcolormesh plot for var2 in axis 3: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax3, cmap=cmap2, norm=norm2, ticks=bounds2, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)
   cbar.set_label(plot.var2_title, fontsize=10)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=10)
   for label in cbar.ax.xaxis.get_ticklabels()[1::2]:
       label.set_visible(False)

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   print( fig_name)
   P.savefig(fig_name, format="png", bbox_inches='tight')

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

   if (blank == 'False'):
      p1.remove()
      p2.remove()
      for coll in p3.collections:
         coll.remove()
      for coll in p4.collections:
         coll.remove()

   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()
   t5.remove()

   if (blank == 'False'):
      t6.remove()
      t7.remove()
      t8.remove()

#######################


#######################################################################################

def ob_cent_plot(map, fig, ax1, ax2, ax3, x, y, vx, vy, plot, var1, var3, t, init_label, valid_label, domain, outdir, t_star
, t_plus):

   #This subroutine is used to create object centroid scatter plots.  If passed an
   #empty basemap and figure object as well as the data and specifications for the plot, it
   #will create a plot over the model domain with appropriate labels, save the image us a
   #.png using the appropriate timestamp, then remove the plots from the underlying basemap
   #template

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, ob_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted
   #y - 2d numpy array of y-axis, basemap-relative locations of data to be plotted
   #vx - 2d numpy array of x-axis, basemap-relative locations of verification data to be plotted
   #vy - 2d numpy array of y-axis, basemap-relative locations of verification data to be plotted
   #plot - web_plot object to be plotted
   #var1 - 2d numpy array of the variable (e.g. centroid locations) to be scattered
   #var3 - 2d numpy array of the verification data to be shaded
   #t - time index associated with the plot (used for website indexing)
   #init_label - string of forecast initialization time
   #valid_label - string of time forecast valid
   #domain - string specifying subdomain over which data is plotted (currently always 'full')
   #outdir - string of directory path to save .png image in
   #t_star - multiplication factor for time index (used for website indexing)
   #t_plus - addition factor for time index (used for website indexing)

   #Returns:

   #Nothing, but all plots are removed from base template of the basemap after plot is saved

#######################

############# Calculate appropriate index for filename (e.g. 10-min forecast with 5-min output ###############
############# would be:  t_star = 5, t_plus = 0, t = 2, outtime = 'f010'                       ###############

   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

   P.sca(ax1) #switch to Axis 1

########## Get colormap information for scatter plot: #####################

   bounds1 = plot.var1_levels
   cmap1 = plot.var1_cmap
   norm1 = matplotlib.colors.BoundaryNorm(bounds1, cmap1.N)

########## Scatter object centroids, with coloring scaled by var1 levels, shade verification field: #####################

   p1 = map.scatter(x, y, c=var1, s=15, linewidth=0.75, marker='+', cmap=cmap1, vmin=np.min(bounds1), vmax=np.max(bounds1), alpha=plot.alpha)
   p3 = map.contourf(vx, vy, var3, colors=plot.var3_color, levels=plot.var3_level, alpha=0.6)

############ Add text labels of plot titles, initialization and valid forecast time, and OTS info: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold', color=plot.var1_tcolor)
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)

########### Switch to axis 2 and set axis formatting: ################

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Plot and format colormap for scatter plot in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap1, norm=norm1, ticks=bounds1, orientation = 'horizontal', alpha=plot.alpha)
   cbar.set_label('Total Interest Score', fontsize=10)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=10)
   for label in cbar.ax.xaxis.get_ticklabels()[1::2]:
       label.set_visible(False)

########### Set Axis 3 to blank: ################

   P.sca(ax3)
   ax3.axis('off')   

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight')

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

   p1.remove()
   for coll in p3.collections:
      coll.remove()

   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()

#######################################################################################

def dot_plot(map, fig, ax1, ax2, ax3, x, y, ob_x, ob_y, ob_var, ob_levels, ob_cmap, ob_label, plot, var1, var2, t, init_label, valid_label, domain, outdir, q_u, q_v, scale, t_star, t_plus, spec, quiv, showmax='False'):

   #This subroutine is the same as env_plot with the exception that scatter plots of OK
   #mesonet observations differences from model fields are overlain in a scatter plot. 

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted
   #y - 2d numpy array of y axis, basemap-relative locations of data to be plotted
   #ob_x - 1d numpy array of the x positions (basemap relative) of OK mesonet stations to be plotted 
   #ob_y - 1d numpy array of the y positions (basemap relative) of OK mesonet stations to be plotted 
   #ob_var - 1d numpy array of OK mesonet data to be plotted
   #ob_levels - intervals for colormap of OK mesonet plot data
   #ob_cmap - colormap for OK mesonet plot data
   #ob_label - string of label for OK mesonet plot data
   #plot - web_plot object to be plotted
   #var1 - 2d numpy array of the data to be plotted as a filled contour
   #var2 - 2d numpy array of the data to be plotted as a contour
   #t - time index associated with the plot (used for website indexing)
   #init_label - string of forecast initialization time
   #valid_label - string of time forecast valid
   #domain - string specifying subdomain over which data is plotted (currently always 'full')
   #outdir - string of directory path to save .png image in
   #q_u - if quiv=True, u values to produce quiver plot
   #q_v - if quiv=True, v values to produce quiver plot
   #scale - if quiv=True, scaling factor for quiver vectors
   #t_star - multiplication factor for time index (used for website indexing)
   #t_plus - addition factor for time index (used for website indexing)
   #spec - boolean specifying if special contour level is to be plotted (e.g. 60 deg. dewpoint)
   #quiv - boolean specifying if quiver plot is to be overlain
   #showmax - boolean specifying if scatter of domain max value and associated label are to be plotted

   #Returns:

   #Nothing, but all plots are removed from base template of the basemap after plot is saved

#######################

############# Calculate appropriate index for filename (e.g. 10-min forecast with 5-min output ###############
############# would be:  t_star = 5, t_plus = 0, t = 2, outtime = 'f010'                       ###############

   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   elif (len(temp_outtime) == 3):
      outtime = 'f' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

   P.sca(ax1) #set plotting axis to 1

############# Create filled contour plot, with intervals highlighted for var1 data:    #################

   p1 = map.contourf(x, y, var1, cmap=plot.cmap, levels = plot.var1_levels, alpha=plot.alpha, extend=plot.extend)
   p2 = map.contour(x, y, var1, colors='k', levels = plot.var1_levels, linewidths=0.4, alpha=0.25, linestyles='solid')

   p1.cmap.set_over(plot.over_color)
   p1.cmap.set_under(plot.under_color)

############# Create contour plot for var2 data:  #################

   p3 = map.contour(x, y, var2, colors=plot.var2_colors, levels=plot.var2_levels, linewidths=[1., 1.5])

############# If special contour is to be plotted, plot it: ###############

   if (spec == 'True'):
      p4 = map.contour(x, y, var1, colors=plot.var1_cont_colors, levels=plot.var1_cont_levels, linewidths=1.25)

############ If quiver plot is to be included, include it: #################

   if (quiv == 'True'):
      quiv_x = x[0:-1:8,0:-1:8]
      quiv_y = y[0:-1:8,0:-1:8]

      q1 = map.quiver(quiv_x, quiv_y, q_u, q_v, color='k', linewidth=0.3, minshaft=2.5, scale=scale, edgecolor='none', pivot='tail', alpha=0.35)

############ Add text labels of plot titles and initialization and valid forecast time: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold')
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)

############ If domain max value is to be plotted, plot it: #################

   if (showmax == 'True'):
      mx, my = np.unravel_index(var1.argmax(), var1.shape)
      p4 = map.scatter(x[mx,my], y[mx,my], s=30, linewidth=0.85, marker='+', color='k', alpha=0.8)

############ Scatter OK mesonet data: #################

   p5 = map.scatter(ob_x, ob_y, c=ob_var, s=300, linewidth=0.75, edgecolors=cb_colors.gray8, cmap=plot.cmap, alpha=1., vmin=np.min(plot.var1_levels), vmax=np.max(plot.var1_levels))
#old way   p5 = map.scatter(ob_x, ob_y, c=ob_var, s=200, linewidth=0.5, edgecolors=cb_colors.gray8, cmap=ob_cmap, alpha=0.8, vmin=np.min(ob_levels), vmax=np.max(ob_levels))

########### Switch to axis 2 and set axis formatting: ################

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Get colormap information: #####################

   cmap = plot.cmap
   bounds = plot.var1_levels
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

########## Plot and format colormap for filled contour plot in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)   
   cbar.set_label(plot.var1_title, fontsize=10)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=10)
   for label in cbar.ax.xaxis.get_ticklabels()[1::2]:
       label.set_visible(False)

########### Switch to axis 3 and set axis to be blank: ################

   P.sca(ax3)
   ax3.axis('off')

########### Plot legend and description of contour plots (and max value if needed): ################

   for i in range(0, len(plot.var2_levels)):
      temp_label = str(plot.var2_levels[i])
      label = temp_label[0:2] + ' dBZ'

      P.axhline(y=0.75, xmin=(0.04 + (0.23 * i)), xmax=(0.17 + (0.23 * i)), linewidth=(1. + i), color=plot.var2_colors[i])
      P.text((0.025 + (0.235 * i)), 0.3, label, fontsize=10, color=plot.var2_colors[i])

   if (showmax == 'True'):
      var1_max = str(np.max(var1))
      var1_label = 'Max Val.: ' + var1_max[0:6]
      P.text(0.55, 0.55, var1_label, fontsize=10)

#
########### Switch to axis 3 and set axis formatting: ################
#
#   P.sca(ax3)
#   ax3.spines["top"].set_alpha(0.9)
#   ax3.spines["right"].set_alpha(0.9)
#   ax3.spines["top"].set_linewidth(0.5)
#   ax3.spines["right"].set_linewidth(0.5)
#   ax3.spines["bottom"].set_alpha(0.9)
#   ax3.spines["left"].set_alpha(0.9)
#   ax3.spines["bottom"].set_linewidth(0.5)
#   ax3.spines["left"].set_linewidth(0.5)

########## Get OK Mesonet colormap information: #####################

#   cmap2 = ob_cmap
#   bounds2 = ob_levels
#   norm2 = matplotlib.colors.BoundaryNorm(bounds2, cmap2.N)

########## Plot and format colormap for OK Mesonet plot in axis 3: #####################

#   cbar2 = matplotlib.colorbar.ColorbarBase(ax3, cmap=cmap2, norm=norm2, ticks=bounds2, orientation = 'horizontal', extend=plot.extend, alpha=0.8)
#   cbar2.set_label(ob_label, fontsize=10)
#   cl2 = P.getp(cbar2.ax, 'xmajorticklabels')
#   P.setp(cl2, fontsize=10)
#   for label in cbar2.ax.xaxis.get_ticklabels()[1::2]:
#       label.set_visible(False)

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight')

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

   p5.remove()

   for coll in p1.collections:
      coll.remove()
   for coll in p2.collections:
      coll.remove()
   for coll in p3.collections:
      coll.remove()

   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()

   if (showmax == 'True'):
      p4.remove()

   if (spec == 'True'):
      for coll in p4.collections:
         coll.remove()
   if (quiv == 'True'):
      q1.remove()

#######################################################################################

#######################################################################################

def mem_plot(map, fig, ax1, ax2, ax3, x, y, plot, var1, var2, t, init_label, valid_label, domain, outdir, q_u, q_v, scale, t_star, t_plus, spec, quiv, showmax='False'):
#### FIX COMMENTS PSS
   #This subroutine is used to create the majority of the NEWS-e plots used on the website
   #(all except the paintball, object-matching, and MRMS verification plots).  If passed an
   #empty basemap and figure object as well as the data and specifications for the plot, it
   #will create a plot over the model domain with appropriate labels, save the image us a 
   #.png using the appropriate timestamp, then remove the plots from the underlying basemap
   #template 

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted 
   #y - 2d numpy array of y axis, basemap-relative locations of data to be plotted 
   #plot - web_plot object to be plotted
   #var1 - 2d numpy array of the data to be plotted as a filled contour
   #var2 - 2d numpy array of the data to be plotted as a contour
   #t - time index associated with the plot (used for website indexing)
   #init_label - string of forecast initialization time
   #valid_label - string of time forecast valid
   #domain - string specifying subdomain over which data is plotted (currently always 'full') 
   #outdir - string of directory path to save .png image in
   #q_u - if quiv=True, u values to produce quiver plot
   #q_v - if quiv=True, v values to produce quiver plot
   #scale - if quiv=True, scaling factor for quiver vectors
   #t_star - multiplication factor for time index (used for website indexing)
   #t_plus - addition factor for time index (used for website indexing)
   #spec - boolean specifying if special contour level is to be plotted (e.g. 60 deg. dewpoint)
   #quiv - boolean specifying if quiver plot is to be overlain
   #showmax - boolean specifying if scatter of domain max value and associated label are to be plotted  

   #Returns:

   #Nothing, but all plots are removed from base template of the basemap after plot is saved

#######################

############# Calculate appropriate index for filename (e.g. 10-min forecast with 5-min output ###############
############# would be:  t_star = 5, t_plus = 0, t = 2, outtime = 'f010'                       ###############

   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   elif (len(temp_outtime) == 3): 
      outtime = 'f' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

   P.sca(ax1) #set plotting axis to 1

#   r1, r2, r3, r4 = pmm_rectangle(map, fig, x, y, plot) #plot pmm rectangle - deprecated

############# Create contour plot for var2 data:  #################

   p3 = map.contourf(x, y, var2, colors=plot.var2_colors, levels=plot.var2_levels)

############# Create filled contour plot, with intervals highlighted for var1 data:    #################

   p1 = map.contourf(x, y, var1, cmap=plot.cmap, levels = plot.var1_levels, alpha=plot.alpha, extend=plot.extend)
   p2 = map.contour(x, y, var1, colors='k', levels = plot.var1_levels, linewidths=0.4, alpha=0.25, linestyles='solid')

   p1.cmap.set_over(plot.over_color)
   p1.cmap.set_under(plot.under_color)

############# If special contour is to be plotted, plot it: ###############

   if (spec == 'True'):
      p4 = map.contour(x, y, var1, colors=plot.var1_cont_colors, levels=plot.var1_cont_levels, linewidths=1.25)

############ If quiver plot is to be included, include it: #################

   if (quiv == 'True'):
      quiv_x = x[0:-1:8,0:-1:8]
      quiv_y = y[0:-1:8,0:-1:8]

      q1 = map.quiver(quiv_x, quiv_y, q_u, q_v, color='k', linewidth=0.3, minshaft=2.5, scale=scale, edgecolor='none', pivot='tail', alpha=0.35)

############ Add text labels of plot titles and initialization and valid forecast time: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold')
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)
   
############ If domain max value is to be plotted, plot it: #################

   if (showmax == 'True'): 
      mx, my = np.unravel_index(var1.argmax(), var1.shape)
      p4 = map.scatter(x[mx,my], y[mx,my], s=30, linewidth=0.85, marker='+', color='k', alpha=0.8)

########### Switch to axis 2 and set axis formatting: ################

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Get colormap information: #####################

   cmap = plot.cmap
   bounds = plot.var1_levels
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

########## Plot and format colormap for filled contour plot in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)
   cbar.set_label(plot.var1_title, fontsize=10)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=10)
   for label in cbar.ax.xaxis.get_ticklabels()[1::2]: #only label every other contour interval 
       label.set_visible(False)

########### Switch to axis 3 and set axis to be blank: ################

   P.sca(ax3)
   ax3.axis('off')   

########### Plot legend and description of contour plots (and max value if needed): ################

   for i in range(0, (len(plot.var2_levels)-1)):
      temp_label = str(plot.var2_levels[i])
      label = temp_label + ' m$^{2}$ s$^{-2}$'

      P.axhline(y=0.75, xmin=(0.04 + (0.5 * i)), xmax=(0.4 + (0.45 * i)), linewidth=(2.5), color=plot.var2_colors[i])
      P.text((0.025 + (0.5 * i)), 0.25, label, fontsize=10, color=plot.var2_colors[i])

   if (showmax == 'True'): 
      var1_max = str(np.max(var1))
      var1_label = 'Max Val.: ' + var1_max[0:6]
      P.text(0.55, 0.55, var1_label, fontsize=10) 

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight') 

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

   for coll in p1.collections:
      coll.remove()
   for coll in p2.collections:
      coll.remove()
   for coll in p3.collections:
      coll.remove()
   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()

#   rem_pmm_rect(r1, r2, r3, r4) #remove pmm rectangle plot - deprecated

   if (showmax == 'True'): 
      p4.remove()

   if (spec == 'True'):
      for coll in p4.collections: 
         coll.remove()
   if (quiv == 'True'):
      q1.remove()

#######################################################################################

def mrms_plot(map, fig, ax1, ax2, ax3, x, y, plot, var1, var2, t, init_label, valid_label, domain, outdir, fnum, showmax='False'):

   #This subroutine will plot MRMS data on the NEWS-e grid for comparison.  Plots of 
   #composite reflectivity, with Az. shear rotation tracks, and NWS warnings and LSRs
   #will be produced. 

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted
   #y - 2d numpy array of y axis, basemap-relative locations of data to be plotted
   #plot - web_plot object to be plotted
   #var1 - 2d numpy array of the data to be plotted as a filled contour (i.e. reflectivity)
   #var2 - 2d numpy array of the data to be plotted as a single filled contour (i.e. az. shear)
   #t - time index associated with the plot (used for website indexing)
   #init_label - string of forecast initialization time
   #valid_label - string of time forecast valid
   #domain - string specifying subdomain over which data is plotted (currently always 'full')
   #outdir - string of directory path to save .png image in
   #fnum - string of time index of forecast (i.e. '000' or '055') (calculated within the subroutine for other plots)
   #showmax - boolean specifying if scatter of domain max value and associated label are to be plotted

   #Returns:

   #Nothing, but all plots are removed from base template of the basemap after plot is saved

#######################

   outtime = 'f' + fnum #create figure index label

   P.sca(ax1) #set current axis to 1

############# Create filled single contour plot for var2 data:    #################

   p3 = map.contourf(x, y, var2, colors=plot.var2_colors, levels=plot.var2_levels, alpha=1.)

############# Create filled contour plot, with intervals highlighted for var1 data:    #################

   p1 = map.contourf(x, y, var1, cmap=plot.cmap, levels = plot.var1_levels, alpha=plot.alpha, extend=plot.extend)
   p2 = map.contour(x, y, var1, colors='k', levels = plot.var1_levels, linewidths=0.4, alpha=0.25, linestyles='solid')

   p1.cmap.set_over(plot.over_color)
   p1.cmap.set_under(plot.under_color)

############ Add text labels of plot titles and initialization and valid forecast time: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold')
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)

############ If domain max value is to be plotted, plot it: #################

   if (showmax == 'True'):
      mx, my = np.unravel_index(var1.argmax(), var1.shape)
      p4 = map.scatter(x[mx,my], y[mx,my], s=30, linewidth=0.85, marker='+', color='k', alpha=0.8)

########### Switch to axis 2 and set axis formatting: ################

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Get colormap information: #####################

   cmap = plot.cmap
   bounds = plot.var1_levels
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

########## Plot and format colormap for filled contour plot in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)
   cbar.set_label(plot.var1_title, fontsize=10)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=10)
   for label in cbar.ax.xaxis.get_ticklabels()[1::2]:
       label.set_visible(False)

############ Axis 3 is not used ... so turn it off: ####################

   P.sca(ax3)
   ax3.axis('off')

########### Add label of max value if showmax is set to True: ################

   if (showmax == 'True'):
      var1_max = str(np.max(var1))
      var1_label = 'Max Val.: ' + var1_max[0:6]
      P.text(0.55, 0.55, var1_label, fontsize=10)

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight')

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

   for coll in p1.collections:
      coll.remove()
   for coll in p2.collections:
      coll.remove()
   for coll in p3.collections:
      coll.remove()
   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()

#######################################################################################
#######################################################################################

#######################################################################################

def plot_warn(map, fig, ax1, ax2, ax3, shapefile, start_time, end_time, svr_color, tor_color, ff='True'):

   #This subroutine will plot NWS warning products from shapefiles downloaded from the 
   #Iowa State Mesonet archive:  https://mesonet.agron.iastate.edu/request/gis/watchwarn.phtml

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #shapefile - string of full path to shapefile to read/plot
   #time - float of time in seconds forecast is valid
   #svr_color - color to plot svr warnings
   #tor_color - color to plot tor warnings
   #ff - boolean for whether or not to plot flash flood warnings instead of svr/tor
   
   #Returns:

   #Either the svr/tor or ff plot objects

#######################

   P.sca(ax1) #set axis to 1

   map.readshapefile(shapefile, 'warnings', drawbounds = False) #read shapefile

########### If plotting svr/tor warnings: ################

#   if (ff == 'False'):
   svr = []
   tor = []
   ff = []

########### Get issuance and expiration times of warnings, convert to seconds: ################
########### If warning is valid at the current forecast time, append plot of polygon to svr/tor variable: ################

   for info, shape in zip(map.warnings_info, map.warnings):
      if (info['PHENOM'] == 'SV'):
         temp_inithr = info['ISSUED'][8:10]
         temp_initmin = info['ISSUED'][10:12]
         init_sec = double(temp_inithr) * 3600. + double(temp_initmin) * 60.

         temp_expirehr = info['EXPIRED'][8:10]
         temp_expiremin = info['EXPIRED'][10:12]
         expire_sec = double(temp_expirehr) * 3600. + double(temp_expiremin) * 60.

         if (init_sec < 35000.):  ########REMOVE THESE LATER - FUDGE FOR TIMES > 00Z
            init_sec = init_sec + 86400.
            expire_sec = expire_sec + 86400.
         elif (expire_sec < 35000.):
            expire_sec = expire_sec + 86400.

         if (start_time <= expire_sec):
            if (end_time >= init_sec):
               wx, wy = zip(*shape)
               svr.append(map.plot(wx, wy, marker=None, color=svr_color, linewidth=1.5, alpha=0.8))
#               poly_out = zip(wx,wy)
#               poly_svr = Polygon(poly_out, facecolor=cb_colors.blue3, edgecolor=cb_colors.blue8, linewidth = 2., alpha = 0.5, zorder = 1)
#               ax1.add_patch(poly_svr)
#                  svr.append(map.plot(wx, wy, marker=None, color=svr_color, linewidth=1.75, alpha=0.6))

      elif (info['PHENOM'] == 'TO'):
         temp_inithr = info['ISSUED'][8:10]
         temp_initmin = info['ISSUED'][10:12]
         init_sec = double(temp_inithr) * 3600. + double(temp_initmin) * 60.

         temp_expirehr = info['EXPIRED'][8:10]
         temp_expiremin = info['EXPIRED'][10:12]
         expire_sec = double(temp_expirehr) * 3600. + double(temp_expiremin) * 60.

         if (init_sec < 35000.):  ########REMOVE THESE LATER - FUDGE FOR TIMES > 00Z
            init_sec = init_sec + 86400.
            expire_sec = expire_sec + 86400.
         elif (expire_sec < 35000.):
            expire_sec = expire_sec + 86400.

         if (start_time <= expire_sec):
            if (end_time >= init_sec):
               wx, wy = zip(*shape)
#               poly_out = zip(wx,wy)
#               poly_tor = Polygon(poly_out, facecolor=cb_colors.red3, edgecolor=cb_colors.red8, linewidth = 2., alpha = 0.5, zorder = 1)
#               ax1.add_patch(poly_tor)
               tor.append(map.plot(wx, wy, marker=None, color=tor_color, linewidth=1.5, alpha=0.8))
#      return svr, tor
#   else:

########### If plotting ff warnings: ################

#      ff = []
#      for info, shape in zip(map.warnings_info, map.warnings):

      elif (info['PHENOM'] == 'FF'):
         temp_inithr = info['ISSUED'][8:10]
         temp_initmin = info['ISSUED'][10:12]
         init_sec = double(temp_inithr) * 3600. + double(temp_initmin) * 60.

         temp_expirehr = info['EXPIRED'][8:10]
         temp_expiremin = info['EXPIRED'][10:12]
         expire_sec = double(temp_expirehr) * 3600. + double(temp_expiremin) * 60.

         if (init_sec < 35000.):  ########REMOVE THESE LATER - FUDGE FOR TIMES > 00Z
            init_sec = init_sec + 86400.
            expire_sec = expire_sec + 86400.
         elif (expire_sec < 35000.):
            expire_sec = expire_sec + 86400.

         if (start_time <= expire_sec):
            if (end_time >= init_sec):
               wx, wy = zip(*shape)
#               poly_out = zip(wx,wy)
#               poly_ff = Polygon(poly_out, facecolor=cb_colors.green3, edgecolor=cb_colors.green8, linewidth = 2., alpha = 0.5, zorder = 1)
#               ax1.add_patch(poly_ff)
               ff.append(map.plot(wx, wy, marker=None, color=cb_colors.green6, linewidth=1.5, alpha=0.8))
   return svr, tor, ff

#######################

#######################################################################################

def plot_lsr(map, fig, ax1, ax2, ax3, shapefile, init_time, time, plot_h='True', plot_w='True', plot_t='True'):

   #This subroutine will plot NWS Local Storm Reports (LSRs) from shapefiles downloaded from the
   #Iowa State Mesonet archive:  https://mesonet.agron.iastate.edu/request/gis/lsrs.phtml

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #shapefile - string of full path to shapefile to read/plot
   #init_time - float of forecast initialization time in seconds
   #time - float of time in seconds forecast is valid
   #plot_h - boolean on whether to plot hail reports
   #plot_w - boolean on whether to plot wind reports
   #plot_t - boolean on whether to plot tornado reports

   #Returns:

   #hail, wind, and tornado variables - lists of scatter plots of LSRs

#######################

   P.sca(ax1) #set current axis to 1

   map.readshapefile(shapefile, 'lsr', drawbounds = False) #read shapefile

   hail = []
   wind = []
   tornado = [] 

########### For each storm report in shapefile, scatter with the appropriate symbol if  ################
########### report occurs between forecast initialization and valid time:   ################
   
   for info, shape in zip(map.lsr_info, map.lsr):
      if ((info['TYPECODE'] == 'H') and (plot_h == 'True')):
         temp_inithr = info['VALID'][8:10]
         temp_initmin = info['VALID'][10:12]
         init_sec = double(temp_inithr) * 3600. + double(temp_initmin) * 60.
         if (init_sec < 25000.):  ########REMOVE THESE LATER - FUDGE FOR TIMES > 00Z
            init_sec = init_sec + 86400.

         if ((init_sec >= init_time) and ((time+150.) >= init_sec)):
            temp_lat = double(info['LAT'])
            temp_lon = double(info['LON'])
            wx, wy = map(temp_lon, temp_lat)
            hail.append(map.scatter(wx, wy, s=25, marker='o', facecolor=cb_colors.green6, edgecolor='k', linewidth=0.4, zorder=10))

      elif ((info['TYPECODE'] == 'D') and (plot_w == 'True')):
         temp_inithr = info['VALID'][8:10]
         temp_initmin = info['VALID'][10:12]
         init_sec = double(temp_inithr) * 3600. + double(temp_initmin) * 60.
         if (init_sec < 25000.):
            init_sec = init_sec + 86400.

         if ((init_sec >= init_time) and ((time+150.) >= init_sec)):
            temp_lat = double(info['LAT'])
            temp_lon = double(info['LON'])
            wx, wy = map(temp_lon, temp_lat)
            wind.append(map.scatter(wx, wy, s=25, marker='s', facecolor=cb_colors.blue6, edgecolor='k', linewidth=0.4, zorder=10))

      if ((info['TYPECODE'] == 'T') and (plot_t == 'True')):
         temp_inithr = info['VALID'][8:10]
         temp_initmin = info['VALID'][10:12]
         init_sec = double(temp_inithr) * 3600. + double(temp_initmin) * 60.
         if (init_sec < 25000.):
            init_sec = init_sec + 86400.

         if ((init_sec >= init_time) and ((time+150.) >= init_sec)):
            temp_lat = double(info['LAT'])
            temp_lon = double(info['LON'])
            wx, wy = map(temp_lon, temp_lat)
            tornado.append(map.scatter(wx, wy, s=25, marker='v', facecolor=cb_colors.red6, edgecolor='k', linewidth=0.4, zorder=10))

   return hail, wind, tornado

#######################

#######################################################################################

def remove_warn(var):

   #This subroutine will remove NWS warning plots from a basemap 

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot, plot_warn

   #Input:

   #var - List of Basemap warning plots to remove

   #Returns:

   #NOTHING!

#######################

   for i in range(0, len(var)):
      var[i][0].remove()

#######################

#######################################################################################

def remove_lsr(hail, wind, tornado):

   #This subroutine will remove NWS LSR plots from a basemap 

   #Dependencies:
   
   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot, plot_LSR

   #Input:
   
   #hail - List of Basemap hail LSR plots to remove
   #wind - List of Basemap wind LSR plots to remove
   #tornado - List of Basemap tornado LSR plots to remove

   #Returns:
   
   #NOTHING!

#######################

   for i in range(0, len(hail)):
      hail[i].remove()
   for i in range(0, len(wind)):
      wind[i].remove()
   for i in range(0, len(tornado)):
      tornado[i].remove()

#######################

#######################################################################################

def pmm_rectangle(map, fig, x, y, plot):

   #NO LONGER USED!
   #This subroutine will draw a rectangle around the domain where probability-matched
   #mean values are available    
   
   #Dependencies:
   
   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot
   
   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted
   #y - 2d numpy array of y axis, basemap-relative locations of data to be plotted
   #plot - web_plot object to be plotted
   
   #Returns:
   
   #r1, r2, r3, r4:  plot variables for each side of the rectangle

#######################

   lx_rect = x[(plot.neighborhood):-(plot.neighborhood-1),(plot.neighborhood)] - 1500.
   rx_rect = x[(plot.neighborhood+1):-(plot.neighborhood-1),-(plot.neighborhood)] - 1500.
   upx_rect = x[-(plot.neighborhood),(plot.neighborhood):-(plot.neighborhood-1)] - 1500.
   lowx_rect = x[(plot.neighborhood),(plot.neighborhood+1):-(plot.neighborhood-1)] - 1500.

   ly_rect = y[(plot.neighborhood):-(plot.neighborhood-1),(plot.neighborhood)] - 1500.
   ry_rect = y[(plot.neighborhood+1):-(plot.neighborhood-1),-(plot.neighborhood)] - 1500.
   upy_rect = y[-(plot.neighborhood),(plot.neighborhood):-(plot.neighborhood-1)] - 1500.
   lowy_rect = y[(plot.neighborhood),(plot.neighborhood+1):-(plot.neighborhood-1)] - 1500.

   r1 = P.plot(lx_rect,ly_rect, linewidth=2.5, color='k', alpha=0.1)
   r2 = P.plot(rx_rect,ry_rect, linewidth=2.5, color='k', alpha=0.1)
   r3 = P.plot(upx_rect,upy_rect, linewidth=2.5, color='k', alpha=0.1)
   r4 = P.plot(lowx_rect,lowy_rect, linewidth=2.5, color='k', alpha=0.1)

   return r1, r2, r3, r4

#######################

#######################################################################################

def rem_pmm_rect(r1, r2, r3, r4):

   #NO LONGER USED!
   #This subroutine will remove a rectangle around the domain where probability-matched
   #mean values are available

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot

   #Input:
 
   #r1 - r4:  Plot variables for each side of the rectangle

   #Returns:

   #NOTHING!

   #Dependencies:

   #basemap, matplotlib, pyplot
   #mymap

   #Input:

   #Returns:
   #NOTHING!

#######################

   for i in range(0, len(r1)):
      r1[i].remove()
   for i in range(0, len(r2)):
      r2[i].remove()
   for i in range(0, len(r3)):
      r3[i].remove()
   for i in range(0, len(r4)):
      r4[i].remove()

#######################################################################################
def sat_plot(map, fig, ax1, ax2, ax3, x, y, plot, var1, var2, t, init_label, valid_label, domain, outdir, q_u, q_v, scale, t_star, t_plus, spec, quiv, showmax='False'):

############# Calculate appropriate index for filename (e.g. 10-min forecast with 5-min output ###############
############# would be:  t_star = 5, t_plus = 0, t = 2, outtime = 'f010'                       ###############

   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   elif (len(temp_outtime) == 3):
      outtime = 'f' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

   P.sca(ax1) #set plotting axis to 1

############# Create filled contour plot, with intervals highlighted for var1 data:    #################
   cmap = plot.cmap
   bounds = plot.var1_levels
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

   p1 = map.contourf(x, y, var1, cmap=plot.cmap, norm=norm, levels = plot.var1_levels, alpha=plot.alpha, extend=plot.extend)
#   p2 = map.contour(x, y, var1, colors='k', levels = plot.var1_levels, linewidths=0.4, alpha=0.25, linestyles='solid')

   p1.cmap.set_over(plot.over_color)
   p1.cmap.set_under(plot.under_color)


############ Add text labels of plot titles and initialization and valid forecast time: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold')
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)

########### Switch to axis 2 and set axis formatting: ################

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Plot and format colormap for filled contour plot in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha)
   cbar.set_label(plot.var1_title, fontsize=10)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=10)
   if (cmap.N < 30 ):
       for label in cbar.ax.xaxis.get_ticklabels()[1::2]: #only label every th contour interval 
          label.set_visible(False)
   if (cmap.N >= 30 ):
       for label in cbar.ax.xaxis.get_ticklabels()[1::2]:
          label.set_visible(False)

########### Switch to axis 3 and set axis to be blank: ################

   P.sca(ax3)
   ax3.axis('off')

   if (showmax == 'True'):
      var1_max = str(np.max(var1))
      var1_label = 'Max Val.: ' + var1_max[0:6]
      P.text(0.55, 0.55, var1_label, fontsize=10)

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight')

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

   for coll in p1.collections:
      coll.remove()
#   for coll in p2.collections:
#      coll.remove()
#   for coll in p3.collections:
#      coll.remove()
   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()

#   rem_pmm_rect(r1, r2, r3, r4) #remove pmm rectangle plot - deprecated

   if (showmax == 'True'):
      p4.remove()

   if (spec == 'True'):
      for coll in p4.collections:
         coll.remove()
   if (quiv == 'True'):
      q1.remove()

#######################################################################################

def env_plot_rain(map, fig, ax1, ax2, ax3, x, y, plot, var1, var2, t, init_label, valid_label, domain, outdir, q_u, q_v, scale, t_star, t_plus, spec, quiv, showmax='False'):

   #This subroutine is used to create the majority of the NEWS-e plots used on the website
   #(all except the paintball, object-matching, and MRMS verification plots).  If passed an
   #empty basemap and figure object as well as the data and specifications for the plot, it
   #will create a plot over the model domain with appropriate labels, save the image us a 
   #.png using the appropriate timestamp, then remove the plots from the underlying basemap
   #template 

   #Dependencies:

   #python:  basemap, matplotlib, pyplot
   #plotting_cbook:  mymap, create_fig, web_plot

   #Input:

   #map - Basemap instance for plotting (from mymap)
   #fig - Figure instance for plotting (from create_fig)
   #ax1 - Axis 1 instance for plotting (from create_fig)
   #ax2 - Axis 2 instance for plotting (from create_fig)
   #ax3 - Axis 3 instance for plotting (from create_fig)
   #x - 2d numpy array of x-axis, basemap-relative locations of data to be plotted 
   #y - 2d numpy array of y axis, basemap-relative locations of data to be plotted 
   #plot - web_plot object to be plotted
   #var1 - 2d numpy array of the data to be plotted as a filled contour
   #var2 - 2d numpy array of the data to be plotted as a contour
   #t - time index associated with the plot (used for website indexing)
   #init_label - string of forecast initialization time
   #valid_label - string of time forecast valid
   #domain - string specifying subdomain over which data is plotted (currently always 'full') 
   #outdir - string of directory path to save .png image in
   #q_u - if quiv=True, u values to produce quiver plot
   #q_v - if quiv=True, v values to produce quiver plot
   #scale - if quiv=True, scaling factor for quiver vectors
   #t_star - multiplication factor for time index (used for website indexing)
   #t_plus - addition factor for time index (used for website indexing)
   #spec - boolean specifying if special contour level is to be plotted (e.g. 60 deg. dewpoint)
   #quiv - boolean specifying if quiver plot is to be overlain
   #showmax - boolean specifying if scatter of domain max value and associated label are to be plotted  

   #Returns:

   #Nothing, but all plots are removed from base template of the basemap after plot is saved

#######################

############# Calculate appropriate index for filename (e.g. 10-min forecast with 5-min output ###############
############# would be:  t_star = 5, t_plus = 0, t = 2, outtime = 'f010'                       ###############

   temp_outtime = str(t * t_star + t_plus)
   if (len(temp_outtime) == 1):
      outtime = 'f00' + temp_outtime
   elif (len(temp_outtime) == 3):
      outtime = 'f' + temp_outtime
   else:
      outtime = 'f0' + temp_outtime

   P.sca(ax1) #set plotting axis to 1

#   r1, r2, r3, r4 = pmm_rectangle(map, fig, x, y, plot) #plot pmm rectangle - deprecated

############# Create filled contour plot, with intervals highlighted for var1 data:    #################
   cmap = plot.cmap
   bounds = plot.var1_levels
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

   temp_cmap = plot.cmap
   temp_cmap.set_under(plot.under_color)
   temp_cmap.set_over(plot.over_color)

   p1 = map.contourf(x, y, var1, cmap=temp_cmap, levels= plot.var1_levels, norm=norm, alpha=plot.alpha, extend=plot.extend)
#   p1 = map.contourf(x, y, var1, cmap=temp_cmap, levels = plot.var1_levels, alpha=plot.alpha, extend=plot.extend)
#   p1 = map.contourf(x, y, var1, cmap=plot.cmap, levels = plot.var1_levels, alpha=plot.alpha, extend=plot.extend)
   p2 = map.contour(x, y, var1, colors='k', levels = plot.var1_levels, linewidths=0.4, alpha=0.25, linestyles='solid')

#   p1.cmap.set_over(plot.over_color)
#   p1.cmap.set_under(plot.under_color)

############# Create contour plot for var2 data:  #################

   p3 = map.contour(x, y, var2, colors=plot.var2_colors, levels=plot.var2_levels, linewidths=[1., 1.5])

############# If special contour is to be plotted, plot it: ###############

   if (spec == 'True'):
      p4 = map.contour(x, y, var1, colors=plot.var1_cont_colors, levels=plot.var1_cont_levels, linewidths=1.25)

############ If quiver plot is to be included, include it: #################

   if (quiv == 'True'):
      quiv_x = x[0:-1:13,0:-1:13]
      quiv_y = y[0:-1:13,0:-1:13]
#      q1 = map.quiver(quiv_x, quiv_y, q_u, q_v, color='k', linewidth=0.3, minshaft=2.5, scale=scale, edgecolor='none', pivot='tail', alpha=0.35)
      barbs1, barbs2 = map.barbs(quiv_x, quiv_y, q_u, q_v, length=5, color='k', alpha=0.35, barb_increments=dict(half=5., full=10., flag=50.), fill_empty="True")

############ Add text labels of plot titles and initialization and valid forecast time: ###############

   t1 = fig.text(0.03, 0.965, plot.var1_title, fontsize=10, fontweight='bold')
   t2 = fig.text(0.03, 0.945, plot.var2_title, fontsize=10, fontweight='bold', color=plot.var2_tcolor)
   t3 = fig.text(0.74, 0.965, init_label, fontsize=9)
   t4 = fig.text(0.74, 0.945, valid_label, fontsize=9)

############ If domain max value is to be plotted, plot it: #################

   if (showmax == 'True'):
      mx, my = np.unravel_index(var1.argmax(), var1.shape)
      p4 = map.scatter(x[mx,my], y[mx,my], s=30, linewidth=0.85, marker='+', color='k', alpha=0.8)

########### Switch to axis 2 and set axis formatting: ################

   P.sca(ax2)
   ax2.spines["top"].set_alpha(0.9)
   ax2.spines["right"].set_alpha(0.9)
   ax2.spines["top"].set_linewidth(0.5)
   ax2.spines["right"].set_linewidth(0.5)
   ax2.spines["bottom"].set_alpha(0.9)
   ax2.spines["left"].set_alpha(0.9)
   ax2.spines["bottom"].set_linewidth(0.5)
   ax2.spines["left"].set_linewidth(0.5)

########## Get colormap information: #####################

   cmap = plot.cmap
   bounds = plot.var1_levels
   norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

########## Plot and format colormap for filled contour plot in axis 2: #####################

   cbar = matplotlib.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, ticks=bounds, orientation = 'horizontal', extend=plot.extend, alpha=plot.alpha, spacing='uniform')
   cbar.set_label(plot.var1_title, fontsize=8)
   cl = P.getp(cbar.ax, 'xmajorticklabels')
   P.setp(cl, fontsize=8)
   for label in cbar.ax.xaxis.get_ticklabels()[1::2]: #label every other contour interval 
       label.set_visible(True)

########### Switch to axis 3 and set axis to be blank: ################

   P.sca(ax3)
   ax3.axis('off')

########### Plot legend and description of contour plots (and max value if needed): ################

   for i in range(0, len(plot.var2_levels)):
      temp_label = str(plot.var2_levels[i])
      label = temp_label[0:2] + ' dBZ'

      P.axhline(y=0.75, xmin=(0.04 + (0.23 * i)), xmax=(0.17 + (0.23 * i)), linewidth=(1. + i), color=plot.var2_colors[i])
      P.text((0.025 + (0.235 * i)), 0.3, label, fontsize=10, color=plot.var2_colors[i])

   if (showmax == 'True'):
      var1_max = str(np.max(var1))
      var1_label = 'Max Val.: ' + var1_max[0:6]
      P.text(0.55, 0.55, var1_label, fontsize=10)

########### Save figure as .png, then clear all plots from basemap template: ################

   fig_name = outdir + plot.name + '_' + outtime + '.png'
   P.savefig(fig_name, format="png", bbox_inches='tight')

   P.cla()
   P.sca(ax2)
   P.cla()
   P.sca(ax1)

   for coll in p1.collections:
      coll.remove()
   for coll in p2.collections:
      coll.remove()
   for coll in p3.collections:
      coll.remove()
   t1.remove()
   t2.remove()
   t3.remove()
   t4.remove()

#   rem_pmm_rect(r1, r2, r3, r4) #remove pmm rectangle plot - deprecated

   if (showmax == 'True'):
      p4.remove()

   if (spec == 'True'):
      for coll in p4.collections:
         coll.remove()
   if (quiv == 'True'):
      barbs1.remove()
      barbs2.remove()
#      q1.remove()

