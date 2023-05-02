## Copyright (C) GNSS ACADEMY 
##
## Name          : SatFunctions.py
## Purpose       : Satellite Analyses functions
## Project       : WP0-JSNP
## Component     : 
## Author        : GNSS Academy
## Creation date : 2021
## File Version  : 1.0
## Version date  : 
##

import sys, os
from pandas import unique
from interfaces import LOS_IDX
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from interfaces import POS_IDX
sys.path.append(os.getcwd() + '/' + \
    os.path.dirname(sys.argv[0]) + '/' + 'COMMON')
from COMMON import GnssConstants
from COMMON.Plots import generatePlot
import numpy as np
from math import pi
# from pyproj import Transformer
from COMMON.Coordinates import xyz2llh


# Plot Satellite Visibility Figures
def plotSatVisibility(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    PlotConf["Title"] = "Satellite Visibility from TLSA on Year 2015"\
        " DoY 006"

    PlotConf["yLabel"] = "GPS-PRN"
    PlotConf["yTicks"] = sorted(unique(LosData[LOS_IDX["PRN"]]))
    PlotConf["yTicksLabels"] = sorted(unique(LosData[LOS_IDX["PRN"]]))
    PlotConf["yLim"] = [0, max(unique(LosData[LOS_IDX["PRN"]])) + 1]

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    for prn in sorted(unique(LosData[LOS_IDX["PRN"]])):
        Label = "G" + ("%02d" % prn)
        FilterCond = LosData[LOS_IDX["PRN"]] == prn
        PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]][FilterCond] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = LosData[LOS_IDX["PRN"]][FilterCond]
        PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]][FilterCond]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/SAT/' + 'SAT_VISIBILITY_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Satellite Geometrical Range Figures
def plotSatGeomRnge(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Satellite Geometical Range from TLSA on Year 2015"\
        " DoY 006"

    PlotConf["yLabel"] = "Range [km]"

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["RANGE[m]"]]/1000
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/SAT/' + 'SAT_GEOMETRICAL_RANGE_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot Satellite Tracks Figures
def plotSatTracks(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (16.8,15.2)
    PlotConf["Title"] = "Satellite Tracks during visibility periods from "\
        "TLSA on Year 2015 DoY 006"

    PlotConf["LonMin"] = -135
    PlotConf["LonMax"] = 135
    PlotConf["LatMin"] = -35
    PlotConf["LatMax"] = 90
    PlotConf["LonStep"] = 15
    PlotConf["LatStep"] = 10

    # PlotConf["yLabel"] = "Latitude [deg]"
    PlotConf["yTicks"] = range(PlotConf["LatMin"],PlotConf["LatMax"]+1,10)
    PlotConf["yLim"] = [PlotConf["LatMin"], PlotConf["LatMax"]]

    # PlotConf["xLabel"] = "Longitude [deg]"
    PlotConf["xTicks"] = range(PlotConf["LonMin"],PlotConf["LonMax"]+1,15)
    PlotConf["xLim"] = [PlotConf["LonMin"], PlotConf["LonMax"]]

    PlotConf["Grid"] = True

    PlotConf["Map"] = True

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    # Transform ECEF to Geodetic
    LosData[LOS_IDX["SAT-X[m]"]].to_numpy()
    LosData[LOS_IDX["SAT-Y[m]"]].to_numpy()
    LosData[LOS_IDX["SAT-Z[m]"]].to_numpy()
    DataLen = len(LosData[LOS_IDX["SAT-X[m]"]])
    Longitude = np.zeros(DataLen)
    Latitude = np.zeros(DataLen)
    # transformer = Transformer.from_crs('epsg:4978', 'epsg:4326')
    for index in range(DataLen):
        x = LosData[LOS_IDX["SAT-X[m]"]][index]
        y = LosData[LOS_IDX["SAT-Y[m]"]][index]
        z = LosData[LOS_IDX["SAT-Z[m]"]][index]
        Longitude[index], Latitude[index], h = xyz2llh(x, y, z)
        # Latitude[index], Longitude[index], h = transformer.transform(x, y, z)

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0
    PlotConf["xData"][Label] = Longitude
    PlotConf["yData"][Label] = Latitude
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/SAT/' + 'SAT_TRACKS_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Satellite Range Velocity Figures
def plotSatVelocity(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Satellite Range Velocity from TLSA on Year 2015"\
        " DoY 006"
    PlotConf["yLabel"] = "Absolute Velocity [km/s]"
    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = np.sqrt(LosData[LOS_IDX["VEL-X[m/s]"]]**2 + LosData[LOS_IDX["VEL-Y[m/s]"]]**2 + LosData[LOS_IDX["VEL-Z[m/s]"]]**2)/1000
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]
    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/SAT/' + 'SAT_VELOCITY_TLSA_D006Y15.png'
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

#Plot Satellite NAV Clock Figures
def plotPRN19NavSatClock(LosData):
    PlotConf = {}
    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,8.8)
    PlotConf["Title"] = "PRN19 NAV CLK from TLSA on Year 2015"\
        " DoY 006"
    PlotConf["yLabel"] = "CLK[km]"

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 11)

    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '|'
    PlotConf["LineWidth"] = 1

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = []
    for prn in sorted(unique(LosData[LOS_IDX["PRN"]])):
        Label = "G" + ("%02d" % prn)
        FilterCond = LosData[LOS_IDX["PRN"]] == 19
        PlotConf["xData"][Label] = ((LosData[LOS_IDX["SOD"]])[FilterCond]) / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = ((LosData[LOS_IDX["SV-CLK[m]"]])[FilterCond])/1000
        PlotConf["zData"][Label] = (LosData[LOS_IDX["PRN"]])[FilterCond]
        PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/SAT/' + 'SAT_CLK_TLSA_D006Y15_PRN19.png' 
    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Satellite Clock Corrected Figures
def plotSatClockCorrected(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    PlotConf["Title"] = "Satellite CLK + DTR - TGD from TLSA on Year 2015"\
        " DoY 006"

    PlotConf["yLabel"] = "CLK[km]"
    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '_'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"
    PlotConf["ColorBarMin"] = 1
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarTicks"] = sorted(unique(LosData[LOS_IDX["PRN"]]))

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    #Construction of the Clock Correction Array for each Epoch.This vector contain the clock correction according the equation Corrected-Clock =CLK + DTR - TGD
    clk_corr_array = (LosData[LOS_IDX["SV-CLK[m]"]] + LosData[LOS_IDX["DTR[m]"]] - LosData[LOS_IDX["TGD[m]"]]) / 1000 
    #This LOOP is to creat a Y-DATA that contain the clock corrected(including all PRNs)
    #CDATA is for Time
    #ZDATA for the Color bar(according to the PRN)
    Label = []
    for prn in sorted(unique(LosData[LOS_IDX["PRN"]])):
        Label =("%02d" % prn)
        FilterCond = LosData[LOS_IDX["PRN"]]== prn
        PlotConf["xData"][Label] = ((LosData[LOS_IDX["SOD"]]))[FilterCond] / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = clk_corr_array[FilterCond]
        PlotConf["zData"][Label] = (LosData[LOS_IDX["PRN"]])[FilterCond]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/SAT/' + 'SAT_CLK_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Satellite Clock Total TGD Figures
def plotSatTGD(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    PlotConf["Title"] = "Satellite TGD (Total Group Delay) from TLSA on Year 2015"\
        " DoY 006"

    PlotConf["yLabel"] = "TGD[m]"
   

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '_'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Satellite PRN"
    PlotConf["ColorBarMin"] = 1
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarTicks"] = sorted(unique(LosData[LOS_IDX["PRN"]]))

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = []
    for prn in sorted(unique(LosData[LOS_IDX["PRN"]])):
        Label =("%02d" % prn)
        FilterCond = LosData[LOS_IDX["PRN"]]== prn 
        PlotConf["xData"][Label] = ((LosData[LOS_IDX["SOD"]])[FilterCond]) / GnssConstants.S_IN_H
        PlotConf["yData"][Label] = ((LosData[LOS_IDX["TGD[m]"]])[FilterCond])
        PlotConf["zData"][Label] = (LosData[LOS_IDX["PRN"]])[FilterCond]

        PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/SAT/' + 'SAT_TGD_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot Satellite clock relativistic effect (DTR) Figures
def plotSatDTR(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Satellite DTR (Clock Relativistic Effect) from TLSA on Year 2015"\
        " DoY 006"
    PlotConf["yLabel"] = "DTR[m]"
    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0

    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["DTR[m]"]]
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/SAT/' + 'SAT_DTR_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    # Plot Slant Ionospheric Delays Figures
def plotIonoDelays(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Ionospheric Klobuchar Delays (STEC) from TLSA on Year 2015"\
        " DoY 006"
    PlotConf["yLabel"] = "STEC[m]"
    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.5
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    Label = 0

    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["STEC[m]"]]
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/ION/' + 'IONO_STEC_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    # Plot Satellite visibility periods Figures
def plotSatVisibilityPeriods(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Satellite Visibility vs STEC from TLSA on Year 2015"\
        " DoY 006"
    PlotConf["yLabel"] = "GPS-PRN"
    PlotConf["yLim"] = [0,max(unique(LosData[LOS_IDX["PRN"]]))+1]
    PlotConf["yTicks"] = sorted(unique(LosData[LOS_IDX["PRN"]]))

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "STEC[m]"
    PlotConf["ColorBarMin"] = min(unique(LosData[LOS_IDX["STEC[m]"]]))
    PlotConf["ColorBarMax"] = max(unique(LosData[LOS_IDX["STEC[m]"]]))

    
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
    

    Label = []
    for prn in sorted(unique(LosData[LOS_IDX["PRN"]])):
        Label ="G" + ("%02d" % prn)
        FilterCond = LosData[LOS_IDX["PRN"]]== prn 
        PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]][FilterCond]/ GnssConstants.S_IN_H
        PlotConf["yData"][Label] = LosData[LOS_IDX["PRN"]][FilterCond]
        PlotConf["zData"][Label] = LosData[LOS_IDX["STEC[m]"]][FilterCond]
        PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/ION/' + 'IONO_STEC_vs_PRN_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Plot VTEC vs Time Figures
def plotVTECvsTime(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Ionospheric Klobuchar Delays (VTEC) from TLSA on Year 2015"\
        " DoY 006"
    PlotConf["yLabel"] = "VTEC[m]"
    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.25
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
 
    Label = 0
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["VTEC[m]"]]
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/ION/' + 'IONO_VTEC_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

#Plot satellite visibility periods Figures
def PlotPRNvsTIME(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Satellite Visibility vs VTEC from TLSA on Year 2015"\
        " DoY 006"
    PlotConf["yLabel"] = "GPS-PRN"
    PlotConf["yLim"] = [0,max(unique(LosData[LOS_IDX["PRN"]]))+1]
    PlotConf["yTicks"] = sorted(unique(LosData[LOS_IDX["PRN"]]))
    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.5

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "VTEC[m]"
    PlotConf["ColorBarMin"] = min(unique(LosData[LOS_IDX["VTEC[m]"]]))
    PlotConf["ColorBarMax"] = max(unique(LosData[LOS_IDX["VTEC[m]"]]))

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = []
    for prn in sorted(unique(LosData[LOS_IDX["PRN"]])):
        Label =("%02d" % prn)
        FilterCond = LosData[LOS_IDX["PRN"]]== prn 
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["PRN"]]
    PlotConf["zData"][Label] = LosData[LOS_IDX["VTEC[m]"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/ION/' + 'IONO_VTEC_vs_PRN_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)


# Plot STD (Slant Tropospheric Delay) vs Time(Elevation) Figures
def plotSTDvsTIME(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Slant Tropospheric Delays (STD) from TLSA on Year 2015"\
        " DoY 006"
    PlotConf["yLabel"] = "STD[m]"
    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.25
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
 
    Label = 0
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["TROPO[m]"]]
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/TRO/' + 'TROPO_STD_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

# Compute and plot the Zenith Tropo Delay (ZTD) vs Time Figures
def plotZTDvsTIME(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Zenith Tropospheric Delays (ZTD) from TLSA on Year 2015"\
        " DoY 006"
    PlotConf["yLabel"] = "ZTD[m]"
    PlotConf["yLim"] = [2.3150, 2.3180]

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.25
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    Label = 0
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["TROPO[m]"]]/(1.001/np.sqrt(0.002001+np.sin(LosData[LOS_IDX["ELEV"]]*pi/180)**2))
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/TRO/' + 'TROPO_ZTD_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    # Plot Pseudo-ranges (Code Measurements C1) Figures
def plotPSRvsTime(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = " Pseudo-range C1C from TLSA "
    PlotConf["yLabel"] = "pseudo-range[km]"
    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.25
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
 
    Label = 0
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["MEAS[m]"]] /1000
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/MSR/' + 'MEAS_CODES_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

        # Plot Tau = C1C/c Figures
def plotTAUvsTime(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = " Tau=Rho/c from TLSA on Year 2015 "\
        " DoY 006"
    PlotConf["yLabel"] = "Tau[ms]"
   

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.5
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
   
    Label = 0
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["MEAS[m]"]]*(1/299792458*1000) 
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/MSR/' + 'TAU_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    # Plot Time of Flight (ToF) Figures
def plotTOFvsTime(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = " Time of Flight (ToF) from TLSA on Year 2015 "\
        " DoY 006"
    PlotConf["yLabel"] = "ToF[ms]"

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.25
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
 
    Label = 0
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = LosData[LOS_IDX["TOF[ms]"]] 
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/MSR/' + 'TOF_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    #Plot Doppler Frequency Figures
def plotDOPPLER_FREQ(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = " Doppler Frequency from TLSA on Year 2015 "\
        " DoY 006"
    PlotConf["yLabel"] = "Doppler Frequency [KHz]"

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    PlotConf["Grid"] = 1
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 0.25
    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "Elevation [deg]"
    PlotConf["ColorBarMin"] = 0.
    PlotConf["ColorBarMax"] = 90.

    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}

    c = 299792458
    f_L1 = 1575.42
    #R_rcv = [4637952.29175333, 121207.93198126, 4362375.76379941] 
    Rx_los = LosData[LOS_IDX["SAT-X[m]"]] - 4637952.29175333
    Ry_los = LosData[LOS_IDX["SAT-Y[m]"]] - 121207.93198126
    Rz_los = LosData[LOS_IDX["SAT-Z[m]"]] - 4362375.76379941
    norm_Rlos = np.sqrt((LosData[LOS_IDX["SAT-X[m]"]] - 4637952.29175333)**2 + (LosData[LOS_IDX["SAT-Y[m]"]] - 121207.93198126 )**2 + (LosData[LOS_IDX["SAT-Z[m]"]] - 4362375.76379941)**2)
    ux_los = np.divide(Rx_los, norm_Rlos)
    uy_los = np.divide(Ry_los, norm_Rlos)
    uz_los = np.divide(Rz_los, norm_Rlos)
    #R_los = np.substract(R_rcv,R_sat)
    v_los = np.multiply(ux_los, LosData[LOS_IDX["VEL-X[m/s]"]]) + np.multiply(uy_los, LosData[LOS_IDX["VEL-Y[m/s]"]]) + np.multiply(uz_los, LosData[LOS_IDX["VEL-Z[m/s]"]])
    
    Label = 0
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]]/ GnssConstants.S_IN_H
    PlotConf["yData"][Label] = - (v_los/ c)*(f_L1 * 1000)
    PlotConf["zData"][Label] = LosData[LOS_IDX["ELEV"]]

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/MSR/' + 'DOPPLER_FREQ_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)     

    # Plot Residuals C1 Figures 
def plotRESC1(LosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,6.6)
    PlotConf["Title"] = "Residuals C1C vs Time for TLSA"
    
    

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    
    PlotConf["yLabel"] = "Residuals [km]"
    PlotConf["yLim"] = [2522.250, 2522.280]

    PlotConf["ColorBarMin"] = 1
    PlotConf["ColorBarMax"] = 32
    PlotConf["ColorBarTicks"] = sorted(unique(LosData[LOS_IDX["PRN"]]))
    
    PlotConf["Grid"] = 1
    
    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "GPS-PRN"   
         
        
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
      
    Label = 0
    
    CLKp1 = LosData[LOS_IDX["SV-CLK[m]"]]+LosData[LOS_IDX["DTR[m]"]]-LosData[LOS_IDX["TGD[m]"]]
    
    RESc1 = LosData[LOS_IDX["MEAS[m]"]]  - (LosData[LOS_IDX["RANGE[m]"]] - CLKp1 + LosData[LOS_IDX["STEC[m]"]] + LosData[LOS_IDX["TROPO[m]"]]) 
    
    PlotConf["xData"][Label] = LosData[LOS_IDX["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = RESc1 / 1000
    PlotConf["zData"][Label] = LosData[LOS_IDX["PRN"]]
    

    PlotConf["Path"] = sys.argv[1] + '/OUT/LOS/MSR/' + 'MEAS_RESIDUALS_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    #Plot instantaneous number of satellites Figures
def plotPOS_SATS(PosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Number of Satellites in PVT vs Time for TLSA on Year 2015"\
        " DoY 006"

    PlotConf["yLabel"] = "Number of Satellites"
    PlotConf["yTicks"] = range (0, max(PosData[POS_IDX["NSATS"]]+1))
    PlotConf["yLim"] = [0, max(PosData[POS_IDX["NSATS"]]+1)]

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]
    
    

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = ''
    PlotConf["LineWidth"] = 1.5
    
      
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
   
    Label = 0
    PlotConf["xData"][Label] = PosData[POS_IDX["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"][Label] = PosData[POS_IDX["NSATS"]]
   
    PlotConf["Path"] = sys.argv[1] + '/OUT/POS/POS/' + 'POS_SATS_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    #Plot the PDOP, GDOP, TDOP Figures
def plotGPTDOP(PosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Dilution Of Precision (DOP) from TLSA on Year 2015"\
        " DoY 006"

    PlotConf["yLabel"] = "DOP"
    PlotConf["yLim"] = [0.0, 4.5]

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = ''
    PlotConf["LineWidth"] = 1.5
    
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["Color"] = {}
    PlotConf["Legend"] = {}
    
    Label = ["GDOP","PDOP","TDOP"]
    Color = ["purple","green","skyblue"]

    for index, label in enumerate(Label):
        PlotConf["xData"][label] = PosData[POS_IDX["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][label] = PosData[POS_IDX[label]]
        PlotConf["Color"][label] = Color[index]
          
    PlotConf["yData"].keys()
    PlotConf["Legend"].keys()
    PlotConf["Color"].keys()
        
    PlotConf["Path"] = sys.argv[1] + '/OUT/POS/POS/' + 'POS_DOP_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    # Plot POS_HVDOP vs TIME Figures
def plotPOSHVDOP(PosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "Dilution Of Precision (DOP) from TLSA on Year 2015"\
        " DoY 006"


    PlotConf["yLabel"] = "DOP , Number of Satellites"


    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]


    PlotConf["Grid"] = 1

    PlotConf["Marker"] = ''
    PlotConf["LineWidth"] = 1.5
    
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["Color"] = {}
    PlotConf["Legend"] = {}
    
    Label = ["HDOP", "VDOP", "NSATS"]
    Color = ["purple","green", "yellow"]

    for index, label in enumerate(Label):
        PlotConf["xData"][label] = PosData[POS_IDX["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][label] = PosData[POS_IDX[label]]
        PlotConf["Color"][label] = Color[index]
          
        PlotConf["yData"].keys()
        PlotConf["Legend"].keys()
        PlotConf["Color"].keys()
    
    PlotConf["Path"] = sys.argv[1] + '/OUT/POS/POS/' + 'POS_HVDOP_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    # Plot ENU_PE vs TIME Figures
def plotENUPETIME(PosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "ENU Position Error from TLSA on Year 2015"\
        " DoY 006"

    PlotConf["yLabel"] = "ENU_PE[m]"

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = ''
    PlotConf["LineWidth"] = 1.5
    
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["Color"] = {}
    PlotConf["Legend"] = {}
    
    Label = ["UPE[m]","EPE[m]","NPE[m]",]
    Color = ["purple","green","skyblue",]

    for index, label in enumerate(Label):
        PlotConf["xData"][label] = PosData[POS_IDX["SOD"]] / GnssConstants.S_IN_H
        PlotConf["yData"][label] = PosData[POS_IDX[label]]
        PlotConf["Color"][label] = Color[index]
          
    PlotConf["yData"].keys()
    PlotConf["Legend"].keys()
    PlotConf["Color"].keys()
    
    PlotConf["Path"] = sys.argv[1] + '/OUT/POS/POS/' + 'POS_ENU_PE_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    # Plot the Horizontal(HPE) and Vertical Position Error(VPE) Figures
def plotHVPETIME(PosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "HPE-VPE Position Error from TLSA on Year 2015"\
        " DoY 006"

    
    PlotConf["yLabel"] = "H/VPE[m]"
    PlotConf["yLim"] = [0, 8]

    PlotConf["xLabel"] = "Hour of DoY 006"
    PlotConf["xTicks"] = range(0, 25)
    PlotConf["xLim"] = [0, 24]

    PlotConf["Grid"] = 1

    PlotConf["Marker"] = ''
    PlotConf["LineWidth"] = 1.5
    
    
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["Color"] = {}
    PlotConf["Legend"] = {}

    HPE = np.sqrt(PosData[POS_IDX["EPE[m]"]]**2+PosData[POS_IDX["NPE[m]"]]**2)
    VPE = abs(PosData[POS_IDX["UPE[m]"]]) 
    
    Label = ["HPE","VPE"]
    Color = ["purple","green"]
    
    PlotConf["xData"]["VPE"] = PosData[POS_IDX["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"]["VPE"] = VPE
    PlotConf["Color"][1] = VPE
    
    PlotConf["xData"]["HPE"] = PosData[POS_IDX["SOD"]] / GnssConstants.S_IN_H
    PlotConf["yData"]["HPE"] = HPE
    PlotConf["Color"][0] = HPE

    PlotConf["yData"].keys()
    PlotConf["Legend"].keys()
    PlotConf["Color"].keys()
    
    PlotConf["Path"] = sys.argv[1] + '/OUT/POS/POS/' + 'POS_HVPE_vs_TIME_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)

    # Plot Satellite EPE vs. NPE Figures 
def plotEPENPE(PosData):
    PlotConf = {}

    PlotConf["Type"] = "Lines"
    PlotConf["FigSize"] = (8.4,7.6)
    PlotConf["Title"] = "EPE vs NPE from TLSA on Year 2015"\
        " DoY 006"

    PlotConf["yLabel"] = "NPE[m]"
    PlotConf["yLim"] = [-3, 5]
    
    PlotConf["xLabel"] = "EPE[m]"
    PlotConf["xLim"] = [-2, 3]

    PlotConf["ColorBar"] = "gnuplot"
    PlotConf["ColorBarLabel"] = "HDOP"
    PlotConf["ColorBarMin"] = 0.6
    PlotConf["ColorBarMax"] = 2.2
     
    PlotConf["Grid"] = 1

    PlotConf["Marker"] = '.'
    PlotConf["LineWidth"] = 1.5
    
      
    PlotConf["xData"] = {}
    PlotConf["yData"] = {}
    PlotConf["zData"] = {}
   
    Label = 0
    PlotConf["xData"][Label] = PosData[POS_IDX["EPE[m]"]]
    PlotConf["yData"][Label] = PosData[POS_IDX["NPE[m]"]]
    PlotConf["zData"][Label] = PosData[POS_IDX["HDOP"]]
   
    PlotConf["Path"] = sys.argv[1] + '/OUT/POS/POS/' + 'POS_NPE_vs_EPE_TLSA_D006Y15.png'

    # Call generatePlot from Plots library
    generatePlot(PlotConf)
