#!/usr/bin/env python

## Copyright (C) GNSS ACADEMY 
##
## Name          : receiver_analysis.py
## Purpose       : WP0 Takss: Plot Receiver SPP Analyses
## Project       : WP0-JSNP
## Component     : 
## Author        : GNSS Academy
## Creation date : 2021
## File Version  : 1.0
##

import sys, os

# Add path to find all modules
Common = os.path.dirname(os.path.dirname(
    os.path.abspath(sys.argv[0]))) + '/COMMON'
sys.path.insert(0, Common)

from collections import OrderedDict
from interfaces import LOS_IDX
from interfaces import POS_IDX
from pandas import read_csv 
from yaml import dump
import SatFunctions

#######################################################
# INTERNAL FUNCTIONS 
#######################################################

def displayUsage():
    sys.stderr.write("ERROR: Please provide path to SCENARIO as a unique \nargument\n")

def readConf(CfgFile):
    Conf = OrderedDict({})
    with open(CfgFile, 'r') as f:
        # Read file
        Lines = f.readlines()

        # Read each configuration parameter which is compound of a key and a value
        for Line in Lines:
            if "#" in Line or Line.isspace(): continue
            LineSplit = Line.split('=')
            try:
                LineSplit = list(filter(None, LineSplit))
                Conf[LineSplit[0].strip()] = LineSplit[1].strip()

            except:
                sys.stderr.write("ERROR: Bad line in conf: %s\n" % Line)

    return Conf

#######################################################
# MAIN PROCESSING
#######################################################

print( '-----------------------------')
print( 'RUNNING RECEIVER ANALYSES ...')
print( '-----------------------------')

if len(sys.argv) != 2:
    displayUsage()
    sys.exit()

# Take the arguments
Scen = sys.argv[1]

# Path to conf
CfgFile = Scen + '/CFG/receiver_analysis.cfg'

# Read conf file
Conf = readConf(CfgFile)

# Print 
print('Reading Configuration file',CfgFile)

#print(dump(Conf))

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>> LOS FILE ANALYSES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Get LOS file full path
LosFile = Scen + '/OUT/LOS/' + Conf["LOS_FILE"]

#-----------------------------------------------------------------------
# PLOT SATELLITE ANALYSES
#-----------------------------------------------------------------------

# Plot Satellite Visibility figures
if(Conf["PLOT_SATVIS"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["PRN"],LOS_IDX["ELEV"]])
    
    print( 'Plot Satellite Visibility Periods ...')

    # Configure plot and call plot generation function
    SatFunctions.plotSatVisibility(LosData)

# Plot Satellite Geometrical Ranges figures
if(Conf["PLOT_SATRNG"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["RANGE[m]"],LOS_IDX["ELEV"]])

    print( 'Plot Satellite Geometrical Ranges ...')
    
    # Configure plot and call plot generation function
    SatFunctions.plotSatGeomRnge(LosData)

# Plot Satellite Tracks figures
if(Conf["PLOT_SATTRK"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],
    LOS_IDX["SAT-X[m]"],
    LOS_IDX["SAT-Y[m]"],
    LOS_IDX["SAT-Z[m]"],
    LOS_IDX["ELEV"]])
    
    print( 'Plot Satellite Tracks ...')

    # Configure plot and call plot generation function
    SatFunctions.plotSatTracks(LosData)

    # Plot Satellite Range Velocity figures
if(Conf["PLOT_SATVEL"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],
    LOS_IDX["VEL-X[m/s]"], 
    LOS_IDX["VEL-Y[m/s]"],
    LOS_IDX["VEL-Z[m/s]"], 
    LOS_IDX["ELEV"]])

    print( 'Plot Satellite Range Velocity ...')
    
    # Configure plot and call plot generation function
    SatFunctions.plotSatVelocity(LosData)

#Plot PRN19 NAV CLK figures
if(Conf["PLOT_NAVSATCLK"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["SV-CLK[m]"],LOS_IDX["PRN"]])
    

    print( 'Plot Satellite NAV Clock ...')
    
    # Configure plot and call plot generation function
    SatFunctions.plotPRN19NavSatClock(LosData)

    #Plot Satellite Clock Corrected figures
if(Conf["PLOT_SATCLK"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["SV-CLK[m]"],LOS_IDX["DTR[m]"],LOS_IDX["TGD[m]"],LOS_IDX["PRN"]])
    

    print( 'Plot Satellite Clock Corrected ...')
    
    # Configure plot and call plot generation function
    SatFunctions.plotSatClockCorrected(LosData)

    #Plot Satellite TGD  figures
if(Conf["PLOT_SATTGD"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["TGD[m]"],LOS_IDX["PRN"]])
    

    print( 'Plot Satellite TGD ...')
    
    # Configure plot and call plot generation function
    SatFunctions.plotSatTGD(LosData)

    #Plot Satellite DTR  figures
if(Conf["PLOT_SATDTR"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["DTR[m]"],LOS_IDX["ELEV"]])
    

    print( 'Plot Satellite DTR ...')
    
    # Configure plot and call plot generation function
    SatFunctions.plotSatDTR(LosData)
  
#-----------------------------------------------------------------------
# PLOT IONOSPHERE ANALYSES
#-----------------------------------------------------------------------


# Plot STEC vs TIME(ELEV) figures
if(Conf["PLOT_IONODELAYS"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["STEC[m]"],LOS_IDX["ELEV"]])
    
    print( 'Plot of Slant Ionospheric Delays ...')

    # Configure plot and call plot generation function
    SatFunctions.plotIonoDelays(LosData)

# Plot PRN vs TIME(STEC) figures
if(Conf["PLOT_SATVISPERIODS"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["PRN"],LOS_IDX["STEC[m]"]])
    
    print( 'Plot Satellite Visibility Periods ...')

    # Configure plot and call plot generation function
    SatFunctions.plotSatVisibilityPeriods(LosData)

# Plot VTEC vs TIME figures
if(Conf["PLOT_VTECvsTime"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["VTEC[m]"],LOS_IDX["STEC[m]"],LOS_IDX["MPP[elev]"],LOS_IDX["ELEV"]])
    
    print( 'Plot Ionospheric Klobuchar Delays VTEC ...')

    # Configure plot and call plot generation function
    SatFunctions.plotVTECvsTime(LosData)

    # Plot PRN vs. TIME(VTEC) figures
if(Conf["PLOT_VTECvsPRN"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["PRN"],LOS_IDX["VTEC[m]"]])
    
    print( 'Plot satellite visibility periods vs VTEC ...')

    # Configure plot and call plot generation function
    SatFunctions.PlotPRNvsTIME(LosData)

    # Plot STD vs Time(Elevation) figures
if(Conf["PLOT_STDvsTime"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["TROPO[m]"],LOS_IDX["ELEV"]])
    
    print( 'Plot Slant Tropospheric Delays (STD)...')

    # Configure plot and call plot generation function
    SatFunctions.plotSTDvsTIME(LosData)

    # Plot ZTD vs Time(Elevation) figures
if(Conf["PLOT_ZTDvsTime"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["TROPO[m]"],LOS_IDX["ELEV"]])
    
    print( 'Plot Zenith Tropospheric Delays (ZTD)...')

    # Configure plot and call plot generation function
    SatFunctions.plotZTDvsTIME(LosData)

    # Plot PSRvsTime Figures
if(Conf["PLOT_PSRvsTime"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["MEAS[m]"],LOS_IDX["ELEV"]])
    
    print( 'Plot Pseudo-ranges (Code Measurements C1)...')

    # Configure plot and call plot generation function
    SatFunctions.plotPSRvsTime(LosData)

    # Plot ToF vs Time Figures
if(Conf["PLOT_TAUvsTime"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["MEAS[m]"],LOS_IDX["ELEV"]])
    
    print( 'Plot Tau = C1C/c...')

    # Configure plot and call plot generation function
    SatFunctions.plotTAUvsTime(LosData)

     # Plot ToF vs Time Figures
if(Conf["PLOT_TOFvsTime"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["TOF[ms]"],LOS_IDX["ELEV"]])
    
    print( 'Plot Time of Flight (ToF)...')

    # Configure plot and call plot generation function
    SatFunctions.plotTOFvsTime(LosData)

    # Plot the Doppler Frequency Figures
if(Conf["plotDOPPLER_FREQ"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["SAT-X[m]"],LOS_IDX["SAT-Y[m]"],LOS_IDX["SAT-Z[m]"],LOS_IDX["VEL-X[m/s]"],LOS_IDX["VEL-Y[m/s]"],LOS_IDX["VEL-Z[m/s]"],LOS_IDX["ELEV"]])
    
    
    print( 'Plot the Doppler Frequency...')

    # Configure plot and call plot generation function
    SatFunctions.plotDOPPLER_FREQ(LosData)

    # Build and Plot the PVT filter residuals Figures
if(Conf["plotMEAS_RESIDUALS_vs_TIME"] == '1'):
    # Read the cols we need from LOS file
    LosData = read_csv(LosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[LOS_IDX["SOD"],LOS_IDX["MEAS[m]"],LOS_IDX["RANGE[m]"],LOS_IDX["SV-CLK[m]"],LOS_IDX["TGD[m]"],LOS_IDX["DTR[m]"],LOS_IDX["TROPO[m]"],LOS_IDX["STEC[m]"],LOS_IDX["PRN"]])
    
    print( 'Plot the PVT filter residuals...')

    # Configure plot and call plot generation function
    SatFunctions.plotRESC1(LosData)
    
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>> POS FILE ANALYSES
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Get POS file full path
PosFile = Scen + '/OUT/POS/' + Conf["POS_FILE"]
#-----------------------------------------------------------------------
# PLOT Position Analyses
#-----------------------------------------------------------------------

    #  Plot the instantaneous number of satellites Figures
if(Conf["plotPOS_SATS"] == '1'):
    # Read the cols we need from POS file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[POS_IDX["SOD"],POS_IDX["NSATS"]])
    
    print( 'Plot Satellites Used in PVT...')

    # Configure plot and call plot generation function
    SatFunctions.plotPOS_SATS(PosData)

    #  Plot the (X)DOPS Figures
if(Conf["plot_GPTDOPS"] == '1'):
    # Read the cols we need from POS file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[POS_IDX["SOD"],POS_IDX["GDOP"],POS_IDX["PDOP"],POS_IDX["TDOP"]])
    
    print( 'Plot Dilution of Precision (DOP)...')

    # Configure plot and call plot generation function
    SatFunctions.plotGPTDOP(PosData)

    #  Plot the H/V-DOPs Figures
if(Conf["plot_HVDOPs"] == '1'):
    # Read the cols we need from POS file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[POS_IDX["SOD"],POS_IDX["HDOP"],POS_IDX["VDOP"],POS_IDX["NSATS"]])
    
    print( 'Plot the HDOP and VDOP...')

    # Configure plot and call plot generation function
    SatFunctions.plotPOSHVDOP(PosData)

    #Plot ENU_PE vs TIME Figures
if(Conf["plot_ENU"] == '1'):
    # Read the cols we need from POS file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[POS_IDX["SOD"],POS_IDX["EPE[m]"],POS_IDX["NPE[m]"],POS_IDX["UPE[m]"]])
    
    print( 'Plot the East/North/Up Position Error (EPE, NPE, UPE) ...')

    # Configure plot and call plot generation function
    SatFunctions.plotENUPETIME(PosData)

    #Plot HPE and VPE vs Time Figures
if(Conf["plot_HVPE"] == '1'):
    # Read the cols we need from POS file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[POS_IDX["SOD"],POS_IDX["EPE[m]"],POS_IDX["NPE[m]"],POS_IDX["UPE[m]"]])
    
    print( 'Plot the Horizontal and Vertical Position Error (HPE) and (VPE) ...')

    # Configure plot and call plot generation function
    SatFunctions.plotHVPETIME(PosData)

    #Plot NPE vs EPE Figures
if(Conf["plot_POS_NPE_vs_EPE"] == '1'):
    # Read the cols we need from POS file
    PosData = read_csv(PosFile, delim_whitespace=True, skiprows=1, header=None,\
    usecols=[POS_IDX["SOD"],POS_IDX["HDOP"],POS_IDX["NPE[m]"],POS_IDX["EPE[m]"]])
    
    print( 'Plot Horizontal Scatter plot with NPE vs. EPE ...')

    # Configure plot and call plot generation function
    SatFunctions.plotEPENPE(PosData)


   