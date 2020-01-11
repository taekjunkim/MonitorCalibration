#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 17:01:08 2019
Minotor Color calibration

Usage: import MonColorCalib_v2 as MC
       e.g., Lum = MC.main(RGB=[50,50,50]); ## Lum = [x,y,Y(cd/m2)];
       e.g., RGBfitted = MC.main(xyY=[0.33,0.33,8]);       
       e.g., MC.makeTable();       


Luminance = cd/m2 = lumen/sr/m2
          = Radiance * CMF(Color Matching Function) * nm
          = (Watt/sr/nm/m2) * CMF * nm
          (%% 1.0 on CMF = 683*lumen/Watt)
          = (Watt/sr/nm/m2) * (683*lumen/Watt) * nm

Radiance = ((r/255)^Rgamma * Rintensity(lambda))
           + ((g/255)^Ggamma * Gintensity(lambda))
           + ((b/255)^Bgamma * Bintensity(lambda))

### In version 2, I'm trying spline fitting instead of gamma fitting, 
because I found undesirable spectral density function shape change 
depending on RGB values. For example, SPD at very low luminance cannot be 
expressed by a gain control of SPD at high luminance. 

### Possible solution: at each nm, spectral outputs are interpolated in 
the range of 0~255.



CMF is given by CIE1931(RGB to XYZ)


target xyY --- order from Bushnell et al. (2011)
           --- see below (def makeTable) to see order from colors.py
1: (0.33,0.33) at 4,8,12,16,20cd/m2; 
2: (0.3499,0.3483) at 4,8,12,16,20cd/m2; 
3: (0.3866,0.3266) at 4,8,12,16,20cd/m2; 
4: (0.3316,0.2883) at 4,8,12,16,20cd/m2; 
5: (0.2766,0.25) at 4,8,12,16,20cd/m2; 
6: (0.2949,0.31) at 4,8,12,16,20cd/m2; 
7: (0.3133,0.37) at 4,8,12,16,20cd/m2;
8: (0.3699,0.3666) at 4,8,12,16,20cd/m2; 
9: (0.4433,0.3233) at 4,8,12,16,20cd/m2; 
10: (0.3333,0.2466) at 4,8,12,16,20cd/m2; 
11: (0.2233,0.17) at 4,8,12,16,20cd/m2; 
12: (0.2599,0.29) at 4,8,12,16,20cd/m2; 
13: (0.2966,0.41) at 4,8,12,16,20cd/m2;
14: (0.39,0.385) at 4,8,12,16,20cd/m2;
15: (0.445,0.3525) at 4,8,12,16,20cd/m2;
16: (0.5,0.32) at 4,8,12,16,20cd/m2;
17: (0.4175,0.2625) at 4,8,12,16,20cd/m2;
18: (0.335,0.205) at 4,8,12,16,20cd/m2;
19: (0.2525,0.1475) at 4,8,12,16,20cd/m2;
20: (0.17,0.09) at 4,8,12,16,20cd/m2;
21: (0.1975,0.18) at 4,8,12,16,20cd/m2;
22: (0.225,0.27) at 4,8,12,16,20cd/m2;
23: (0.2525,0.36) at 4,8,12,16,20cd/m2;
24: (0.28,0.45) at 4,8,12,16,20cd/m2;
25: (0.335,0.4175) at 4,8,12,16,20cd/m2;


@author: taekjunkim
"""

import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt;
from scipy.interpolate import interp2d
from scipy.optimize import least_squares;

def main(RGB=None,xyY=None,fig=0,out=1):

    ### read photometer measurement
    ### First goal is to generate spectral density function matrix (256x401x3)
    ### for intensity (0~255), for wavelength 380~780 (401), for R, G, B (3)
    
    ### To do this, 
    ### 1. make an empty ndarray (256x401x3)
    ### 2. fill the ndarray with the measured values
    ### 3. interpolate the missing values with spline fitting
    
    spdMtx = np.empty((256,401,3));
    spdMtx[:] = np.nan;
    wavelengths = np.arange(380,781,1);
    intensity = np.arange(0,256,1);
    
    intensitySampled = np.arange(0,306,51);
    cName = ['red','green','blue'];

    Direc = './Dec2019_BenQ_Rig3/';    
    
    fit_RGB = [None]*3;    
    for i in range(0,3): # R,G,B
        for j in range(1,len(intensitySampled)): # intensitySampled
            stepNum = str(intensitySampled[j]);
            txtName = Direc+stepNum+'-'+cName[i]+'.txt';
            dNow = pd.read_csv(txtName);
            dNow.columns = ['nm',stepNum];
            dNow = dNow.drop(index=0);
            
            dataNow = dNow[stepNum].get_values().reshape(1,-1);
            if j==1:
                measured = dataNow;
            else:
                measured = np.vstack((measured,dataNow));

        measured = np.vstack((np.zeros((1,np.shape(measured)[1])),measured));
        x_val = dNow['nm'].get_values();
        y_val = intensitySampled;
        fit_RGB[i] = interp2d(x_val, y_val, measured, kind='cubic')
        fitted = fit_RGB[i](wavelengths,intensity);

        spdMtx[:,:,i] = fitted;

    if fig == 1:
        plt.figure(figsize=(6,3));        
        for i in range(3):
            plt.subplot(2,3,i+1);
            plt.imshow(spdMtx[:,:,i],vmin=0,vmax=np.max(spdMtx));
            
            plt.subplot(2,3,i+4);
            for j in range(0,275,25):
                plt.plot(wavelengths,spdMtx[j,:,i]);

    ### Color matching function: CIE1931 - RGB to XYZ    
    CMFraw = pd.read_csv('cie1931_XYZ_CMF.csv',header=None);
    CMF = pd.DataFrame()
    CMF['X'] = CMFraw.iloc[:,1];
    CMF['Y'] = CMFraw.iloc[:,2];
    CMF['Z'] = CMFraw.iloc[:,3];
    CMF.index = range(360,360+471);
    CMF = CMF.loc[range(380,781,1)];   ### range matching

    if xyY != None:
        ### Spectral power distribution
        RGBfit = least_squares(getRGBAtxyY,np.array([100.0,100.0,100.0]),
                               bounds=([0,0,0],[255,255,255]),
                               args=(xyY,fit_RGB,CMF));
        if out==1:                       
            print(RGBfit);                       
        return RGBfit;
    if RGB != None:
        Lum = getLumAtRGB(RGB,spdMtx,CMF);
        if out==1:
            print(Lum);
        return Lum;

def getRGBAtxyY(RGB,xyY,fit_RGB,CMF):

    wavelengths = np.arange(380,781,1);    
    
    XYZ = np.array([None]*3,dtype=float);
    XYZ[0] = xyY[0]*xyY[2]/xyY[1];
    XYZ[1] = xyY[2];
    XYZ[2] = (1-xyY[0]-xyY[1])*xyY[2]/xyY[1];

    SPD_CMF = fit_RGB[0](wavelengths,RGB[0]) + \
              fit_RGB[1](wavelengths,RGB[1]) + \
              fit_RGB[2](wavelengths,RGB[2]);
                     
    XYZnow = [np.dot(SPD_CMF,CMF['X'])*683,
              np.dot(SPD_CMF,CMF['Y'])*683,
              np.dot(SPD_CMF,CMF['Z'])*683];
    XYZnow = np.array(XYZnow);
    
    return XYZnow - XYZ;


def getLumAtRGB(RGB,spdMtx,CMF):
    ### spectral power distribution
    RGB = np.round(RGB).astype('int');
    
    SPD_CMF = spdMtx[int(np.round(RGB[0])),:,0] + \
              spdMtx[int(np.round(RGB[1])),:,1] + \
              spdMtx[int(np.round(RGB[2])),:,2];
                     
    XYZ = [np.dot(SPD_CMF,CMF['X'])*683,
              np.dot(SPD_CMF,CMF['Y'])*683,
              np.dot(SPD_CMF,CMF['Z'])*683];
    x = XYZ[0]/np.sum(XYZ);       
    y = XYZ[1]/np.sum(XYZ);           
    
    ### xyY is provided
    return [x,y,XYZ[1]];

def makeTable():
    ### color order from colors.py
    xyTable = np.array(([0.17, 0.09], [0.198, 0.18], [0.223, 0.17],
                        [0.225, 0.27], [0.253, 0.148], [0.253, 0.36],
                        [0.26, 0.29], [0.277, 0.25], [0.28, 0.45],
                        [0.295, 0.31], [0.297, 0.41], [0.313, 0.37],
                        [0.33, 0.33], [0.335, 0.205], [0.333, 0.247],
                        [0.332, 0.288], [0.335, 0.418], [0.35, 0.348],
                        [0.37, 0.367], [0.387, 0.327], [0.39, 0.385],
                        [0.418, 0.263], [0.443, 0.323], [0.445, 0.353],
                        [0.5, 0.32]));
    LumTable = np.array([4,8,12,16,20]);
    LumTable = np.reshape(LumTable,(5,1));
    
    for i in range(0,len(LumTable)):
        print('Luminance = ',LumTable[i]);
        for j in range(0,xyTable.shape[0]):
            xyYnow = np.hstack((xyTable[j,:],LumTable[i]));
            xyYnow = list(xyYnow);
            RGBfitted = main(xyY=xyYnow,out=0);
            print(np.round(RGBfitted.x));