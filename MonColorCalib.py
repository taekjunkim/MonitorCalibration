#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 16:09:20 2018
Minotor Color calibration

Usage: import MonColorCalib as MC
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

CMF is given by CIE1931(RGB to XYZ)


target xyY
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
from scipy.optimize import least_squares;

def main(RGB=None,xyY=None,fig=0,out=1):
    

    ### read photometer measurement
    ### R,G,B at 51,102,153,204,255 
    ### We can calculate gamma & R,G,B intensity (as a function of lambda) from here
    cMtx = dict();
    cName = ['red','green','blue'];
    #Direc = './colorCalRig1Sept2016/';    
    #Direc = './Dec2018_BenQ_Rig1/';        
    #Direc = './Dec2018_BenQ_Rig1_Std_G5/';
    #Direc = './Jan2019_BenQ_Rig1_Std_G5_T1/';
    #Direc = './Jan2019_BenQ_Rig1_Std_G5_T2/';
    Direc = './Feb2019_BenQ/';    
    for i in range(0,3):
        cMtx[cName[i]] = pd.DataFrame();    
        cMtx[cName[i]]['nm'] = np.arange(380,781,1);
        for j in range(1,5+1):
            stepNum = str(51*j);
            txtName = Direc+stepNum+'-'+cName[i]+'.txt';
            dNow = pd.read_csv(txtName);
            dNow.columns = ['nm',stepNum];
            dNow = dNow.drop(index=0);
            cMtx[cName[i]] = pd.merge(cMtx[cName[i]],dNow,on='nm',how='left');   
        cMtx[cName[i]] = cMtx[cName[i]].interpolate();   ### interpolate to have 1nm resolution        
    
    cdMtx = pd.DataFrame(index=['51','102','153','204','255'],columns=cName);
    cdMtx['red'] = np.sum(cMtx['red'].iloc[:,1:]); 
    cdMtx['green'] = np.sum(cMtx['green'].iloc[:,1:]); 
    cdMtx['blue'] = np.sum(cMtx['blue'].iloc[:,1:]); 
    
    pLevel = np.arange(51,51*(5+1),51)/255;
    fit_RGB = [None]*3;
    #fit_RGB[0] = least_squares(powerFit,[2,0.1],args=(pLevel,cdMtx['red']));
    #fit_RGB[1] = least_squares(powerFit,[2,0.1],args=(pLevel,cdMtx['green']));
    #fit_RGB[2] = least_squares(powerFit,[2,0.1],args=(pLevel,cdMtx['blue']));    
    fit_RGB[0] = least_squares(powerFit,[2,0.1],args=(pLevel[0:4],cdMtx['red'][0:4]));
    fit_RGB[1] = least_squares(powerFit,[2,0.1],args=(pLevel[0:4],cdMtx['green'][0:4]));
    fit_RGB[2] = least_squares(powerFit,[2,0.1],args=(pLevel[0:4],cdMtx['blue'][0:4]));    

    if fig == 1:
        xValues = np.arange(0,256,1)/255;
        plt.figure(figsize=(10,3));
        plt.subplot(1,3,1);
        plt.plot(np.arange(51,306,51),cdMtx['red'][0:5],'ro');
        plt.plot(xValues*255,xValues**fit_RGB[0].x[0]*fit_RGB[0].x[1],'r');     
        plt.xlim((-5,260));
        plt.ylim((-0.01,0.5));
        plt.subplot(1,3,2);
        plt.plot(np.arange(51,306,51),cdMtx['green'][0:5],'go');
        plt.plot(xValues*255,xValues**fit_RGB[1].x[0]*fit_RGB[1].x[1],'g');  
        plt.xlim((-5,260));       
        plt.ylim((-0.01,0.5));
        plt.subplot(1,3,3);
        plt.plot(np.arange(51,306,51),cdMtx['blue'][0:5],'bo');
        plt.plot(xValues*255,xValues**fit_RGB[2].x[0]*fit_RGB[2].x[1],'b');  
        plt.xlim((-5,260));
        plt.ylim((-0.01,0.5));
    
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
                               args=(xyY,cMtx,CMF,fit_RGB));
        if out==1:                       
            print(RGBfit);                       
        return RGBfit;
    if RGB != None:
        Lum = getLumAtRGB(RGB,cMtx,CMF,fit_RGB);
        if out==1:
            print(Lum);
        return Lum;

def getRGBAtxyY(RGB,xyY,cMtx,CMF,fit_RGB):

    XYZ = np.array([None]*3,dtype=float);
    XYZ[0] = xyY[0]*xyY[2]/xyY[1];
    XYZ[1] = xyY[2];
    XYZ[2] = (1-xyY[0]-xyY[1])*xyY[2]/xyY[1];
    
    SPD_CMF = pd.Series();
    """
    SPD_CMF = (RGB[0]/255)**fit_RGB[0].x[0]*cMtx['red']['255'] + \
                     (RGB[1]/255)**fit_RGB[1].x[0]*cMtx['green']['255'] + \
                     (RGB[2]/255)**fit_RGB[2].x[0]*cMtx['blue']['255'];
    """                 
    SPD_CMF = (RGB[0]/204)**fit_RGB[0].x[0]*cMtx['red']['204'] + \
                     (RGB[1]/204)**fit_RGB[1].x[0]*cMtx['green']['204'] + \
                     (RGB[2]/204)**fit_RGB[2].x[0]*cMtx['blue']['204'];

    SPD_CMF.index = cMtx['red']['nm'];   
    XYZnow = [np.dot(SPD_CMF,CMF['X'])*683,
              np.dot(SPD_CMF,CMF['Y'])*683,
              np.dot(SPD_CMF,CMF['Z'])*683];
    XYZnow = np.array(XYZnow);
    #import pdb; pdb.set_trace();
    
    return XYZnow - XYZ;


def getLumAtRGB(RGB,cMtx,CMF,fit_RGB):
    ### spectral power distribution
    SPD_CMF = pd.Series();
    """
    SPD_CMF = (RGB[0]/255)**fit_RGB[0].x[0]*cMtx['red']['255'] + \
                     (RGB[1]/255)**fit_RGB[1].x[0]*cMtx['green']['255'] + \
                     (RGB[2]/255)**fit_RGB[2].x[0]*cMtx['blue']['255'];
    """                 
    SPD_CMF = (RGB[0]/204)**fit_RGB[0].x[0]*cMtx['red']['204'] + \
                     (RGB[1]/204)**fit_RGB[1].x[0]*cMtx['green']['204'] + \
                     (RGB[2]/204)**fit_RGB[2].x[0]*cMtx['blue']['204'];
    SPD_CMF.index = cMtx['red']['nm'];         

    XYZ = [np.dot(SPD_CMF,CMF['X'])*683,
              np.dot(SPD_CMF,CMF['Y'])*683,
              np.dot(SPD_CMF,CMF['Z'])*683];
    x = XYZ[0]/np.sum(XYZ);       
    y = XYZ[1]/np.sum(XYZ);           
    
    ### xyY is provided
    return [x,y,XYZ[1]];
    
def powerFit(X,pLevel,Lum):
    y = (pLevel**X[0])*X[1];
    return y-Lum;

def makeTable():
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
    
