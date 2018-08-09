#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 16:09:20 2018
Minotor Color calibration

Usage: import MonColorCalib as MC
       e.g., Lum(cd/m2) = MC.main(RGB=[50,50,50]);
       e.g., RGBfitted = MC.main(xyY=[0.33,0.33,8]);       


Luminance = cd/m2 = lumen/sr/m2
          = Radiance * CMF(Color Matching Function) * nm
          = (Watt/sr/nm/m2) * CMF * nm
          (%% 1.0 on CMF = 683*lumen/Watt)
          = (Watt/sr/nm/m2) * (683*lumen/Watt) * nm

Radiance = ((r/255)^Rgamma * Rintensity(lambda))
           + ((g/255)^Ggamma * Gintensity(lambda))
           + ((b/255)^Bgamma * Bintensity(lambda))

CMF is given by CIE1931(RGB to XYZ)

@author: taekjunkim
"""

import pandas as pd;
import numpy as np;
from scipy.optimize import least_squares;

def main(RGB=None,xyY=None):
    

    ### read photometer measurement
    ### R,G,B at 51,102,153,204,255 
    ### We can calculate gamma & R,G,B intensity (as a function of lambda) from here
    cMtx = pd.DataFrame();
    cName = ['red','green','blue'];
    Direc = './colorCalRig1Sept2016/';    
    for i in range(0,3):
        cMtx[cName[i]]="";    
        for j in range(1,5+1):
            stepNum = str(51*j);
            txtName = Direc+stepNum+'-'+cName[i]+'.txt';
            dNow = pd.read_csv(txtName);
            dNow.columns = ['nm',stepNum];
            dNow = dNow.drop(index=0);
            if j==1:
                cMtx[cName[i]] = [dNow];
            else:
                cMtx[cName[i]][0][stepNum] = dNow[stepNum];
    
    cdMtx = pd.DataFrame(index=['51','102','153','204','255'],columns=cName);
    cdMtx['red'] = np.sum(cMtx['red'][0])[1:]; 
    cdMtx['green'] = np.sum(cMtx['green'][0])[1:]; 
    cdMtx['blue'] = np.sum(cMtx['blue'][0])[1:]; 
    
    pLevel = np.arange(51,51*(5+1),51)/255;
    fit_RGB = [None]*3;
    fit_RGB[0] = least_squares(powerFit,[2,1],args=(pLevel,cdMtx['red']));
    fit_RGB[1] = least_squares(powerFit,[2,1],args=(pLevel,cdMtx['green']));
    fit_RGB[2] = least_squares(powerFit,[2,1],args=(pLevel,cdMtx['blue']));    
 
    ### Color matching function: CIE1931 - RGB to XYZ    
    CMFraw = pd.read_csv('cie1931_XYZ_CMF.csv',header=None);
    CMF = pd.DataFrame()
    CMF['r'] = CMFraw.iloc[:,1];
    CMF['g'] = CMFraw.iloc[:,2];
    CMF['b'] = CMFraw.iloc[:,3];
    CMF.index = range(360,360+471);
    CMF = CMF.loc[range(380,784,4)];   ### CMF sampled every 4nm (4 will be multiplied later)

    if xyY != None:
        ### Spectral power distribution
        RGBfit = least_squares(getRGBAtxyY,np.array([100.0,100.0,100.0]),
                               args=(xyY,cMtx,CMF,fit_RGB));
        return RGBfit;
    if RGB != None:
        Lum = getLumAtRGB(RGB,cMtx,CMF,fit_RGB);
        return Lum;

def getRGBAtxyY(RGB,xyY,cMtx,CMF,fit_RGB):

    XYZ = np.array([None]*3,dtype=float);
    XYZ[0] = xyY[0]*xyY[2]/xyY[1];
    XYZ[1] = xyY[2];
    XYZ[2] = (1-xyY[0]-xyY[1])*xyY[2]/xyY[1];
    
    SPD_CMF = pd.Series();
    SPD_CMF = (RGB[0]/255)**fit_RGB[0].x[0]*cMtx['red'][0]['255'] + \
                     (RGB[1]/255)**fit_RGB[1].x[0]*cMtx['green'][0]['255'] + \
                     (RGB[2]/255)**fit_RGB[2].x[0]*cMtx['blue'][0]['255'];
    SPD_CMF.index = cMtx['red'][0]['nm'];   
    XYZnow = [np.dot(SPD_CMF,CMF['r'])*683*4,
              np.dot(SPD_CMF,CMF['g'])*683*4,
              np.dot(SPD_CMF,CMF['b'])*683*4];
    XYZnow = np.array(XYZnow);
    #import pdb; pdb.set_trace();
    
    return XYZnow - XYZ;

    


def getLumAtRGB(RGB,cMtx,CMF,fit_RGB):
    ### spectral power distribution
    SPD_CMF = pd.Series();
    SPD_CMF = (RGB[0]/255)**fit_RGB[0].x[0]*cMtx['red'][0]['255'] + \
                     (RGB[1]/255)**fit_RGB[1].x[0]*cMtx['green'][0]['255'] + \
                     (RGB[2]/255)**fit_RGB[2].x[0]*cMtx['blue'][0]['255'];
    SPD_CMF.index = cMtx['red'][0]['nm'];         
    
    ### Y of XYZ is computed
    return np.dot(SPD_CMF,CMF['g'])*683*4;
    
def powerFit(X,pLevel,Lum):
    y = (pLevel**X[0])*X[1];
    return y-Lum;
    
    