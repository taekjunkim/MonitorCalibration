# MonitorCalibration 

## I already have photometer measurements for R,G,B
- e.g., 51-red.txt, 102-red.txt, 153-red.txt, 204-red.txt, 255-red.txt
- Each txt file contains "Radiance = (Watt/sr/nm/m2)" values measured between 380 and 780nm (4nm step)

## In this code
- Estimate Gamma for R,G,B based on photometer measurements (version 1)
- Calculate SPD (Spectral Power Distribution) for a given RGB
- Calculate XYZ from SPD
   - dot product of SPD and CMF (Color Matching Function)
   - For CMF, CIE1931 is used
- Estimate the nearest RGB for a given xyY

## Version1 vs. Version2
- Version1 assumes that R,G,B gun outputs follow a gamma function
- In Version2, intermediate R,G,B values were estimated from spline functions
