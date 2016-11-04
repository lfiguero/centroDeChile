# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 23:02:49 2016

@author: Sebastián
"""
# coding: utf-8
import shapefile
import numpy
import matplotlib.pyplot as plt
import os
import ogr, osr # Paquetes que vienen con GDAL
import os, sys
import scipy as sp
from scipy import integrate

# Este código debe ejecutarse en el mismo directorio donde estan
# division_comunal.shp, division_comunap.prj, etc. Este código dibuja cada
# parte (incluyendo islas, etc.) de cada comuna por separado.
#cambiamos el directorio
os.chdir("C:\\proyectotaller")
# Función auxiliar que separa arreglo de puntos en las distintas partes
def separarPartes(puntos, inicioPartes):
    """ puntos es un arreglo de numpy; inicioPartes una lista de primer índice
    de cada parte """
    salida = []
    npartes = len(inicioPartes)
    for k in range(0, npartes):
        if k < npartes-1:
            puntos_parte_k = puntos[inicioPartes[k]:inicioPartes[k+1], :]
        else:
            puntos_parte_k = puntos[inicioPartes[k]:, :]
        salida.append(puntos_parte_k)
    return salida
#Funciones integrales
def areaBruta(puntos_parte):
    salida = 0.0
    for i in range(len(puntos_parte)-1):
        salida = salida - (puntos_parte[i+1][1]-puntos_parte[i][1])*(puntos_parte[i][0]+puntos_parte[i+1][0])/2.0
    return salida
def Ix(Pi,Pi1,Ti,Ti1):
    F= lambda s: (5.0/4)*((sp.cos(Pi - s*(Pi - Pi1))*sp.sin(2.0*Ti - 2.0*s*(Ti - Ti1))*(Ti - Ti1))/5.0 - sp.sin(Pi - s*(Pi - Pi1))*(sp.cos(2*Ti + 2*s)/10.0 - 1.0/2)*(Pi - Pi1))
    I=integrate.quad(F,0,1)
    I=-I[0]
    return I
def Iy(Pi,Pi1,Ti,Ti1):
    F= lambda s: (sp.cos(Pi - s*(Pi - Pi1))*(sp.sin(Ti - s*(Ti - Ti1))*sp.sin(Ti - s*(Ti - Ti1)))*(Ti - Ti1))/2.0 - (sp.sin(Pi - s*(Pi - Pi1))*(Pi - Pi1)*(sp.sin(2.0*Ti - 2.0*s*(Ti - Ti1))/4.0 - Ti/2.0 + (s*(Ti - Ti1))/2.0))/2.0
    I=integrate.quad(F,0,1)
    I=-I[0]
    return I
def Iz(Pi,Pi1,Ti,Ti1):
    F= lambda s: -(sp.cos(2.0*Ti - 2.0*s*(Ti - Ti1))*(Pi - Pi1))/4.0
    I=integrate.quad(F,0,1)
    I=-I[0]
    return I
def Areap(Pi,Pi1,Ti,Ti1):
    F= lambda s: -sp.cos(Ti+(Ti1-Ti)*s)*(Pi-Pi1)
    I=integrate.quad(F,0,1)
    I=-I[0]
    return I
# Código principal
sf = shapefile.Reader("division_comunal")
i=0
k=0
j=0
cod=['']*346
comunas=['']*346 
R=6371000   
IX=0
IY=0
IZ=0
POBLA=0   
for comuna in sf.shapeRecords():
    nombre = comuna.record[2]
    cods = comuna.record[6]
    inicioPartes = comuna.shape.parts
    puntos = numpy.array(comuna.shape.points)
    psp = separarPartes(puntos, inicioPartes)   
    psp2=psp
    psp3=['']*len(psp) 
    INTpspx=0
    INTpspy=0
    INTpspz=0
    AREApsp =0      
    for k in range(0,len(psp)):        
        for i in range(0,len(psp[k])):
             pointX = psp[k][i,0]
             pointY = psp[k][i,1] 
             inputEPSG = 32719
             outputEPSG = 4326
             inSpatialRef = osr.SpatialReference()
             inSpatialRef.ImportFromEPSG(inputEPSG)
             outSpatialRef = osr.SpatialReference()
             outSpatialRef.ImportFromEPSG(outputEPSG)
             point = ogr.Geometry(ogr.wkbPoint)
             coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)        
             # transform point
             point.Transform(coordTransform)
             # create a geometry from coordinates           
             point.AddPoint(pointX, pointY)
             # transform point
             point.Transform(coordTransform)        
             # print point in EPSG 4326
             psp[k][i,0]=point.GetX()
             psp[k][i,1]=point.GetY()
        #zona poligonal en coordenadas lat/long           
        AREAi=0
        INTxi=0
        INTyi=0
        INTzi=0
        for i in range(0,len(psp[k])-1):
            Pi=psp[k][i,1]
            Pi1=psp[k][i+1,1]
            Ti=psp[k][i,0]
            Ti1=psp[k][i+1,0]
            IIX=Ix(Pi,Pi1,Ti,Ti1)
            IIY=Iy(Pi,Pi1,Ti,Ti1)
            IIZ=Iz(Pi,Pi1,Ti,Ti1)
            A=Areap(Pi,Pi1,Ti,Ti1)
            AREAi=AREAi+A
            INTxi=INTxi+IIX
            INTyi=INTyi+IIY
            INTzi=INTzi+IIZ
        INTpspx=INTpspx+INTxi
        INTpspy=INTpspy+INTyi
        INTpspz=INTpspz+INTzi
        AREApsp=AREApsp+AREAi
    Densicomunal=H[j]/AREApsp
    ICX=Densicomunal*R*INTpspx
    ICy=Densicomunal*R*INTpspy
    ICz=Densicomunal*R*INTpspz
    IX=IX+ICX
    IY=IY+ICy
    IZ=IZ+ICz
    POBLA=POBLA+H[j]
    j=j+1
#COORDENADAS EN X Y Z
X=IX/POBLA
Y=IY/POBLA
Z=IZ/POBLA
            
            
