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

# Este código debe ejecutarse en el mismo directorio donde están
# division_comunal.shp, division_comunal.prj, etc.

# Este arreglo contiene la población de cada comuna en el orden en que vienen en division_comunal.*
H = [9148, 11341, 6490, 10479, 11789, 9028, 18357, 33127, 50764, 14639, 3884, 14548, 4097, 5319, 30964, 9191, 9043, 34902, 4326, 12454, 17756, 18151, 16522, 13914, 248945, 39414, 153797, 9480, 11370, 11370, 13425, 7997, 9150, 7122, 1836, 8426, 1665, 5084, 15811, 700, 2208, 348669, 9752, 1206, 13493, 147886, 332, 5605, 4593, 20091, 14146, 16452, 16150, 158261, 13912, 5488, 6531, 9015, 52099, 12725, 7248, 17620, 23857, 44566, 14972, 14158, 18637, 12144, 250638, 56188, 26235, 23430, 47716, 4821, 5591, 7557, 5650, 6802, 13916, 15820, 35664, 6630, 18921, 12373, 10933, 9624, 3458, 75221, 35255, 9182, 30598, 18483, 25671, 11116, 202441, 4331, 211275, 4252, 26029, 13818, 30137, 10418, 21023, 13430, 112956, 4149, 25718, 10825, 19024, 7540, 27616, 30372, 11324, 15481, 35451, 22191, 5451, 32109, 72135, 6911, 23776, 29987, 13481, 269992, 15793, 10403, 23823, 56178, 51136, 21705, 15728, 7569, 28517, 26145, 8933, 3722, 10396, 3627, 14612, 13442, 20572, 227768, 107508, 9867, 86176, 21934, 48045, 52389, 13331, 171673, 55764, 21542, 175585, 29199, 5173, 15461, 115366, 94827, 11830, 5278, 10063, 9120, 10883, 5061, 15525, 11851, 4983, 51335, 3503, 15770, 10051, 4995, 18234, 9571, 6636, 12501, 8062, 32109, 397497, 97230, 99527, 94766, 104018, 225509, 94255, 85195, 32468, 33723, 78887, 177766, 24714, 40245, 9977, 6689, 40329, 9878, 14460, 17839, 28228, 4335, 19361, 107311, 28968, 39404, 18760, 40881, 15428, 17158, 41416, 9714, 4250, 21957, 7959, 9187, 13020, 39554, 8512, 264842, 113340, 79421, 16405, 20732, 757721, 14455, 74232, 66512, 277802, 4646, 28439, 38544, 9380, 16663, 18733, 7006, 19823, 19823, 35185, 31343, 163148, 631, 731, 0, 210936, 679, 739, 1462, 110866, 8097, 17371, 21094, 38524, 5643, 10074, 6726, 19132, 17346, 26385, 88803, 13472, 17029, 11861, 87697, 9205, 13902, 23680, 55207, 17369, 72121, 15829, 14504, 25165, 63210, 16371, 28066, 308137, 311399, 41153, 33535, 5188, 33855, 25790, 20723, 10859, 10039, 3774, 29532, 20964, 195813, 184953, 27781, 107698, 79195, 128090, 121118, 41960, 162671, 249621, 126487, 203946, 84195, 50696, 101737, 152985, 142136, 90846, 237369, 87667, 80910, 121214, 202146, 91927, 289949, 111436, 808000, 19541, 87741, 94455, 873, 5650, 5334, 925, 531, 2759, 803, 1677, 1862, 363, 1163, 59221, 27187, 6166, 21327, 125483, 11110, 119292, 5761, 11519, 1156, 16243, 1384, 162320, 4194, 2360, 135368, 45864, 792]

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
cod=['']*346
comunas=['']*346
R=6371000
IX=0
IY=0
IZ=0
POBLA=0
# Cómputos geográficos preliminares
inputEPSG = 32719
outputEPSG = 4326
inSpatialRef = osr.SpatialReference()
inSpatialRef.ImportFromEPSG(inputEPSG)
outSpatialRef = osr.SpatialReference()
outSpatialRef.ImportFromEPSG(outputEPSG)
coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
# Ciclo por comunas
for j, comuna in enumerate(sf.shapeRecords()):
    nombre = comuna.record[2]
    cods = comuna.record[6]
    inicioPartes = comuna.shape.parts
    puntos = numpy.array(comuna.shape.points)
    psp = separarPartes(puntos, inicioPartes)
    INTpspx=0
    INTpspy=0
    INTpspz=0
    AREApsp =0
    for k in range(0,len(psp)):
        for i in range(0,len(psp[k])):
             pointX = psp[k][i,0]
             pointY = psp[k][i,1]
             point = ogr.Geometry(ogr.wkbPoint)
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
#COORDENADAS EN X Y Z
X=IX/POBLA
Y=IY/POBLA
Z=IZ/POBLA
