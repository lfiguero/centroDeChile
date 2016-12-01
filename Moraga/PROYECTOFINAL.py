# coding: utf-8
import shapefile
import numpy
import matplotlib.pyplot as plt
import ogr, osr # Paquetes que vienen con GDAL
from scipy import integrate, cos, sin, pi, arctan, sqrt
from tqdm import tqdm # Barra de progreso para ciclos largos
from H import H # Esto está en H.py

# Este código debe ejecutarse en el mismo directorio donde están
# division_comunal.shp, division_comunal.prj, etc.

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
def Ix(Pi,Pi1,Ti,Ti1):
    F= lambda s: -(cos(Ti - s*(Ti - Ti1)))*(cos(Ti - s*(Ti - Ti1)))*sin(Pi - s*(Pi - Pi1))*(Ti - Ti1)
    I=integrate.quad(F,0,1)
    I=-I[0]
    return I
def Iy(Pi,Pi1,Ti,Ti1):
    F= lambda s: cos(Pi - s*(Pi - Pi1))*sin(Ti - s*(Ti - Ti1))*(sin(Ti - s*(Ti - Ti1)))*(Ti - Ti1)
    I=integrate.quad(F,0,1)
    I=-I[0]
    return I
def Iz(Pi,Pi1,Ti,Ti1):
    F= lambda s: (sin(2.0*Ti - 2.0*s*(Ti - Ti1))*(Pi - Pi1))/(8.0*(Ti - Ti1))
    I=F(1)-F(0)
    I=-I
    return I 
def Areap(Pi,Pi1,Ti,Ti1):
    if Ti == Ti1:
        aux = cos(Ti)
    else:
        aux = (sin(Ti1) - sin(Ti)) / (Ti1 - Ti)
    I = (Pi1-Pi)*aux
    I=-I
    return I
# Código principal
sf = shapefile.Reader("division_comunal")
ncomunas = 346
R = 6371000.0
IX = 0.0
IY = 0.0
IZ = 0.0
POBLA = 0.0
# Cómputos geográficos preliminares
inputEPSG = 32719
outputEPSG = 4326
inSpatialRef = osr.SpatialReference()
inSpatialRef.ImportFromEPSG(inputEPSG)
outSpatialRef = osr.SpatialReference()
outSpatialRef.ImportFromEPSG(outputEPSG)
coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)
# Alocación de memoria para integrales sobre comunas
INTpspx = numpy.zeros(ncomunas)
INTpspy = numpy.zeros(ncomunas)
INTpspz = numpy.zeros(ncomunas)
AREApsp = numpy.zeros(ncomunas)
# Ciclo por comunas
for j, comuna in enumerate(tqdm(sf.shapeRecords())):
    nombre = comuna.record[2]
    cods = comuna.record[6]
    inicioPartes = comuna.shape.parts
    puntos = numpy.array(comuna.shape.points)
    psp = separarPartes(puntos, inicioPartes)
    # Ciclo por polígonos
    for k in range(0,len(psp)):
        for i in range(0,len(psp[k])):
             pointX = psp[k][i,0]
             pointY = psp[k][i,1]
             point = ogr.Geometry(ogr.wkbPoint)
             # create a geometry from coordinates
             point.AddPoint(pointX, pointY)
             # transform point
             point.Transform(coordTransform)
             # Sobreescribimos con nuevas coordenadas
             psp[k][i,0] = (point.GetX())*pi/180 # Ángulo azimutal (φ)
             psp[k][i,1] = -((point.GetY()*pi/180)+pi/2) # Ángulo polar/colatitud (θ)
        #zona poligonal en coordenadas lat/long
        AREAi = 0.0
        INTxi = 0.0
        INTyi = 0.0
        INTzi = 0.0
        # Ciclo por segmentos
        for i in range(0,len(psp[k])-1):
            #                φ[i]            φ[i+1]        θ[i]          θ[i+1]
            #                  ↓               ↓            ↓              ↓ 
            AREAi += Areap(psp[k][i,0], psp[k][i+1,0], psp[k][i,1], psp[k][i+1,1])
            INTxi += Ix(psp[k][i,0], psp[k][i+1,0], psp[k][i,1], psp[k][i+1,1])
            INTyi += Iy(psp[k][i,0], psp[k][i+1,0], psp[k][i,1], psp[k][i+1,1])
            INTzi += Iz(psp[k][i,0], psp[k][i+1,0], psp[k][i,1], psp[k][i+1,1])
        INTpspx[j] += INTxi
        INTpspy[j] += INTyi
        INTpspz[j] += INTzi
        AREApsp[j] += AREAi
    Densicomunal=H[j]/AREApsp[j]
    IX += Densicomunal*R*INTpspx[j]
    IY += Densicomunal*R*INTpspy[j]
    IZ += Densicomunal*R*INTpspz[j]
    POBLA += H[j]
#COORDENADAS EN X Y Z
X = IX/POBLA
Y = IY/POBLA
Z = IZ/POBLA
THETA = (arctan(sqrt(X*X+Y*Y)/Z))*180/pi-90 # Latitud en grados
PHI = ((arctan(Y/X))*180/pi) # Longitud en grados
print("LON",PHI)
print("LAT",THETA)
