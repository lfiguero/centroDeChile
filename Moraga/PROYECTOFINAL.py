# coding: utf-8
import shapefile
import numpy
import matplotlib.pyplot as plt
import ogr, osr # Paquetes que vienen con GDAL
from scipy import integrate, cos, sin, pi, arctan, arctan2, sqrt
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

# Funciones que transforman entre longitud-latitud (grados) y
# ángulo azimutal-ángulo polar/colatitud (radianes)
def lonlat2phitheta(lon, lat):
    phi = lon*pi/180.0
    theta = pi/2.0 - lat*pi/180.0
    return phi, theta
def phitheta2lonlat(phi, theta):
    lon = phi*180.0/pi
    lat = 90.0 - theta*180.0/pi
    return lon, lat
# Funciones que transforman entre coordenadas cartesianas y coordenadas
# esféricas (en el orden radio, ángulo azimutal, ángulo polar/colatitud)
def cart2sph(x, y, z):
    r = sqrt(x**2 + y**2 + z**2)
    theta = arctan2(sqrt(x**2+y**2), z)
    phi = arctan2(y, x)
    return r, theta, phi
def sph2cart(r, theta, phi):
    x = r*cos(phi)*sin(theta)
    y = r*sin(phi)*sin(theta)
    z = r*cos(theta)
    return x, y, z


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
intComX = numpy.zeros(ncomunas)
intComY = numpy.zeros(ncomunas)
intComZ = numpy.zeros(ncomunas)
intCom1 = numpy.zeros(ncomunas)
# Variable para contener las geometrías de las comunas
psp = [None]*ncomunas
# Ciclo por comunas
for j, comuna in enumerate(tqdm(sf.shapeRecords())):
    nombre = comuna.record[2]
    cods = comuna.record[6]
    inicioPartes = comuna.shape.parts
    puntos = numpy.array(comuna.shape.points)
    psp[j] = separarPartes(puntos, inicioPartes)
    # Ciclo por polígonos
    for k in range(0,len(psp[j])):
        # Ciclo por vértices
        for i in range(0,len(psp[j][k])):
             pointX = psp[j][k][i,0]
             pointY = psp[j][k][i,1]
             point = ogr.Geometry(ogr.wkbPoint)
             # create a geometry from coordinates
             point.AddPoint(pointX, pointY)
             # transform point
             point.Transform(coordTransform)
             # Sobreescribimos con nuevas coordenadas: De longitud-latitud a
             # ángulo azimutal (φ)-ángulo polar/colatitud (θ)
             psp[j][k][i,0], psp[j][k][i,1] = \
                     lonlat2phitheta(point.GetX(), point.GetY())
        AREAk = 0.0
        INTxk = 0.0
        INTyk = 0.0
        INTzk = 0.0
        # Ciclo por segmentos
        for i in range(0,len(psp[j][k])-1):
            #                  φ[i]            φ[i+1]             θ[i]           θ[i+1]
            #                    ↓               ↓                 ↓               ↓ 
            AREAk += Areap(psp[j][k][i,0], psp[j][k][i+1,0], psp[j][k][i,1], psp[j][k][i+1,1])
            INTxk += Ix(psp[j][k][i,0], psp[j][k][i+1,0], psp[j][k][i,1], psp[j][k][i+1,1])
            INTyk += Iy(psp[j][k][i,0], psp[j][k][i+1,0], psp[j][k][i,1], psp[j][k][i+1,1])
            INTzk += Iz(psp[j][k][i,0], psp[j][k][i+1,0], psp[j][k][i,1], psp[j][k][i+1,1])
        intComX[j] += INTxk
        intComY[j] += INTyk
        intComZ[j] += INTzk
        intCom1[j] += AREAk
    Densicomunal=H[j]/intCom1[j]
    IX += Densicomunal*R*intComX[j]
    IY += Densicomunal*R*intComY[j]
    IZ += Densicomunal*R*intComZ[j]
    POBLA += H[j]
#COORDENADAS EN X Y Z
X = IX/POBLA
Y = IY/POBLA
Z = IZ/POBLA
R, THETA, PHI = cart2sph(X, Y, Z)
LON, LAT = phitheta2lonlat(PHI, THETA)
print("LON", LON)
print("LAT", LAT)
