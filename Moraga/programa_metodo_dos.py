# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 17:06:52 2016

@author: Sebastián
"""

import shapefile
import numpy 
import matplotlib.pyplot as plt
import triangle

# FUNCIONES
def separarPartes(puntos, inicioPartes):
    """ puntos es un arreglo de np; inicioPartes una lista de primer índice
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
#AREA DEL TRIANGULO
def areatriangulo(a,b,c):
    s=(a+b+c)/2
    area=numpy.sqrt((s-a)*(s)*(s-b)*(s-c))
    return area
# Código principal
sf = shapefile.Reader("division_comunal")
cod=['']*364
ACHILE=0
AreaC=['']*346
SumA=0;
Icomunax=0
Icomunay=0
Intsub=0
ICx=0
ICy=0
i=0

for i, comuna in enumerate(sf.shapeRecords()):
    if i !=321 and i !=324 and i!=331 and i!=332 and i !=345 and i!=335:    
        nombre = comuna.record[2]
        cod = comuna.record[6]
        inicioPartes = comuna.shape.parts
        puntos = numpy.array(comuna.shape.points)
        psp = separarPartes(puntos, inicioPartes)
        # Puntos separados por parte   
        IC=0
        Acomuna=0
        ICC=0 
        ICCx=0
        ICCy=0
        Area=0
        R=6371000
        for j, p in enumerate(psp):
            PP=p
            TRI=triangle.delaunay(p)
            T=numpy.empty(len(TRI)-1)
            Atri=0       
            Icomunax=0
            Icomunay=0                          
            for k in range(0,len(TRI)-1):                
                xx=numpy.abs(p[TRI[k,0],0]-p[TRI[k,1],0])
                yy=numpy.abs(p[TRI[k,0],1]-p[TRI[k,1],1])
                a=numpy.sqrt(xx*xx+yy*yy)            
                xx=numpy.abs(p[TRI[k,0],0]-p[TRI[k,2],0])
                yy=numpy.abs(p[TRI[k,0],1]-p[TRI[k,2],1])
                b=numpy.sqrt(xx*xx+yy*yy)            
                xx=numpy.abs(p[TRI[k,1],0]-p[TRI[k,2],0])
                yy=numpy.abs(p[TRI[k,1],1]-p[TRI[k,2],1])
                c=numpy.sqrt(xx*xx+yy*yy)            
                T[k]=areatriangulo(a,b,c)
                pcx=((p[TRI[k,0],0]+p[TRI[k,1],0]+p[TRI[k,2],0]))*T[k]/3
                pcy=((p[TRI[k,0],1]+p[TRI[k,1],1]+p[TRI[k,2],1]))*T[k]/3
                Atri=Atri+T[k]                
                ICCx=ICCx+pcx
                ICCy=ICCy+pcy
        AreaC[i]=Atri*H[i]
        SumA=Atri*H[i]+SumA
        ACHILE=ACHILE+Atri 
        Intsub=Atri*H[i]
        Icomunax=ICCx*H[i]
        Icomunay=ICCy*H[i]
        ICx=ICx+Icomunax
        ICy=ICy+Icomunay
        X=ICx/SumA
        Y=ICy/SumA





