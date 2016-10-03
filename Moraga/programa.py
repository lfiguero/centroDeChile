# coding: utf-8

# PREÁMBULO
import shapefile
import numpy
import matplotlib.pyplot as plt
import triangle
import os
import utm

# FUNCIÓN QUE SEPARA PARTES DEL ARCHIVO SHP
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

# FUNCIONES PARA INTEGRAR POR REGLA DE SIMPSON ANIDADA
def I3(x,c,d,a):
    x=x
    R=6371000
    h=((x-c)/6)*((R*R)*numpy.cos(c)*(numpy.cos(x)-numpy.sin(x))+4*(numpy.cos(x)-numpy.sin(x))+((R*R)*numpy.cos((x-a)*d+c)*(numpy.cos(x)-numpy.sin(x))))
    return h
def I4(a,b,c,d):
    a=a
    b=b
    c=c
    I=((b-a)/6)*(I3(a,c,d,a)+4*I3((a+b)/2,c,d,a)+I3(b,c,d,a))
    return I
    ########################
def I1(x,c,d,a):
    R=6371000
    h=((x-c)/6)*((R*R)*numpy.cos(c)*(numpy.cos(c)*numpy.cos(x)+numpy.cos(c)*numpy.sin(x)-numpy.sin(c))+4*((R*R)*numpy.cos((c+(x-a)*d+c)/2)*(numpy.cos((c+(x-a)*d+c)/2)*numpy.cos(x)+numpy.cos((c+(x-a)*d+c)/2)*numpy.sin(x)-numpy.sin((c+(x-a)*d+c)/2)))+((R*R)*numpy.cos((x-a)*d+c)*(numpy.cos((x-a)*d+c)*numpy.cos(x)+numpy.cos((x-a)*d+c)*numpy.sin(x)-numpy.sin((x-a)*d+c))))
    return h

def I2(a,b,c,d):
    a=a
    b=b
    c=c
    I=((b-a)/6)*(I1(a,c,d,a)+4*I1((a+b)/2,c,d,a)+I1(b,c,d,a))
    return I
    ###########################
def II1(x,c,d,b):
    x=x
    R=6371000
    h=((x-c)/6)*((R*R)*numpy.cos(c)*(numpy.cos(c)*numpy.cos(x)+numpy.cos(c)*numpy.sin(x)-numpy.sin(c))+4*((R*R)*numpy.cos((c+(b-x)*d+c)/2)*(numpy.cos((c+(b-x)*d+c)/2)*numpy.cos(x)+numpy.cos((c+(b-x)*d+c)/2)*numpy.sin(x)-numpy.sin((c+(b-x)*d+c)/2)))+((R*R)*numpy.cos((b-x)*d+c)*(numpy.cos((b-x)*d+c)*numpy.cos(x)+numpy.cos((b-x)*d+c)*numpy.sin(x)-numpy.sin((b-x)*d+c))))
    return h


def II2(a,b,c,d):
    a=a
    b=b
    c=c
    I=((b-a)/6)*(II1(a,c,d,b)+4*II1((a+b)/2,c,d,b)+II1(b,c,d,b))
    return I
    ##########################

def II3(x,c,d,b):
    x=x
    R=6371000
    h=((x-c)/6)*((R*R)*numpy.cos(c)*(numpy.cos(x)-numpy.sin(x))+4*(numpy.cos(x)-numpy.sin(x))+((R*R)*numpy.cos((b-x)*d+c)*(numpy.cos(x)-numpy.sin(x))))
    return h

def II4(a,b,c,d):
    I=((b-a)/6)*(II3(a,c,d,b)+4*II3((a+b)/2,c,d,b)+II3(b,c,d,b))
    return I
    #########################
def IP(a,b,c,d):
    R=6371000
    I=((b-a)/6)*((R*R)*(numpy.sin((c-a)*d+c)-numpy.sin(c))+4*((R*R)*(numpy.sin(((a+b)*0.5-a)*d+c)-numpy.sin(c)))+((R*R)*(numpy.sin((b-a)*d+c)-numpy.sin(c))))
    return I
def IP2(a,b,c,d):
    R=6371000
    I=((b-a)/6)*((R*R)*(numpy.sin((c-a)*d+c)-numpy.sin(c))+4*((R*R)*(numpy.sin(((a+b)*0.5-a)*d+c)-numpy.sin(c)))+((R*R)*(numpy.sin((b-a)*d+c)-numpy.sin(c))))
    return I

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
SUBCHILE=0
comunasNaN = [321, 324, 331, 332, 345, 335]
comunasFueraDeRango = [37, 8, 9, 11, 14, 19, 22, 23, 28, 30, 31, 37, 142, 229, 280, 294, 325, 329, 330, 276, 281]
for i, comuna in enumerate(sf.shapeRecords()):
    #COMUNAS
    if i not in comunasNaN and i not in comunasFueraDeRango:
        nombre = comuna.record[2]
        cod = comuna.record[6]
        inicioPartes = comuna.shape.parts
        puntos = numpy.array(comuna.shape.points)
        psp = separarPartes(puntos, inicioPartes)
        # Puntos separados por parte
        Acomuna=0
        ICCx=0
        ICCy=0
        Icomunax=0
        Icomunay=0
        SUBCOMUNA=0
        for j, p in enumerate(psp):
            for l in range(0,len(p)):
                if abs(p[l,1])>100000 and abs(p[l,1])< 9999999 and abs(p[l,0])>100000 and abs(p[l,0])< 9999999:
                    #TRANSFORMACION A UTM
                    [p[l,0],p[l,1]]=utm.to_latlon(abs(p[l,0]),abs(p[l,1]), 19, 'k')
                else:
                    #ALARMA POR SI ALGUNA COMUNA SE SALE DEL RANGO DE UTM
                    print("----------------------------------")
                    print("ALARMA1")
                    print([i, j])
                    print("----------------------------------")
             #TRIANGULARIZACION
            TRI=triangle.delaunay(p)
            T=numpy.empty(len(TRI)-1)
            for k in range(0,len(TRI)-1):
                #CALCULO DE INTEGRAL POR MEDIO DE SEPARACION DE PUNTOS EN TRIANGULOS
                x1=p[TRI[k,0],0]
                y1=p[TRI[k,0],1]
                x2=p[TRI[k,1],0]
                y2=p[TRI[k,1],1]
                x3=p[TRI[k,2],0]
                y3=p[TRI[k,2],1]
                x4=p[TRI[k,1],0]
                y4=p[TRI[k,1],1]
                a=x1
                b=x4
                c=y4
                d=y3
                a2=x2
                b2=x4
                c2=y4
                d2=y3
                pcx=I4(a,b,c,d)
                pcy=I2(a,b,c,d)
                pcx2=-II4(a2,b2,c2,d2)
                pcy2=-II2(a2,b2,c2,d2)
                AAA=IP(a,b,c,d)
                AAA2=-IP2(a2,b2,c2,d2)
                SUBCOMUNA=SUBCOMUNA+AAA+AAA2
                ICCx=ICCx+pcx+pcx2
                ICCy=ICCy+pcy+pcy2
        #RESULTADOS
        SumA=SUBCOMUNA*H[i]+SumA
        Icomunax=ICCx*H[i]
        Icomunay=ICCy*H[i]
        ICx=ICx+Icomunax
        ICy=ICy+Icomunay
#PARA OBTENER EL RESULTADO
# \theta=ICx/SumA
# \Phi  =ICy/SumA






