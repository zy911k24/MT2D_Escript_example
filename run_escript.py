##############################################################################
#
# Copyright (c) 2003-2015 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development since 2012 by School of Earth Sciences
#
##############################################################################
#from __future__ import print_function

__copyright__="""Copyright (c) 2003-2015 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import csv
import codecs
import sys

import re
import esys.escript            as escript
import esys.escript.pdetools   as pdetools
import pandas as pd
from esys.downunder import *
from esys.escript import *
from esys.finley import MakeDomain
from esys.finley import ReadMesh
from esys.weipa import saveSilo
from esys.weipa import saveVTK
from esys.escript.linearPDEs import LinearPDE, SolverOptions, LinearSinglePDE, LinearPDESystem
import esys.escript.unitsSI as U
from esys.escript.pdetools import Locator
from esys.escript.pdetools import Projector
from esys.pycad.gmsh import *

import cmath
import logging
import cPickle as pickle
import numpy as np
import os
import matplotlib.pyplot as plt
import json, time
import esys.pycad              as pycad

try:
    import esys.finley         as finley
    HAVE_FINLEY = True
except ImportError:
    HAVE_FINLEY = False

HAVE_DIRECT = escript.hasFeature("PASO_DIRECT") or escript.hasFeature('trilinos')
#os.mkdir("data")
#os.mkdir("output")
#save_path= os.path.join("data") 

############################
###Parameters Setup###
############################

frequency = np.logspace(-1, -1, 1)  #frequencies
#frequency = (0.001, 0.1)
print ("frequency =", frequency)
PERIODES = 1./frequency
df = pd.DataFrame(frequency)   
df.to_csv('frequency.csv', index=False, header=False)
SIGMA_BASE=1/100. 

sigma_bg = 0.002  # background conductivity
MUE=4*np.pi*1e-7

############################
###Mesh Design###
############################

#mesh_file = "data/smaller_air_layered_horizontal_anisotropy.msh"
mesh_file = "flip_anisotropy.msh"

domain=finley.ReadGmsh(mesh_file, 2)

#############

def dphase(z16):

  z = np.array(z16)
  pom = np.arctan(np.imag(z)/np.real(z))
  pom = np.degrees(pom)  
  dphase = pom #np.angle(z, deg=True)
  return dphase

def app_res(z16, period):
  z = np.array(z16)
  MUE=4*np.pi*1e-7
  o=2*np.pi/period
  omega_mue = o*MUE
  app_res = abs(z)**2/omega_mue
  return app_res

filepath = 'MT_TAB_ROT.DAT'
app_res_xx = []
app_res_xy = []
app_res_yx = []
app_res_yy = []
phase_xx = []
phase_xy = []
phase_yx = []
phase_yy = []
a = []
stations = []
result = []
with open(filepath) as fp:
  for line in fp:
    result.append(line.strip().split())

for i in range( len(result) ):
  a.append(result[i][4:12])
 
df = pd.DataFrame(a)   
df.to_csv('impedance.csv', index=False, header=False) 

for i in range( len(result) ):
  stations.append(result[i][1])
 
df1 = pd.DataFrame(stations)   
df1.to_csv('stations.csv', index=False, header=False) 


Zxx_real = np.array(list(df[0]))
Zxx_imag = np.array(list(df[1]))
Zxx = [ complex(0.0,00) ]*len(result)
for n in range( len(result) ): 
  Zxx[n] = complex(float(Zxx_real[n]),float(Zxx_imag[n]))
  
df1 = pd.DataFrame(Zxx)   
df1.to_csv('Zxx.csv', index=False, header=False) 

Zxy_real = np.array(list(df[2]))
Zxy_imag = np.array(list(df[3]))
Zxy = [ complex(0.0,00) ]*len(result)
for n in range( len(result) ): 
  Zxy[n] = complex(float(Zxy_real[n]),float(Zxy_imag[n]))

df1 = pd.DataFrame(Zxy)   
df1.to_csv('Zxy.csv', index=False, header=False) 

Zyx_real = np.array(list(df[4]))
Zyx_imag = np.array(list(df[5]))
Zyx = [ complex(0.0,00) ]*len(result)
for n in range( len(result) ): 
  Zyx[n] = complex(float(Zyx_real[n]),float(Zyx_imag[n]))

df1 = pd.DataFrame(Zyx)   
df1.to_csv('Zyx.csv', index=False, header=False) 

Zyy_real = np.array(list(df[6]))
Zyy_imag = np.array(list(df[7]))
Zyy = [ complex(0.0,00) ]*len(result)
for n in range( len(result) ): 
  Zyy[n] = complex(float(Zyy_real[n]),float(Zyy_imag[n]))
  
df1 = pd.DataFrame(Zyy)   
df1.to_csv('Zyy.csv', index=False, header=False)   

for i in range (0,len(result)):
  o = app_res(Zxx[i], 10.)
  app_res_xx.append(o)
a=np.array(app_res_xx, dtype=float)
df = pd.DataFrame(a)   
df.to_csv('app_res_xx.csv', index=False, header=False) 

for i in range (0,len(result)):
  o = app_res(Zxy[i], 10.)
  app_res_xy.append(o)
b=np.array(app_res_xy, dtype=float)
df = pd.DataFrame(b)   
df.to_csv('app_res_xy.csv', index=False, header=False) 

for i in range (0,len(result)):
  o = app_res(Zyx[i], 10.)
  app_res_yx.append(o)
c=np.array(app_res_yx, dtype=float)
df = pd.DataFrame(c)   
df.to_csv('app_res_yx.csv', index=False, header=False) 

for i in range (0,len(result)):
  o = app_res(Zyy[i], 10.)
  app_res_yy.append(o)
d=np.array(app_res_yy, dtype=float)
df = pd.DataFrame(d)   
df.to_csv('app_res_yy.csv', index=False, header=False) 


for i in range (0,len(result)):
  o = dphase(Zxx[i])
  phase_xx.append(o)
e=np.array(phase_xx, dtype=float)
df = pd.DataFrame(e)   
df.to_csv('phase_xx.csv', index=False, header=False) 

for i in range (0,len(result)):
  o = dphase(Zxy[i])
  phase_xy.append(o)
f1=np.array(phase_xy, dtype=float)
df = pd.DataFrame(f1)   
df.to_csv('phase_xy.csv', index=False, header=False) 

for i in range (0,len(result)):
  o = dphase(Zyx[i])
  phase_yx.append(o)
g=np.array(phase_yx, dtype=float)
df = pd.DataFrame(g)   
df.to_csv('phase_yx.csv', index=False, header=False) 

for i in range (0,len(result)):
  o = dphase(Zyy[i])
  phase_yy.append(o)
h=np.array(phase_yy, dtype=float)
df = pd.DataFrame(h)   
df.to_csv('phase_yy.csv', index=False, header=False) 

X = np.array(stations, dtype=float)-238.
XX = 1000.*X
print "measurement", XX
print XX.tolist()
#########################
###Model Setup###
#########################
DIR='output'
##horizontal Anisotropy definition
rhox_layers = np.r_[1.0e+14, 80., 30., 80.]    #principal resistivity x
rhoy_layers = np.r_[1.0e+14, 100., 60., 100.]   #principal resistivity y
rhoz_layers = np.r_[1.0e+14, 70., 20., 70.]  #principal resistivity z
sigmax_layers = 1./rhox_layers
sigmay_layers = 1./rhoy_layers
sigmaz_layers = 1./rhoz_layers
##########
sigma_base = 100.
angle_base = 0.
ropbase = 100.
ubase = 0.

WIDTH=20000.

DEPTH=1.2

loc=Locator(Function(domain), [ (s,DEPTH) for s in XX.tolist()])
stations=np.array(loc.getX())

dfstations = pd.DataFrame(stations)
dfstations.to_csv('stations_escript.csv', index=False, header=False)

tags = escript.getTagNames(domain)
print ("tags:", tags)
##'tags:', ['air_layer', 'left', 'middle', 'right'])
#Anisotropy definition
def aniso_conductivity(rop1=ropbase, rop2=ropbase, rop3=ropbase, ustr=ubase, udip=ubase, usla=ubase):
  
  sgp1 = 1./rop1
  sgp2 = 1./rop2
  sgp3 = 1./rop3
  
  rstr = np.pi*ustr/180.0 
  rdip = np.pi*udip/180.0
  rsla = np.pi*usla/180.0

  sps = np.sin(rstr)
  cps = np.cos(rstr)

  sth = np.sin(rdip)
  cth = np.cos(rdip)

  sfi = np.sin(rsla)
  cfi = np.cos(rsla)
 
  pom1 = sgp1*np.square(cfi) + sgp2*np.square(sfi)
  pom2 = sgp1*np.square(sfi) + sgp2*np.square(cfi)
  pom3 = (sgp1-sgp2)*sfi*cfi

  s2ps = np.square(sps)
  c2ps = np.square(cps)
 
  s2th = np.square(sth) 
  c2th = np.square(cth)

  csps = cps*sps
  csth = cth*sth

  sigmaxx = pom1*c2ps+pom2*s2ps*c2th-2.*pom3*cth*csps+sgp3*s2th*s2ps
  sigmaxy = pom1*csps-pom2*c2th*csps+pom3*cth*(c2ps-s2ps)-sgp3*s2th*csps
  sigmaxz = -pom2*csth*sps+pom3*sth*cps+sgp3*csth*sps
 
  sigmayx = sigmaxy
  sigmayy = pom1*s2ps+pom2*c2ps*c2th+2.*pom3*cth*csps+sgp3*s2th*c2ps
  sigmayz = pom2*csth*cps+pom3*sth*sps-sgp3*csth*cps

  sigmazx = sigmaxz
  sigmazy = sigmayz
  sigmazz = pom2*s2th+sgp3*c2th

  print ("sigmaxx=", sigmaxx)
  print ("sigmaxy=", sigmaxy)
  print ("sigmaxz=", sigmaxz)

  print ("sigmayx=", sigmayx)
  print ("sigmayy=", sigmayy)
  print ("sigmayz=", sigmayz)

  print ("sigmazx=", sigmazx)
  print ("sigmazy=", sigmazy)
  print ("sigmazz=", sigmazz)

  return sigmaxx, sigmaxy, sigmaxz, sigmayx, sigmayy, sigmayz, sigmazx, sigmazy, sigmazz

########
sigma_fulltensor_medium1=aniso_conductivity(rop1=rhox_layers[1], rop2=rhoy_layers[1], rop3=rhoz_layers[1], ustr=55, udip=25., usla=30.)
sigma_fulltensor_medium2=aniso_conductivity(rop1=rhox_layers[2], rop2=rhoy_layers[2], rop3=rhoz_layers[2], ustr=30, udip=10., usla=20.)
sigma_fulltensor_medium3=aniso_conductivity(rop1=rhox_layers[3], rop2=rhoy_layers[3], rop3=rhoz_layers[3], ustr=30, udip=10., usla=20.)

sigma_xx_list=np.r_[1./1.0e+14, sigma_fulltensor_medium1[0], sigma_fulltensor_medium2[0], sigma_fulltensor_medium3[0]]
sigma_xy_list=np.r_[0., sigma_fulltensor_medium1[1], sigma_fulltensor_medium2[1], sigma_fulltensor_medium3[1]]
sigma_xz_list=np.r_[0., sigma_fulltensor_medium1[2], sigma_fulltensor_medium2[2], sigma_fulltensor_medium3[2]]
print "xx, xy, xz", sigma_xx_list, sigma_xy_list, sigma_xz_list
sigma_yx_list=np.r_[0., sigma_fulltensor_medium1[3], sigma_fulltensor_medium2[3], sigma_fulltensor_medium3[3]]
sigma_yy_list=np.r_[1./1.0e+14, sigma_fulltensor_medium1[4], sigma_fulltensor_medium2[4], sigma_fulltensor_medium3[4]]
sigma_yz_list=np.r_[0., sigma_fulltensor_medium1[5], sigma_fulltensor_medium2[5], sigma_fulltensor_medium3[5]]
print "yx, yy, yz", sigma_yx_list, sigma_yy_list, sigma_yz_list
sigma_zx_list=np.r_[0., sigma_fulltensor_medium1[6], sigma_fulltensor_medium2[6], sigma_fulltensor_medium3[6]]
sigma_zy_list=np.r_[0., sigma_fulltensor_medium1[7], sigma_fulltensor_medium2[7], sigma_fulltensor_medium3[7]]
sigma_zz_list=np.r_[1./1.0e+14, sigma_fulltensor_medium1[8], sigma_fulltensor_medium2[8], sigma_fulltensor_medium3[8]]
print "zx, zy, zz", sigma_zx_list, sigma_zy_list, sigma_zz_list

sigma_xx = Scalar(0, Function(domain))
for i in range( len(tags) ):
  
  sig = sigma_xx_list[i]
  sigma_xx += sig * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})

sigma_xx.expand()
#saveSilo("x", sigma_xx=sigma_xx)

sigma_xy = Scalar(0, Function(domain))
for i in range( len(tags) ):
  
  sig = sigma_xy_list[i]
  sigma_xy += sig * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})

sigma_xy.expand()

sigma_xz = Scalar(0, Function(domain))
for i in range( len(tags) ):
  
  sig = sigma_xz_list[i]
  sigma_xz += sig * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})

sigma_xz.expand()

sigma_yx = Scalar(0, Function(domain))
for i in range( len(tags) ):
  
  sig = sigma_yx_list[i]
  sigma_yx += sig * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})

sigma_yx.expand()

sigma_yy = Scalar(0, Function(domain))
for i in range( len(tags) ):
  
  sig = sigma_yy_list[i]
  sigma_yy += sig * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})

sigma_yy.expand()

sigma_yz = Scalar(0, Function(domain))
for i in range( len(tags) ):
  
  sig = sigma_yz_list[i]
  sigma_yz += sig * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})

sigma_yz.expand()

sigma_zx = Scalar(0, Function(domain))
for i in range( len(tags) ):
  
  sig = sigma_zx_list[i]
  sigma_zx += sig * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})

sigma_zx.expand()

sigma_zy = Scalar(0, Function(domain))
for i in range( len(tags) ):
  
  sig = sigma_zy_list[i]
  sigma_zy += sig * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})

sigma_zy.expand()

sigma_zz = Scalar(0, Function(domain))
for i in range( len(tags) ):
  
  sig = sigma_zz_list[i]
  sigma_zz += sig * escript.insertTaggedValues(escript.Scalar(0,escript.Function(domain)),**{ tags[i] : 1})

sigma_zz.expand()

#########################
###Forward Simulation###
#########################

def generate(fn, DOMAIN=domain, SIGMA_XX=SIGMA_BASE, SIGMA_XY=SIGMA_BASE, SIGMA_XZ=SIGMA_BASE, SIGMA_YX=SIGMA_BASE, SIGMA_YY=SIGMA_BASE, SIGMA_YZ=SIGMA_BASE, SIGMA_ZX=SIGMA_BASE, SIGMA_ZY=SIGMA_BASE, SIGMA_ZZ=SIGMA_BASE, locator=loc):
  
  domain=DOMAIN
  x=domain.getX()[0]     
  z=domain.getX()[1] 

  
  z=Solution(domain).getX()[1]
  x=Solution(domain).getX()[0]
  print ("x", x)
  print ("z", z)
  tags = escript.getTagNames(domain)
  print ("tags:", tags)
  
  proj = Projector(domain, reduce=False, fast=False)
  loc = locator
  
  print "Start  ..."

  sigma_xx=SIGMA_XX
  sigma_xy=SIGMA_XY
  sigma_xz=SIGMA_XZ
  sigma_yx=SIGMA_YX
  sigma_yy=SIGMA_YY
  sigma_yz=SIGMA_YZ
  sigma_zx=SIGMA_ZX
  sigma_zy=SIGMA_ZY
  sigma_zz=SIGMA_ZZ
   
  zb=FunctionOnBoundary(domain).getX()[1]
  print zb
  En0=whereZero(zb-inf(zb))  
  En1=whereZero(zb-sup(zb))
  TEdata={}
  rho_app_xydata={}
  xyphasedata={}
  TMdata={}
  rho_app_yxdata={}
  yxphasedata={}
  xxdata={}
  rho_app_xxdata={}
  xxphasedata={}
  yydata={}
  rho_app_yydata={}
  yyphasedata={}
  yyphasdata={}    
  n=0
  for p in PERIODES:
    o=2*np.pi/p
    print "frequency = ", 1./p 
    pde_tol=1.e-8
    #pde1=LinearPDESystem(domain, isComplex=True)   ####PDEs follows Guo 2018
    pde1=LinearPDE(domain, numEquations=2, numSolutions=2, isComplex=True)
    pde1.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
    pde1.getSolverOptions().setVerbosityOn()
    
    para_U = sigma_yy*sigma_zz-sigma_yz**2
    U_data = np.array(loc.getValue(para_U))   
    print "U_data", U_data
    para_K = (sigma_xy*sigma_yz-sigma_xz*sigma_yy)/para_U
    K_data = np.array(loc.getValue(para_K))   
    print "K_data", K_data
    para_B = (sigma_xz*sigma_yz-sigma_xy*sigma_zz)/para_U
    B_data = np.array(loc.getValue(para_B))   
    print "B_data", B_data
    para_C = sigma_xx+sigma_xy*para_B+sigma_zx*para_K
    C_data = np.array(loc.getValue(para_C))   
    print "C_data", C_data
    q=pde1.createCoefficient("q")
    r=pde1.createCoefficient("r")
    y=pde1.createCoefficient("y")
    A=pde1.createCoefficient('A') 
    B=pde1.createCoefficient('B') 
    D=pde1.createCoefficient('D') 
    C=pde1.createCoefficient('C') 

    q[0]=whereZero(z-inf(z))  #Ex  at top of the domain
    q[1]=whereZero(z)    #Hx
    r[0]=1.*whereZero(z-inf(z)) 
    r[1]=0.             #Hx
    

    A[0,0,0,0]=1./(1j*o*MUE)  #Ex
    A[0,1,0,1]=1./(1j*o*MUE)  #Ex
    A[1,0,1,0]=sigma_yy/para_U
    A[1,0,1,1]=sigma_zy/para_U
    A[1,1,1,0]=sigma_yz/para_U
    A[1,1,1,1]=sigma_zz/para_U  #Hx

    D[0,0]=para_C  #Ex
    D[1,1]=1j*o*MUE                                    #Hx   
    
    C[0,1,0]=para_K
    C[0,1,1]=-para_B
    B[1,0,0]=-para_K
    B[1,1,0]=para_B            
    
    pde1.setValue(A=-A, B=-B, C=C, D=D, q=q, r=r)#, y=y)
    mtfields1=pde1.getSolution()
    
    Ex1=mtfields1[0]
    Hx1=mtfields1[1]
    gEx1=grad(Ex1, where=loc.getFunctionSpace())
    gHx1=grad(Hx1, where=loc.getFunctionSpace())
    Ex1data = np.array(loc.getValue(Ex1))
    print "Ex1", Ex1data
    Hx1data = np.array(loc.getValue(Hx1))    
    Hy1=1.*gEx1[1]/(1j*o*MUE)
    Ey1=(sigma_yz/para_U)*gHx1[0]+(sigma_zz/para_U)*gHx1[1]+para_B*Ex1

    Ey1data = np.array(loc.getValue(Ey1))
    Hy1data = np.array(loc.getValue(Hy1)) 
    ###pde1 end
    pde2=LinearPDE(domain, numEquations=2, numSolutions=2, isComplex=True)
    pde2.getSolverOptions().setSolverMethod(SolverOptions.DIRECT)
    pde2.getSolverOptions().setVerbosityOn()

    para_U = sigma_yy*sigma_zz-sigma_yz**2
    U_data = np.array(loc.getValue(para_U))   
    print "U_data", U_data
    para_K = (sigma_xy*sigma_yz-sigma_xz*sigma_yy)/para_U
    K_data = np.array(loc.getValue(para_K))   
    print "K_data", K_data
    para_B = (sigma_xz*sigma_yz-sigma_xy*sigma_zz)/para_U
    B_data = np.array(loc.getValue(para_B))   
    print "B_data", B_data
    para_C = sigma_xx+sigma_xy*para_B+sigma_zx*para_K
    C_data = np.array(loc.getValue(para_C))   
    print "C_data", C_data
    q=pde2.createCoefficient("q")
    r=pde2.createCoefficient("r")
    y=pde2.createCoefficient("y")
    A=pde2.createCoefficient('A') 
    B=pde2.createCoefficient('B') 
    D=pde2.createCoefficient('D') 
    C=pde2.createCoefficient('C') 

    q[0]=whereZero(z-inf(z)) 
    q[1]=whereZero(z)   #Hx   should be at air and erath interface
    r[0]=0.*whereZero(z-inf(z))                 #Ex
    r[1]=1.*np.cos(np.pi)             #Hx   Bx(y,0)=B0         #Hx   Bx(y,0)=B0
    #y[0]=En0    
    #y[1]=En1

    A[0,0,0,0]=1./(1j*o*MUE)  #Ex
    A[0,1,0,1]=1./(1j*o*MUE)  #Ex
    A[1,0,1,0]=sigma_yy/para_U
    A[1,0,1,1]=sigma_zy/para_U
    A[1,1,1,0]=sigma_yz/para_U
    A[1,1,1,1]=sigma_zz/para_U  #Hx

    D[0,0]=para_C  #Ex
    D[1,1]=1j*o*MUE                                    #Hx   
    
    C[0,1,0]=para_K
    C[0,1,1]=-para_B
    B[1,0,0]=-para_K
    B[1,1,0]=para_B         
    
    pde2.setValue(A=-A, B=-B, C=C, D=D, q=q, r=r)#, y=y)
    mtfields2=pde2.getSolution()
    
    Ex2=mtfields2[0]
    Hx2=mtfields2[1]
    gEx2=grad(Ex2, where=loc.getFunctionSpace())
    gHx2=grad(Hx2, where=loc.getFunctionSpace())
    Ex2data = np.array(loc.getValue(Ex2))
    Hx2data = np.array(loc.getValue(Hx2))    
    Hy2=1.*gEx2[1]/(1j*o*MUE)
    Ey2=(sigma_yz/para_U)*gHx2[0]+(sigma_zz/para_U)*gHx2[1]+para_B*Ex2

    Ey2data = np.array(loc.getValue(Ey2))
    Hy2data = np.array(loc.getValue(Hy2))
    det = Hx1*Hy2-Hx2*Hy1
    Zxx= (Ex1*Hy2-Ex2*Hy1)/det
    Zxy= (Ex2*Hx1-Ex1*Hx2)/det
    Zyx= (Ey1*Hy2-Ey2*Hy1)/det
    Zyy= (Ey2*Hx1-Ey1*Hx2)/det

    TEdata[p]=np.array(loc.getValue(Zxy))    
    TMdata[p]=np.array(loc.getValue(Zyx))
    print ("Zxy %s of periods %s"%(TEdata[p], p))
    print ("Zyx %s of periods %s"%(TMdata[p], p))
    xxdata[p]=np.array(loc.getValue(Zxx))
    yydata[p]=np.array(loc.getValue(Zyy))
    omega_mue=o*MUE
    TEdata_abs_square = abs(TEdata[p])**2
    TMdata_abs_square = abs(TMdata[p])**2
    xxdata_abs_square = abs(xxdata[p])**2
    yydata_abs_square = abs(yydata[p])**2

    rho_app_xydata[p]=TEdata_abs_square/omega_mue
    #xyphase=atan2(Zxy.imag(),Zxy.real())/np.pi*180.
    #xyphasedata[p]=np.array(loc(xyphase))
    xyphasedata[p]=np.angle(TEdata[p], deg=True)
    rho_app_yxdata[p]=TMdata_abs_square/omega_mue
    yxphase=atan2(Zyx.imag(),Zyx.real())/np.pi*180.
    yxphasedata[p]=np.array(loc(yxphase))

    rho_app_xxdata[p]=xxdata_abs_square/omega_mue
    xxphase=atan2(Zxx.imag(),Zxx.real())/np.pi*180.
    xxphasedata[p]=np.array(loc(xxphase))
    
    rho_app_yydata[p]=yydata_abs_square/omega_mue
    yyphase=atan2(Zyy.imag(),Zyy.real())/np.pi*180.
    yyphasedata[p]=np.array(loc(yyphase))
    yyphasdata[p]=np.array(np.delete(yyphasedata[p], [0, 1]))
    n+=1  

  df = pd.DataFrame(Ex1data)   
  df.to_csv('Ex1data.csv', index=False, header=False) 
  df = pd.DataFrame(Hx1data)   
  df.to_csv('Hx1data.csv', index=False, header=False) 
  df = pd.DataFrame(Ey1data)   
  df.to_csv('Ey1data.csv', index=False, header=False) 
  df = pd.DataFrame(Hy1data)   
  df.to_csv('Hy1data.csv', index=False, header=False) 
  df = pd.DataFrame(Ex2data)   
  df.to_csv('Ex2data.csv', index=False, header=False) 
  df = pd.DataFrame(Hx2data)   
  df.to_csv('Hx2data.csv', index=False, header=False) 
  df = pd.DataFrame(Ey2data)   
  df.to_csv('Ey2data.csv', index=False, header=False) 
  df = pd.DataFrame(Hy2data)   
  df.to_csv('Hy2data.csv', index=False, header=False) 

  df = pd.DataFrame(TEdata)   
  df.to_csv('TEdata.csv', index=False, header=False) 
  df = pd.DataFrame(TMdata)   
  df.to_csv('TMdata.csv', index=False, header=False)
  df = pd.DataFrame(rho_app_xydata)
  df.to_csv('rho_app_xydata.csv', index=False, header=False) 
  df = pd.DataFrame(xyphasedata)
  df.to_csv('xyphasedata.csv', index=False, header=False)
  df = pd.DataFrame(rho_app_yxdata)
  df.to_csv('rho_app_yxdata.csv', index=False, header=False) 
  df = pd.DataFrame(yxphasedata)
  df.to_csv('yxphasedata.csv', index=False, header=False) 

  df = pd.DataFrame(xxdata)   
  df.to_csv('xxdata.csv', index=False, header=False) 
  df = pd.DataFrame(yydata)   
  df.to_csv('yydata.csv', index=False, header=False) 
  df = pd.DataFrame(rho_app_xxdata)
  df.to_csv('rho_app_xxdata.csv', index=False, header=False) 

  df = pd.DataFrame(rho_app_yydata)
  df.to_csv('rho_app_yydata.csv', index=False, header=False) 
  df = pd.DataFrame(xxphasedata)
  df.to_csv('xxphasedata.csv', index=False, header=False) 
  df = pd.DataFrame(yyphasedata)
  df.to_csv('yyphasedata.csv', index=False, header=False) 

  return PERIODES, rho_app_xydata[p], xyphasedata[p], rho_app_yxdata[p], yxphasedata[p], rho_app_xxdata[p], xxphasedata[p], rho_app_yydata[p], yyphasdata[p]
###plot


A = generate("data/horizontal_anisotropy", DOMAIN=domain, SIGMA_XX=sigma_xx, SIGMA_XY=sigma_xy, SIGMA_XZ=sigma_xz, SIGMA_YX=sigma_yx, SIGMA_YY=sigma_yy, SIGMA_YZ=sigma_yz, SIGMA_ZX=sigma_zx, SIGMA_ZY=sigma_zy, SIGMA_ZZ=sigma_zz, locator=loc)

X1=np.array( loc.getX() )[:,0]/1000.	
print "X1", X1
X2=np.array(np.delete(X1, [0, 1]))
print "X2", X2

ylbl0xx = r'Apparent Resistivity $\rho_{xx}$ $(\Omega\cdot\,m)$'
ylbl1xx = r'Phase $\phi_{xx}$ $(^{\circ})$'
ylbl0xy = r'Apparent Resistivity $\rho_{xy}$ $(\Omega\cdot\,m)$'
ylbl1xy = r'Phase $\phi_{xy}$ $(^{\circ})$'
ylbl0yx = r'Apparent Resistivity $\rho_{yx}$ $(\Omega\cdot\,m)$'
ylbl1yx = r'Phase $\phi_{yx}$ $(^{\circ})$'
ylbl0yy = r'Apparent Resistivity $\rho_{yy}$ $(\Omega\cdot\,m)$'
ylbl1yy = r'Phase $\phi_{yy}$ $(^{\circ})$'
xlbl1 = 'Y(km)'
fsize = 24
f,ax = plt.subplots(2, figsize=(6,14), dpi=300)
f.subplots_adjust(top=0.9)  # Little extra space for 'suptitle'
f.suptitle('')           # This is actually the plot-title
my_x_ticks = np.arange(-16,18,4)
my_y_ticks = np.arange(-60,-25,5)
my_y_ticks1 = np.arange(-90,120,30)
my_y_ticks2 = np.arange(-90,30,15)
ax[0].scatter(X, a, marker='D', c='b', label="FD")
ax[0].plot(X1, A[5], color="red", linewidth=1.5, linestyle="-", label="Escript")
ax[0].legend(loc='upper right', fontsize=22)
ax[0].set_xlabel(xlbl1, fontsize=fsize)
ax[0].set_ylabel(ylbl0xx, fontsize=fsize)
ax[0].set_xlim([-12,12])
ax[0].set_yscale('log')
ax[0].set_ylim([0.01,100])
ax[0].set_xticks(my_x_ticks)
ax[0].xaxis.set_tick_params(labelsize=22)
ax[0].yaxis.set_tick_params(labelsize=22)

ax[1].scatter(X, e, marker='D', c='b', label="FD")
ax[1].plot(X1, A[6], color="red", linewidth=1.5, linestyle="-", label="Escript")
ax[1].set_xlabel(xlbl1, fontsize=fsize) 
ax[1].set_ylabel(ylbl1xx, fontsize=fsize)
ax[1].set_xlim([-12,12])
ax[1].set_ylim([-90, 0])
ax[1].set_xticks(my_x_ticks)
ax[1].set_yticks(my_y_ticks2)
ax[1].xaxis.set_tick_params(labelsize=22)
ax[1].yaxis.set_tick_params(labelsize=22)
plt.savefig("new_Horizontalxx.png", bbox_inches='tight')

f,ax = plt.subplots(2, figsize=(6,14), dpi=300)#A[1],[2], b, f1 

ax[0].scatter(X, b, marker='D', c='b', label="FD")
ax[0].plot(X1, A[1], color="red", linewidth=1.5, linestyle="-", label="Escript")
ax[0].legend(loc='upper right', fontsize=22)
ax[0].set_xlabel(xlbl1, fontsize=fsize)
ax[0].set_ylabel(ylbl0xy, fontsize=fsize)
ax[0].set_xlim([-12,12])
ax[0].set_yscale('log')
ax[0].set_ylim([1,1000])
ax[0].set_xticks(my_x_ticks)
ax[0].xaxis.set_tick_params(labelsize=22)
ax[0].yaxis.set_tick_params(labelsize=22)

ax[1].scatter(X, f1, marker='D', c='b', label="FD")
ax[1].plot(X1, A[2], color="red", linewidth=1.5, linestyle="-", label="Escript")
ax[1].set_xlabel(xlbl1, fontsize=fsize) 
ax[1].set_ylabel(ylbl1xy, fontsize=fsize)
ax[1].set_xlim([-12,12])
ax[1].set_ylim([-60, -30])
ax[1].set_xticks(my_x_ticks)
ax[1].set_yticks(my_y_ticks)
ax[1].xaxis.set_tick_params(labelsize=22)
ax[1].yaxis.set_tick_params(labelsize=22)
plt.savefig("new_Horizontalxy.png", bbox_inches='tight')

f,ax = plt.subplots(2, figsize=(6,14), dpi=300)#c, #g

ax[0].scatter(X, c, marker='D', c='b', label="FD")
ax[0].plot(X1, A[3], color="red", linewidth=1.5, linestyle="-", label="Escript")
ax[0].legend(loc='upper right', fontsize=22)
ax[0].set_xlabel(xlbl1, fontsize=fsize)
ax[0].set_ylabel(ylbl0yx, fontsize=fsize)
ax[0].set_xlim([-12,12])
ax[0].set_yscale('log')
ax[0].set_ylim([1,1000])
ax[0].set_xticks(my_x_ticks)
ax[0].xaxis.set_tick_params(labelsize=22)
ax[0].yaxis.set_tick_params(labelsize=22)

ax[1].scatter(X, g, marker='D', c='b', label="FD")
ax[1].plot(X1, A[4], color="red", linewidth=1.5, linestyle="-", label="Escript")
ax[1].set_xlabel(xlbl1, fontsize=fsize) 
ax[1].set_ylabel(ylbl1yx, fontsize=fsize)
ax[1].set_xlim([-12,12])
ax[1].set_ylim([-60, -30])
ax[1].set_xticks(my_x_ticks)
ax[1].set_yticks(my_y_ticks)
ax[1].xaxis.set_tick_params(labelsize=22)
ax[1].yaxis.set_tick_params(labelsize=22)
plt.savefig("new_Horizontalyx.png", bbox_inches='tight')

f,ax = plt.subplots(2, figsize=(6,14), dpi=300)#d,h

ax[0].scatter(X, d, marker='D', c='b', label="FD")
ax[0].plot(X1, A[7], color="red", linewidth=1.5, linestyle="-", label="Escript")
ax[0].legend(loc='upper right', fontsize=22)
ax[0].set_xlabel(xlbl1, fontsize=fsize)
ax[0].set_ylabel(ylbl0yy, fontsize=fsize)
ax[0].set_xlim([-12,12])
ax[0].set_yscale('log')
ax[0].set_ylim([0.01,100])
ax[0].set_xticks(my_x_ticks)
ax[0].xaxis.set_tick_params(labelsize=22)
ax[0].yaxis.set_tick_params(labelsize=22)

ax[1].scatter(X, h, marker='D', c='b', label="FD")
ax[1].plot(X2, A[8], color="red", linewidth=1.5, linestyle="-", label="Escript")
ax[1].set_xlabel(xlbl1, fontsize=fsize) 
ax[1].set_ylabel(ylbl1yy, fontsize=fsize) 
ax[1].set_xlim([-12,12])
ax[1].set_ylim([-90, 90])
ax[1].set_xticks(my_x_ticks)
ax[1].set_yticks(my_y_ticks1)
ax[1].xaxis.set_tick_params(labelsize=22)
ax[1].yaxis.set_tick_params(labelsize=22)
plt.savefig("new_Horizontalyy.png", bbox_inches='tight')

