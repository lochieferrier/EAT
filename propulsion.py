# Propulsion simulation for MIT Electric Aircraft Team
# At each point in the flow, mach number M, velocity V, static and dynamic temperature T_0 and T_t, static and dynamic temperature P_0 and P
# as well as density rho are computed.
# input parameters are compressor face pressure ratios, areas, density and velocity ranges of interest
# output parameters are V-rho-Thrust curves
# different types of transitions include front normal shock, isentropic transition, pressure ratio non-isentropic, mach de laval re expansion
#-----------
# Stations
# 0 (infinity/free stream)
# 1 (post mach diffuser)
# 2 (first axial compressor face)
# 3 (post first compressor (constant A with face)
# 4 (Face of centrifugal compressor)
# 5 (Post centrifugal compressor)
# 6 (nozzle exit)

#All units are in SI, mks standard with T in K and P in Pa

#nomenclature
#gamma is specific heat ratio, which is 1.4 for air
#mdot is mass flow, which is set by inlet density, velocity and area

import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams


#normal shock equations
def calcNormalShockMach(mach,gamma):
	postShockMach = math.sqrt(((gamma-1)*mach**2 + 2)/((2*gamma*mach**2)-(gamma-1)))
	return postShockMach

def calcNormalShockStaticPressure(staticpressure,mach,gamma):
	postShockStaticPressure = staticpressure*(((2*gamma*mach**2)-(gamma-1))/(gamma+1))
	return postShockStaticPressure

def calcNormalShockDynamicPressure(dynamicPressure,mach,gamma):
	postShockDynamicPressure = dynamicPressure*((((gamma+1)*mach**2)/((gamma-1)*mach**2 + 2))**(gamma/(gamma-1)) * ((gamma+1)/(2*gamma*mach**2-(gamma-1)))**(1/gamma-1))
	return postShockDynamicPressure

def calcNormalShockDensity(density,mach,gamma):
	postShockDensity = density*(((gamma+1)*mach**2)/((gamma-1)*mach**2 + 2))
	return postShockDensity

def calcNormalShockStaticTemperature(staticTemperature,mach,gamma):
	postShockStaticTemperature = staticTemperature * (((2*gamma*mach**2-(gamma-1))*((gamma-1)*mach**2 + 2))/(((gamma+1)**2)*mach**2))
	return postShockStaticTemperature

def calcIsentropicAreaVelocity(areaInitial,areaFinal,machInitial,vInitial):
	dA = areaFinal - areaInitial
	V = vInitial + (dA * vInitial)/((machInitial**2 - 1)*areaInitial)
	if (machInitial == 1):
		if (dA < 0):
			V=V
		if (dA > 0):
			V = vInitial + (dA * vInitial)/((machInitial**2 - 1)*areaInitial)
	return V

def calcIsentropicTotalPressure(refTotalPressure,mach,gamma):
	isoTotalPressure = refTotalPressure/((1+((gamma-1)/2)*mach**2)**(gamma/(gamma-1)))
	return isoTotalPressure
def calcIsentropicTotalTemperature(refTotalTemperature,mach,gamma):
	isoTotalTemperature = refTotalTemperature/(1+((gamma-1)/2)*M**2)
	return isoTotalTemperature
def calcIsentropicDensity(refDensity,mach,gamma):
	isoDensity = refDensity/((1+((gamma-1)/2)*mach**2)**(1/(gamma-1)))
	return isoDensity

#reference for International Standard Atmosphere (ISA) block calculation: https://www.grc.nasa.gov/www/K-12/airplane/atmos.html
#corrected to K and Pa units, from original degC and K-Pa units

def calcISAStaticTemperature(altitude):
	if altitude > 25000:
		staticTemp = -131.21 + 0.00299 * altitude
	if 11000 < altitude <= 25000:
		staticTemp = -56.46
	if altitude <= 11000:
		staticTemp = 15.04 - 0.00649*altitude
	#convert to Kelvin from Celcius
	staticTemp = staticTemp + 273.15
	return staticTemp

def calcISAStaticPressure(altitude):
	staticTemp = calcISAStaticTemperature(altitude)
	if altitude > 25000:
		staticPressure = 2.488 * ((staticTemp/216.6)**(-11.388))
	if 11000 < altitude <= 25000:
		staticPressure = 22.65 * math.exp(1.73 - 0.000157*altitude) 
	if altitude <= 11000:
		staticPressure = 101.29 * (staticTemp/288.08)**5.256
	#convert to Pa from k-Pa
	staticPressure = staticPressure * 1000
	return staticPressure

def calcISADensity(altitude):
	staticPressure = calcISAStaticPressure(altitude)
	#Convert from Pa to K-Pa
	staticPressure = staticPressure / 1000
	rhoISA = staticPressure/(0.2869*(calcISAStaticTemperature(altitude)))
	return rhoISA

def calcSpeedofSound(gamma,gasConstantR,staticTemperature):
	speedOfSound = math.sqrt(gamma*gasConstantR*staticTemperature)
	return speedOfSound

#compressor calculations
def calcAxialTempRatio(axialDFactor, axialSolidity,axialAlpha1,machInlet,gamma):
	#first need to determine alpha2
	alphaTwo = calcAxialAlphaTwo(axialDFactor,axialSolidity,axialAlpha1)
	tempRatio = ((((gamma-1)*mach**2)/(1+0.5*(gamma-1)*mach**2)) * ((math.cos(axialAlpha1)**2)/(math.cos(calcAxialAlphaTwo)**2) - 1)) + 1
	return tempRatio
def calcAxialPressureRatio(axialEc,gamma,tempRatio):
	#first need to obtain temperature ratio, from which we can derive pressure ratio
	pressureRatio = (tempRatio)**((gamma*axialEc)/(gamma-1))
	return pressureRatio
def calcAxialMachRatio(tempRatio,gamma,mach):
	machRatio = math.sqrt(1/(tempRatio*(1+(gamma-1)*0.5*mach**2) - ((gamma-1)*0.5)*mach**2))
	return machRatio
def calcAxialAlphaTwo(D, s, alphaOne):
	#G is capital gamma
	G = (2*solidity + sin(alphaOne))/math.cos(alphaOne)
	alphaTwo = math.acos((2*s*(1-D)*G + math.sqrt(G**2 + 1 - (4s**2)*(1-D)**2))/(G**2 + 1))
	return alphaTwo

def calcCentrifugalPressureRatio(Cp,centrifugalU,centrifugalSlipFactor,inletT_t,centrifugalEc):
	pressureRatio = (1+ (centrifugalSlipFactor*centrifugalU**2)/(Cp*inletT_t))**((gamma*centrifugalEc)/(gamma-1))
	return pressureRatio

def calcCentrifugalTempRatio(Cp,centrifugalU,centrifugalSlipFactor,inletT_t):
	tempRatio = (1+(centrifugalSlipFactor*centrifugalU**2)/(Cp*inletT_t))
	return tempRatio

def display():
	global M
	print("mach number: ",M)

def model():
	if station == 1:
		#here we need to determine whether a shock wave is present at the front face
		if M > 1:
			mach = calcNormalShockMach(mach,gamma)
			rho = calcNormalShockDensity(rho,M,gamma)
			T = calcNormalShockStaticTemperature(staticTemperature,mach,gamma)
			T_t = T_t #normal shock has no effect on total temperature, ideally
			P = calcNormalShockStaticPressure(P,mach,gamma)
		if M < 1:
			#continue as normal as it is subsonic?
	if station

#define inlet conditions (single case test is at 10,000 m @M1.1, ISA)
gamma = 1.4
gasConstantR = 286
Cp = 1004 #J/(kg*K) specific heat

T_0 = calcISAStaticTemperature(10000) #static temperature
P_0 = calcISAStaticPressure(10000) #static pressure
a = calcSpeedofSound(gamma,gasConstantR,T_0)	#speed of sound
M = 1.1
T_t = 291

P_t = 0.34
A = 0.2
V = a*M
rho = 0.2
mdot = rho*A*V

station = 1

#compressor 1 setup (axial)
axialDFactor = 0.6 #diffusion factor
axialSolidity = 1
axialEc = 0.88 #polytropic efficiency, baselined from 1980s level
axialAlpha1 = 4
axialAlpha1 = math.radians(axialAlpha1)

#compressor 2 setup (centrifugal)
centrifugalU = 0.6 #rotor speed
centrifugalSlipFactor = 0.9 #slip factor with 20 vanes, as approximated by 1-(2/n) where n is #vanes
centrifugalEc = 0.8

while station < 0:
	model()
	store()
	station = station + 1
	print(station)
display()