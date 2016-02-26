#Goal is to simulate the entire flight using input paramaters, in order to get the most useful data possible
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
rcParams['axes.labelsize'] = 12
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 12
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
#environment vars
g = 9.8 #m/s/s
#Record characteristics
target_v = 360 #m/s

#runway characteristics
runway_alt = 2000 #m

#aircraft design characteristics
mass = 100 #kg
Cd = 0.02 #unitless
Cl = 0.2 #unitless
Ref_A = 3.5 #m^2
E_capacity = 32400000 #J
E_endMachFraction = 0.2 #unitless
E_remaining = E_capacity #J
max_power = 200000 #W
max_power_duration = 120 #s
motorThrustToPowerCurve = 0.010 #N/W
v_climb = 100
climbAngle = 15 #degrees
glideAngle = 15 #degrees

#aircraft flight variables
aircraftOutsideTemp = 15
v_airspeed = 0 #m/s
pitchAngle = 0 #radians
v_y = 0 #m/s
v_x = 0 #m/s
final_v = 330 #m/s
thrust_max = 2300 #N
timestep = 0.2 #s
t_current = 0 #s
aircraft_rho = 1.225 #kg/m^3
aircraft_alt = runway_alt
power = max_power
aircraft_lift = 0 #N

#staging
climbPower = 80000 #W

flightStage = 1
levelRho = 0
#data record
v_airspeedList = []
timeList = []
aircraft_altList = []
aircraftOutsideTempList = []

def model():

    global v_airspeed, t_current, aircraft_alt, aircraft_lift, aircraft_rho, v_y,E_remaining,aircraftOutsideTemp
    #takeoff
    if flightStage == 1:
        power = max_power
    #climb
    if flightStage == 2:
        power = climbPower
    #level    
    if flightStage == 3:
        power = max_power
    #glide
    if flightStage == 4:
        power = 0
        #hard set glide airspeed
        #v_airspeed = math.sqrt((math.sin(math.atan(Cd/Cl))*mass*g)/(0.5*aircraft_rho*Cd*Ref_A))

    #force calculation
    aircraft_rho = calculateRho(aircraft_alt)
    #aircraftOutsideTemp = calculateOutsideTemp(aircraft_alt)
    
    thrust = power * motorThrustToPowerCurve;
    if flightStage == 4:
        thrust = math.sin(glideAngle*0.0174533)*mass*g

    drag = 0.5 * Cd * Ref_A * aircraft_rho * pow(v_airspeed,2)
    lift = 0.5 * Cl * Ref_A * aircraft_rho * pow(v_airspeed,2)
    aircraft_lift = lift
    #velocity update for accelerating case
    v_airspeed = v_airspeed + ((thrust - drag)/mass)*timestep
    if flightStage ==2:
        #v_y = ((thrust-drag)/(mass*g) )*v_airspeed
        v_y = math.sin(climbAngle*0.0174533)*v_airspeed
        print(v_y)
    else:
        v_y = 0
        if flightStage ==4:
            v_y = -v_airspeed * (Cd/Cl)
    
    
    #position update
    aircraft_alt = aircraft_alt + v_y*timestep
    print('alt: ',aircraft_alt)
    print('v_airspeed: ',v_airspeed)
    print('v_y: ',v_y)
    print('levelRho: ',levelRho)

    #Energy usage
    E_remaining = E_remaining - power*timestep

    #timestep
    t_current = t_current + timestep
    
def store():
    global v_airspeedList,timeList
    v_airspeedList.append(v_airspeed)
    timeList.append(t_current)
    aircraft_altList.append(aircraft_alt)
    
def check():
    global flightStage
    #takoff roll
    if flightStage == 1:
        if aircraft_lift >= mass*g:
            flightStage = 2
    #climb     
    if flightStage == 2:
        if aircraft_rho <= levelRho:
            flightStage = 3
    if flightStage == 3:
        if E_remaining/E_capacity <= E_endMachFraction:
            flightStage = 4
    if flightStage == 4:
        if aircraft_alt<=runway_alt:
            flightStage = 5
    
def calculateLevelRho():
    global thrust_max, Ref_A, Cd, target_v 
    rho_min = thrust_max/(0.5*Cd*Ref_A*pow(target_v,2))
    return rho_min
    
def calculateRho(altAboveSeaLevel):
    #only valid up to 11000 m
    #src: https://www.grc.nasa.gov/www/K-12/airplane/atmosmet.html
    T = 15.04 - 0.00649*altAboveSeaLevel
    p = 101.29 * pow((T+273.1)/(288.08),5.526)
    rho = p/((0.2869*(T+273.1)))
    return rho

def display():
    print(t_current)
    print('Remaining capacity:' ,E_remaining/E_capacity)
    plt.plot(timeList, aircraft_altList)
    plt.ylabel('alt (m)')
    plt.xlabel('Time (s)')
    plt.grid(True)
    plt.show()

    
levelRho = calculateLevelRho()
while flightStage < 5:
    model()
    store()
    check()
    print(flightStage)
display()