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

#select single case or multi case mode
mode = 2

startAltList = []
v_achievableList = []
v_achievable = 0
def model():

    global Cd, peakCd, baseCd, v_achievable, localMachNumber, v_airspeed, t_current, aircraft_alt, aircraft_lift, aircraft_rho, v_y,E_remaining,aircraftOutsideTemp
    
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
    Cd = calculateCd(v_airspeed,baseCd,localMachNumber)
    drag = 0.5 * Cd * Ref_A * aircraft_rho * pow(v_airspeed,2)
    lift = 0.5 * Cl * Ref_A * aircraft_rho * pow(v_airspeed,2)
    aircraft_lift = lift
    #velocity update for accelerating case
    v_airspeed = v_airspeed + ((thrust - drag)/mass)*timestep
    if flightStage ==2:
        #v_y = ((thrust-drag)/(mass*g) )*v_airspeed
        v_airspeed = v_airspeed - math.sin(climbAngle*0.0174533)*g*timestep
        #print v_airspeed
        v_y = math.sin(climbAngle*0.0174533)*v_airspeed
        #print(v_y)
    else:
        
        v_y = 0
        if flightStage ==4:
            v_y = -v_airspeed * (Cd/Cl)
    
    #position update
    aircraft_alt = aircraft_alt + v_y*timestep
    # print('alt: ',aircraft_alt)
    # print('v_airspeed: ',v_airspeed)
    # print('v_y: ',v_y)
    # print('levelRho: ',levelRho)

    #mach update
    localMachNumber = calculateMach(aircraft_alt)
    #Energy usage
    E_remaining = E_remaining - power*timestep

    #timestep
    t_current = t_current + timestep
    
    if v_airspeed > v_achievable:
        v_achievable = v_airspeed
    if Cd > peakCd:
        peakCd = Cd

def store():
    global v_airspeedList,timeList
    v_airspeedList.append(v_airspeed)
    timeList.append(t_current)
    aircraft_altList.append(aircraft_alt)
    machList.append(localMachNumber)
    energyList.append(E_remaining)
    T_list.append(calculateT(aircraft_alt))
def check():
    global flightStage,v_achievable
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
            flightStage = 5

    # if flightStage == 4:
    #     if aircraft_alt<=runway_alt:
    #         flightStage = 5

    if E_remaining/E_capacity <= E_endMachFraction:
            #print('triggered')
            flightStage = 5
    #print('alt: ', aircraft_alt, 'stg: ',flightStage, 'E_remaining/E_capacity', E_remaining/E_capacity)
    
def calculateLevelRho():
    global thrust_max, Ref_A, Cd, target_v 
    rho_min = thrust_max/(0.5*baseCd*1.5*Ref_A*pow(target_v,2))
    #rho_min = (mass*g)/(0.5*Cl*Ref_A*pow(target_v,2))
    return rho_min
    
def calculateRho(altAboveSeaLevel):
    #src: https://www.grc.nasa.gov/www/K-12/airplane/atmosmet.html
    #metric units
    if altAboveSeaLevel <= 11000:
        T = 15.04 - 0.00649*altAboveSeaLevel
        p = 101.29 * pow((T+273.1)/(288.08),5.526)
        rho = p/((0.2869*(T+273.1)))
        return rho
    if 11000 < altAboveSeaLevel <= 25000:
        T= -56.46
        p = 22.65 * math.exp(1.73-0.000157*altAboveSeaLevel)
        rho = p/((0.2869*(T+273.1)))
        return rho
    if altAboveSeaLevel > 25000:
        T = -131.21 + 0.00299*altAboveSeaLevel
        p = 2.488 * pow((T+273.1)/(216.1),-11.388)
        return rho

def calculateT(altAboveSeaLevel):
    #src: https://www.grc.nasa.gov/www/K-12/airplane/atmosmet.html
    if altAboveSeaLevel <= 11000:
        T = 15.04 - 0.00649*altAboveSeaLevel
        return T
    if 11000 < altAboveSeaLevel <= 25000:
        T= -56.46
        return T
    if altAboveSeaLevel > 25000:
        T = -131.21 + 0.00299*altAboveSeaLevel
        return T
def calculateMach(altAboveSeaLevel):
    #src: http://www.tscm.com/mach-as.pdf
    altInFt = altAboveSeaLevel * 3.28084
    machInKnots = 29.06 * math.sqrt(518.7-3.57*(altInFt/1000))
    return 0.514444*machInKnots

def calculateCd(airspeed,baseCd,localMachNumber):
    machNumber = airspeed/localMachNumber

    adjustedCd = 0
    if machNumber< 0.7:
        return baseCd
    if 0.7 <= machNumber < 1.05:
        adjustedCd = 1.5*baseCd
        #adjustedCd = baseCd * ((5/3)*machNumber + (1/6))
    if 1.05 <= machNumber:
        adjustedCd = 1.2*baseCd
        #adjustedCd = baseCd * (-1.5*machNumber + 3)
    return adjustedCd

def display():
    #print('achievable v:', v_achievable, 'at Cd: ', Cd)

    if mode == 1:
        plt.figure(1)
        plt.subplot(211)
        plt.ylabel('airspeed (m/s)')
        plt.xlabel('Time (s)')
        plt.plot(timeList,v_airspeedList)
        plt.plot(timeList,machList)

        plt.subplot(212)
        plt.ylabel('alt (m)')
        plt.xlabel('Time (s)')
        plt.plot(timeList,aircraft_altList)
        
        plt.figure(2)
        plt.subplot(111)
        plt.ylabel('Temp (degC)')
        plt.xlabel('Time(s)')
        plt.plot(timeList,T_list)
        plt.show()
    if mode == 2:
        #multi case
        fig, ax = plt.subplots()
        ax.scatter(resultY, resultZ, alpha=0.5)
        plt.xlabel('subsonic Cd')
        plt.ylabel('battery mass (kg)')
        plt.figtext(.02,.02,'')
        plt.show()



if mode == 1:
    #single case mode
    CdList = [0.036]
    mass_battList = [44]

if mode == 2:
    CdList = [x / 1000.0 for x in range(20, 80, 1)]
    mass_battList = [x for x in range(10, 50, 1)]

resultX = []
resultY = []
resultZ = []
print CdList
for drag_coeff in CdList:
    for mass_batt in mass_battList:
        E_capacity = 360000 * mass_batt #J @100 Wh/kg or 360000 J/kg
        #print (E_capacity)
        E_remaining = E_capacity #J
        
        baseCd = drag_coeff
        v_achievable = 0

        #reset aircraft

        #environment vars
        g = 9.8 #m/s/s
        #Record characteristics
        target_v = 340 #m/s

        #runway characteristics
        runway_alt = 9000 #m 7012 for C-130, 0 for ground

        #aircraft design characteristics
        mass = 80 + mass_batt #kg
        Cl = 0.2 #unitless
        Ref_A = 2.251 #m^2
        E_endMachFraction = 0.1 #unitless
        max_power = 200000 #W
        max_power_duration = 120 #s
        motorThrustToPowerCurve = 0.010 #N/W max power checks with Cabo figure of 444lb of thrust at max power
        climbAngle = 10 #degrees
        glideAngle = 15 #degrees
        climbPower = 0.4 * max_power

        #aircraft flight variables
        aircraftOutsideTemp = 15
        v_airspeed = 50 #m/s Release velocity from a C-130 in cruise is 150 m/s relative to air, de-rate to 100
        v_y = 0 #m/s
        v_x = 0 #m/s
        thrust_max = max_power*motorThrustToPowerCurve #N
        timestep = 0.1 #s
        t_current = 0 #s
        aircraft_rho = calculateRho(runway_alt) #kg/m^3
        aircraft_alt = runway_alt
        power = max_power
        aircraft_lift = 0 #N
        localMachNumber = 340 #m/s
        Cd = baseCd
        #lists reset
        v_airspeedList = []
        rho_list = []
        timeList = []
        aircraft_altList = []
        aircraftOutsideTempList = []
        machList = []
        energyList = []
        T_list = []
        flightStage = 1

        levelRho = calculateLevelRho()
        peakCd = 0
        while flightStage < 5:
            
            model()
            store()
            check()
            #print(flightStage)
        
        v_achievableList.append(v_achievable)
        #check for those that make it significantly above
        
        if((v_achievable/localMachNumber) > 1.2):

            print(v_achievable,peakCd,mass_batt)
            resultX.append(v_achievable/localMachNumber)
            resultY.append(baseCd)
            resultZ.append(mass_batt)
display()