import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import Boltzmann

DISTANCE_BETWEEN_PLATES = 50e-4 #cm

c_axis = False

BIAS_VOLTAGE = -200

high_field_limit = 2e5 #[V/cm]


#Low field mobility
e_lowfield_mobility = 950 #[cm2/Vs]
h_lowfield_mobility = 125 #[cm2/Vs]

#High field mobility calculation parameters
alpha_sat = 1.2
beta_sat = 1
delta_sat = -0.44

#Saturation velocity
e_v_sat = 2.2e7 #[cm/s]
h_v_sat = 2.2e7 #[cm/s]

#hall factors
e_h_factor = 1
h_h_factor = 1

#low field parameters
E_MU_MAX_300 = 950
H_MU_MAX_300 = 125

E_MU_MIN_300 = 40
H_MU_MIN_300 = 15.9

E_N_REF = 1.94e17
H_N_REF = 1.76e19

E_LOW_ALPHA = 0.61
H_LOW_ALPHA = 0.34

E_LOW_BETA = -0.5
H_LOW_BETA = -0.5

E_LOW_GAMMA = -2.40
H_LOW_GAMMA = -2.15

#impact ionisation parameters
e_alpha_impact = 1.4686e6 #[1/cm]
h_alpha_impact = 5.5222e6 #[1/cm]
e_beta_impact = 1.2075e+07 #[1/cm]
h_beta_impact = 1.2724e+07  #[1/cm]

e_c_Mobratio = 0.83 #mu//=m_perp/0.83
h_c_Mobratio = 1

e_c_Velratio = 1.16

def electric_field(distance_between_plates, bias_voltage):
    return bias_voltage/distance_between_plates

def weighting_field(distance_between_plates):
    return 1/distance_between_plates

def calculate_electron_mobility(E,T):
    #high field model
    T0 = 300
    delta = -0.44
    beta = 1
    alpha = 1.2
    v_sat = e_v_sat*(T/300)**delta

    alpha_sat = alpha*(T/300)**beta
    #alpha_sat = alpha
    mu_high = e_lowfield_mobility/(1+((e_lowfield_mobility*E)/v_sat)**alpha_sat)**(1/alpha_sat)

    return mu_high

def calculate_electron_mobility_aniso(E,T):
    #high field model along the c axis
    T0 = 300
    delta = -0.44
    beta = 1
    alpha = 1.2
    
    v_aniso = e_v_sat/e_c_Velratio
    mu_aniso = e_lowfield_mobility/e_c_Mobratio
    
    v_sat = v_aniso*(T/T0)**delta

    alpha_sat = alpha*(T/T0)**beta
    #alpha_sat = alpha
    mu_high = mu_aniso/(1+((mu_aniso*E)/v_sat)**alpha_sat)**(1/alpha_sat)

    return mu_high


def calculate_hole_mobility(E,T):
     #high field model
    T0 = 300
    delta = -0.44
    beta = 1
    alpha = 1.2
    
    v_sat = h_v_sat*(T/300)**delta

    alpha_sat = alpha*(T/300)**beta
    #alpha_sat = alpha
    mu_high = h_lowfield_mobility/(1+((h_lowfield_mobility*E)/v_sat)**alpha_sat)**(1/alpha_sat)

    return mu_high

def calculate_hole_velocity(E_field, c_axis, temperature):
    #the c axis is the z axis
    E_mag = E_field
    E_field_x= E_field_z = 0
    if c_axis:
        E_field_x = 0
        E_field_z = E_mag
        mu_x = calculate_hole_mobility(E_field_x,temperature)
        mu_z = calculate_hole_mobility(E_field_z,temperature)
    else:
        E_field_x = E_mag
        E_field_z = 0
        mu_x = calculate_hole_mobility(E_field_x,temperature)
        mu_z = calculate_hole_mobility(E_field_z,temperature)
    
    vx = mu_x * E_field_x
    vz = mu_z * E_field_z
    
    return vx, vz

def calculate_electron_velocity(E_field, c_axis, temperature):
    #the c axis is the z axis
    E_mag = E_field
    E_field_x = E_field_z = 0
    
    if c_axis:
        E_field_x = 0
        E_field_z = E_mag

    else:
        E_field_x = E_mag
        E_field_z = 0
    
    mu_x = calculate_electron_mobility(E_mag,temperature)
    mu_z = calculate_electron_mobility_aniso(E_mag,temperature)

    vx = mu_x * E_field_x
    vz = mu_z * E_field_z

    return vx, vz

def calculate_current(vx, vz, distance_between_plates, charge, c_axis):
    print("Compiting current")
    if c_axis:
        print("CURRENT: ", charge*weighting_field(distance_between_plates)*vz)
        return charge*weighting_field(distance_between_plates)*vz
    else:
        print("CURRENT: ", charge*weighting_field(distance_between_plates)*vx)
        return charge*weighting_field(distance_between_plates)*vx

def e_low_field_mobility(T,Nd, Na):
    """The lattice scattering (acoustic phonons) and ionized impurity scattering, together with
    piezoelectric scattering are the most relevant mechanisms which limit the mean free path
    of carriers at low electric fields in SiC [125, 126]"""

    mu_max = E_MU_MAX_300*(T/300)**E_LOW_GAMMA
    mu_min = E_MU_MIN_300*(T/300)**E_LOW_BETA

    mu_low_perp = mu_min + (mu_max-mu_min)/(1+((Nd+Na)/E_N_REF)**E_LOW_ALPHA)
    mu_low_par = mu_low_perp*1.2

    return mu_low_perp, mu_low_par

def h_low_field_mobility(T,Nd, Na):
    """The lattice scattering (acoustic phonons) and ionized impurity scattering, together with
    piezoelectric scattering are the most relevant mechanisms which limit the mean free path
    of carriers at low electric fields in SiC [125, 126]"""

    mu_max = H_MU_MAX_300*(T/300)**H_LOW_GAMMA
    mu_min = H_MU_MIN_300*(T/300)**H_LOW_BETA

    mu_low_perp = mu_min + (mu_max-mu_min)/(1+((Nd+Na)/H_N_REF)**H_LOW_ALPHA)
    mu_low_par = mu_low_perp

    return mu_low_perp, mu_low_par

def impact_ionisation(gamma,e_b_low,e_a_low, e_b_high, e_a_high,h_b_low,h_a_low, h_b_high, h_a_high,E, e_density, h_density, drift_e, drift_h):
    """describes the generation of free carriers, resulting in the additional generation of carriers due to impact ionisation
    
    everything atm is in cm and eV"""
    hbaromega = 1e3
    T0 = 300
    kb=1.380649e-23
    kb_eV = kb/1.602e-19

    if E < 4e-5:
        e_a,e_b = e_a_low, e_b_low
        h_a,h_b = h_a_low, h_b_low
    else:
        e_a,e_b = e_a_high, e_b_high
        h_a,h_b = h_a_high, h_b_high

    gamma = np.tanh(hbaromega/(2*kb_eV*T0))/np.tanh(hbaromega/(2*kb_eV*T0))

    alpha_e = gamma*e_a*np.exp(-e_b*gamma/E)
    alpha_h = gamma*h_a*np.exp(-h_b*gamma/E)

    generation_rate = alpha_e*e_density*drift_e + alpha_h*h_density*drift_h

    return generation_rate




def drift(step_size, distance_between_plates, temperature, e_z_start, e_x_start, h_z_start, h_x_start, charge, bias, c_axis):
    e_mag = electric_field(distance_between_plates, bias)
    
    hole_current = electron_current = 0

    electron_x_position = e_x_start
    electron_z_position = e_z_start

    hole_x_position = h_x_start
    hole_z_position = h_z_start

    scaled_time_step = step_size * 1e-9
    """print("electric field",e_mag)"""

    h_vx , h_vz = calculate_hole_velocity(np.abs(e_mag), c_axis,temperature)
    e_vx, e_vz = calculate_electron_velocity(np.abs(e_mag), c_axis,temperature)


    electron_z_displacement = e_vz*scaled_time_step
    hole_z_displacement = h_vz*scaled_time_step

    electron_x_displacement = e_vx*scaled_time_step
    hole_x_displacement = h_vx*scaled_time_step

    electron_z_position += electron_z_displacement
    hole_z_position -= hole_z_displacement

    electron_x_position += electron_x_displacement
    hole_x_position -= hole_x_displacement



    if c_axis and 0 < electron_z_position < distance_between_plates:
        print("electron in region", electron_z_position)
        electron_current = calculate_current(e_vx, e_vz, distance_between_plates, -charge, c_axis)
        #print(electron_current)
    
    elif c_axis and( (0 > electron_z_position) or  (electron_z_position > distance_between_plates)):
        print("electron not in region", electron_z_position)
        electron_z_position = distance_between_plates
        electron_current = 0
    
    if c_axis and 0 < hole_z_position < distance_between_plates:
        print("hole in region", hole_z_position)
        hole_current = calculate_current(h_vx, h_vz, distance_between_plates, charge, c_axis)
        #print(hole_current)
    
    elif c_axis and ( (0 > hole_z_position) or ( hole_z_position > distance_between_plates)):
        print("hole not in region", hole_z_position)
        hole_z_position = 0
        hole_current = 0
    
    if (not c_axis) and (0 < electron_x_position < distance_between_plates):
        electron_current = calculate_current(e_vx, e_vz, distance_between_plates, -charge, c_axis)

    elif (not c_axis) and ( (0 > electron_x_position) or  (electron_x_position > distance_between_plates)):
        print("electron out of bounds")
        electron_x_position = distance_between_plates
        electron_current = 0
    
    if (not c_axis) and 0 < hole_x_position < distance_between_plates:
        print("hole in bounds")
        hole_current = calculate_current(h_vx, h_vz, distance_between_plates, charge, c_axis)
       
    
    elif (not c_axis) and ( (0 > hole_x_position) or  (hole_x_position > distance_between_plates)):
        print("holes out of bounds")
        hole_x_position = 0
        hole_current = 0

    
    return electron_current, hole_current, electron_x_position, electron_z_position, hole_x_position, hole_z_position


def theoretical_signal():
    e_field = BIAS_VOLTAGE/DISTANCE_BETWEEN_PLATES
    weighting_field = -1/DISTANCE_BETWEEN_PLATES

    electron_mobility = 950
    hole_mobility = 125

    electron_drift = electron_mobility*e_field
    hole_drift = hole_mobility*e_field

    if abs(electron_drift) > 2.2e7:
        electron_drift = -2.2e7
    if abs(hole_drift) > 2.2e7:
        hole_drift=-2.2e7

    electron_last = -DISTANCE_BETWEEN_PLATES/electron_drift/2
    hole_last = -DISTANCE_BETWEEN_PLATES/hole_drift/2

    time_space = np.linspace(0, 10, 1000)

    electron_signal = np.zeros(1000)
    hole_signal = np.zeros(1000)

    electron_peak = weighting_field*1.602e-19*electron_drift
    hole_peak = weighting_field*1.602e-19*hole_drift
    print(electron_last)
    electron_signal[np.where(time_space<electron_last*1e9)] = electron_peak
    hole_signal[np.where(time_space<hole_last*1e9)] = hole_peak

    total_signal = electron_signal + hole_signal
    print("theoretical velocity set to :", electron_drift)
    #plt.plot(time_space, total_signal)
    #plt.show()
    return time_space, total_signal

def theoretical_signal_c():
    e_field = BIAS_VOLTAGE/DISTANCE_BETWEEN_PLATES
    weighting_field = -1/DISTANCE_BETWEEN_PLATES

    electron_mobility = 950
    hole_mobility = 125

    electron_drift = electron_mobility*e_field/e_c_Mobratio
    hole_drift = hole_mobility*e_field

    if abs(electron_drift) > 2.2e7/1.16:
        electron_drift = -2.2e7/1.16
    if abs(hole_drift) > 2.2e7/1.16:
        hole_drift=-2.2e7/1.16

    electron_last = -DISTANCE_BETWEEN_PLATES/electron_drift/2
    hole_last = -DISTANCE_BETWEEN_PLATES/hole_drift/2

    time_space = np.linspace(0, 10, 1000)

    electron_signal = np.zeros(1000)
    hole_signal = np.zeros(1000)

    electron_peak = weighting_field*1.602e-19*electron_drift
    hole_peak = weighting_field*1.602e-19*hole_drift
    print(electron_last)
    electron_signal[np.where(time_space<electron_last*1e9)] = electron_peak
    hole_signal[np.where(time_space<hole_last*1e9)] = hole_peak

    total_signal = electron_signal + hole_signal

    #plt.plot(time_space, total_signal)
    #plt.show()
    return time_space, total_signal


def main():


    """field_vals = np.logspace(2,7,100)
    h_velocity_vals = []
    e_x_velocity_vals = []
    e_z_velocity_vals = []

    for field in field_vals:
        h_velocity_vals=np.append(h_velocity_vals, calculate_hole_velocity(field,True, 300)[1])
        e_x_velocity_vals = np.append(e_x_velocity_vals, calculate_electron_velocity(field,False, 300)[0])
        e_z_velocity_vals=np.append(e_z_velocity_vals, calculate_electron_velocity(field,True, 300)[1])
    
    plt.plot(np.log10(field_vals), h_velocity_vals[:], label="holes")
    plt.plot(np.log10(field_vals), e_x_velocity_vals[:], label = "electron velocity")
    plt.plot(np.log10(field_vals), e_z_velocity_vals[:], label = "electron velocity, field along c axis}")
    plt.axhline(y=(calculate_electron_velocity(abs(electric_field(DISTANCE_BETWEEN_PLATES, BIAS_VOLTAGE)), False, 300)[0]))
    plt.axvline(x=(np.log10(abs(electric_field(DISTANCE_BETWEEN_PLATES, BIAS_VOLTAGE)))))
    plt.xlabel("electric field (V/cm)[log scale]")
    plt.ylabel("velocity (cm/s)")
    plt.legend()
    plt.show()"""
    
    
    print("E field: ", electric_field(DISTANCE_BETWEEN_PLATES,BIAS_VOLTAGE))
    hole_mobility = h_lowfield_mobility
    electron_mobility = e_lowfield_mobility
    bias = BIAS_VOLTAGE
    n_steps = 10000
    t_max = 20
    t_min  = 0
    step_size = (t_max-t_min)/n_steps
    time = t_min 
    start_position_x = DISTANCE_BETWEEN_PLATES/2
    start_position_y = 0
    if c_axis:
        hole_position_x = 0
        hole_position_z = DISTANCE_BETWEEN_PLATES/2
        electron_position_x = 0
        electron_position_z = DISTANCE_BETWEEN_PLATES/2
    else:
        hole_position_z = 0
        hole_position_x = DISTANCE_BETWEEN_PLATES/2
        electron_position_z = 0
        electron_position_x = DISTANCE_BETWEEN_PLATES/2

    electron_current_vector = []
    electron_x_vector = [electron_position_x]
    electron_z_vector = [electron_position_z]
    
    hole_current_vector  = []
    hole_x_vector = [hole_position_x]
    hole_z_vector = [hole_position_z]

    simulation_time = np.linspace(t_min, t_max, n_steps)

    while time < t_max:
        time += step_size
        electron_current, hole_current, electron_position_x, electron_position_z, hole_position_x, hole_position_z = drift(step_size, DISTANCE_BETWEEN_PLATES, 300, 
                                                                                                                           electron_position_z, electron_position_x,
                                                                                                                           hole_position_z, hole_position_x, 1.602e-19, bias, c_axis)
        electron_x_vector = np.append(electron_x_vector, electron_position_x)
        electron_z_vector = np.append(electron_z_vector, electron_position_z)

        hole_x_vector = np.append(hole_x_vector, hole_position_x)
        hole_z_vector = np.append(hole_z_vector, hole_position_z)

        electron_current_vector = np.append(electron_current_vector, electron_current)
        hole_current_vector = np.append(hole_current_vector, hole_current)
    
    theory_time, theory_signal= theoretical_signal()
    theory_time_c, theory_signal_c= theoretical_signal_c()
    
    plt.scatter(start_position_x, start_position_y)
    plt.plot(electron_x_vector, electron_z_vector, label = "electron")
    plt.plot(hole_x_vector, hole_z_vector, label = "hole")
    plt.legend()
    #plt.xlim(-150e-15,150e-15)
    plt.show()
    total_current = electron_current_vector - hole_current_vector
    plt.plot(theory_time, theory_signal, label = "simple theory")
    #plt.plot(theory_time_c, theory_signal_c, label = "simple theory, c axis values")
    plt.plot(simulation_time, -electron_current_vector[:], linewidth = 4,label = "electron")
    plt.plot(simulation_time, hole_current_vector[:], linewidth = 4, label = "hole")
    plt.plot(simulation_time, -total_current[:], linestyle = 'dashdot',color = 'black',label = "total")
    
    plt.xlabel("Time (ns)")
    plt.ylabel("Current (A)")
    plt.legend()
    plt.show()

"""  temperature_values = np.linspace(10,500,10)
    h_mobility_low_field = []
    e_mobility_low_field = []
    e_c_mobility_low_field = []
    for temperature in temperature_values:
        h_mobility_low_field = np.append(h_mobility_low_field, h_low_field_mobility(temperature,0, 3e15)[0])
        e_mobility_low_field = np.append(e_mobility_low_field, e_low_field_mobility(temperature,0, 3e15)[0])
        e_c_mobility_low_field = np.append(e_c_mobility_low_field, e_low_field_mobility(temperature,0, 3e15)[1])
    plt.plot(temperature_values, h_mobility_low_field, label = "holes")
    plt.plot(temperature_values, e_mobility_low_field, label = "electrons")
    plt.plot(temperature_values, e_c_mobility_low_field, label = "electrons, c axis")
    plt.show()
"""





"""    plt.loglog((field_vals), h_velocity_vals[:], label="holes")
    plt.loglog((field_vals), e_x_velocity_vals[:], label = "electron velocity")
    plt.loglog((field_vals), e_z_velocity_vals[:], label = "electron velocity, field along c axis}")
    plt.xlabel("electric field (V/cm)[log scale]")
    plt.ylabel("velocity (cm/s)")
    plt.legend()
    plt.show()
    
    field_vals = np.logspace(2,7,100)
    h_mobility_vals_300 = []
    e_c_mobility_vals_300 = []
    e_mobility_vals_300 = []
    h_mobility_vals_400 = []
    e_c_mobility_vals_400 = []
    e_mobility_vals_400 = []
    e_c_mobility_vals_600 = []
    for field in field_vals:
        h_mobility_vals_300=np.append(h_mobility_vals_300, calculate_hole_mobility(field, 300))
        e_c_mobility_vals_300 = np.append(e_c_mobility_vals_300, calculate_electron_mobility_aniso(field, 300))
        e_mobility_vals_300=np.append(e_mobility_vals_300, calculate_electron_mobility(field, 300))
        e_c_mobility_vals_600=np.append(e_c_mobility_vals_600, calculate_electron_mobility(field, 600))
        h_mobility_vals_400=np.append(h_mobility_vals_400, calculate_hole_mobility(field, 400))
        e_c_mobility_vals_400 = np.append(e_c_mobility_vals_400, calculate_electron_mobility_aniso(field, 400))
        e_mobility_vals_400=np.append(e_mobility_vals_400, calculate_electron_mobility(field, 400))
    #plt.plot(np.log10(field_vals), h_mobility_vals_300[:], label="holes 300K")
    plt.plot(np.log10(field_vals), e_c_mobility_vals_300[:], label = "electron c-axis mobility 300K")
    #plt.plot(np.log10(field_vals), e_mobility_vals_300, label = "electron mobility 300K")
    #plt.plot(np.log10(field_vals), h_mobility_vals_400[:], label="holes 400K")
    plt.plot(np.log10(field_vals), e_c_mobility_vals_400[:], label = "electron c-axis mobility 400K")
    plt.plot(np.log10(field_vals), e_c_mobility_vals_600[:], label = "electron c-axis mobility 600K")
    #plt.plot(np.log10(field_vals), e_mobility_vals_400, label = "electron mobility 400K")
    plt.xlabel("electric field (V/cm)[log scale]")
    plt.ylabel("mobility (cm^2/Vs)")
    plt.legend()
    plt.show()"""


main()