import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import integrate
from scipy.signal import butter, lfilter, freqz

TEMP = 300
THICKNESS = 163E-4
VOLTAGE = -200
N_SIM_TIMES = 10000
DETECTOR_THICKNESS = 162.77

DATA_START_INDICE = 100
DATA_END_INDICE = 350
DATA_START_TIME = 4750
DATA_END_TIME = 5100

AMPLIFIER_RESISTANCE = 50
AMPLIFIER_GAIN = 100


DATA_CONTOUR_START_INDICE = 200
DATA_CONTOUR_END_INDICE = 280
DATA_CONTOUR_START_TIME = 4910
DATA_CONTOUR_END_TIME = 4970

REFRACTIVE_INDEX = 2.6154
EFFECTIVE_NUMERICAL_APERTURE = 0.5
REFRACTIVE_INDEX_CORRECTION = 1/np.sqrt((REFRACTIVE_INDEX**2-EFFECTIVE_NUMERICAL_APERTURE**2)/(1-EFFECTIVE_NUMERICAL_APERTURE**2))

SILICON_MAT = "/Users/celyn1/Documents/Sem2MPhys/SiBeta/mat_files/silicon_depth_data.mat"
DATA_FILE_NAME = "/Users/celyn1/Documents/Sem2MPhys/TCT_BetaFit/200V_Depth/TCT_x3y3.mat"

def filter_signal(signal, resistance, capacitance, dt):
    cut_off = 1/(2*np.pi*capacitance*resistance)
    pad=np.zeros(50000) 
    
    padded_signal = np.concatenate([signal, pad])
    sampling_rate = 1/(dt)
    RC = resistance*capacitance
    
    b,a = butter(1, cut_off, fs=sampling_rate, btype='low', analog = False)

    filtered = lfilter(b, a, padded_signal)
    
    return filtered, len(pad)

def rc_impulse_response(RC, fs, duration=1.0):
    t = np.arange(0, duration, 1/fs)
    h = (1/RC) * np.exp(-t / RC)
    h /= np.sum(h) 
    
    return h

def apply_rc_smear(signal, R, C, dt, gain=100):
    RC = R * C
    fs = 1/dt
    pad=np.zeros(50000) 

    padded_signal = np.concatenate([signal, pad])
    h = rc_impulse_response(RC, fs, duration=len(padded_signal)/fs)
    smeared = np.convolve(padded_signal, h, mode='full')[:len(padded_signal)]
    return smeared, len(pad)

def plt_signal(comparison_depth, file, no_sim_trials=1, sim_t_min = -16.5, sim_t_max = 16.5):
    f = h5py.File(file, 'r')

    untransferred_file = "build/data/signal_unt_SiC_0000_300_0.0050_-200V_0.0008m.txt"
    #signal time window indexes
    
    START_POSITION = 0
    END_POSITION = 81
    DATA_START_TIME = 1500
    DATA_END_TIME = 3500
    
    itime0=DATA_START_TIME #4910
    itimef=DATA_END_TIME #4970
    #simulated results
    averaged_files = []
    experiment_time_space = np.linspace(0,2000, 20000)*1e-9
    
    experiment_time_segment = experiment_time_space[DATA_START_TIME:DATA_END_TIME]

    sim_files = []

    #filename = "/Users/celyn1/Documents/Sem2MPhys/GarfieldSim_SiC/build/data/signal_0000_300_0.0050_-200V_0.0008m.txt"
    filename_sic = "build/data/signal_SiC_0000_300_0.0050_-200V_0.0008m.txt"
    #filename = "build/temp_data/raw_voxel.txt"
    #sim_results = np.genfromtxt(filename, dtype='float', delimiter='  ', skip_header=False)
    sim_results_sic = np.genfromtxt(filename_sic, dtype='float', delimiter='  ', skip_header=False)
    unt_results_sic = np.genfromtxt(untransferred_file, dtype='float', delimiter='  ', skip_header=False)

    averaged_signal = []

    #signals = sim_results[:,1]*1e-6
    signals_sic = sim_results_sic[:,1]*1e-6
    unt_signal = unt_results_sic[:,1]*1e-6
   

    averaged_signal = np.array(averaged_signal)
    
    sim_areas = []

    simulated_time = np.linspace(sim_t_min, sim_t_max, 10000)*1e-9
    dt = simulated_time[1]-simulated_time[0]
    print("DT=",dt)

    sim_segment = simulated_time[4000:6000]
    total_charge_sim = integrate.simpson(signals_sic, x = simulated_time)
    unt_total_charge = integrate.simpson(unt_signal, x = simulated_time)
    total_charge_exp = integrate.simpson(f["data"][:,0,0,0,0,0,0,40][DATA_START_TIME:DATA_END_TIME]/AMPLIFIER_RESISTANCE, experiment_time_segment)
    print("rescaled: ",-total_charge_sim/1.602e-19*4.242425911557220854e+06)
    print("unscaled: ",-total_charge_sim/1.602e-19)
    print(total_charge_exp/1.602e-19)
   
    
    gain_corrected_charge = total_charge_exp #as we are not considering amplifier effects here for the RC behaviour


    signal_segment = -1*signals_sic[4000:6000]*4.242425911557220854e+06
    unt_segment = -1*unt_signal[4000:6000]*4.242425911557220854e+06
    capacitance = 19.45e-12
    capacitance = total_charge_exp/200 
    #print("CAPACITANCE = ", capacitance)


    
    
    resistance = 200/np.max(f["data"][:,0,0,0,0,0,0,40][DATA_START_TIME:DATA_END_TIME]/AMPLIFIER_RESISTANCE) #divide by gain  here as we are looking purely at the diode behaviour
    #resistance=AMPLIFIER_RESISTANCE+resistance
    #resistance=200/(10**-11)
    print("RESISTANCE:",resistance)    
    print("CAPACITANCE:",capacitance)
    filtered_signal_full, pad = filter_signal(-1*signals_sic*4.242425911557220854e+06, resistance, capacitance,dt)
    
    filtered_signal,pad = filter_signal(signal_segment, resistance, capacitance, dt)

    new_time = np.linspace(simulated_time[4000],simulated_time[6000]+dt*(pad), (6000-4000)+(pad))
    new_time_2 = np.linspace(sim_t_min, sim_t_max+pad*dt, 10000+pad)*1e-9
    
    total_filtered_charge = integrate.simpson(filtered_signal, x = new_time)
    full_filtered_charge = integrate.simpson(filtered_signal_full, x = new_time_2)
    
    print("filtered:",total_filtered_charge/1.602e-19)
    print(full_filtered_charge/1.602e-19)
   
    print(np.shape(filtered_signal))

    plt.plot(experiment_time_segment-DATA_START_TIME*10**-10-0.25e-7,f["data"][:,0,0,0,0,0,0,40][DATA_START_TIME:DATA_END_TIME]/AMPLIFIER_RESISTANCE, label = "Experiment")
    print(7945226/4.242425911557220854e+06)
    """plt.plot(simulated_time[4000:6000], signal_segment, label = "Transferred")
    plt.plot(simulated_time[4000:6000], unt_segment, label = "Raw signal")"""


    print("Ratio = ", -total_charge_sim*4.242425911557220854e+06/total_charge_exp)
    print("The simulation integrates to:" , -total_charge_sim/1.602e-19)
    print("The untransferred signal integrates to:" , -unt_total_charge/1.602e-19)
    print("The ratio between the transferred and untransferred signals is: ", total_charge_sim/unt_total_charge)
    #new_time_2 = np.linspace(simulated_time[0],simulated_time[-1]+0.1e-9*(pad), len(simulated_time)+(pad))
    
    
    plt.plot(new_time,filtered_signal   , label = "Transferred signal")

    #plt.plot(new_time_2,filtered_signal_full   , label = "Transferred signal full")
    plt.xlabel("Time (s)")
    plt.ylabel("Signal (A)")
    plt.legend()
    plt.show()





def plt_signal_silicon(silicon_file, sic_file, no_sim_trials=1, sim_t_min = -16.5, sim_t_max = 16.5):
    f = h5py.File(silicon_file, 'r')

    p = h5py.File(sic_file, 'r')
    #signal time window indexes
    
    SI_START_POSITION = 100
    SI_END_POSITION = 350
    SI_DATA_START_TIME = 4750
    SI_DATA_END_TIME = 5100

    SIC_START_POSITION = 0
    SIC_END_POSITION = 81
    SIC_DATA_START_TIME = 1500
    SIC_DATA_END_TIME = 3500
   
    averaged_files = []
    silicon_time_space = np.linspace(0,10000, 100000)*10**-9
    silicon_time_segment = silicon_time_space[SI_DATA_START_TIME:SI_DATA_END_TIME]
    
    sic_time_space = np.linspace(0,2000, 20000)*10**-9
    
    sic_time_segment = sic_time_space[SIC_DATA_START_TIME:SIC_DATA_END_TIME]

    total_charge_sic = integrate.simpson(p["data"][:,0,0,0,0,0,0,40][SIC_DATA_START_TIME:SIC_DATA_END_TIME]/AMPLIFIER_RESISTANCE, sic_time_segment)
    total_charge_si = integrate.simpson(f["data"][:,0,0,0,0,0,0,260][SI_DATA_START_TIME:SI_DATA_END_TIME]/AMPLIFIER_RESISTANCE, silicon_time_segment)


    print("silicon total charge: ", total_charge_si)
    print("silicon carbide total charge: ", total_charge_sic)

    plt.plot(sic_time_segment-SIC_DATA_START_TIME*10**-10, p["data"][:,0,0,0,0,0,0,40][SIC_DATA_START_TIME:SIC_DATA_END_TIME]/AMPLIFIER_RESISTANCE, label = "SiC")
    plt.plot(silicon_time_segment-SI_DATA_START_TIME*10**-10,f["data"][:,0,0,0,0,0,0,260][SI_DATA_START_TIME:SI_DATA_END_TIME]/AMPLIFIER_RESISTANCE, label = "Si")


    plt.legend()
    plt.show()

#plt_signal_silicon(SILICON_MAT, "/Users/celyn1/Documents/Sem2MPhys/TCT_BetaFit/200V_Depth/TCT_x3y3.mat")
plt_signal(5,"/Users/celyn1/Documents/Sem2MPhys/TCT_BetaFit/200V_Depth/TCT_x3y3.mat",1)