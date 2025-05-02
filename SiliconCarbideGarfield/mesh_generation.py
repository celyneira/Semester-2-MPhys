import numpy as np

import sympy as sp
from sympy import *

import math as math

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import Normalize

from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
from scipy import integrate
from scipy.signal import find_peaks

from matplotlib import rcParams

import progressbar

rcParams['font.family'] = 'DejaVu Serif'
R_LIMIT = 25*10**(-6)
Z_LIMIT = 50*10**(-6)

WAVELENGTH = 720*10**(-9)

EXTRAORDINARY_REFRACTIVE_INDEX = np.sqrt(10.52+0.1701/(WAVELENGTH**2-0.0258)+729.2/(WAVELENGTH**2-194.72))
ORDINARY_REFRACTIVE_INDEX = np.sqrt(9.0+0.1364/(WAVELENGTH**2-0.0334)+545/(WAVELENGTH**2-163.69))

PULSE_DURATION = 160*10**(-15) #s
BETA_2 = 6e-14 #m/W
REFRACTIVE_INDEX = 2.6154
H = 6.62607015*10**(-34)

PULSE_ENERGY = 1000e-12
EFFECTIVE_NUMERICAL_APERTURE = 0.5

REF_INDEX_AIR = 1.0003
REF_INDEX_ALUMINIUM = 1.51
DELAY = (2*REFRACTIVE_INDEX)/(3*10**8) #*VOXEL_POSITION
REFLECTANCE = np.abs((REFRACTIVE_INDEX-REF_INDEX_AIR)/(REFRACTIVE_INDEX+REF_INDEX_AIR))**2
REFLECTANCE = np.abs((REFRACTIVE_INDEX-REF_INDEX_ALUMINIUM)/(REFRACTIVE_INDEX+REF_INDEX_ALUMINIUM))**2

REF_INDEX_CORRECTION = 1/np.sqrt((REFRACTIVE_INDEX**2-EFFECTIVE_NUMERICAL_APERTURE**2)/(1-EFFECTIVE_NUMERICAL_APERTURE**2))

print(REF_INDEX_CORRECTION)
print(1/REF_INDEX_CORRECTION)


def gaussian_fit(x, a, mu, sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2))

def lorentzian_fit(x, A, x0, gamma):
    return A * (1 / np.pi) * (0.5 * gamma) / ((x - x0)**2 + (0.5 * gamma)**2)

def exponential_fit(x, a, b, c):
    return a*np.exp(-b*x)+c

def direct_charge_density_integrand(z, r, V, sigma):
    l = WAVELENGTH
    NAv = EFFECTIVE_NUMERICAL_APERTURE
    nr = REFRACTIVE_INDEX
    Ep = PULSE_ENERGY
    tau = PULSE_DURATION
    beta2 = BETA_2
    h = H
    freq = (3*10**8)/WAVELENGTH
    deltat = DELAY*V
    R = REFLECTANCE

    #z = z*REF_INDEX_CORRECTION

    ndirect = (4 * beta2 * Ep**2 * NAv**4 * nr**4 * np.pi**1.5 * np.log(2)) / (
    np.exp((4 * NAv**2 * nr**2 * np.pi**2 * r**2) /
           (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V - z)**2))) *
    freq * h * tau *
    (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V - z)**2)) *
    (l**2 * nr**2 + NAv**4 * np.pi**2 * (V - z)**2) *
    np.sqrt(np.log(4)))
    
    return ndirect*2*np.pi*r

def reflected_charge_density_integrand(z, r, V, sigma):
    l = WAVELENGTH
    NAv = EFFECTIVE_NUMERICAL_APERTURE
    nr = REFRACTIVE_INDEX
    Ep = PULSE_ENERGY
    tau = PULSE_DURATION
    beta2 = BETA_2
    h = H
    freq = (3*10**8)/WAVELENGTH
    deltat = DELAY*V
    R = REFLECTANCE

    #z = z*REF_INDEX_CORRECTION

    nreflected = (4 * beta2 * Ep**2 * NAv**4 * nr**4 * np.pi**1.5 * R**2 * np.log(2)) / (
    np.exp((4 * NAv**2 * nr**2 * np.pi**2 * r**2) /
           (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V + z)**2))) *
    freq * h * tau *
    (l**2 * nr**2 + NAv**4 * np.pi**2 * (V + z)**2) *
    (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V + z)**2)) *
    np.sqrt(np.log(4)))

    return nreflected*2*np.pi*r

def interference_charge_density_integrand(z, r, V, sigma):
    l = WAVELENGTH
    NAv = EFFECTIVE_NUMERICAL_APERTURE
    nr = REFRACTIVE_INDEX
    Ep = PULSE_ENERGY
    tau = PULSE_DURATION
    beta2 = BETA_2
    h = H
    freq = (3*10**8)/WAVELENGTH
    deltat = DELAY*V
    R = REFLECTANCE

    #z = z*REF_INDEX_CORRECTION

    ninterference = (
    (2**(2.5 - (2 * deltat**2) / tau**2) * beta2 * np.exp(
        (-4 * NAv**2 * nr**2 * np.pi**2 * r**2) /
        (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V - z)**2)) -
        (4 * NAv**2 * nr**2 * np.pi**2 * r**2) /
        (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V + z)**2))
    ) / 2) *
    Ep**2 * NAv**4 * nr**4 * np.pi**1.5 * R**2 * np.sqrt(np.log(2))) / (
    freq * h * tau *
    np.sqrt(l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V - z)**2)) *
    np.sqrt(l**2 * nr**2 + NAv**4 * np.pi**2 * (V - z)**2) *
    np.sqrt(l**2 * nr**2 + NAv**4 * np.pi**2 * (V + z)**2) *
    np.sqrt(l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V + z)**2)))

    return ninterference*2*np.pi*r

def charge_density_with_reflection_and_smearing(z, r, V, sigma):
    l = WAVELENGTH
    NAv = EFFECTIVE_NUMERICAL_APERTURE
    nr = REFRACTIVE_INDEX
    Ep = PULSE_ENERGY
    tau = PULSE_DURATION
    beta2 = BETA_2
    h = H
    freq = (3*10**8)/WAVELENGTH
    deltat = DELAY*V
    R = REFLECTANCE

    z=z*REF_INDEX_CORRECTION

    ndirect = (4 * beta2 * Ep**2 * NAv**4 * nr**4 * np.pi**1.5 * np.log(2)) / (
    np.exp((4 * NAv**2 * nr**2 * np.pi**2 * r**2) /
           (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V - z)**2))) *
    freq * h * tau *
    (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V - z)**2)) *
    (l**2 * nr**2 + NAv**4 * np.pi**2 * (V - z)**2) *
    np.sqrt(np.log(4)))

    nreflected = (4 * beta2 * Ep**2 * NAv**4 * nr**4 * np.pi**1.5 * R**2 * np.log(2)) / (
    np.exp((4 * NAv**2 * nr**2 * np.pi**2 * r**2) /
           (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V + z)**2))) *
    freq * h * tau *
    (l**2 * nr**2 + NAv**4 * np.pi**2 * (V + z)**2) *
    (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V + z)**2)) *
    np.sqrt(np.log(4)))

    ninterference = (
    (2**(2.5 - (2 * deltat**2) / tau**2) * beta2 * np.exp(
        (-4 * NAv**2 * nr**2 * np.pi**2 * r**2) /
        (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V - z)**2)) -
        (4 * NAv**2 * nr**2 * np.pi**2 * r**2) /
        (l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V + z)**2))
    ) / 2) *
    Ep**2 * NAv**4 * nr**4 * np.pi**1.5 * R**2 * np.sqrt(np.log(2))) / (
    freq * h * tau *
    np.sqrt(l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V - z)**2)) *
    np.sqrt(l**2 * nr**2 + NAv**4 * np.pi**2 * (V - z)**2) *
    np.sqrt(l**2 * nr**2 + NAv**4 * np.pi**2 * (V + z)**2) *
    np.sqrt(l**2 * nr**2 + NAv**2 * np.pi**2 * (4 * nr**2 * sigma**2 + NAv**2 * (V + z)**2)))

    result = ndirect + nreflected + ninterference
    
    return result, ndirect

def charge_density_normalisation_with_reflection_and_smearing(r_lim, z_lim, V, sigma):

    z_lim = z_lim*REF_INDEX_CORRECTION


    direct, error = integrate.dblquad(direct_charge_density_integrand, 0, r_lim, 0, z_lim, args=(V, sigma))

    reflected, error = integrate.dblquad(reflected_charge_density_integrand, 0, r_lim, 0, z_lim, args=(V, sigma))

    interference, error = integrate.dblquad(interference_charge_density_integrand, 0, r_lim, 0, z_lim, args=(V, sigma))

    total_normalisation = direct + reflected + interference

    return total_normalisation, direct

def make_mesh_reflection_smearing(r_lim, r_points, z_lim, z_points, V, sigma, unit, return_normalisation, sim_type):
    r_mesh = np.linspace(-r_lim, r_lim, r_points)
    #z_mesh = np.linspace(-z_lim, z_lim, z_points)
    z_mesh = np.linspace(0, z_lim, z_points)

    charge_density_mesh_r, charge_density_mesh_z = np.meshgrid(r_mesh, z_mesh)

    charge_density_mesh = []

    if sim_type == 'direct':
        charge_density_mesh = charge_density_with_reflection_and_smearing(charge_density_mesh_z, charge_density_mesh_r, V, sigma)[1]
        normalisation_constant = charge_density_normalisation_with_reflection_and_smearing(r_lim, z_lim, V, sigma)[1]

    if sim_type == 'reflect':
        charge_density_mesh = charge_density_with_reflection_and_smearing(charge_density_mesh_z, charge_density_mesh_r, V, 0)[0]
        normalisation_constant = charge_density_normalisation_with_reflection_and_smearing(r_lim, z_lim, V, 0)[0]

    if sim_type == 'reflect and smear':
        charge_density_mesh = charge_density_with_reflection_and_smearing(charge_density_mesh_z, charge_density_mesh_r, V, sigma)[0]
        normalisation_constant = charge_density_normalisation_with_reflection_and_smearing(r_lim, z_lim, V, sigma)[0]

    # divide by normalisatiuon constant if required

    if unit:
        print('Total charge carriers before scaling: {:e}'.format(np.sum(charge_density_mesh)))
        charge_density_mesh_normalised = charge_density_mesh/normalisation_constant
        print('Total charge carriers after scaling: {:e}'.format(np.sum(charge_density_mesh_normalised)))
        if return_normalisation:
            return charge_density_mesh_r, charge_density_mesh_z, charge_density_mesh_normalised, normalisation_constant
        
        else:
            return charge_density_mesh_r, charge_density_mesh_z, charge_density_mesh_normalised
        
    else:

        if return_normalisation:
            return charge_density_mesh_r, charge_density_mesh_z, charge_density_mesh, normalisation_constant 
        
        else: 
            return charge_density_mesh_r, charge_density_mesh_z, charge_density_mesh

def scale_mesh_deleted(r_lim, r_points, z_lim, z_points, charge_density_mesh, V, make_plots, replace_deleted = False):
    r_space = np.linspace(-r_lim, r_lim, r_points)
    z_space = np.linspace(0, z_lim, z_points)
    
    r_mesh, z_mesh = np.meshgrid(r_space, z_space)
    
    dr = 2*r_lim/r_points
    dz = z_lim/z_points

    r_space_corrected = np.linspace(-r_lim, r_lim, r_points)
    z_space_corrected = np.linspace(0, z_lim*REF_INDEX_CORRECTION, z_points)    

    r_mesh_corrected, z_mesh_corrected = np.meshgrid(r_space_corrected, z_space_corrected)

    dr_corrected = dr
    dz_corrected = dz*REF_INDEX_CORRECTION

    charge_density_mesh_scaled = charge_density_mesh*dr_corrected*dz_corrected*np.pi*np.abs(r_space_corrected)*10000

    if make_plots:
        plt.clf()
        plt.pcolor(r_mesh, z_mesh, charge_density_mesh_scaled, cmap='BuPu')
        cbar = plt.colorbar()
        cbar.set_label(label= "Number of Charge Carriers", size=14)
        plt.contour(r_mesh, z_mesh, charge_density_mesh_scaled, cmap='BuPu')
        plt.xlabel("Transverse Distance, r [m]", fontsize=15)
        plt.ylabel("Depth, z [m]", fontsize=15)
        plt.ticklabel_format(style='sci', scilimits=(0.01,100))
        plt.axhline(0, color='white', linestyle='dashed')
        plt.show()
        plt.clf()

    deleted_mesh = np.zeros((z_points, r_points))

    number_deleted = 0
    for i in range(z_points):
        for j in range(r_points):
            if charge_density_mesh_scaled[i,j] < 1:
                number_deleted+=1*charge_density_mesh_scaled[i,j]
                #append the deleted mesh with the carriers that are manually removed
                deleted_mesh[i,j] = charge_density_mesh_scaled[i,j]
                #set the carriers at these points to zero
                charge_density_mesh_scaled[i,j] = 0
    print(number_deleted, " carriers were initially deleted.")
    if replace_deleted:
        positive_counts = []
        negative_counts = []

        #distributing the deleted carriers using a normal distribution
        for i in range(z_points):
            j = i
            while True:
                try:
                    positive_data = deleted_mesh[i,int(r_points/2):]
                    negative_data = deleted_mesh[i,:int(r_points/2)]
                    
                    positive_r = r_mesh[i,int(r_points/2):]
                    negative_r = r_mesh[i,:int(r_points/2)]
                    """
                    if V == 0:
                        return 0, charge_density_mesh_scaled
                    """
                    popt_pos,_ = curve_fit(gaussian_fit, positive_r, positive_data,
                                        p0=[max(positive_data), positive_r[np.where(positive_data == np.max(positive_data))][0], 1e-6])
                    popt_neg,_ = curve_fit(gaussian_fit, negative_r, negative_data,
                                        p0=[max(negative_data), negative_r[np.where(negative_data == np.max(negative_data))][0], 1e-6])
                    break
                except RuntimeError:
                    i-=1
                    if i<0:
                        print('Curvefit failed for all indices. Using an approximate Gaussian.')
                        popt_pos = [max(positive_data), positive_r[np.where(positive_data == np.max(positive_data))][0], 1e-5]
                        popt_neg = [max(negative_data), negative_r[np.where(negative_data == np.max(negative_data))][0], 1e-5]
                        break
            
            mean_pos = popt_pos[1]
            std_dev_pos = popt_pos[2]
            mean_neg = popt_neg[1]
            std_dev_neg = popt_neg[2]

            positive_counts_temp, _ = np.histogram(np.random.normal(mean_pos, np.abs(std_dev_pos), int(np.round(np.sum(positive_data)))), bins = r_mesh[j,int(r_points/2)-1:])
            negative_counts_temp, _ = np.histogram(np.random.normal(mean_neg, np.abs(std_dev_neg), int(np.round(np.sum(negative_data)))), bins = r_mesh[j,:int(r_points/2)+1])
            positive_counts.append(positive_counts_temp)
            negative_counts.append(negative_counts_temp)
        
        full_deleted_mesh = np.hstack([negative_counts, positive_counts])

        #redistributing the remaining carriers along the z distribution
        number_remaining = number_deleted-np.sum(full_deleted_mesh)

        peaks_in_z = find_peaks(deleted_mesh[:,int(r_points/2)+1])[0]

        if make_plots:

            plt.plot(z_mesh[:,int(r_points/2)+1], deleted_mesh[:,int(r_points/2)+1])
            plt.xlabel('Depth into detector, z [m]', fontsize=15)
            plt.ylabel('Numer of redistributed charge carriers', fontsize=15)
            plt.ticklabel_format(style='sci', scilimits=(0.1,1))
            plt.show()         
        

        if len(peaks_in_z) > 1:
            print('Two peaks found.')
            try:
                distibution_size = abs(peaks_in_z[0]-peaks_in_z[1])
                positive_z = z_mesh[(peaks_in_z[0]+int(distibution_size/2))-1:,0]
                negative_z = z_mesh[:(peaks_in_z[1]-int(distibution_size/2)),0]

                positive_z_data = deleted_mesh[(peaks_in_z[0]+int(distibution_size/2))-1:,int(r_points/2)+1]
                negative_z_data = deleted_mesh[:(peaks_in_z[1]-int(distibution_size/2)),int(r_points/2)+1]

                if positive_z.size < 3 or negative_z.size<3:
                    raise RuntimeError


                popt_pos,_ = curve_fit(gaussian_fit, positive_z, positive_z_data,
                                                p0=[max(positive_z_data), positive_z[np.where(positive_z_data == np.max(positive_z_data))][0], 1e-6])
                popt_neg,_ = curve_fit(gaussian_fit, negative_z, negative_z_data,
                                                p0=[max(negative_z_data), negative_z[np.where(negative_z_data == np.max(negative_z_data))][0], 1e-6])

                mean_pos = popt_pos[1]
                std_dev_pos = popt_pos[2]
                mean_neg = popt_neg[1]
                std_dev_neg = popt_neg[2]

                positive_counts_temp, pos_bin_edges = np.histogram(np.random.normal(mean_pos, np.abs(std_dev_pos), int((number_remaining/2))), bins = positive_z[:])
                negative_counts_temp, neg_bin_edges = np.histogram(np.random.normal(mean_neg, np.abs(std_dev_neg), int((number_remaining/2))), bins = negative_z[:])
                
                if make_plots:
                    plt.bar(pos_bin_edges[:-1], positive_counts_temp, width=np.diff(pos_bin_edges), edgecolor="black", align="edge")
                    plt.bar(neg_bin_edges[:-1], negative_counts_temp, width=np.diff(neg_bin_edges), edgecolor="black", align="edge")
                    plt.xlabel('Depth into detector, z [m]', fontsize=15)
                    plt.ylabel('Numer of redistributed charge carriers', fontsize=15)
                    plt.ticklabel_format(style='sci', scilimits=(0.1,1))
                    plt.show()
                
            
            except RuntimeError:
                print("Curvefit failed")
                distibution_size = abs(peaks_in_z[0]-peaks_in_z[1])
                positive_z = z_mesh[(peaks_in_z[0]+int(distibution_size/2))-1:,0]
                negative_z = z_mesh[:(peaks_in_z[1]-int(distibution_size/2)),0]

                positive_z_data = deleted_mesh[(peaks_in_z[0]+int(distibution_size/2))-1:,int(r_points/2)+1]
                negative_z_data = deleted_mesh[:(peaks_in_z[1]-int(distibution_size/2)),int(r_points/2)+1]

                positive_counts_temp, _p = np.histogram(np.random.normal(positive_z[np.where(positive_z_data == np.max(positive_z_data))][0], np.abs(10e-6), int((number_remaining/2))), bins = positive_z)
                negative_counts_temp, _n = np.histogram(np.random.normal(negative_z[np.where(negative_z_data == np.max(negative_z_data))][0], np.abs(10e-6), int((number_remaining/2))), bins = negative_z)

            deleted_counts_z=np.append(negative_counts_temp, positive_counts_temp)

            if np.shape(deleted_counts_z) != np.shape(z_mesh[:,int(z_points/2)]):
                while True:
                    deleted_counts_z=np.append(deleted_counts_z,0)
                    break

        elif len(peaks_in_z) == 1:
            print('1 peak found')

            try:
                deleted_data = deleted_mesh[:,int(r_points/2)+1]
                z = z_mesh[:,int(r_points/2)+1]

                popt, pcov = curve_fit(gaussian_fit, z, deleted_data,
                                                p0=[max(deleted_data), z[np.where(deleted_data == np.max(deleted_data))][0], 1e-7])
                mean = popt[1]
                std_dev = popt[2]

                z_bins = np.append(z, z_lim+z_lim/z_points)
                
                deleted_counts_z, bin_edges = np.histogram(np.random.normal(mean, np.abs(std_dev), int((number_remaining))), bins = z_bins)
                
                if make_plots:
                    plt.bar(bin_edges[:-1], deleted_counts_z, width=np.diff(bin_edges), edgecolor="black", align="edge")
                    plt.show()
                
                
            except RuntimeError:
                try:
                    deleted_data = deleted_mesh[:,int(r_points/2)+1]
                    z = z_mesh[:,int(r_points/2)+1]
                    
                    popt, pcov = curve_fit(lorentzian_fit, z, deleted_data, p0 = [max(deleted_data), z[np.where(deleted_data == np.max(deleted_data))][0], 1e-7])

                    print('Gaussian curvefit in z unsuccessful. Lorentzian curvefit was successful.')
                    A = popt[0]
                    z_0 = popt[1] 
                    gamma = popt[2] 

                    z_bins = np.append(z, z_lim+z_lim/z_points)

                    lorentzian_values = lorentzian_fit(z, A, z_0, gamma)
                    
                    #decrease tail
                    alpha = 1
                    prob_weights = (lorentzian_values ** alpha) / np.sum(lorentzian_values ** alpha)

                    # Sample from z values according to power-law probabilities
                    sampled_z = np.random.choice(z, size=int(number_remaining), p=prob_weights)

                    # Create histogram
                    deleted_counts_z, bin_edges = np.histogram(sampled_z, bins=z_bins)                
                    
                    if make_plots:
                        plt.bar(bin_edges[:-1], deleted_counts_z, width=np.diff(bin_edges), edgecolor="black", align="edge")
                        plt.show()

                
                except RuntimeError:

                    print('Curvefit in z unsuccessful.')
                    
                    deleted_data = deleted_mesh[:,deleted_mesh[:,int(r_points/2)+1]]
                    mean = np.where(deleted_data==np.max(deleted_data))[0]
                    z_bins = np.append(z_mesh[:,int(r_points/2)+1], z_lim+z_lim/z_points)
                    std_dev = 1e-6
                    deleted_counts_z, _ = np.histogram(np.random.normal(mean, np.abs(std_dev), int((number_remaining))*2), bins = z_bins)

        else:

            print('Peaks cannot be found')

            deleted_data = deleted_mesh[:,int(r_points/2)+1]
            point = np.where(deleted_data==np.max(deleted_data))[0]
            deleted_counts_z = np.zeros(z_points)
            deleted_counts_z[point] = number_remaining

        full_deleted_mesh[:,int(r_points/2)] = full_deleted_mesh[:,int(r_points/2)] + deleted_counts_z
        ""
        #reintroducing the deleted carriers
        total_charge_carrier_mesh = charge_density_mesh_scaled + full_deleted_mesh

        return full_deleted_mesh, total_charge_carrier_mesh
    else:
        return charge_density_mesh_scaled
    
def plot_mesh(charge_density_mesh_r, charge_density_mesh_z, charge_density_mesh, contour_label, save_title):
    plt.clf()
    plt.pcolor(charge_density_mesh_r, charge_density_mesh_z, charge_density_mesh, cmap='BuPu')
    cbar = plt.colorbar()
    cbar.set_label(label=contour_label, size=14)
    plt.contour(charge_density_mesh_r, charge_density_mesh_z, charge_density_mesh, cmap='BuPu')
    plt.xlabel("Transverse Distance, r [m]", fontsize=15)
    plt.ylabel("Depth, z [m]", fontsize=15)
    plt.ticklabel_format(style='sci', scilimits=(0.01,100))
    plt.axhline(0, color='white', linestyle='dashed')
    plt.savefig('plots/{}.png'.format(save_title), bbox_inches='tight')
    #plt.show()
    plt.clf()

    return None

def plot_scatter(x_values, y_values, x_label, y_label, save_title):
    #plt.scatter(x_values, y_values, s=1)
    plt.plot(x_values, y_values)
    plt.xlabel(x_label, fontsize=15)
    plt.ylabel(y_label, fontsize=15)
    plt.ticklabel_format(style='sci', scilimits=(0.1,1))
    #plt.yscale('log')
    plt.savefig('plots/{}.png'.format(save_title), bbox_inches='tight')
    plt.show()
    
    return None

def plot_two(x_values, y_1, y_2, x_label, y_label, y_label_1, y_label_2, save_title):
    plt.plot(x_values, y_1, label=y_label_1, color = 'r')
    plt.plot(x_values, y_2, label=y_label_2, linestyle = "dashed", color = 'b')
    plt.xlabel(x_label, fontsize=15)
    plt.ylabel(y_label, fontsize=15)
    plt.ticklabel_format(style='sci', scilimits=(0.1,1))
    plt.legend()

    plt.savefig('plots/{}.png'.format(save_title), bbox_inches='tight')
    plt.show()
    
    return None

def main(r_lim, r_points, z_lim, z_points, positions, sim_type, make_plots, replace_deleted = True):

    #This section is used to define the number of voxel positions, and any quantities we want to examine at different positions
    

    position_array = np.linspace(0, positions, positions+1)
    #positions_array = [8.464040e-05]
    voxel_positions = np.zeros(len(position_array))
    
    normalisation_values = np.zeros(len(position_array))
    normalisation_values_no_smear = np.zeros(len(position_array))
    
    #amount of smearing
    sigma_array = [0]

    #This section is the loop. It returns voxels of charge density or number of charge carriers for a range of depths into the detector.

    for j in range(len(sigma_array)):
        #if sim_type == 'reflect and smear':
        if sim_type == 'reflect and smear':
            sigma = sigma_array[j]
            print('Sigma value:', sigma_array[j])

        pbar = progressbar.ProgressBar(widgets=[progressbar.Percentage(), progressbar.Bar()], maxval=positions+10)
        pbar.start()

        for i in range(len(position_array)):
            #V = (position_array[i]/positions)*z_lim*REF_INDEX_CORRECTION
            V = z_lim/2*REF_INDEX_CORRECTION
            voxel_positions[i] = V
            #V = positions_array[i]
            print("Voxel position: {:e}".format(V))

            #reflection and smearing
            chargeDensityMeshR, chargeDensityMeshZ, chargeDensityMesh = make_mesh_reflection_smearing(r_lim, r_points, z_lim, z_points, V, sigma, unit = True, return_normalisation=False, sim_type=sim_type)
            print('Total charge carriers before scaling: {:e}'.format(np.sum(chargeDensityMesh)))
            
            if replace_deleted:
                chargeCarrierMeshDeleted, chargeCarrierMeshScaled = scale_mesh_deleted(r_lim, r_points, z_lim, z_points, chargeDensityMesh, V, make_plots=False, replace_deleted=True)
            else:
                chargeCarrierMeshScaled = scale_mesh_deleted(r_lim, r_points, z_lim, z_points, chargeDensityMesh, V, make_plots=False, replace_deleted=True)
            
            print('Total charge carriers: {:e}'.format(np.sum(chargeCarrierMeshScaled)))
            #print('Redistributed charge carriers: {:e}'.format(np.sum(chargeCarrierMeshDeleted)))
            
            np.savetxt('voxels/charge_carrier_mesh_{:e}_{:e}.txt'.format(V, sigma), chargeCarrierMeshScaled)
            pbar.update(i+10)

            if make_plots:
                plot_mesh(chargeDensityMeshR, chargeDensityMeshZ, chargeDensityMesh, "Density of charge carriers", "active/smeared_charge_density_mesh_{:e}_{:e}".format(V, sigma))
                plot_mesh(chargeDensityMeshR, chargeDensityMeshZ, chargeCarrierMeshScaled, "Number of charge carriers", "active/smeared_charge_carrier_mesh_{:e}_{:e}".format(V, sigma))
                if replace_deleted:
                    plot_mesh(chargeDensityMeshR, chargeDensityMeshZ, chargeCarrierMeshDeleted, "Number of deleted charge carriers", "active/deleted_charge_carrier_mesh_{:e}_{}".format(V, sigma))
            

            normalisation_values[i] \
                = charge_density_normalisation_with_reflection_and_smearing(r_lim, z_lim, V, sigma)[0]
            print('Normalisation constant: {:e}'.format(normalisation_values[i]))
            
            #normalisation_values_no_smear[i] \
            #    = charge_density_normalisation_with_reflection_and_smearing(r_lim, z_lim, V, 0)[0]
            print()
    
        pbar.finish()

        np.savetxt('voxels/normalisation/normalisation_constants_{:e}.txt'.format(sigma), normalisation_values)
        np.savetxt('voxels/positions/voxel_positions_{:e}.txt'.format(sigma), voxel_positions)        

    #plot_scatter((position_array/positions)*z_lim, voxel_sizes, "Depth into detector, V [m]", "Size of voxel [indices]", "normalisation/active_sizes")
    plot_scatter((position_array/positions)*z_lim, normalisation_values, "Voxel position, V [m]", "Normalisation constant [C]", "normalisation/active_normalisation_constants")
    #plot_two((position_array/positions)*z_lim, normalisation_values_no_smear, normalisation_values, "Voxel position, V [m]", "Normalisation constant [C]", 'No smearing', 'Smearing (1e-6)', "normalisation/normalisation_constants_compare_{:e}".format(sigma))
    #plot_scatter(position_array, max_values_reflect, "Depth into detector index", "Maximum charge density", "maximum_values_reflection")

    #plot_two((position_array/positions)*z_lim, max_values_direct, max_values_reflect, "Depth into detector, V [m]", "Maximum charge density", "Direct component only", "Including reflection", "maximum_values_both")

    #plot_scatter((position_array/positions)*z_lim, charge_density_with_reflection(0, (position_array/positions)*z_lim, 160*10**(-6)), "Depth into detector, z [m]", "Charge density", "charge_density_plot")
    return None

main(25*10**(-6), 251, 50*10**(-6), 500,1, sim_type='reflect and smear', make_plots=True)