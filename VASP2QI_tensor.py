#!/usr/bin/env python3
import os
import sys
import numpy as np
from multiprocessing import Pool, Value, Lock
from time import time
sys.path.append(os.path.abspath("/home/piratmori28/Desktop/Thesis/VASPQI"))
from VASP2QI_parser import get_args
try:
    #A library to extract bands, kpoints and momentum matrices from WAVECAR, see: https://github.com/QijingZheng/VaspBandUnfolding/blob/master/vaspwfc.py
    import vaspwfc as vp
except ImportError:
    print("Error: vaspwfc not found! Try installing VaspBandUnfolding with pip install git+https://github.com/QijingZheng/VaspBandUnfolding.", file = sys.stderr)
    sys.exit(1)

#Physical constants
H_PLANC = 6.582119569E-16 #Reduced Planc constant in (eV/s)
M_ELECTRON = 9.1093837E-31 #Electron mass in (kg)
E_CHARGE = 1.602176634E-19 #Elementary charge in (C)

def Read_input():
    '''
    Read name of a directory containing WAVECAR, weights of individual kpoints from OUTCAR (which cannot be read using vaspwfc, so they are extracted using a bash script) and a number of kpoints.
    '''

    file = arguments.source + '/WAVECAR'
    wavecar_data = vp.vaspwfc(fnm = file, lsorbit = True)
    num_kpoints = len(wavecar_data._kvecs)
    k_weights = Make_weights(input_weights, num_kpoints)
    return k_weights, wavecar_data, num_kpoints

def Correct_moment_mat(ks_i, ks_j, vasp_wfc):
    # ks_i and ks_j are list containing spin-, kpoint- and band-index of the

        # k-points in direct coordinate
        if not arguments.exclude_k:
            k0 = vasp_wfc._kvecs[ks_i[1] - 1]
        else:
            k0 = 0
        # plane-waves in direct coordinates
        G0 = vasp_wfc.gvectors(ikpt=ks_i[1])
        # G + k in Cartesian coordinates
        Gk = np.dot(
            G0 + k0,                            # G in direct coordinates
            vasp_wfc._Bcell * 2 * np.pi            # reciprocal basis x 2pi
        )

        # plane-wave coefficients for initial (mk) and final (nk) states
        CG_mk = vasp_wfc.readBandCoeff(*ks_i)
        CG_nk = vasp_wfc.readBandCoeff(*ks_j)
        ovlap = CG_nk.conj() * CG_mk

        ################################################################################
        # Momentum operator matrix element between pseudo-wavefunctions
        ################################################################################
        if vasp_wfc._lgam:
            # for gamma-only, only half the plane-wave coefficients are stored.
            # Moreover, the coefficients are multiplied by a factor of sqrt2

            # G > 0 part
            moment_mat_ps = np.sum(ovlap[:,None] * Gk, axis=0)

            # For gamma-only version, add the other half plane-waves, G_ = -G
            # G < 0 part, C(G) = C(-G).conj()
            moment_mat_ps -= np.sum(
                    ovlap[:,None].conj() * Gk,
                    axis=0)

            # remove the sqrt2 factor added by VASP
            moment_mat_ps /= 2.0

        elif vasp_wfc._lsoc:
            moment_mat_ps = np.sum(
                ovlap[:, None] * np.r_[Gk, Gk],
                axis=0)
        else:
            moment_mat_ps = np.sum(
                ovlap[:,None] * Gk, axis=0
            )

        return moment_mat_ps

def Make_weights(input_data, size):
    '''
    Process the input string of the weights of individual kpoints.
    '''
    weights = np.zeros(size,int)
    data_split = input_data.split()
    for ii in range(0,size):
        try:
            weights[ii] = int(float(data_split[ii]))
        except ValueError:
            print("Error: Cannot interpret ", data_split[ii], " as an integer. OUTCAR corrupted?", file = sys.stderr)
            sys.exit(1)
    return weights

def Progres_bar(act_index, num_kpoints):
    '''
    Draw a progres bar to the stdout stating the percentage of sum over k-space computed
    '''
    progres_made = int(act_index/num_kpoints * 100)
    print("\033[F", end = "", file = sys.stderr)
    print("\033[F", end = "", file = sys.stderr)
    print("Sum in kspace started: ", progres_made, " % done.", file = sys.stderr)
    print("|", end = "", file = sys.stderr)
    for ii in range(100):
        if ii < progres_made:
            print("\u2588", end = "", file = sys.stderr)
        else:
            print(" ", end = "", file = sys.stderr)
    print("|", end = "\n", file = sys.stderr)

def Get_gamma(kvecs):
    dist = abs(kvecs[1]) - abs(kvecs[0])
    omega_dist = np.linalg.norm(dist) * H_PLANC / 2 / M_ELECTRON * E_CHARGE * 10**20
    return omega_dist

def Enter_Sum_Wrapper(arguments_in):
    '''
    Wrapper around the sum, which allows for a parallel run utilizing the multiprocessing library
    Two global variables declared here:
    lock: Allows for an update of the progres bar one worker at a time.
    num: number of k-points already computed, shared variable between all workers
    '''
    global arguments
    global input_weights
    global lock
    global num
    arguments = arguments_in[1]
    input_weights = arguments_in[0]
    lock = Lock()
    num = Value('i', 0)
    NP = arguments.number_of_processes
    #Does not overwrite the user input in the shell
    #print("\n") 
    print("\n", file = sys.stderr) 
    with Pool(NP) as p:
        return sum(p.map(Enter_Sum, range(0,NP)))

def Find_valence_cond(bands, num_bands, efermi):
    '''
    Find indices of valnece and condution bands by comparing their energies in gamma point to the fermi energy
    Also takes only as many bands as specified (if specified number is smaller 
    than the number of bands from the VASP calculation)
    '''
    conduction = []
    valence = []
    for ii in range(0, num_bands):
        if (bands[0, 0, ii] < efermi):
            valence.append(ii)
        else:
            conduction.append(ii)
    if len(conduction) > arguments.conduction_bands:
        conduction = conduction[0 : arguments.conduction_bands]
    if len(valence) > arguments.valence_bands:
        num_valence = len(valence)
        index = num_valence - arguments.valence_bands
        valence = valence[index : num_valence]
    return valence, conduction

def Enter_Sum(index):
    '''
    Implementation of https://journals.aps.org/prb/abstract/10.1103/PhysRevB.68.085208
    Start the sum over kspace. The actual sum contains 7 for loops, so they are broken into smaller functions.
    First read the WAVECAR to wavecar_data and get band energies, k-points for the calculation and number of bands.
    Frequency of the incident light omega is read from stdin. 
    The delta function in the actual sum is replaced by Lorentzian centered at omega. Its width is a fraction of omega.
    '''
    #Read WAVECAR 
    [k_weights, wavecar_data, num_kpoints] = Read_input()
    [ _ , energies] = wavecar_data.readWFBand()
    k_vects = wavecar_data._kvecs
    num_bands = wavecar_data._nbands
    [valence_states, conduction_states] = Find_valence_cond(energies, num_bands, wavecar_data._efermi)
    gamma = Get_gamma(k_vects)
    qi_tensor = np.zeros([3,3,3,3], complex)
    omega = arguments.omega
    for k in range(int(index), num_kpoints, arguments.number_of_processes):
        if k >= arguments.number_of_processes:
            num.value += 1
        lock.acquire()
        try: 
            pass
            Progres_bar(num.value, num_kpoints)
        finally:
            lock.release()
        gap_at_k = (energies[0,k,min(conduction_states)] - energies[0,k,max(valence_states)]) / H_PLANC
       # if (gap_at_k  > 2 * omega):
       #     continue
        qi_tensor += Band_Sum(k, energies[0,k,:], num_bands, k_weights[k], gamma, wavecar_data, valence_states, conduction_states)
    num.value += 1
    lock.acquire()
    try:    
        pass
        Progres_bar(num.value, num_kpoints)
    finally:
        lock.release()

    volume = wavecar_data._Omega * 10**(-30) #Volume of the supercell used in VASP calculations in m^3
    prefactor = (2 * np.pi) **3 / (sum(k_weights) * volume)  / 4 / np.pi ** 2 * 1j * (E_CHARGE / M_ELECTRON)**4 * H_PLANC * 10**(40) * E_CHARGE 
    qi_tensor *= prefactor
    return qi_tensor

def Band_Sum(k_index, bands_energies, n_bands, weight, gamma, wf_obj, valence_states, conduction_states):
    '''
    Sum over all valence, conduction and intemediate sattes
    '''
    #delta = max(np.matmul(wf_obj._kvecs[1], wf_obj._Bcell) - np.matmul(wf_obj._kvecs[0], wf_obj._Bcell))
    #gamma = 2 * np.linalg.norm(np.matmul(wf_obj._kvecs[k_index], wf_obj._Bcell)) * delta
    #gamma += delta ** 2
    #gamma *= (H_PLANC / 2 / M_ELECTRON * E_CHARGE * 10**20)
    inner_tensor_output = np.zeros([3,3,3,3], complex)
    omega = arguments.omega
    all_bnds = valence_states + conduction_states
    for initial in valence_states:
        for final in conduction_states:
            omega_fv = (bands_energies[final] - bands_energies[initial]) / H_PLANC
            for inter in all_bnds:
                omega_jv = (bands_energies[inter] - bands_energies[initial]) / H_PLANC
                inner_tensor = Get_all_elements(wf_obj, initial, final, inter, omega_fv, k_index, wf_obj)
                inner_tensor *= Lorentzian(omega_fv, 2 * omega, gamma)
                inner_tensor /= (omega_fv / 2 - omega_jv) * ( omega_fv / 2 )**3
                inner_tensor *= weight
                inner_tensor_output += inner_tensor
    return inner_tensor_output

def Get_all_elements(wf_obj, init, fin, iner, omega_fv, k_index, wf):
    '''
    Compute all 81 elements of the tensor
    '''
    inner_inner_tensor = np.zeros([3,3,3,3], complex)
    #Indices in vaspwfc start with 1
    k_index += 1
    fin += 1
    init += 1
    iner += 1
    '''
    Momentum matrices in vaspwfc are calculated as ~ <W_n| k + G| W_m>, where |W_n> is the pseudo wavefunction form VASP, G is a point in reciprocal lattice
    and k is wavevector. We dont want to consider the part with k in QI calculations, so the p matrices have to be recalculated
    '''
    p_vf = Correct_moment_mat([1, k_index, init], [1, k_index, fin], wf)
    p_fj = Correct_moment_mat([1, k_index, fin], [1, k_index, iner], wf)
    p_jv = Correct_moment_mat([1, k_index, iner], [1, k_index, init], wf)
    p_vv = Correct_moment_mat([1, k_index, init], [1, k_index, init], wf)
    p_ff = Correct_moment_mat([1, k_index, fin], [1, k_index, fin], wf)
    #Sum over x,y,z coordinates
    for i1 in range(0,3):
        for i2 in range(0,3):
            for i3 in range(0,3):
                for i4 in range(0,3):
                    inner_inner_tensor[i1][i2][i3][i4] = (p_vv[i1] - p_ff[i2]) * p_vf[i2] * (p_fj[i3] * p_jv[i4] + p_fj[i4] * p_jv[i3])/  2
    return inner_inner_tensor

def Lorentzian(x, x0, gamma):
    return 1 / np.pi * gamma / ( (x - x0)**2 + gamma**2 )
