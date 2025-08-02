#!/usr/bin/env python3
import sys
import os
import numpy as np
#sys.path.append(os.path.abspath("/home/piratmori28/Desktop/Thesis/VASPQI"))
from VASP2QI_parser import get_args
from VASP2QI_tensor import Enter_Sum_Wrapper

H_PLANC = 6.582119569E-16 #Reduced Planc constant in (eV/s)

def Write_tensor(arguments):
    destination = arguments[1].directory_name
    file_name = arguments[1].file_name
    file_name = Check_for_dir_file(file_name, destination)
    signs_1 = ['X', 'Y', 'Z']
    signs_2 = ['XX', 'YY', 'ZZ', 'YZ', 'XZ', 'XY']
    line = 'Omega' + ' '*15
    with open(file_name, "a") as file:
        for a1 in range(0,3):
            for a2 in range(0,3):
                for a3 in range(0,6):
                    element = ' '*24 + signs_1[a1] + signs_1[a1] + signs_2[a2] + ' '*24
                    line += element
        file.write(line)

        omega_run = arguments[1].omega_run
        if (omega_run):
            if (arguments[1].omega_step == None) or (arguments[1].omega_max == None):
                Write_omega_max_step_missing(file, arguments)
            elif (arguments[1].omega_max < 1):
                Write_wrong_omega_max(file, arguments)
            else:
                Write_multiple_omega(file, arguments)
        else:
            Write_single_omega(file, arguments)

def Write_omega_max_step_missing(f, arg):
    print('Warning: When --omega_run is on, --omega_step and --omega_max must be set! \nAbonding the frequency loop!', file = sys.stderr)
    tensor = Enter_Sum_Wrapper(arg)
    Write_to_file(f, arg[1].omega, tensor)

def Write_wrong_omega_max(f, arg):
    print('Warning: --omega must not be greater than --omega_max! \nAbonding the frequency loop!', file = sys.stderr)
    tensor = Enter_Sum_Wrapper(arg)
    Write_to_file(f, arg[1].omega, tensor)

def Write_single_omega(f, arg):
    tensor = Enter_Sum_Wrapper(arg)
    Write_to_file(f, arg[1].omega, tensor)

def Write_multiple_omega(f, arg):
    print('Entering a loop over frequencies:',file = sys.stderr)
    omega = arg[1].omega
    omega_step = arg[1].omega_step * omega
    omega_max = arg[1].omega_max * omega
    print('Going from ', '%.2E' % omega, 'to ', '%.2E' % omega_max, 'with a step of ', '%.2E' % omega_step, file = sys.stderr)
    freq = np.arange(omega, omega_max + omega_step, omega_step)
    ii = 1
    for om in freq:
        print('', file = sys.stderr)
        print('Calculation: ', ii, '/', len(freq),file = sys.stderr)
        print('Omega: ', '%.2E' % om,file = sys.stderr)
        tensor = Enter_Sum_Wrapper(arg)
        Write_to_file(f, arg[1].omega, tensor)
        arg[1].omega = om
        ii += 1

def Write_to_file(file, omega, tensor):
    line = f"{omega * H_PLANC:<20.6f}"
    for a1 in range(0,3):
        for a2 in range(0,3):
            for a3 in range(0,6):
                element = tensor[a1][a2][a3]
                print(element.real)
                line += f"{element.real:>20.5e}{'+' if element.imag >= 0 else '-'}i{abs(element.imag):<20.5e}"
    print(line)
    file.write(line)


def Check_for_dir_file(file, directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    if os.path.exists(directory+'/'+file + '.dat'):
        ii = 1
        while True:
            newfile = file + "_" + str(ii)
            if not os.path.exists(directory + '/'+ newfile + '.dat'):
                file = newfile
                break
            ii += 1
    file += '.dat'
    return directory + '/' + file


