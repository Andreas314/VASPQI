#!/usr/bin/env python3
import sys
import os
import numpy as np
sys.path.append(os.path.abspath("/home/piratmori28/Desktop/Thesis/VASPQI"))
from VASP2QI_parser import get_args
from VASP2QI_tensor import Enter_Sum_Wrapper

if __name__ == '__main__':
    arguments = get_args()
    omega_run = arguments[1].omega_run
    if (omega_run):
        if (arguments[1].omega_step == None) or (arguments[1].omega_max == None):
            print('Warning: When --omega_run is on, --omega_step and --omega_max must be set! \nAbonding the frequency loop!', file = sys.stderr)
            tensor = Enter_Sum_Wrapper(arguments)
            print(arguments[1].omega, abs(tensor[0,0,0,0]), abs(tensor[0,0,1,1]), abs(tensor[0,1,1,0]))
        elif (arguments[1].omega_max < 1):
            print('Warning: --omega must not be greater than --omega_max! \nAbonding the frequency loop!', file = sys.stderr)
            tensor = Enter_Sum_Wrapper(arguments)
            print(arguments[1]/omega, abs(tensor[0,0,0,0]), abs(tensor[0,0,1,1]), abs(tensor[0,1,1,0]))
        else:
            print('Entering a loop over frequencies:',file = sys.stderr)
            omega = arguments[1].omega
            omega_step = arguments[1].omega_step * omega
            omega_max = arguments[1].omega_max * omega
            print('Going from ', '%.2E' % omega, 'to ', '%.2E' % omega_max, 'with a step of ', '%.2E' % omega_step, file = sys.stderr)
            freq = np.arange(omega, omega_max + omega_step, omega_step)
            ii = 1
            for om in freq:
                print('', file = sys.stderr)
                print('Calculation: ', ii, '/', len(freq),file = sys.stderr)
                print('Omega: ', '%.2E' % om,file = sys.stderr)
                tensor = Enter_Sum_Wrapper(arguments)
                print(om, abs(tensor[0,0,0,0]), abs(tensor[0,0,1,1]), abs(tensor[0,1,1,0]))
                arguments[1].omega = om
                ii += 1
    else:
        tensor = Enter_Sum_Wrapper(arguments)
        print(arguments[1].omega, abs(tensor[0,0,0,0]), abs(tensor[0,0,1,1]), abs(tensor[0,1,1,0]))
