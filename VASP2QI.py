#!/usr/bin/env python3

import argparse
import sys
def BZ_RESTRICTION(value):
        value = float(value)
        if (value <= 0) or (value > 1):
            raise argparse.ArgumentTypeError('Value for the restriction of the BZ must be:  0 < k <= 1!')
        return value
def BAND_RESTRICTION(value):
    value = int(value)
    if value < 0:
        raise argparse.ArgumentTypeError('Value for the restriction of number of band must positive!')
    return value
def OMEGA_RESTRICTION(value):
    value = float(value)
    if value < 0:
        raise argparse.ArgumentTypeError('Frequency must be positive!')
    return value
def NP_RESTRICTION(value):
    value = int(value)
    if value < 0:
        raise argparse.ArgumentTypeError('Number of processes must be positive!')
    return value
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description = 'Input arguments for VASPQI:',
            usage = 'Compute quantum interference tensor to get current density from VASP\'s WAVECAR and OUTCAR. '
            )
    parser.add_argument('-o', '--omega',
                       type = OMEGA_RESTRICTION, 
                       required = True,
                       help = 'Frequency of the incident light (eg. 10E16 s^-1)',
                       )
    parser.add_argument('-f', '--file_name',
                       type = str, 
                       required = False,
                        help = 'Name of the file to which the output is written (DEFAULT: QI_TENSOR)',
                        default = 'QI_TENSOR'
                       )
    parser.add_argument('-d', '--directory_name',
                       type = str, 
                       required = False,
                        help = 'Name of the directory to which the output files are saved (DEFAULT: VASP2QI_RESULTS)',
                        default = 'VASP2QI_RESULTS'
                       )
    parser.add_argument('-s', '--source',
                       type = str, 
                       required = False,
                        help = 'Name of the directory which the inputs are stored in (DEFAULT: .)',
                        default = '.'
                        )
    parser.add_argument('-c', '--conduction_bands',
                       type = BAND_RESTRICTION, 
                       required = False,
                        help = 'How many bands above the fermi level are taken into the calculations (DEFAULT: 100)',
                        default = 100
                        )
    parser.add_argument('-v', '--valence_bands',
                       type = BAND_RESTRICTION, 
                       required = False,
                        help = 'How many bands below the fermi level are taken into the calculations (DEFAULT: 100)',
                        default = 100
                        )
    parser.add_argument('-k', '--kpoints_restriction',
                       type = BZ_RESTRICTION, 
                       required = False,
                        help = 'What portion of the Brillouin zone is taken into the calculation in each direction (DEFAULT: 1 1 1)',
                        default = [1,1,1],
                        nargs = 3
                        )
    parser.add_argument('-np', '--number_of_processes',
                        type = NP_RESTRICTION,
                        required = False,
                        help = 'Number of processes (parallel run) to which the kpoint loop is divided (DEFAULT: 1)',
                        default = 1
                        )
    arguments =  parser.parse_args()
    print(arguments)


