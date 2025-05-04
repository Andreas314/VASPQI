#!/usr/bin/env python3

import argparse
import sys
import subprocess
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
def get_args():
    parser = argparse.ArgumentParser(
            description = 'Input arguments for VASPQI:',
            usage = 'Compute quantum interference tensor to get current density from VASP\'s WAVECAR and OUTCAR. '
            )
    parser.add_argument( '--omega',
                       type = OMEGA_RESTRICTION, 
                       required = True,
                       help = 'Frequency of the incident light (eg. 10E16 s^-1)',
                       )
    parser.add_argument( '--file_name',
                       type = str, 
                       required = False,
                        help = 'Name of the file to which the output is written (DEFAULT: QI_TENSOR)',
                        default = 'QI_TENSOR'
                       )
    parser.add_argument( '--directory_name',
                       type = str, 
                       required = False,
                        help = 'Name of the directory to which the output files are saved (DEFAULT: VASP2QI_RESULTS)',
                        default = 'VASP2QI_RESULTS'
                       )
    parser.add_argument( '--source',
                       type = str, 
                       required = False,
                        help = 'Name of the directory which the inputs are stored in (DEFAULT: .)',
                        default = '.'
                        )
    parser.add_argument( '--conduction_bands',
                       type = BAND_RESTRICTION, 
                       required = False,
                        help = 'How many bands above the fermi level are taken into the calculations (DEFAULT: 100)',
                        default = 100
                        )
    parser.add_argument( '--valence_bands',
                       type = BAND_RESTRICTION, 
                       required = False,
                        help = 'How many bands below the fermi level are taken into the calculations (DEFAULT: 100)',
                        default = 100
                        )
    parser.add_argument( '--number_of_processes',
                        type = NP_RESTRICTION,
                        required = False,
                        help = 'Number of processes (parallel run) to which the kpoint loop is divided (DEFAULT: 1)',
                        default = 1
                        )
    parser.add_argument('--omega_run',
                        action = 'store_true',
                        help = 'Switch for a loop over multiple frequencies. Starts at --omega. When on, --omega_step and --omega_max must be set as well')
    parser.add_argument( '--omega_step',
                        type = OMEGA_RESTRICTION,
                        required = False,
                        help = 'Size of one step between --omega and --omega_max in --omega_run (relative to --omega)',
                        )
    parser.add_argument( '--omega_max',
                        type = OMEGA_RESTRICTION,
                        required = False,
                        help = 'Maximal frequency for which to compute the interference tensor (relative to --omega)',
                        )
    parser.add_argument('--exclude_k',
                        action = 'store_true',
                        help = 'Whether to exclude the k terms in the momentum matrices')
    
    arguments =  parser.parse_args()
    Process = subprocess.run(args=['./VASP2QI_kparsing.sh', arguments.source], capture_output = True)
    if Process.returncode != 0:
        print(Process.stderr)
        sys.exit(Process.returncode)
    return Process.stdout, arguments

