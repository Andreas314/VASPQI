#!/usr/bin/env python3
import sys
import os
import numpy as np
#sys.path.append(os.path.abspath("/home/piratmori28/Desktop/Thesis/VASPQI"))
from VASP2QI_parser import get_args
from VASP2QI_tensor import Enter_Sum_Wrapper
from VASP2QI_filewriter import Write_tensor
if __name__ == '__main__':
    arguments = get_args()
    Write_tensor(arguments)
