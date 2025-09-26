# VASP2QI

VASP2QI is a Python tool for computing the **Quantum Interference (QI) tensor** from VASP output files. The code processes `KPOINTS`, `WAVECAR`, and `OUTCAR` files to evaluate photocurrent generation in bulk semiconductors as described in *Coherent Control of Photocurrent Generation in Bulk Semiconductors*.

---

## Features
- Reads **VASP output files** (`KPOINTS`, `WAVECAR`, `OUTCAR`).
- Computes the **QI tensor** for bulk crystals.
- Supports calculation at a single frequency or over a range of frequencies.
- Parallelizable over k-points.
- Flexible options for band and frequency selection.

---

## Installation

Clone the repository and install the dependencies:

```bash
git clone https://github.com/Andreas314/VASPQI.git
cd VASP2QI
pip install -r requirements.txt
```

---

## Usage

```bash
./VASP2QI.py [OPTIONS]
```

### Required Inputs
Place the following files in the input directory (default: `.`):
- `KPOINTS`
- `WAVECAR`
- `OUTCAR`

### Options
```
  -h, --help            Show this help message and exit

  --omega OMEGA         Frequency of the incident light (in eV)

  --file_name FILE_NAME
                        Name of the file to which the output is written
                        (DEFAULT: QI_TENSOR)

  --directory_name DIRECTORY_NAME
                        Directory where output files are saved
                        (DEFAULT: VASP2QI_RESULTS)

  --source SOURCE       Directory containing input files
                        (DEFAULT: .)

  --conduction_bands CONDUCTION_BANDS
                        Number of conduction bands included above the
                        Fermi level (DEFAULT: 100)

  --valence_bands VALENCE_BANDS
                        Number of valence bands included below the
                        Fermi level (DEFAULT: 100)

  --np NUMBER_OF_PROCESSES
                        Number of processes for parallel k-point loops
                        (DEFAULT: 1)

  --omega_run           Switch for running over multiple frequencies.
                        Requires --omega_step and --omega_max.

  --omega_step OMEGA_STEP
                        Frequency step size for --omega_run (in eV)

  --omega_max OMEGA_MAX
                        Maximum frequency for --omega_run (in eV)

  --exclude_k           Exclude k terms in momentum matrices

  --silence             Suppress progress bar output
```

## Output
- The program writes the QI tensor to the specified output file (`QI_TENSOR.dat` by default).
- If `--omega_run` is enabled, a series of tensors are generated for frequencies between `--omega` and `--omega_max` with steps of `--omega_step`.
- All results are saved in the chosen output directory (`VASP2QI_RESULTS` by default).

---

## Notes
- This implementation is currently valid **only for bulk crystals**.
- The QI tensor is computed as described in *Coherent Control of Photocurrent Generation in Bulk Semiconductors*.

---

## License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Citation
If you use this tool in your research, please cite the original theoretical framework:

> Atanasov, R.; Haché, A.; Hughes, J. L. P.; Driel, H. M. van; Sipe,
J. E. Coherent Control of Photocurrent Generation in Bulk Semiconductors.
Phys. Rev. Lett. 1996, vol. 76, pp. 1703–1706. Available from doi: 10.1103/
PhysRevLett.76.1703.
