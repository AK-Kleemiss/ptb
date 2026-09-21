# PTB
A density matrix (P) tight-binding (TB) method based on a polarized valence double-zeta basis set.
It is available for all elements and structures until Z = 86.

In this development version, before execution, two files have to be given:
The `.atompara` file containing all empirical parameters. If no specific location is given via a `-par` flag, it is assumed to be located in the `$HOME` directory (e.g., `~/.atompara`).
An individual location can be defined via the `-par <path of .atompara>` command.
The `.basis_vDZP` file contains the vDZP basis set in the correct format. Its file location can be defined via the `-bas` flag, otherwise it is assumed to be in 
`$HOME` (e.g., `~/.basis_vDZP`).

Molecular (total) charge can be incorporated via the presence of a `.CHRG` file. A similar procedure follows with the number of unpaired electrons (`.UHF`),
eventhough PTB is mainly developed for closed-shell systems.

## Deviation from upstream vDZP: Yb second d shell (2026-09-21)

This fork's `.basis_vDZP` (and the copy embedded in `source/default_files.f90`) differs from the published vDZP in one place: the second d shell of Yb (Z=70). Upstream vDZP (verified against ORCA 6.1.1's built-in copy) uses exponents 0.0487 / 0.00207; the 0.00207 primitive is ~25x more diffuse than the corresponding Tm or Lu functions and cannot be handled by PTB's fixed tight-binding Hamiltonian (it produced positive-energy occupied MOs and densities not converged within 30 A). It was replaced by the geometric mean of the Tm and Lu exponents, 0.1269 / 0.0511, and the Yb `.atompara` block was refitted against wB97X-3c reference densities on that basis. The startup banner prints a note about this. Results for Yb are therefore not comparable to upstream PTB; all other elements are unchanged.


## Building
You can use a statically-linked release binary (recommended), but you can also build it with the source code.
```
git clone https://github.com/grimme-lab/ptb.git
cd source
```
You can the build the project via `make`.
