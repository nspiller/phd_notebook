# Generate files for chemshell from protonated PDB file

When following the steps below, one can create the files necessary to run a chemshell calculation from a PDB file. 
This way one can manually mutate residues in a converged QMMM calculation and feed the modified protein into chemshell.
The procedure follows largely the [QMMM setup](https://sites.google.com/site/ragnarbjornsson/mm-and-qm-mm-setup), 
but has been modified to work with the file format generated _from_ a QMMM calculation.

Required:
- `protein.pdb`: e.g. from a converged calculation with the mutation
- `psfgen`: standalone executable for [psfgen plugin](https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/) 
(unfortunately not available on their homepage anymore)
- CHARMM topology file `top_all36_prot.rtf` with all residues in the PDB file
- `qmatoms` and `act` files or knowledge how to generate those 
- `chemshell` to create the fragment file (not necessary for PSF file)

Creates:
- `newxplor.psf`: new protein structure file
- `save-new.chm`: lists for chemshell including same information in the CHARMM format 
- `system.c`: fragment file with accurate coordinates

1. check format of PDB file
	- The PDB file requiresuse a very specific format. Check with`chemshell pdb4psfgen.chm FRAGMENT PSF` if the PDB file has the exact same column widths.
2. Cut PDB file 
	- Modify PDB file name in `cut_pdb.zsh`.
	- Define `cofactornames` and `extraspecies` in `cut_pdb.zsh` according to the PDB file.
		- note: Order is important here if you want to reuse the `act` and `qmatoms` files.
		- note: Script will not work if deprotonated CYS are called CSD in PDB file
	- Copy `psfgen-template.tlc` to WDIR. Run `./cut_pdb.zsh`.
3. Create new PSF and PDB
	- Modify `psfgen.tcl` to point to the correct topology file.
	- Run `psfgen psfgen.tcl >& psfgen.out`.
	- note: The created files are: `new.pdb`, `new.psf` (CHARMM format), and `newxplor.psf` (x-plor format)
	- note: Initial (`result.pdb`) and final PDB files (`new.pdb`) should only differ only in the alignment of the atom name name and length of segment ID
4. Create `save-new.chm`
	- Run `./chemshell-listprep.py` 
5. Create FRAG file
	- Run `chemshell create_fragment.chm` 
	- Run `./correct_fragment.py system.c new.pdb` to correct element names in fragment file
6. update `qmatoms` and `act`
	- option 1: redefine QM and active using e.g. `QMregiondefine.chm` or `actregiondefine.chm`
	- option 2: update old files directly
		- determine how old system differs from the new. 
			- e.g. the mutation deleted atoms 100 and 101: N=99 X=-2
			- or the mutation introduced 3 additional atoms after position 100: N=100 X=3
			- default numbering starts with 1 which is used by chemshell and PDB file format
		- run `update_region.py qmatoms X N` and `update_region.py act X N`
		- note: for multiple changes run several times with N1 > N2 > N3 etc.
