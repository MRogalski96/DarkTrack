# DarkTrack
Repository for simulating and reconstructing lensless digital in-line holographic microscopy (DIHM) data
## Contents
- MSHoloSim.m - engine for simulating a set of lensless digital in-line holographic microscopy (DIHM) holograms with respect to multiple scattering for a group of microbeads with given 4D (space + time) locations and diameters
- DarkTrack.m - algorithm for reconstructing the 4D (space + time) location of objects/particles present in input set of DIHM holograms. Additionally it performs the extended depth of focus (EDOF) reconstruction (with all objects moved to the focus plane and minimized twin-image effect)
- UsageExample.m - script that demonstrates how to use the MSHoloSim and DarkTrack algorithms. Requires exemplary data to run (see below)

## Exemplary data
Our exemplary datasets may be downloaded here:<br>
https://drive.google.com/drive/folders/1UNtZ3IeEX5ms_Vx85b0S685D-uipLe_l?usp=sharing

## How to cite
M. Rogalski, J. A. Picazo-Bueno, J. Winnik, P. Zdańkowski, V. Micó, and M. Trusiak, “Accurate automatic object 4D tracking in digital in-line holographic microscopy based on computationally rendered dark fields,” Scientific Reports, vol. 12, no. 1, p. 12909, Dec. 2022, doi: 10.1038/s41598-022-17176-1.

## Created by
Mikołaj Rogalski, <br>
mikolaj.rogalski.dokt@pw.edu.pl <br>
Institute of Micromechanics and Photonics, <br>
Warsaw University of Technology, 02-525 Warsaw, Poland
