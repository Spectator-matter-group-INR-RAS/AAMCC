# Abrasion-Ablation Monte Carlo for Colliders or AAMCC

Abrasion-Ablation Monte Carlo for Colliders (AAMCC) is Monte-Carlo model specially designed to describe spectator matter production in a wide range 
of colliding nuclei and energy of the collision on the event-by-event basis

## Main assumptions of AAMCC
 
 - Nucleus-nucleus collisions are simulated by means of the Glauber Monte Carlo model. Non-participated nucleons form spectator matter (prefragment).
 - Excitation energy of prefragment can be calculated via one of 4 parametrization: 
   - Ericson formula
   - Gamaird-Schmidt formula
   - ALADIN collaboration parametrization
   - Hybrid approach. Excitation energy is calculated as follows:
     - in peripheral collisions with less then ~15% of removed nucleons the particle-hole model is used 2) (Ericson formula);
     - otherwise a parabolic ALADIN approximation 3) is applied with parameters tuned to data obtained in nuclear emulsions.
 - Decays of prefragments are simulated as follows:
   - pre-equilibrium decays modeled with MST-clustering algorithm;
   - Fermi break-up model from Geant4 v9.2;
   - Multifragmentation is calculated via SMM in realization of Geant4 v10.4.
   - Evaporation-fission calculation is carried out via Weisskopf-Ewing model adopted from Geant4 v10.4.
   
## How to use AAMCC

- To compile AAMCC one will need **Geant4 v10.4** and **ROOT v6.*** installed with CMake or as a part of FairROOT/AliROOT/etc. and [**MCini** library](https://github.com/eugene274/mcini)
- To compile create a build dir, move to it and run `cmake /path/to/CMakeList.txt` and `make`
- After it executable file *GRATE* will be created 
- Run executable, e.g. `./GRATE` 
- AAMCC can be used from terminal or using bash pipe 

### Abrasion-Ablation mode

- Answer **0** to question that means **no**
  > Do you want to read Abrasion stage from file?
- Enter the projectile (or side A) nucleus name, list is given by the output
- Enter the target (or side B) nucleus name, list is given by the output
- Input the lower limit for the impact parameter (minimum bias is negative)
- Input the upper limit for the impact parameter (will be requested only if lower limit is non-negative)
- Choose the geometry of the collision (**1** for collider, **0** for fixed target)
- Input the sqrt S_nn (collider mode) or kinetic energy of projectile nucleus (fixed target mode) in GeV per nucleon
- Input number of events to be generated
- Choose the excitation energy parametrization (preferable option is **4**)
- Enter the file name to write histograms (.root will be supplied)
- Wait for the run to finish

### As afterburner to precalculated UniGen format files

- You will need a file compatible UniGen format with *UParticle* field *Status* that is set to 0 for spectator nucleons and >0 otherwise
  - You may use any desired version you want, for UrQMD you may use [our fork of UniGen](https://github.com/alexsvetlichnyy/UniGen)
- Answer **1** to question that means **yes**
  > Do you want to read Abrasion stage from file?
- Choose a primary generator (currently only option is UrQMD, so anwser **1**)
- Enter full path to appropriate file
- Choose the excitation energy parametrisation (preferable option is **4**)
- Enter the file name to write histograms (.root will be supplied)
- Wait for the run to finish

### Using with bash pipe

Just pipe the text files with lines, containing all the answers in sequential lines

### Output format 

Output is written in both AAMCC internal format and [MCini](https://github.com/eugene274/mcini) format that is basically UniGen format with initial state.

- In internal format all is stored in `std::vector<double>`, or just `double` leafs. Results of the modelling will contain only spectators along with possible initial stage information.
- In MCini format...
  - in Abrasion-Ablation mode only spectator fragments are stored along with provided by Glauber model initial stage information
  - in afterburner mode participant and produced particles will be stored along with spectator fragments

## Selected papers published using AAMCC

I. A. Pshenichnov, U. A. Dmitrieva and A. O. Svetlichnyi  Secondary Nuclei from Peripheral and Ultraperipheral Collisions of Relativistic Heavy Ions, http://dx.doi.org/10.3103/S1062873820080249

A. O. Svetlichnyi and I. A. Pshenichnov, Formation of Free and Bound Spectator Nucleons in Hadronic Interactions between Relativistic Nuclei, http://dx.doi.org/10.3103/S1062873820080110

R. Nepeivoda, A. Svetlichnyi, N. Kozyrev, I. Pshenichnov, Pre-Equilibrium Clustering in Production of Spectator Fragments in Collisions of Relativistic Nuclei, http://dx.doi.org/10.3390/particles5010004

N. Kozyrev, A. Svetlichnyi, R. Nepeivoda, I. Pshenichnov, Peeling away neutron skin in ultracentral collisions of relativistic nuclei, http://dx.doi.org/10.1140/epja/s10050-022-00832-5

A. Svetlichnyi, S. Savenkov, R.Nepeivoda, I. Pshenichnov, Clustering in Oxygen Nuclei and Spectator Fragments in 16O–16O Collisions at the LHC, http://dx.doi.org/10.3390/physics5020027

E. Vasyagina, S. Savenkov, A. Svetlichnyi, Simulation of the Momentum Distributions of the Spectator Fragments in 124Xe + CsI Collisions at the BM@N with Accounting for Pre-Equilibrium Clusterization, http://dx.doi.org/10.1134/S1547477125701110

## Selected conference papers published using AAMCC

http://dx.doi.org/10.22323/1.398.0310

http://dx.doi.org/10.22323/1.397.0223

http://dx.doi.org/10.1063/5.0063284
