import numpy as np
#from openmmplumed import PlumedForce
from openmm.app import *
from openmm import *
from openmm.unit import *
import pdbfixer
import US_utils as op            # openmm functions US version
import sys 
import os
import pickle as pckl
import mdtraj as md
import kinaseCVs as kcv
from kinaseCVs import resids_DDR1
from kinaseCVs import resids_Abl1

oss=os.system



if __name__=='__main__':
    '''
    Python file to perform Molecular dynamics using openmm - basic example to use openmm_utils package
    
        (A better script will be updated)
    
    Args:
    
    Input:
     -file            : PDB file to start
     -cpt             : path to cpt file
     -nsteps          : Number simulation steps (default=25000000)
     -Abl1            : Add this flag for Abl1 simulation
     -DDR1            : Add this flag for DDR1 simulation
     -temp            : Temperature of simulation (default=300K)
     -pressure        : Pressure of simulation (default=1bar)
     -plumed          : Plumed file to bias or calculate CVs on the fly
    
    Output:
      fixed_1.pdb     : PDB after ficing the missing residues and terminal atoms
      fixedH_1.pdb    : PDB after adding hydrogen
      solvated_1.pdb  : After solvating adding ions
      minim_1.pdb     : After running a small energy minimization step
      mdout_1.pdb     : final structure file after runnung 'run_NVT()' once
      mdlog_1.txt     : Log file after runnung 'run_NVT()' once
      mdout_1.dcd     : trajectory file after runnung 'run_NVT()' once
    '''
    
    if '-source' in sys.argv:
        source = sys.argv[sys.argv.index('-source')+1]
        flag=1
    else:
        flag=0

    us_bias = sys.argv[sys.argv.index('-usbias')+1]
    us_eps_grid = np.loadtxt(sys.argv[sys.argv.index('-useps')+1])
    us_eps_id = int(sys.argv[sys.argv.index('-usid')+1])
    us_eps = us_eps_grid[us_eps_id]

    if '-nsteps' in sys.argv:
        nsteps = int(sys.argv[sys.argv.index('-nsteps')+1])
    else:
        nsteps = 5000000000        # 100ns simulation                                
    
    if '-temp' in sys.argv:
        temp = int(sys.argv[sys.argv.index('-temp')+1])   #Integer value
    else:
        temp = 300              # 300K temperature
    
    if '-press' in sys.argv:
        pressure = float(sys.argv[sys.argv.index('-press')+1])   #float value
    else:
        pressure = 1.0              # 1.0 bar pressure
    
    #Should add an argument for pressure example as 1bar
     
    
    if '-restart' in sys.argv:
        restart = sys.argv[sys.argv.index('-restart')+1]
    else:
        restart = False

    if '-Abl1' in sys.argv:
        resids = resids_Abl1
    if '-DDR1' in sys.argv:
        resids = resids_DDR1

    if '-diheds' in sys.argv:
        helixk =  float(sys.argv[sys.argv.index('-helixk')+1])
        diheds_file = sys.argv[sys.argv.index('-diheds')+1]
        with open(diheds_file, 'rb') as file:
            diheds = pckl.load(file)
    else:
        diheds = False
        helixk = False

    if flag:
            
            oss("pwd")           
            print(us_eps)
            with open(f"{us_bias}", "rb") as file:
                us_biasp=pckl.load(file)
                
            pdb=PDBFile("%s.pdb"%source)
          
            topology=pdb.topology
            positions=pdb.positions
            
            ff = ForceField('amber99sbildn.xml', 'tip3p.xml')
            simulation=op.get_LangevinM_system(topology,ff,temp,0.002)

            #simulation.loadState("%s.state"%source)
            #velocities=simulation.context.getState(getVelocities=True,enforcePeriodicBox=True).getVelocities()
            top=md.Topology.from_openmm(simulation.topology)
            traj=md.load("%s.pdb"%source)
            
            classk=kcv.kinase_cvs("WT",traj.top,resids,traj)
            atoms=classk.getatoms()
            # Restart Simulations - NPT production run
            
            positions,velocities,simulation=op.run_MD(positions,simulation,nsteps,ens='NPT',run_type='prod',\
                                                   writeXTC=True,usparams=us_biasp,us_eps=us_eps,printdists=atoms,diheds=diheds,helixk=helixk)           #No '_' in run_type
