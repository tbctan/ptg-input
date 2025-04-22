from ase.build import molecule
from ase import Atoms
from ase.io import read
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import QuasiNewton,LBFGS
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo,IdealGasThermo

atoms = read('ptg-oh.xsf',index=-1)
atoms[-2].magmom=0.8 #Assign magnetic moment to O
nPt=atoms.get_chemical_symbols().index("Pt")
atoms[nPt].magmom=0.8


pseudodir = '/home/bctan/QE/pseudoatsume/GBRV_pbe_UPF'
pseudopotentials = {'H': 'h_pbe_v1.4.uspp.F.UPF', 'O':'o_pbe_v1.2.uspp.F.UPF',
                    'C':'c_pbe_v1.2.uspp.F.UPF', 'Pt':'pt_pbe_v1.4.uspp.F.UPF', 
                    'C1':'c_pbe_v1.2.uspp.F.UPF'}
input_data = {
    'control': {'etot_conv_thr':0.0000001,'forc_conv_thr':0.0001, 'disk_io':'low',
        'outdir':'./tmp', 'pseudo_dir':pseudodir, },
    'system': {'ecutwfc':36, 'ecutrho':400,'occupations':'smearing',
        'degauss':0.007,'smearing':'m-p','input_dft':'rpbe'},
    'electrons': {'electron_maxstep':1001,'conv_thr':1.0e-08 , ##34e-10,
        'mixing_mode':'plain','mixing_beta':0.3,
        'diagonalization':'rmm-davidson'}
    }

profile = EspressoProfile(command='mpirun pw.x', pseudo_dir = pseudodir)

calc = Espresso(profile = profile, pseudopotentials=pseudopotentials, input_data=input_data,
                      tprnfor=True,kpts=(8,2,1))

atoms.calc = calc

potentialenergy = atoms.get_potential_energy()
nat=atoms.get_global_number_of_atoms()
vib = Vibrations(atoms,indices=(nat-1,nat-2), nfree=2)
vib.run()
vib.write_mode()
vib.summary(log='out_2.txt')
vib_energies = vib.get_energies()


thermo = HarmonicThermo(vib_energies=vib_energies, ignore_imag_modes=True,
                        potentialenergy=potentialenergy)

G=thermo.get_entropy(temperature=298.15,verbose=True)
G.summary(log='out_3.txt')