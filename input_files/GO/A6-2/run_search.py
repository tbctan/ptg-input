from ase import Atoms
from ase.io import Trajectory
from ase.constraints import FixAtoms

from gpaw import GPAW, FermiDirac, Mixer
from gpaw import extra_parameters

from gofee.candidates.candidate_generation import CandidateGenerator, StartGenerator
from gofee.utils import OperationConstraint
from gofee.candidates.basic_mutations import RattleMutation, PermutationMutation
from gofee.gofee_modified import GOFEE

import numpy as np
import math

extra_parameters['blacs'] = True

a = 2.47011790882142
c = 4.0
len_cc = a/math.sqrt(3.0)
padding = 3.0

scaff = Atoms('C28H4', cell = [2*math.sqrt(3.0)*a, 5.5*a+2*padding, c+2*padding], pbc=(1, 0, 0),
              positions = [(len_cc,     padding,       0.5*c+padding),
                           (2*len_cc,   padding,       0.5*c+padding),
                           (0.5*len_cc, padding+0.5*a, 0.5*c+padding), 
                           (2.5*len_cc, padding+0.5*a, 0.5*c+padding),
                           (len_cc,     padding+a,     0.5*c+padding), 
                           (2*len_cc,   padding+a,     0.5*c+padding),
                           (0.5*len_cc, padding+1.5*a, 0.5*c+padding), 
                           (2.5*len_cc, padding+1.5*a, 0.5*c+padding),
                           (len_cc,     padding+2*a,   0.5*c+padding), 
                           (2*len_cc,   padding+2*a,   0.5*c+padding),
                           (0.5*len_cc, padding+2.5*a, 0.5*c+padding),
                           (2.5*len_cc, padding+2.5*a, 0.5*c+padding),
                           (len_cc,     padding+3*a,   0.5*c+padding),
                           (2*len_cc,   padding+3*a,   0.5*c+padding),
                           (4*len_cc,   padding,       0.5*c+padding),
                           (5*len_cc,   padding,       0.5*c+padding),
                           (3.5*len_cc, padding+0.5*a, 0.5*c+padding),
                           (5.5*len_cc, padding+0.5*a, 0.5*c+padding),
                           (4*len_cc,   padding+a,     0.5*c+padding),
                           (5*len_cc,   padding+a,     0.5*c+padding),
                           (3.5*len_cc, padding+1.5*a, 0.5*c+padding),
                           (5.5*len_cc, padding+1.5*a, 0.5*c+padding),
                           (4*len_cc,   padding+2*a,   0.5*c+padding),
                           (5*len_cc,   padding+2*a,   0.5*c+padding),
                           (3.5*len_cc, padding+2.5*a, 0.5*c+padding),
                           (5.5*len_cc, padding+2.5*a, 0.5*c+padding),
                           (4*len_cc,   padding+3*a,   0.5*c+padding),
                           (5*len_cc,   padding+3*a,   0.5*c+padding),
                           (0.90662637, 2.03260543,    0.5*c+padding),
                           (3.37174335, 2.03260543,    0.5*c+padding),
                           (5.18499609, 2.03260543,    0.5*c+padding),
                           (7.65011306, 2.03260543,    0.5*c+padding)])

constraints = FixAtoms(indices=range(0, 32))
scaff.set_constraint(constraints)

stoichiometry = 6*[6] + 2*[1] + [78] + [7] 

cell = scaff.get_cell()

orig = np.array((0., a+padding, padding))
vecs = np.copy(cell)
vecs[1,1] -= 2*padding + a
vecs[2,2] -= 2*padding
box = [orig, vecs]

sg = StartGenerator(scaff, stoichiometry, box)

box_constraint = OperationConstraint(box=box)

n_to_optimize = len(stoichiometry)

candidate_generator = CandidateGenerator([0.2, 0.2, 0.6],
                                        [sg,
                                         PermutationMutation(n_to_optimize, Npermute=2),
                                         RattleMutation(n_to_optimize, Nrattle=3, rattle_range=4)])



calc = GPAW(mode='lcao', basis='dzp', xc='PBE',
            occupations=FermiDirac(0.05), maxiter=999,
            mixer=Mixer(nmaxold=5, beta=0.05, weight=75),
            nbands='110%', kpts=(8, 2, 1), txt = 'eval.txt')

search = GOFEE(structures=None, calc=calc, gpr=None, startgenerator=sg,
               candidate_generator=candidate_generator, max_steps=500,
               population_size=10, position_constraint=box_constraint,
               dualpoint=True, restart='restart',
               kappa='decay',similarity_thr=0.9999)

search.run()

