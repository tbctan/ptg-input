from ase import Atoms
from ase.io import Trajectory
from ase.constraints import FixAtoms

from gpaw import GPAW, FermiDirac, Mixer
from gpaw import extra_parameters

from gofee.candidates.candidate_generation import CandidateGenerator, StartGenerator
from gofee.utils import OperationConstraint
from gofee.candidates.basic_mutations import RattleMutation, PermutationMutation
from gofee_modified.gofee import GOFEE

import numpy as np
import math

extra_parameters['blacs'] = True

a = 2.47011790882142
c = 4.0
len_cc = a/math.sqrt(3.0)
len_ch = 1.09
padding = 3.0

scaff = Atoms('C32H4', cell = [4*a, 5*a+2*padding, c+2*padding], pbc=(1, 0, 0),
              positions = [(len_cc*0.5,          0.5*padding+len_cc,        0.5*c+padding),
                           (len_cc*0.5,          0.5*padding+2*len_cc,      0.5*c+padding),
                           (len_cc*0.5+0.5*a,    0.5*padding+0.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+0.5*a,    0.5*padding+2.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+a,        0.5*padding+len_cc,        0.5*c+padding),
                           (len_cc*0.5+a,        0.5*padding+2*len_cc,      0.5*c+padding),
                           (len_cc*0.5+1.5*a,    0.5*padding+0.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+1.5*a,    0.5*padding+2.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+2*a,      0.5*padding+len_cc,        0.5*c+padding),
                           (len_cc*0.5+2*a,      0.5*padding+2*len_cc,      0.5*c+padding),
                           (len_cc*0.5+2.5*a,    0.5*padding+0.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+2.5*a,    0.5*padding+2.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+3*a,      0.5*padding+len_cc,        0.5*c+padding),
                           (len_cc*0.5+3*a,      0.5*padding+2*len_cc,      0.5*c+padding),
                           (len_cc*0.5+3.5*a,    0.5*padding+0.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+3.5*a,    0.5*padding+2.5*len_cc,    0.5*c+padding),
                       #    (len_cc*0.5+4*a,      0.5*padding+len_cc,        0.5*c+padding),
                       #    (len_cc*0.5+4*a,      0.5*padding+2*len_cc,      0.5*c+padding),
                           (len_cc*0.5,          0.5*padding+4*len_cc,      0.5*c+padding),
                           (len_cc*0.5,          0.5*padding+5*len_cc,      0.5*c+padding),
                           (len_cc*0.5+0.5*a,    0.5*padding+3.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+0.5*a,    0.5*padding+5.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+a,        0.5*padding+4*len_cc,      0.5*c+padding),
                           (len_cc*0.5+a,        0.5*padding+5*len_cc,      0.5*c+padding),
                           (len_cc*0.5+1.5*a,    0.5*padding+3.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+1.5*a,    0.5*padding+5.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+2*a,      0.5*padding+4*len_cc,      0.5*c+padding),
                           (len_cc*0.5+2*a,      0.5*padding+5*len_cc,      0.5*c+padding),
                           (len_cc*0.5+2.5*a,    0.5*padding+3.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+2.5*a,    0.5*padding+5.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+3*a,      0.5*padding+4*len_cc,      0.5*c+padding),
                           (len_cc*0.5+3*a,      0.5*padding+5*len_cc,      0.5*c+padding),
                           (len_cc*0.5+3.5*a,    0.5*padding+3.5*len_cc,    0.5*c+padding),
                           (len_cc*0.5+3.5*a,    0.5*padding+5.5*len_cc,    0.5*c+padding),
                      #     (len_cc*0.5+4*a,      0.5*padding+4*len_cc,      0.5*c+padding),
                      #     (len_cc*0.5+4*a,      0.5*padding+5*len_cc,      0.5*c+padding),
              #
                      #     (len_cc*0.5,          0.5*padding+7*len_cc,      0.5*c+padding),
                      #     (len_cc*0.5+0.5*a,    0.5*padding+6.5*len_cc,    0.5*c+padding),
                      #     (len_cc*0.5+a,        0.5*padding+7*len_cc,      0.5*c+padding),
                      #     (len_cc*0.5+1.5*a,    0.5*padding+6.5*len_cc,    0.5*c+padding),
                      #     (len_cc*0.5+2*a,      0.5*padding+7*len_cc,      0.5*c+padding),
                      #     (len_cc*0.5+2.5*a,    0.5*padding+6.5*len_cc,    0.5*c+padding),
                      #     (len_cc*0.5+3*a,      0.5*padding+7*len_cc,      0.5*c+padding),
                      #     (len_cc*0.5+3.5*a,    0.5*padding+6.5*len_cc,    0.5*c+padding),
                      #     (len_cc*0.5+4*a,      0.5*padding+7*len_cc,      0.5*c+padding),
                           (1.94775322,  1.11137831, 0.5*c+padding),
                           (4.41787113,  1.11137831, 0.5*c+padding),
                           (6.88798904,  1.11137831, 0.5*c+padding),
                           (9.35810695,  1.11137831, 0.5*c+padding)])

constraints = FixAtoms(indices=range(0, 36))
scaff.set_constraint(constraints)
NCAR = 6
NHYD = 2
stoichiometry = NCAR*[6] + NHYD*[1] + [78]

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

