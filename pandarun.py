#!/usr/local/bin/python

import os, sys
#
# os.system('rsync -avz /Users/andcj/PD_scripts/*.py jakobsen@panda.hpc.dtu.dk:dfxm/.')
# os.system('rsync -avz /Users/andcj/PD_scripts/lib jakobsen@panda.hpc.dtu.dk:dfxm/.')
# os.system('rsync -avz /Users/andcj/PD_scripts/resolution_paper/*.py jakobsen@panda.hpc.dtu.dk:dfxm/.')
#
# os.system('ssh jakobsen@panda.hpc.dtu.dk \"cd dfxm; mpiexec -n %s python %s\"' % (sys.argv[1], sys.argv[2]))
#
# # os.system('ssh jakobsen@panda.hpc.dtu.dk \"cd dfxm; sh make_rec_data.sh\"')
#
# os.system('rsync -avz jakobsen@panda.hpc.dtu.dk:dfxm/output/ /Users/andcj/PD_scripts/output/. ')

os.system('rsync -avz /Users/andcj/PD_scripts/*.py andcj@panda2.fysik.dtu.dk:dfxm/.')
os.system('rsync -avz /Users/andcj/PD_scripts/lib andcj@panda2.fysik.dtu.dk:dfxm/.')
os.system('rsync -avz /Users/andcj/PD_scripts/resolution_paper/*.py andcj@panda2.fysik.dtu.dk:dfxm/.')
os.system('rsync -avz /Users/andcj/PD_scripts/reconstruct andcj@panda2.fysik.dtu.dk:dfxm/.')

os.system('ssh andcj@panda2.fysik.dtu.dk \"cd dfxm; mpiexec -n %s python %s\"' % (sys.argv[1], sys.argv[2]))

# os.system('ssh andcj@panda2.fysik.dtu.dk \"cd dfxm; sh make_rec_data.sh\"')

os.system('rsync -avz andcj@panda2.fysik.dtu.dk:dfxm/output/ /Users/andcj/PD_scripts/output/. ')
