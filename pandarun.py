#!/usr/local/bin/python

import os, sys

os.system('rsync -avz /Users/andcj/PD_scripts/*.py jakobsen@panda.hpc.dtu.dk:dfxm/.')
os.system('rsync -avz /Users/andcj/PD_scripts/lib jakobsen@panda.hpc.dtu.dk:dfxm/.')
os.system('rsync -avz /Users/andcj/PD_scripts/resolution_paper jakobsen@panda.hpc.dtu.dk:dfxm/.')

os.system('ssh jakobsen@panda.hpc.dtu.dk \"cd dfxm; mpiexec -n %s python %s\"' % (sys.argv[1], sys.argv[2]))

# os.system('ssh jakobsen@panda.hpc.dtu.dk \"cd dfxm; sh make_rec_data.sh\"')

os.system('rsync -avz jakobsen@panda.hpc.dtu.dk:dfxm/output/ /Users/andcj/PD_scripts/output/. ')
