#!/usr/local/bin/python

import os, sys

os.system('rsync -avz -e \'ssh -p 5022\' /Users/andcj/PD_scripts/*.py jakobsen@firewall.esrf.eu:dfxm/.')
os.system('rsync -avz -e \'ssh -p 5022\' /Users/andcj/PD_scripts/lib jakobsen@firewall.esrf.eu:dfxm/.')

#os.system('ssh -p 5022 jakobsen@firewall.esrf.eu \"cd dfxm; python %s\"' % (sys.argv[1]))

os.system('rsync -avz -e \'ssh -p 5022\' jakobsen@firewall.esrf.eu:dfxm/output/ /Users/andcj/PD_scripts/output/. ')
