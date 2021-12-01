
from calcMK import *

if __name__ == '__main__':

### uncomment for unpolar MK test ###
#    #print(calcMK.__doc__)
#    #fasta_list = 'test_list.txt'
    fasta_list = 'alignment_list.txt'
    window     = 90
    polar      = 'unpolar'
    ingroup_name = 'dmel'
    outgroup_name = 'dsim'
    method     = 'MK.pl'
    ncpu       = 64
    cutoff     = 0.05
    screenout  = 'yes'
    #for cutoff in [0.05,0.10,0.15,0.20,0.30,0.50]:
    calcMK(fasta_list,cutoff,method,window,polar,ncpu,screenout,ingroup_name,outgroup_name)

### uncomment for polar MK test ###
#    fasta_list = 'alignment_list.txt'
#    #fasta_list = 'test_list.txt'
#    window     = 0
#    polar      = 'polar'
#    ingroup_name = 'dsim'
#    outgroup_name = 'outgroup'
#    fixation_name = 'dmel'
#    method     = 'MK.pl'
#    ncpu       = 64
#    cutoff     = 0.05
#    screenout  = 'yes'
#    #for cutoff in [0.05,0.10,0.15,0.20,0.30,0.50]:
#    calcMK(fasta_list,cutoff,method,window,polar,ncpu,screenout,ingroup_name,outgroup_name,fixation_name)
