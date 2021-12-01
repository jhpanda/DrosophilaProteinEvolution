"""
    To generate 1000 random sequences for PPISites or NonSites, each sequence contain 1000 randomly selected PPISites or NonSites
    12/05/2019, Junhui Peng
"""

import random
import numpy as np

class Site(object):
    def __init__(self,fbid,site):
        self.fbid  = fbid
        self.index = site

class PPISites():
    def __init__(self,fsites,fgroups='../groups_gene.dat',seqnum=1000,seqlen=1000):
        self.seqnum = seqnum
        self.seqlen = seqlen
        self.fsites = fsites
        nsiteclass  = len(fsites)
        self.nsiteclass = nsiteclass

        with open(fgroups,'r') as p:
            lines = p.readlines()
            m = 0
            ngroup = 0
            genegroups = {}
            for line in lines:
                genes = line.split()
                for gene in genes:
                    genegroups[gene] = ngroup
                m += 1
                #if m==8:
                ngroup += 1
            #ngroup += 1
            self.ngroup = ngroup
            self.genegroups = genegroups

        self.sites = [[[] for s in range(ngroup)] for s in range(nsiteclass+1)]
        for i in range(self.nsiteclass):
            fsite = fsites[i]
            with open(fsite,'r') as f:
                lines = f.readlines()
                for line in lines[1:]:
                    linesplit = line.split()
                    fbid      = linesplit[0]
                    seqlen    = int(linesplit[1]) 
                    numsites  = int(linesplit[2])
                    if numsites == 0:
                        sites = []
                    else:
                        sites = linesplit[3]
                        sites_ = [int(s) for s in sites.split(',')]
                    try:
                        igroup = self.genegroups[fbid]
                        for j in sites_:
                            self.sites[i][igroup] += [Site(fbid,j)]
                        ## genome wide, random expectation ##
                        for j in [s for s in range(seqlen)]:
                            self.sites[-1][igroup] += [Site(fbid,j)]
                    except KeyError:
                        print('%s do not have high quaility alignments'%fbid)
            print(len(self.sites[i]))

    def random_sequences(self):
        seqnum = self.seqnum
        seqlen = self.seqlen
        aa_each_group = int(seqlen/self.ngroup)
        aa_groups     = [aa_each_group for s in range(self.ngroup)]
        aa_groups[-1] = seqlen-aa_each_group*(self.ngroup-1)
        for i in range(self.nsiteclass+1):
            if i<self.nsiteclass:
                fsite = self.fsites[i]
                fout = open('selected_%s'%fsite,'w')
            else:
                ## genome wide, random expectation ##
                fout = open('selected_random_expectation.dat','w')
            for n in range(seqnum):
                fout.write('sequence %d\n'%n)
                random_sites = []
                for m in range(self.ngroup):
                    mlen = aa_groups[m]
                    random_sites += random.choices(self.sites[i][m],k=mlen)
                for site in random_sites:
                    fbid  = site.fbid
                    index = site.index
                    fout.write('%s %d\n'%(fbid,index))
            fout.close()

if __name__ == '__main__':
    fsites ="""
                nonsites.dat
                ppisites.dat
            """.split()

    ppisite = PPISites(fsites,seqnum=100,seqlen=10000)
    ppisite.random_sequences()
