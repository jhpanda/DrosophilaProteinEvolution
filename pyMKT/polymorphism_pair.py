
from codon_pair import *
from mut_code import mut_code

class PolymorphismPair(object):
    """ calculate Pn,Ps,Dn,Ds for a gene
        this is the MK.pl method, i.e. each ingroup/outgroup codon pair
    """
    def __init__(self,name,debug=False):
        """ self.name: name for the MKTest, should be helpful for debug
        """
        self.name = name
        self.debug = debug

    def remove_gaps(self,sequences,reference):
        npopulation = len(sequences)
        self.npopulation = npopulation
        new_ref = ''
        new_seq = ['' for s in sequences]
        for i in range(len(reference)):
            if reference[i] not in gaps and sequences[0][i] not in gaps:
                new_ref += reference[i]
                for j in range(npopulation):
                    new_seq[j] += sequences[j][i]
        return new_seq,new_ref

    def divergent(self,codon1,codon0):
        """ determine if a site is polymorphism or fixed
            I dont want to consider only two codons here as in MK.pl
            This might be complicated, but could help if I want to add more 
            details in.
        """
        if all(codon_table[s]==codon_table[codon0] for s in codon1):
            if codon0 not in codon1:
                ## fixed synonymous site
                ## All codons code for the same amino acid
                ## All codons are different with the outgroup codon
                return 'fixed_syn'
            else:
                ## polymorphism synonymous site
                ## All codons code for the same amino acid
                ## At least one of the codons are the same as outgroup codon
                if all(s!=codon0 for s in codon1):
                    return 'polymorph_syn_codon_not_in_outgroup'
                else:
                    return 'polymorph_syn'
                
        elif all(codon_table[s]!=codon_table[codon0] for s in codon1):
            c10 = codon1[0]
            if all(codon_table[s]==codon_table[c10] for s in codon1):
                ## fixed non-synonymous site
                ## All ingroup codons code for the same amino acid
                ## The amino acid is different with outgroup one
                return 'fixed_nonsyn'
            else:
                ## polymorphism non-synonymous site
                ## Ingroup codons code the different amino acid, not fixed
                ## There might exist fixation and polymorph in the population
                return 'polymorph_nonsyn_codon_not_in_outgroup'
        else:
            ## polymorphism site
            ## some codons code for the same amino acid as the outgroup codon
            if all(s!=codon0 for s in codon1):
                ## The ingroup codons are all different to the outgroup codon
                ## There might exist fixation and polymorph in the population
                return 'polymorph_codon_not_in_outgroup'
            else:
                ## Some of the ingroup codons are the same as the outgroup codon
                ## This situation suggests polymorph in the population
                return 'polymorph_nonsyn'
    
    def bool_stop_codon(self,x,codon1,codon0,stop_codon_type='any'):
        """ This bool function returns:
            if stop_codon_type=='all':
                True if all ingroup/outgroup codons are stop codon else False
                print(warning info)
            if stop_codon_type=='any':
                True if any ingroup/outgroup codons are stop codon else False
        """
        if stop_codon_type=='all':
            if all(s in stop_codon for s in codon1) and codon0 in stop_codon:
                return True
            else:
                if any(s in stop_codon for s in codon1+[codon0]):
                    print("Warn! stop codon detected for %s"%self.name, codon1,x)
                #print(codon1,codon0)
            return False
        
        if stop_codon_type=='any':
            if any(s in stop_codon for s in codon1+[codon0]):
                return True
            else:
                return False

    def site_changes(self,path,codon1_nr,codon0):
        """ Calculate pairwise Pn,Ps,Dn,Ds with 
        """
        divergent = self.divergent(codon1_nr,codon0)
        ps = pn = ds = dn = 0

        if divergent == 'polymorph_syn':
            ## MK.pl dont consider polymorph sites with >2 unique codons  
            n,s = mut_code[path[0]+path[1]]
            pn += n
            ps += s
            ## here, I will add more pn,ps to this site 
            #for i in range(len(path)-1):
            #    n,s = mut_code[path[i]+path[i+1]]
            #    pn += n
            #    ps += s

        if divergent == 'polymorph_nonsyn':
            ## MK.pl dont consider polymorph sites with >2 unique codons  
            n,s = mut_code[path[0]+path[1]]
            pn += n
            ps += s
            ## here, I will add more pn,ps to this site 
            #for i in range(len(path)-1):
            #    n,s = mut_code[path[i]+path[i+1]]
            #    pn += n
            #    ps += s

        if divergent == 'fixed_syn':
            n,s = mut_code[codon0+path[1]]
            dn += n
            ds += s
            for i in range(1,len(path)-1):
                #n,s = mut_code[(path[i],path[i+1])]
                n,s = mut_code[path[i]+path[i+1]]
                pn += n
                ps += s
                
        if divergent == 'fixed_nonsyn':
            #n,s = mut_code[(path[0],path[1])]
            n,s = mut_code[path[0]+path[1]]
            dn += n
            ds += s
            for i in range(1,len(path)-1):
                #n,s = mut_code[(path[i],path[i+1])]
                n,s = mut_code[path[i]+path[i+1]]
                pn += n
                ps += s

        if divergent == 'polymorph_nonsyn_codon_not_in_outgroup':
            #n,s = mut_code[(codon0,path[1])]
            n,s = mut_code[codon0+path[1]]
            dn += n
            ds += s
            for i in range(1,len(path)-1):
                #n,s = mut_code[(path[i],path[i+1])]
                n,s = mut_code[path[i]+path[i+1]]
                pn += n
                ps += s

        if divergent == 'polymorph_codon_not_in_outgroup' or divergent == 'polymorph_syn_codon_not_in_outgroup':
            ##n,s = mut_code[(codon0,path[1])]
            n,s = mut_code[codon0+path[1]]
            dn += n
            ds += s
            for i in range(1,len(path)-1):
                #n,s = mut_code[(path[i],path[i+1])]
                n,s = mut_code[path[i]+path[i+1]]
                pn += n
                ps += s

        return ps,pn,ds,dn

    def window_changes(self,sequences,reference_outgroup,stop_codon_type='any'):
        new_seq,new_ref = self.remove_gaps(sequences,reference_outgroup)
        nsite = len(new_ref)//3 # floor division
        Ds = Dn = 0
        Ps = Pn = 0
        x = 0
        for n in range(nsite):
            codon1 = [seq[n*3:n*3+3] for seq in new_seq]
            codon1_nr = list(set(codon1)) ## unique codons (remove redundant)
            codon0 = new_ref[n*3:n*3+3]
            ps = pn = ds = dn = 0
            ## MK.pl do not consider more than 2 codon states
            codon1_nr_ = [s for s in codon1_nr]
            if len(list(set(codon1_nr+[codon0])))>3:
                codon1_nr_ = [codon0]
            ##
            if self.bool_stop_codon(x,codon1_nr_,codon0,stop_codon_type):
                x += 1
            elif len(codon1_nr_) == 1 and codon1_nr_[0]==codon0:
                if self.debug:
                    divergent = 'None'
                    print(x,codon1_nr_,codon0,pn,dn,ps,ds,divergent)
                x += 1
            else:
                min_path = pair.shortest_path(codon1_nr_,codon0)
                if len(min_path)==1:
                    if self.debug:
                        min_path = []
                        divergent = self.divergent(codon1_nr_,codon0)
                        print(x,codon1_nr,codon0,min_path,pn,dn,ps,ds,divergent)
                else:
                    ps,pn,ds,dn = self.site_changes(min_path,codon1_nr_,codon0)
                    if self.debug:
                        divergent = self.divergent(codon1_nr_,codon0)
                        print(x,codon1_nr,codon0,min_path,pn,dn,ps,ds,divergent)
                Ps += ps
                Pn += pn
                Ds += ds
                Dn += dn
                x += 1
        #print("Ps: %5d Pn: %5d Ds: %5d Dn: %5d"%(Ps,Pn,Ds,Dn))
        return Ps,Pn,Ds,Dn

if __name__ == '__main__':
    from fasta2seq import *
    from fisher import pvalue
    #alignments = fasta2seq('../fasta-aligned/FBgn0015790_FBgn0190898_FBgn0242765.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0000352_FBgn0182039_FBgn0241241.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0038959_FBgn0192402_FBgn0228168.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0011244_FBgn0194432_FBgn0234487.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0000109_FBgn0185140_FBgn0238192.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0283741_FBgn0188176_FBgn0234208.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0000216_FBgn0186193_FBgn0239268.fasta')
    alignments = fasta2seq('../fasta-aligned/FBgn0085328_FBgn0192936_FBgn0277279.fasta')
    alignments = fasta2seq('../fasta-aligned/FBgn0031574_FBgn0194649_FBgn0235692.fasta')
    #alignments = fasta2seq('test.fasta')
    #alignments = fasta2seq('../fasta/mafft.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0036514_FBgn0186237_FBgn0237163.fasta')
    
    #dyak = 'outgroup'
    #namelist1 = ['dmel_%d'%s for s in range(1,5)]
    namelist1 = ['dmel_%d'%s for s in range(1,411)]
    namelist2 = ['dsim_%d'%s for s in range(1,391)]
    refname1  = 'dsim_0'
    refname2  = 'dmel_0'
    
    #namelist1 = ['sp1_%d'%s for s in range(0,2)]
    #namelist2 = ['sp2_%d'%s for s in range(0,2)]
    #refname1  = 'sp2_0'
    #refname2  = 'sp1_0'

    #sequences = [alignments[s] for s in namelist2]
    #reference = alignments[refname2] 

    umktest = PolymorphismPair('test',debug=True)
    #alignments = fasta2seq('test.fasta')
    #sequences = [alignments[s] for s in 'dmel_1 dmel_2 dmel_3 dmel_4'.split()]
    #reference = alignments['dsim_0']
    #umktest._test(sequences,reference)
    sequences = [alignments[s] for s in namelist1]
    reference = alignments[refname1] 
    #umktest._test(sequences,reference)
    #for i in range(100):
    contigency = umktest.window_changes(sequences,reference)
    #contigency = umktest.window_changes(sequences,reference)
    #contigency = umktest.window_changes(sequences,reference)
    #contigency = umktest.window_changes(sequences,reference)
    #contigency = umktest.window_changes(sequences,reference)
    Ps,Pn,Ds,Dn = contigency
    if Ps and Dn:
        alpha = 1-Ds*Pn/(Ps*Dn)
    else:
        alpha = -99.999
    #from scipy.stats import fisher_exact
    ##_,pval = fisher_exact([[Ps,Ds],[Pn,Dn]])
    pval = pvalue(Ps,Ds,Pn,Dn).two_tail
    print(Pn,Dn,Ps,Ds,alpha,pval)
    #main()
    #test_all_genes()
    #divergent = umktest.divergent(['AGT','TGT'],'TGC')
    #path = graph.shortest_path(['AGT'],'TGC')
    #print(path)
    #codon1_nr = ['AGT','AGC']
    #codon0 = 'TGC'
    #ps,pn,ds,dn = umktest.site_changes(path,codon1_nr,codon0)
    #print(ps,pn,ds,dn)
    #print(divergent)

