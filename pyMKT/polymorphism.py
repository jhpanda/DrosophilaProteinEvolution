
from codon_graph import *

class Polymorphism(object):
    """ calculate Pn,Ps,Dn,Ds for a gene
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
                return 'polymorph'
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
                return 'polymorph'
        else:
            ## polymorphism site
            ## these codons do not code for the same amino acid
            return 'polymorph'
    
    def bool_stop_codon(self,x,codon1,codon0):
        """ If all codons are stop codon, then stop.
            If some codons are stop codon, then raise a warning
        """
        if all(s in stop_codon for s in codon1) and codon0 in stop_codon:
            return True
        else:
            if any(s in stop_codon for s in codon1) and not codon0 in stop_codon:
                print("Warn! stop codon detected for %s"%self.name, codon1,x)
                #print(codon1,codon0)
            return False

    def site_changes(self,path,codon1_nr,codon0):
        """ Calculate Pn,Ps,Dn,Ds along the minimum path
        """
        divergent = self.divergent(codon1_nr,codon0)
        #codons_ = []
        #for codon in path:
        #    if codon!=codon0:
        #        codons_ += [codon]
        #divergent = self.divergent(codons_ ,codon0)
        ps = pn = ds = dn = 0
        if path[-1] == codon0:
            ## if outgroup codon is in the end of the path
            #for i in range(len(path)-2):
            #    if codon_table[path[i]] == codon_table[path[i+1]]:
            #        ps += 1
            #    else:
            #        pn += 1
            #if divergent == 'polymorph':
            #    if codon_table[path[-2]] == codon_table[path[-1]]:
            #        ps += 1
            #    else:
            #        pn += 1
            #if divergent == 'fixed_syn':
            #    ds += 1
            #if divergent == 'fixed_nonsyn':
            #    dn += 1
            if divergent == 'polymorph':
                for i in range(len(path)-1):
                    if codon_table[path[i]] == codon_table[path[i+1]]:
                        ps += 1
                    else:
                        pn += 1
            if divergent == 'fixed_syn':
                for i in range(len(path)-1):
                    if path[i] in codon1_nr and path[i+1] in codon1_nr:
                        ps += 1
                    else:
                        ds += 1
            if divergent == 'fixed_nonsyn':
                for i in range(len(path)-1):
                    if codon_table[path[i]] == codon_table[path[i+1]]:
                        if path[i] in codon1_nr and path[i+1] in codon1_nr:
                            ps += 1
                        else:
                            ds += 1
                    else:
                        dn += 1
        else:
            ## if outgroup codon is in middle of the path
            index = path.index(codon0)

            ## Pn,Ps,Dn,Ds before the outgroup codon index
            for i in range(index-1):
                if codon_table[path[i]] == codon_table[path[i+1]]:
                    ps += 1
                else:
                    pn += 1
            if divergent == 'polymorph':
                if codon_table[path[index-1]] == codon_table[path[index]]:
                    ps += 1
                else:
                    pn += 1
            if divergent == 'fixed_syn':
                ds += 1
            if divergent == 'fixed_nonsyn':
                dn += 1

            ## Pn,Ps,Dn,Ds after the outgroup codon index
            if divergent == 'polymorph':
                if codon_table[path[index]] == codon_table[path[index+1]]:
                    ps += 1
                else:
                    pn += 1
            if divergent == 'fixed_syn':
                ds += 1
            if divergent == 'fixed_nonsyn':
                dn += 1
            for i in range(index+1,len(path)-1):
                if codon_table[path[i]] == codon_table[path[i+1]]:
                    ps += 1
                else:
                    pn += 1

        return ps,pn,ds,dn

    def window_changes(self,sequences,reference_outgroup):
        new_seq,new_ref = self.remove_gaps(sequences,reference_outgroup)
        nsite = len(new_ref)//3 # floor division
        Ds = Dn = 0
        Ps = Pn = 0
        x = 0
        for n in range(nsite):
            codon1 = [seq[n*3:n*3+3] for seq in new_seq]
            codon1_nr = list(set(codon1)) ## unique codons (remove redundant)
            codon0 = new_ref[n*3:n*3+3]
            if self.bool_stop_codon(x,codon1_nr,codon0):
                break
            else:
                ps = pn = ds = dn = 0
                if len(codon1_nr)==1:
                    if codon1_nr[0]!=codon0:
                        min_path = graph.shortest_path(codon1_nr,codon0)
                        ps,pn,ds,dn = self.site_changes(min_path,codon1_nr,codon0)
                        ### uncomment for debug ###
                        if self.debug:
                            divergent = self.divergent(codon1_nr,codon0)
                            print(x,codon1_nr,codon0,min_path,[codon_table[s] for s in min_path],ps,pn,ds,dn,divergent)
                    else:
                        if self.debug:
                            min_path = []
                            divergent = self.divergent(codon1_nr,codon0)
                            print(x,codon1_nr,codon0,min_path,[codon_table[s] for s in min_path],ps,pn,ds,dn,divergent)
            
                if len(codon1_nr)>1:
                    min_path = graph.shortest_path(codon1_nr,codon0)
                    ps,pn,ds,dn = self.site_changes(min_path,codon1_nr,codon0)
                    ### uncomment for debug ###
                    if self.debug:
                        divergent = self.divergent(codon1_nr,codon0)
                        print(x,codon1_nr,codon0,min_path,[codon_table[s] for s in min_path],ps,pn,ds,dn,divergent)
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
    alignments = fasta2seq('../fasta-aligned/FBgn0011244_FBgn0194432_FBgn0234487.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0000109_FBgn0185140_FBgn0238192.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0031653_FBgn0184718_FBgn0068255.fasta')
    #alignments = fasta2seq('../fasta-aligned/FBgn0267823_FBgn0196410_FBgn0229423.fasta')
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

    #umktest = Polymorphism('test',debug=True)
    umktest = Polymorphism('test',debug=False)
    #alignments = fasta2seq('test.fasta')
    #sequences = [alignments[s] for s in 'dmel_1 dmel_2 dmel_3 dmel_4'.split()]
    #reference = alignments['dsim_0']
    #umktest._test(sequences,reference)
    sequences = [alignments[s] for s in namelist1]
    reference = alignments[refname1] 
    #umktest._test(sequences,reference)
    for i in range(100):
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
    pval = pvalue(Pn,Dn,Ps,Ds).two_tail
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

