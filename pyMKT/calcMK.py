
import sys,argparse
from pathlib import Path
from scipy.stats import fisher_exact
from polymorphism import *
from polymorphism_pair import *
from fisher import pvalue ## Pypi fast fisher implementation
from fasta2seq import *
from multiprocessing import Pool

def rm_lowfreq(sequences,reference_ingroup,cutoff):
    npopulation  = len(sequences)
    ref_ingroup  = reference_ingroup
    if cutoff > 0.05:
        sequences_   = [list(sequences[s]) for s in range(npopulation)]
        for i in range(len(sequences[0])):
            nucleotides = [sequences_[s][i] for s in range(npopulation)]
            ### this is to remove SNPs with low frequencies
            snp_num = [0 if s==ref_ingroup[i] else 1 for s in nucleotides]
            if sum(snp_num)/npopulation <= cutoff:
                for s in range(npopulation):
                    sequences_[s][i] = reference_ingroup[i]
            ### following removes minor allele SNPs appear less than 3 times 
            ### don't know if this needs to be considered
            unique_nuc = list(set(nucleotides))
            if len(unique_nuc)>1:
                num = [0 for s in unique_nuc]
                for m in range(len(num)):
                    if m in [0,1]:
                        for n in nucleotides:
                            if unique_nuc[m]==n:
                                num[m]+=1
                    else:
                        pass
                flag_nuc = []
                for m in range(len(num)):
                    #if num[m]<0.05*npopulation:
                    if num[m]<3:
                        flag_nuc += [unique_nuc[m]]
                #print(i,flag_nuc)
                for s in range(npopulation):
                    if sequences_[s][i] in flag_nuc:
                        sequences_[s][i] = reference_ingroup[i]
        ## final sequences
        newsequences = []
        for i in range(npopulation):
            newsequences += [''.join(sequences_[i])]
        return newsequences
    else:
        return sequences

def fixation(sequences,reference_outgroup,sequences_fixation=None):
    """ See Begun et al for details about fixed differences
    """
    npopulation  = len(sequences)
    if sequences_fixation:
        nfixation = len(sequences_fixation)
        sequences_   = [list(sequences[s]) for s in range(npopulation)]
        for i in range(0,len(sequences[0]),3):
            codon0 = reference_outgroup[i:i+3] 
            codons = [sequences_fixation[s][i:i+3] for s in range(nfixation)]
            if any(codon!=codon0 for codon in codons):
                for s in range(npopulation):
                    sequences_[s][i]   = '-'
                    sequences_[s][i+1] = '-'
                    sequences_[s][i+2] = '-'
        newsequences = []
        for i in range(npopulation):
            newsequences += [''.join(sequences_[i])]
        return newsequences
    else:
        return sequences

def remove_gaps_in_ref(sequences,reference_ingroup,reference_outgroup):
    npopulation = len(sequences)
    ref_ingrp  = ''
    ref_outgrp = ''
    new_seq    = ['' for s in sequences]
    for i in range(len(reference_ingroup)):
        if reference_ingroup[i] not in gaps:
            ref_ingrp  += reference_ingroup[i]
            ref_outgrp += reference_outgroup[i]
            for j in range(npopulation):
                new_seq[j] += sequences[j][i]
    return new_seq,ref_ingrp,ref_outgrp
   
def window_MKT(sequences,reference,name,method,window_size=0):
    if method == 'shortest_path':
        polymorph = Polymorphism(name)
    elif method == 'MK.pl':
        polymorph = PolymorphismPair(name)
    else:
        print("method is not 'shortest_path' or 'MK.pl'")
        sys.exit(0)
    Ps,Pn,Ds,Dn = polymorph.window_changes(sequences,reference)
    if Ps and Dn:
        alpha = 1-Pn*Ds/(Ps*Dn)
        pval = pvalue(Ps,Ds,Pn,Dn).two_tail ## Pypi fast fisher implementation
    else:
        alpha = -99.99
        pval  = 1.00
    if window_size>0:
        if pval<=0.5 and Ps and Dn:
            return [name,Ps,Pn,Ds,Dn,alpha,pval]
        else:
            return []
    else:
        return [name,Ps,Pn,Ds,Dn,alpha,pval]

def sliding_window_MKT(sequences,reference_ingroup,reference_outgroup,name,method,window_size):
    """ if window_size>0, then a sliding window of MKTest will be performed.
    """
    
    ### to better track current window, we only removed gaps in ingroup ref
    ### other gaps will be removed in polymorphism.py or polymorphism_pair.py

    sequences_ ,ref_ingrp,ref_outgrp = remove_gaps_in_ref(sequences,reference_ingroup,reference_outgroup)

    sliding_step = 9
    sliding_results = []
    seqlen = len(sequences_[0])
    if seqlen<window_size+sliding_step:
        name_  = "%s_%s_%d"%(name,ref_ingrp,0)
        result = window_MKT(sequences_,ref_outgrp,name_,method,window_size)
        sliding_results = [[name_,Ps,Pn,Ds,Dn,alpha,pval]]
        
    else:
        max_ = seqlen-window_size
        for i in range(0,max_,sliding_step):
            #start = i*window_size
            start = i
            ref_ingrp_i = ref_ingrp[start:start+window_size]
            seq_i = [s[start:start+window_size] for s in sequences_]
            ref_i = ref_outgrp[start:start+window_size]
            name_  = "%s_%s_%d"%(name,ref_ingrp_i,i)
            result = window_MKT(seq_i,ref_i,name_,method,window_size)
            if result:
                name_i,Ps,Pn,Ds,Dn,alpha,pval = result
                sliding_results += [[name_i,Ps,Pn,Ds,Dn,alpha,pval]]
    return sliding_results

def test_one_gene(fasta,cutoff,polar,method,window_size,screenout,ingroup_name,outgroup_name,fixation_name=None):
    """ window_size=0 for whole gene; >0 for sliding window MKtest
        polar: "polar"/"unpolar_dsim"/"unpolar_dmel"
    """
    name = Path(fasta).resolve().stem ## basename of fasta file
    #name = Path(fasta).stem() ## basename of fasta file

    alignments = fasta2seq(fasta)
    keys = list(alignments.keys())

    if keys:
        # if ingroup_name == 'dmel' and outgroup_name == 'dsim':
        # outgroup_name = 'dsim' or outgroup_name = ['dyak','dsim']
        # if outgroup_name is string, it's unpolarized, otherwise polarized
        if polar == 'unpolar':
            names_ingrp = []
            for key in keys:
                if ingroup_name in key and '%s_0'%ingroup_name not in key:
                    names_ingrp += [key]
            refnm_ingrp  = '%s_0'%ingroup_name
            refnm_outgrp = '%s_0'%outgroup_name
            sequences  = [alignments[s] for s in names_ingrp]
            ref_ingrp  = alignments[refnm_ingrp]
            ref_outgrp = alignments[refnm_outgrp]
            sequences = rm_lowfreq(sequences,ref_ingrp,cutoff)
        if polar == 'polar':
            if fixation_name == None:
                print('You are running polarized MK test, you need to provide an outgroup name for fixations (fixation_name), see Begun et al, 2007 for detail')
                sys.exit(0)
            names_ingrp = []
            names_fixation = []
            for key in keys:
                if ingroup_name in key and '%s_0'%ingroup_name not in key:
                    names_ingrp += [key]
                if fixation_name in key and '%s_0'%fixation_name not in key:
                    names_fixation += [key]
            sequences  = [alignments[s] for s in names_ingrp]
            seqs_fixation = [alignments[s] for s in names_fixation]
            refnm_ingrp  = '%s_0'%ingroup_name
            ref_ingrp  = alignments[refnm_ingrp]
            ref_outgrp  = alignments['outgroup']
            sequences = rm_lowfreq(sequences,ref_ingrp,cutoff)
            sequences = fixation(sequences,ref_outgrp,seqs_fixation)

        if window_size > 0:
            results = sliding_window_MKT(sequences,ref_ingrp,
                                         ref_outgrp,
                                         name,method,window_size)
        else:
            results = [window_MKT(sequences,ref_outgrp,name,method)]
            
        #print(len(results))
        for result in results:
            name_,ps,pn,ds,dn,a,pval = result
            if screenout == 'yes':
                sys.stdout.write("%s %4d %4d %4d %4d %8.3f %8.2e\n"%(name_,
                                                            pn,dn,ps,ds,a,pval))
        return results
    else:
        #print("Warn! empty alignments for %s"%name)
        return None

def test_all_genes(fasta_list,cutoff,foutput,method='shortest_path',window_size=0,polar='unpolar_dmel',ncpu=1,screenout='yes',ingroup_name='dmel',outgroup_name='dsim',fixation_name=None):
    help_message = """ 
Perform MKTest
Inputs:
 fasta_list  str           fasta files for alignments of genes
 cutoff      float         cutoff to remove low SNPs
 window_size int           perform sliding window test, 0: no, >0: window_size
 polar       polar         ingroup-"dsim/dmel"|outgroup-"dyak"|fixation-"dsim"
             unpolar       ingroup-"dmel"|outgroup-"dsim"
    """
    with open(fasta_list,'r') as p:
        lines = p.readlines()
        fastas = [line.strip() for line in lines]

    fout = open(foutput,'w')
    if ncpu==1:
        for f in fastas:
            results = test_one_gene(f,cutoff,polar,method,window_size,screenout,ingroup_name,outgroup_name,fixation_name)
            #print(results)
            if results:
                for result in results:
                    name,ps,pn,ds,dn,a,pval = result
                    fout.write("%s %4d %4d %4d %4d %8.3f %8.2e\n"%(name,pn,dn,ps,ds,a,pval))

    else:
        pool = Pool(ncpu)
        #arglist = []
        results_all = []
        for f in fastas:
        #    arglist += [(f,cutoff,polar,window_size)]
            #r = pool.apply_async(test_one_gene,
            pool.apply_async(test_one_gene,
                                 args=(f,cutoff,polar,method,window_size,
                                    screenout,ingroup_name,outgroup_name,
                                    fixation_name),
                                    callback=results_all.append) 
                                 #callback=mycallback).get()
        #    results.append(r)
        #for r in results:
        #    r.wait()
        #results = pool.starmap(test_one_gene,arglist) 
        pool.close()
        pool.join()
        for results in results_all:
            if results:
            #print(results)
                for result in results:
                    name,ps,pn,ds,dn,a,pval = result
                    fout.write("%s %4d %4d %4d %4d %8.3f %8.2e\n"%(name,pn,dn,ps,ds,a,pval))
    fout.close()

def calcMK(fasta_list,cutoff,method,window,polar,ncpu,screenout='yes',ingroup_name='dmel',outgroup_name='dsim',fixation_name=None):
    """ Python Package for McDonald–Kreitman test (pyMKT)

 -Following inputs are required
    fasta_list:     filename containing a list of alignments
    cutoff:         a cutoff value or a list of cutoff values
    method:         "shortest_path"; or "MK.pl"
    window_size:    0 for whole gene; >0 for sliding window MKTest
    method:         "polar"; "unpolar_dmel"; or "unpolar_dsim"
    ncpu:           1 for serial; > 1 for parallel
    screenout:      if 'yes': print output to screen; else: no screen output

 -Output format in "{polar}-cutoff{cutoff}-window{window_size}.out"
    mel_sim_yak Pn Dn Ps Ds alpha Pval (window=0)
    mel_sim_yak_{DNAseq_window}_{Seqindex} Pn Dn Ps Ds alpha Pval (window>0)
    """
    if type(cutoff)==float:
        cutoff = [cutoff]
    for c in cutoff:
        if window>0:
            foutput = '%s_%s-cutoff%s-window%s.1.out'%(polar,ingroup_name,c,window)
        else:
            foutput = '%s_%s-cutoff%s.1.out'%(polar,ingroup_name,c)
        test_all_genes(fasta_list,c,foutput,method,window_size=window,
                        polar=polar,ncpu=ncpu,screenout=screenout,
                        ingroup_name=ingroup_name,outgroup_name=outgroup_name,
                        fixation_name=fixation_name)

if __name__ == '__main__':
    fasta_list = 'alignment_list.txt'
    window     = 0
    polar      = 'unpolar'
    ingroup_name = 'dmel'
    outgroup_name = 'dsim'
    #polar      = 'polar'
    method     = 'MK.pl'
    ncpu       = 64
    #for cutoff in [0.05,0.10,0.15,0.20,0.30,0.50]:
    cutoff     = 0.05
    calcMK(fasta_list,cutoff,method,window,polar,ncpu)
#    for cutoff in [0.05]:
#        foutput    = '%s-cutoff%s.out'%(polar,cutoff)
#        test_all_genes(fasta_list,cutoff,foutput,method,
#                       window_size=window,polar=polar,ncpu=ncpu)

#    fasta_list = 'test_list.txt'
#    window     = 0
#    polar      = 'unpolar_dmel'
#    #polar      = 'polar'
#    method     = 'MK.pl'
#    ncpu       = 1
#    #for cutoff in [0.05,0.10,0.15,0.20,0.30,0.50]:
#    for cutoff in [0.05]:
#        foutput    = 'test_%s-cutoff%s.out'%(polar,cutoff)
#        test_all_genes(fasta_list,cutoff,foutput,method,
#                       window_size=window,polar=polar,ncpu=ncpu)

    #fasta_list = 'test_list.txt'
#    fasta_list = 'alignment_list.txt'
#    window     = 120
#    polar      = 'unpolar_dmel'
#    method     = 'MK.pl'
#    cutoff     = 0.05
#    foutput    = '%s-cutoff%s-window%s.out'%(polar,cutoff,window)
#    ncpu = 64
#    test_all_genes(fasta_list,cutoff,foutput,method,
#                        window_size=window,polar=polar,ncpu=ncpu)
