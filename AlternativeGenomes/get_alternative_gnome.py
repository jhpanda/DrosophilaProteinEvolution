"""

    The script is to get alternative genome from VCF. 

    The Dmel VCF file is downloaded from dgrp2.gnets.ncsu.edu/data.html, 
    containing genotypes of 205 Dmel individuals

    The Dsim VCF file is downloaded from zenodo.org/record/154261#.XY0Q-MDPxUQ 
    containing genotypes of 195 Dsim individuals

    The script generates 2*individuals alternative genomes excluding indels.

    Created by jpeng; 09/11/2019

"""

import sys,copy,subprocess,gzip
import numpy as np
from multiprocessing import Pool
from fasta2seq import *
#from tqdm import tqdm


def update(pbar,*a):
    pbar.update()

def readlines(f):
    try:
        #with open(f,'r') as p:
        #    lines = p.readlines()
        lines  = open(f,'r')
        nlines = int(subprocess.check_output("wc -l %s"%f,shell=True).split()[0])
    except IOError:
        print("cannot read file %s"%f)
        sys.exit(0)
    return lines,nlines

class Map():
    """
        Class objects for single site v5-v6 mapping
        - self.v5coor: store V5 coordinates of the mapping, e,g. "2L:1989"
        - self.v6coor: store V6 coordinates of the mapping, e,g. "2L:1989"
        - self.chrom: store the chromosome information
        - self.v6pos: store the V6 sites on the self.chrom chromosome.
    """
    def __init__(self,coor_str5,coor_str6):
        coor5 = coor_str5.split(':')
        coor6 = coor_str6.split(':')
        chrom = coor6[0]
        v5pos = int(coor5[1])
        v6pos = int(coor6[1])
        self.v5coor = coor_str5
        self.v6coor = coor_str5
        self.chrom  = chrom
        self.v6pos  = v6pos

class SNP():
    """
        Class objects for single site SNP
        - self.chrom: store which chromosome the SNP site locates on.
        - self.v6pos: store the specific site of the SNP on self.chrom.
        - self.variants: store all the genotype information of the SNP site.
    """
    def __init__(self,v6map,variants):
        self.chrom    = v6map.chrom
        self.v6pos    = v6map.v6pos
        self.variants = variants

class VCFReader():
    """
        Class objects for all sites v5-v6 mapping and genotypes
        - self.coormap: all the SNP sites v5-v6 mapping: {'V5coor':Map, ...}
        - self.snp_coordinates: all the SNP sites v6coor
        - self.snp_variants: all the SNP sites genotype (the keys are v6coor)
    """
    coormap = {} # this attribute is to store v5-v6 mapping
    snp_coordinates = [] # this attribute is to store coor of all SNP sites
    snp_variants = {} # this attribute is to store v6-variants mapping
    
    def process_map_line(self,line):
        if line.startswith('#'):
            pass
        else:
            coor5,coor6 = line.split()[0:2]
            v6map = Map(coor5,coor6)
            self.coormap[coor5] = v6map

    def read_mapping(self,fmap,ncpu=1):
        """ - Read in a Flybase V5-V6 coordinate mapping file. 
            - Add the mapping to self.coormap

Mapping file format:
#Number of properly formatted coordinates found = 4438427
#Number of skipped lines = 0
#Number of coordinates converted = 4438327
#Number of failed conversions = 100
#Original   Converted   Notes
2L:4998 2L:4998
        """

        lines,nlines = readlines(fmap)
        #plines  = tqdm(lines,ascii=True, unit_scale=True, unit_divisor=1000)
        #plines.set_description("ReadMap")
        #for line in plines:
        freq = max(10**(int(np.log(nlines)/np.log(10))-3),1)
        i_now = 0
        for line in lines:
            self.process_map_line(line)
            if i_now==0 or i_now%freq==0 or i_now==nlines-1:
                percent = i_now/(nlines-1) * 100
                sys.stdout.write("\rReadMap... %5.1f%%\r" % percent)
            i_now += 1
        sys.stdout.write('\nFinished v5-v6 mapping\n')
        ## avoid using multiprocessing here, as it takes more time ?##
        #pool = Pool(ncpu)
        #with tqdm(total=len(lines),ascii=True, unit_scale=True, unit_divisor=1000) as pbar:
        #    for line in lines:
        #        pool.apply_async(self.process_map_line,args=(line,),callback=update(pbar))
        #    pool.close()
        #    pool.join()
        #pool.starmap(self.process_map_line,lines)
        #pool.close()
        #pool.join()
        #print('Finished v5-v6 mapping')

    def add_variants(self,v6map,variants):
        v6coor = v6map.v6coor
        snp    = SNP(v6map,variants)
        self.snp_variants[v6coor] = snp

    def variantsMoreThan5percent(self,snps,ref_nt,var_nt,delimiter):
        """ Returns variants if we find more than 5% of all allels have SNPs.
            In Dmel case, if a site has less than (5%*410=) 20 SNPs, this SNP 
            is disregarded in future analysis.

            delimiter can be '/' or '|'
        """
        ntot = 0
        nvar = 0
        variants = []
        for snpi in snps:
            #vari = snpi.split('/')
            vari = snpi.split(delimiter)
            i = 0 if vari[0]=='1' else 1
            j = 0 if vari[1]=='1' else 1
            a0 = var_nt if vari[0]=='1' else ref_nt
            a1 = var_nt if vari[1]=='1' else ref_nt
            ntot += 2
            nvar += i+j
            variants += [a0,a1]
        if nvar/ntot >= 0.05:
            return variants
        if nvar/ntot <  0.05:
            return False

    def process_vcf_line(self,line,delimiter):
        """specify delimiter that seperate the genotype ('/' or '|')
        """
        if not line.startswith("#"): 
            info  = line.split()
            chrom = info[0] 
            v5pos = info[1]
            ref_nt= info[3]
            var_nt= info[4]
            snps  = info[9:]
            if len(ref_nt)==1 and len(var_nt)==1:
                ## to determin if a SNP is significant ##
                variants = self.variantsMoreThan5percent(snps,ref_nt,var_nt,delimiter)
                if variants:
                    coor5 = "%s:%s"%(chrom,v5pos)
                    try:
                        v6map = self.coormap[coor5]
                        self.snp_coordinates += [v6map.v6coor]
                        self.add_variants(v6map,variants)
                    except KeyError:
                        pass
                    #print(">Warn! Unmapped key: %s"%coor5)
        else:
            pass
         

    def read_vcf(self,fvcf,coor_version='V5',ncpu=1,delimiter='/'):
        """ - Reads in a VCF with 'V5' version coordinates by default. 
            - Please specify delimiter that seperate the genotype ('/' or '|')
            - Reutrns:
              - the list of 'v5coor', self.snp_coordinates
              - the dictionary self.snp_variants: {'v6coor': [variants], ...}
        """
        print('Opening %s to read'%fvcf)
        lines,nlines = readlines(fvcf)
        print('%s opened'%fvcf)
        #plines  = tqdm(lines,ascii=True,unit='B', unit_scale=True, unit_divisor=1024)
        #plines  = tqdm(lines,ascii=True, unit_scale=True, unit_divisor=1000)
        #plines.set_description("ReadVCF")
        #for line in plines:
        freq = max(10**(int(np.log(nlines)/np.log(10))-3),1)
        i_now = 0
        for line in lines:
            self.process_vcf_line(line,delimiter)
            if i_now==0 or i_now%freq==0 or i_now==nlines-1:
                percent = (i_now/(nlines-1))*100
                sys.stdout.write("\rReadVCF... %5.2f%%\r" % percent)
            i_now += 1
        sys.stdout.write('\nFinished VCF file\n')
        #print('Processing vcf')
        #pool = Pool(ncpu)
        #pool.starmap(self.process_map_line,lines)
        #pool.close()
        #pool.join()
        #print('Finished vcf')


    def update_v5(self,fmap,fvcf,ncpu=1,delimiter='/'):
        """ If vcf uses V5 coordinates and Genome is V6 coordinates
            Please specify delimiter that seperate the genotype ('/' or '|')
        """
        self.read_mapping(fmap,ncpu=ncpu) 
        self.read_vcf(fvcf,coor_version='V5',ncpu=ncpu,delimiter=delimiter)
    def update_v6(self,fmap,fvcf,ncpu=1,delimiter='/'):
        """ If vcf and Genome use same coordinates 
            But we can assume they are different by providing a pseudo-map,
            where v5 and v6 are the same.
            So just use update_v5 instead.

            Please specify delimiter that seperate the genotype ('/' or '|')
        """
        self.read_mapping(fmap,ncpu=ncpu) 
        self.read_vcf(fvcf,coor_version='V5',ncpu=ncpu,delimiter=delimiter)
        

class GenomeReader():
    def read_genome(self,fasta):
        """ Read the genome file
            To save time & memory, it's better to remove line breaks of 
            the genome fasta file. 
            *One line per chromosome*
             >2L
             ATCGATCG....//....ATCGATCG
             >X
             ATCGATCG....//....ATCGATCG
        """
        print('Opening %s to read'%fasta)
        self.sequences = fasta2seq(fasta,progress_bar='on')
        print('Finished reading %s'%fasta)

    def get_snp(self,vcf):
        """ Read the genome file
            Get SNP information from vcf file.
            - initialize self.alternative_genome
            - initialize self.snp_coordinates from vcf.snp_coordinates
            - initialize self.snp_variants from vcf.snp_variants
            - Returns:
              - self.snp_num: total number of SNPs (excluding indels)
              - self.npopulations: total number of alter genomes (410 for Dmel)
        """
        self.snp_coordinates = vcf.snp_coordinates
        self.snp_variants    = vcf.snp_variants
        self.snp_num = len(self.snp_coordinates)

        ## populations ##
        npopulations = len(vcf.snp_variants[vcf.snp_coordinates[0]].variants)
        self.npopulations = npopulations

    def mutate_one(self,ialt,outdir):
        """ Mutate genome sequence of the ialt_th genome 
        """

        ## initialize the alternative genome ##
        alternative_genome = dict()
        for key in self.sequences:
            alternative_genome[key] = list(self.sequences[key])

        ## mutate and write this alternative genome ##
        i = 1
        freq = max(10**(int(np.log(self.snp_num)/np.log(10))-3),1)
        for v6coor in self.snp_coordinates:
            snp     = self.snp_variants[v6coor]
            variant = snp.variants[ialt]
            chrom   = snp.chrom
            v6pos   = snp.v6pos
            try:
                alternative_genome[chrom][v6pos-1] = variant
            except (IndexError,KeyError):
                pass
            if i==0 or i%freq==0 or i==self.snp_num:
                percent = i/self.snp_num * 100
                sys.stdout.write("\rMutating %d... %5.1f%%\r"%(ialt,percent))
            i+=1

        ## now write the mutated genome ##
        fasta = '%s/DGRP2-%03d.fasta'%(outdir,ialt+1)
        for key in alternative_genome:
            alternative_genome[key] = ''.join(alternative_genome[key])
        seq2fasta(alternative_genome,fasta,linebreak=False,compress=True)
        
        print('Finished mutating individual %d'%(ialt+1))

    def mutate_and_write(self,outdir,ncpu = 1):
        """ Final main() 
            Note that the multiprocessing is not optimized, 
            it may not save time
        """
        if ncpu==1:
            for i in range(self.npopulations):
                self.mutate_one(i,outdir)        

        elif ncpu > 1:
            #for ialt in range(10):
            pool = Pool(ncpu)
            for i in range(self.npopulations):
                pool.apply_async(self.mutate_one,args=(i,outdir,))
            pool.close()
            pool.join()

if __name__ == '__main__':

    fmap = 'positions_latest.txt'
    #fvcf = 'ALL.clean.vcf'
    #fvcf = 'ALL.valid.vcf'
    fvcf = 'test.vcf'
    #fmap = 'test.txt'
    #fvcf = 'dgrp2-test.vcf'
    vcf = VCFReader()

    #vcf.read_mapping(fmap,ncpu=64)
    #vcf.update_v5(fmap,fvcf,ncpu=64)
    vcf.update_v5(fmap,fvcf,ncpu=1,delimiter='|')

    #fasta  = 'fasta-ref/dsim-all-chromosome-r2.01.fasta.nobreak'
    fasta  = 'GRCh37.primary_assembly.genome.fa.nobreak'
    #fasta  = 'fasta-ref/dmel-test.fasta'
    outdir = 'fasta-AltGenomeM5pSNPs'
    genome = GenomeReader()
    genome.read_genome(fasta)
    genome.get_snp(vcf)
    #genome.mutate_and_write(outdir,ncpu = 4)
    genome.mutate_and_write(outdir,ncpu = 16)
