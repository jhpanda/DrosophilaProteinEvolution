
codon_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X', 'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

stop_codon = ['TAG','TAA','TGA']
gaps = ['-','N']

from itertools import permutations

class Graph():
    """
        A codon graph is constructed to search for minimum paths from outgroup
        codon to all ingrp codons.
    """
    def __init__(self,graph_dict,vertices,edges):
        """ all vertices and edges are pre-computed as above
        """
        self.graph_dict = graph_dict
        self.vertices = vertices
        self.edges = edges


class CodonPair():
    """
        Minimum paths for ingroup codons to outgroup codon is constructed
        no intermediate codons even if there are multiple changes
    """

    def dist_BASE(self,codon1,codon2):
        """ minmum subsitution from codon1 to codon2
        """
        d = [0 if codon1[m]==codon2[m] else 1 for m in range(3)]
        return sum(d) 

    def dist_AA(self,codon1,codon2):
        """ distance from codon1 to codon2 (base + AA substitutions)
        """
        d_  = 0 if codon_table[codon1]==codon_table[codon2] else 1
        d_ += 0 if codon1==codon2 else self.dist_BASE(codon1,codon2)
        return d_ 

    def dist_path(self,path):
        path_len = len(path)
        d = 0
        for i in range(path_len-1):
            d += self.dist_AA(path[i],path[i+1])
        return d

    def shortest_path(self,codons,outgroup_codon):
        """ Shortest path that include all unique codons.
            The shortest path have the minimum AA and BASE substitutions.
        """
        ncodon = len(codons)
        codon0 = outgroup_codon
        ## we only find path from ingroup codons to outgroup codons
        ## distance of each codon to outgroup codon are calculated
        ## the path is like this: outgroup -> closest -> distant 
        unique_ingrp_codons = []
        for codon in codons:
            if codon!=codon0:
                unique_ingrp_codons += [codon]

        min_d = 1e10
        codon_dist = {}
        for codon in unique_ingrp_codons:
            d = self.dist_AA(codon,codon0)
            codon_dist[codon] = d
            if d<min_d:
                min_d = d
        min_dist_codons = []
        for codon in unique_ingrp_codons:
            d = codon_dist[codon]
            if d == min_d:
                min_dist_codons += [codon]
        #print(min_dist_codons)

        raw_paths = []
        for codon1 in min_dist_codons:
            codons_to_permute = []
            for codon in unique_ingrp_codons:
                if codon!=codon1 and codon!=codon0:
                    codons_to_permute += [codon]
            permuts = permutations(codons_to_permute)
            raw_paths += [[codon0]+[codon1]+list(s) for s in permuts]
        #print(raw_paths)
        
        #print(max_d,min_d,max_dist_codons)
        min_d = 1e10
        raw_path_dist = {}
        for raw_path in raw_paths:
            d = self.dist_path(raw_path)
            raw_path_dist[tuple(raw_path)] = d
            if d<min_d:
                min_d = d
                min_path = raw_path

        #min_paths = []
        #for raw_path in raw_paths:
        #    if raw_path_dist[tuple(raw_path)] == min_d:
        #        min_paths += [raw_path]

        return min_path

pair = CodonPair()

if __name__ == '__main__':
    #paths = graph.find_path('AAA','AAT')
    #print(paths)
    #print(paths)
    #paths = graph.path2_dfs('CCC','CAG',2)
    #print(paths)
    #mpath = graph.shortest_path(['ACG','AGG','AAT'],'AGG')
    #mpath = graph.shortest_path(['CCC','CAG'],'CAC')
    #mpath = graph.shortest_path(['TTG', 'CTC'], 'TTC')
    #mpath = graph.shortest_path(['TTC', 'TGC'], 'AGC')
    #mpath = graph.shortest_path(['AGT'], 'TGC')
    #mpath = graph.shortest_path(['GCC'], 'ACT')
    #pair  = CodonPair()
    #mpaths = pair.shortest_path(['GCC'], 'ACG')
    #mpaths = pair.shortest_path(['AGG', 'AGC', 'AGA', 'AGT'],'AGG')
    #print(mpath)
    #for i in range(1):
    mpaths = pair.shortest_path(['AGG', 'AGC', 'AGA', 'AGT','CGG'],'AGG')
    #mpaths = pair.shortest_path(['AAT', 'AGT', 'AAA','AGA'], 'AAT')
    #mpaths = pair.shortest_path(['CGG','CGG','CAA','GAA','CGA','CAG','GAG','AGC'],'ATG')
    #mpath = graph.shortest_path(['CGG','CAA','GAA','CGA','CAG','GAG','AGC'],'ATG')
    #mpaths = pair.shortest_path(['GTC'],'GCT')
    #mpaths = pair.shortest_path(['GTC'],'GTC')
    print(mpaths)
