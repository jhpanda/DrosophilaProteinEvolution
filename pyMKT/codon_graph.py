""" 
    Codon graph to search all possible paths.
    Created by Junhui Peng, 10/16/2019
"""

codon_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X', 'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

stop_codon = ['TAG','TAA']
gaps = ['-','N']

graph_dict = {
    'TTT':['ATT', 'CTT', 'GTT', 'TAT', 'TCT', 'TGT', 'TTA', 'TTC', 'TTG'],
    'TTC':['TTT', 'ATC', 'CTC', 'GTC', 'TAC', 'TCC', 'TGC', 'TTA', 'TTG'],
    'TTA':['TTT', 'TTC', 'ATA', 'CTA', 'GTA', 'TAA', 'TCA', 'TGA', 'TTG'],
    'TTG':['TTT', 'TTC', 'TTA', 'ATG', 'CTG', 'GTG', 'TAG', 'TCG', 'TGG'],
    'TCT':['TTT', 'ACT', 'CCT', 'GCT', 'TAT', 'TGT', 'TCA', 'TCC', 'TCG'],
    'TCC':['TTC', 'TCT', 'ACC', 'CCC', 'GCC', 'TAC', 'TGC', 'TCA', 'TCG'],
    'TCA':['TTA', 'TCT', 'TCC', 'ACA', 'CCA', 'GCA', 'TAA', 'TGA', 'TCG'],
    'TCG':['TTG', 'TCT', 'TCC', 'TCA', 'ACG', 'CCG', 'GCG', 'TAG', 'TGG'],
    'TAT':['TTT', 'TCT', 'AAT', 'CAT', 'GAT', 'TGT', 'TAA', 'TAC', 'TAG'],
    'TAC':['TTC', 'TCC', 'TAT', 'AAC', 'CAC', 'GAC', 'TGC', 'TAA', 'TAG'],
    'TAA':['TTA', 'TCA', 'TAT', 'TAC', 'AAA', 'CAA', 'GAA', 'TGA', 'TAG'],
    'TAG':['TTG', 'TCG', 'TAT', 'TAC', 'TAA', 'AAG', 'CAG', 'GAG', 'TGG'],
    'TGT':['TTT', 'TCT', 'TAT', 'AGT', 'CGT', 'GGT', 'TGA', 'TGC', 'TGG'],
    'TGC':['TTC', 'TCC', 'TAC', 'TGT', 'AGC', 'CGC', 'GGC', 'TGA', 'TGG'],
    'TGA':['TTA', 'TCA', 'TAA', 'TGT', 'TGC', 'AGA', 'CGA', 'GGA', 'TGG'],
    'TGG':['TTG', 'TCG', 'TAG', 'TGT', 'TGC', 'TGA', 'AGG', 'CGG', 'GGG'],
    'CTT':['TTT', 'ATT', 'GTT', 'CAT', 'CCT', 'CGT', 'CTA', 'CTC', 'CTG'],
    'CTC':['TTC', 'CTT', 'ATC', 'GTC', 'CAC', 'CCC', 'CGC', 'CTA', 'CTG'],
    'CTA':['TTA', 'CTT', 'CTC', 'ATA', 'GTA', 'CAA', 'CCA', 'CGA', 'CTG'],
    'CTG':['TTG', 'CTT', 'CTC', 'CTA', 'ATG', 'GTG', 'CAG', 'CCG', 'CGG'],
    'CCT':['TCT', 'CTT', 'ACT', 'GCT', 'CAT', 'CGT', 'CCA', 'CCC', 'CCG'],
    'CCC':['TCC', 'CTC', 'CCT', 'ACC', 'GCC', 'CAC', 'CGC', 'CCA', 'CCG'],
    'CCA':['TCA', 'CTA', 'CCT', 'CCC', 'ACA', 'GCA', 'CAA', 'CGA', 'CCG'],
    'CCG':['TCG', 'CTG', 'CCT', 'CCC', 'CCA', 'ACG', 'GCG', 'CAG', 'CGG'],
    'CAT':['TAT', 'CTT', 'CCT', 'AAT', 'GAT', 'CGT', 'CAA', 'CAC', 'CAG'],
    'CAC':['TAC', 'CTC', 'CCC', 'CAT', 'AAC', 'GAC', 'CGC', 'CAA', 'CAG'],
    'CAA':['TAA', 'CTA', 'CCA', 'CAT', 'CAC', 'AAA', 'GAA', 'CGA', 'CAG'],
    'CAG':['TAG', 'CTG', 'CCG', 'CAT', 'CAC', 'CAA', 'AAG', 'GAG', 'CGG'],
    'CGT':['TGT', 'CTT', 'CCT', 'CAT', 'AGT', 'GGT', 'CGA', 'CGC', 'CGG'],
    'CGC':['TGC', 'CTC', 'CCC', 'CAC', 'CGT', 'AGC', 'GGC', 'CGA', 'CGG'],
    'CGA':['TGA', 'CTA', 'CCA', 'CAA', 'CGT', 'CGC', 'AGA', 'GGA', 'CGG'],
    'CGG':['TGG', 'CTG', 'CCG', 'CAG', 'CGT', 'CGC', 'CGA', 'AGG', 'GGG'],
    'ATT':['TTT', 'CTT', 'GTT', 'AAT', 'ACT', 'AGT', 'ATA', 'ATC', 'ATG'],
    'ATC':['TTC', 'CTC', 'ATT', 'GTC', 'AAC', 'ACC', 'AGC', 'ATA', 'ATG'],
    'ATA':['TTA', 'CTA', 'ATT', 'ATC', 'GTA', 'AAA', 'ACA', 'AGA', 'ATG'],
    'ATG':['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'GTG', 'AAG', 'ACG', 'AGG'],
    'ACT':['TCT', 'CCT', 'ATT', 'GCT', 'AAT', 'AGT', 'ACA', 'ACC', 'ACG'],
    'ACC':['TCC', 'CCC', 'ATC', 'ACT', 'GCC', 'AAC', 'AGC', 'ACA', 'ACG'],
    'ACA':['TCA', 'CCA', 'ATA', 'ACT', 'ACC', 'GCA', 'AAA', 'AGA', 'ACG'],
    'ACG':['TCG', 'CCG', 'ATG', 'ACT', 'ACC', 'ACA', 'GCG', 'AAG', 'AGG'],
    'AAT':['TAT', 'CAT', 'ATT', 'ACT', 'GAT', 'AGT', 'AAA', 'AAC', 'AAG'],
    'AAC':['TAC', 'CAC', 'ATC', 'ACC', 'AAT', 'GAC', 'AGC', 'AAA', 'AAG'],
    'AAA':['TAA', 'CAA', 'ATA', 'ACA', 'AAT', 'AAC', 'GAA', 'AGA', 'AAG'],
    'AAG':['TAG', 'CAG', 'ATG', 'ACG', 'AAT', 'AAC', 'AAA', 'GAG', 'AGG'],
    'AGT':['TGT', 'CGT', 'ATT', 'ACT', 'AAT', 'GGT', 'AGA', 'AGC', 'AGG'],
    'AGC':['TGC', 'CGC', 'ATC', 'ACC', 'AAC', 'AGT', 'GGC', 'AGA', 'AGG'],
    'AGA':['TGA', 'CGA', 'ATA', 'ACA', 'AAA', 'AGT', 'AGC', 'GGA', 'AGG'],
    'AGG':['TGG', 'CGG', 'ATG', 'ACG', 'AAG', 'AGT', 'AGC', 'AGA', 'GGG'],
    'GTT':['TTT', 'CTT', 'ATT', 'GAT', 'GCT', 'GGT', 'GTA', 'GTC', 'GTG'],
    'GTC':['TTC', 'CTC', 'ATC', 'GTT', 'GAC', 'GCC', 'GGC', 'GTA', 'GTG'],
    'GTA':['TTA', 'CTA', 'ATA', 'GTT', 'GTC', 'GAA', 'GCA', 'GGA', 'GTG'],
    'GTG':['TTG', 'CTG', 'ATG', 'GTT', 'GTC', 'GTA', 'GAG', 'GCG', 'GGG'],
    'GCT':['TCT', 'CCT', 'ACT', 'GTT', 'GAT', 'GGT', 'GCA', 'GCC', 'GCG'],
    'GCC':['TCC', 'CCC', 'ACC', 'GTC', 'GCT', 'GAC', 'GGC', 'GCA', 'GCG'],
    'GCA':['TCA', 'CCA', 'ACA', 'GTA', 'GCT', 'GCC', 'GAA', 'GGA', 'GCG'],
    'GCG':['TCG', 'CCG', 'ACG', 'GTG', 'GCT', 'GCC', 'GCA', 'GAG', 'GGG'],
    'GAT':['TAT', 'CAT', 'AAT', 'GTT', 'GCT', 'GGT', 'GAA', 'GAC', 'GAG'],
    'GAC':['TAC', 'CAC', 'AAC', 'GTC', 'GCC', 'GAT', 'GGC', 'GAA', 'GAG'],
    'GAA':['TAA', 'CAA', 'AAA', 'GTA', 'GCA', 'GAT', 'GAC', 'GGA', 'GAG'],
    'GAG':['TAG', 'CAG', 'AAG', 'GTG', 'GCG', 'GAT', 'GAC', 'GAA', 'GGG'],
    'GGT':['TGT', 'CGT', 'AGT', 'GTT', 'GCT', 'GAT', 'GGA', 'GGC', 'GGG'],
    'GGC':['TGC', 'CGC', 'AGC', 'GTC', 'GCC', 'GAC', 'GGT', 'GGA', 'GGG'],
    'GGA':['TGA', 'CGA', 'AGA', 'GTA', 'GCA', 'GAA', 'GGT', 'GGC', 'GGG'],
    'GGG':['TGG', 'CGG', 'AGG', 'GTG', 'GCG', 'GAG', 'GGT', 'GGC', 'GGA'],
}

vertices = ['TTT', 'TTC', 'TTA', 'TTG', 'TCT', 'TCC', 'TCA', 'TCG', 
            'TAT', 'TAC', 'TAA', 'TAG', 'TGT', 'TGC', 'TGA', 'TGG', 
            'CTT', 'CTC', 'CTA', 'CTG', 'CCT', 'CCC', 'CCA', 'CCG', 
            'CAT', 'CAC', 'CAA', 'CAG', 'CGT', 'CGC', 'CGA', 'CGG', 
            'ATT', 'ATC', 'ATA', 'ATG', 'ACT', 'ACC', 'ACA', 'ACG', 
            'AAT', 'AAC', 'AAA', 'AAG', 'AGT', 'AGC', 'AGA', 'AGG', 
            'GTT', 'GTC', 'GTA', 'GTG', 'GCT', 'GCC', 'GCA', 'GCG', 
            'GAT', 'GAC', 'GAA', 'GAG', 'GGT', 'GGC', 'GGA', 'GGG']

edges = [{'TTT','ATT'}, {'CTT','TTT'}, {'GTT','TTT'}, {'TAT','TTT'}, 
         {'TCT','TTT'}, {'TGT','TTT'}, {'TTA','TTT'}, {'TTT','TTC'}, 
         {'TTG','TTT'}, {'ATC','TTC'}, {'CTC','TTC'}, {'GTC','TTC'}, 
         {'TTC','TAC'}, {'TCC','TTC'}, {'TGC','TTC'}, {'TTA','TTC'}, 
         {'TTG','TTC'}, {'TTA','ATA'}, {'TTA','CTA'}, {'GTA','TTA'}, 
         {'TAA','TTA'}, {'TTA','TCA'}, {'TGA','TTA'}, {'TTA','TTG'}, 
         {'ATG','TTG'}, {'TTG','CTG'}, {'GTG','TTG'}, {'TAG','TTG'}, 
         {'TTG','TCG'}, {'TGG','TTG'}, {'TCT','ACT'}, {'TCT','CCT'}, 
         {'TCT','GCT'}, {'TCT','TAT'}, {'TCT','TGT'}, {'TCT','TCA'}, 
         {'TCT','TCC'}, {'TCT','TCG'}, {'TCC','ACC'}, {'CCC','TCC'}, 
         {'TCC','GCC'}, {'TCC','TAC'}, {'TGC','TCC'}, {'TCC','TCA'}, 
         {'TCC','TCG'}, {'TCA','ACA'}, {'CCA','TCA'}, {'GCA','TCA'}, 
         {'TAA','TCA'}, {'TGA','TCA'}, {'TCG','TCA'}, {'TCG','ACG'}, 
         {'TCG','CCG'}, {'TCG','GCG'}, {'TAG','TCG'}, {'TGG','TCG'}, 
         {'AAT','TAT'}, {'CAT','TAT'}, {'GAT','TAT'}, {'TGT','TAT'}, 
         {'TAA','TAT'}, {'TAT','TAC'}, {'TAG','TAT'}, {'AAC','TAC'}, 
         {'CAC','TAC'}, {'TAC','GAC'}, {'TGC','TAC'}, {'TAA','TAC'}, 
         {'TAG','TAC'}, {'TAA','AAA'}, {'TAA','CAA'}, {'GAA','TAA'}, 
         {'TGA','TAA'}, {'TAG','TAA'}, {'TAG','AAG'}, {'TAG','CAG'}, 
         {'TAG','GAG'}, {'TGG','TAG'}, {'AGT','TGT'}, {'CGT','TGT'}, 
         {'TGT','GGT'}, {'TGA','TGT'}, {'TGC','TGT'}, {'TGG','TGT'}, 
         {'AGC','TGC'}, {'TGC','CGC'}, {'GGC','TGC'}, {'TGA','TGC'}, 
         {'TGG','TGC'}, {'TGA','AGA'}, {'TGA','CGA'}, {'TGA','GGA'}, 
         {'TGG','TGA'}, {'TGG','AGG'}, {'TGG','CGG'}, {'TGG','GGG'}, 
         {'CTT','ATT'}, {'GTT','CTT'}, {'CAT','CTT'}, {'CTT','CCT'}, 
         {'CGT','CTT'}, {'CTA','CTT'}, {'CTC','CTT'}, {'CTT','CTG'}, 
         {'CTC','ATC'}, {'GTC','CTC'}, {'CTC','CAC'}, {'CCC','CTC'}, 
         {'CTC','CGC'}, {'CTC','CTA'}, {'CTC','CTG'}, {'CTA','ATA'}, 
         {'GTA','CTA'}, {'CAA','CTA'}, {'CTA','CCA'}, {'CTA','CGA'}, 
         {'CTA','CTG'}, {'ATG','CTG'}, {'GTG','CTG'}, {'CAG','CTG'}, 
         {'CTG','CCG'}, {'CGG','CTG'}, {'ACT','CCT'}, {'GCT','CCT'}, 
         {'CAT','CCT'}, {'CGT','CCT'}, {'CCA','CCT'}, {'CCC','CCT'}, 
         {'CCT','CCG'}, {'CCC','ACC'}, {'CCC','GCC'}, {'CCC','CAC'}, 
         {'CCC','CGC'}, {'CCC','CCA'}, {'CCC','CCG'}, {'CCA','ACA'}, 
         {'GCA','CCA'}, {'CAA','CCA'}, {'CGA','CCA'}, {'CCA','CCG'}, 
         {'ACG','CCG'}, {'GCG','CCG'}, {'CAG','CCG'}, {'CGG','CCG'}, 
         {'CAT','AAT'}, {'GAT','CAT'}, {'CGT','CAT'}, {'CAA','CAT'}, 
         {'CAC','CAT'}, {'CAT','CAG'}, {'AAC','CAC'}, {'CAC','GAC'}, 
         {'CAC','CGC'}, {'CAA','CAC'}, {'CAC','CAG'}, {'CAA','AAA'}, 
         {'GAA','CAA'}, {'CAA','CGA'}, {'CAA','CAG'}, {'AAG','CAG'}, 
         {'GAG','CAG'}, {'CGG','CAG'}, {'CGT','AGT'}, {'CGT','GGT'}, 
         {'CGT','CGA'}, {'CGT','CGC'}, {'CGT','CGG'}, {'AGC','CGC'}, 
         {'GGC','CGC'}, {'CGA','CGC'}, {'CGG','CGC'}, {'AGA','CGA'}, 
         {'GGA','CGA'}, {'CGG','CGA'}, {'CGG','AGG'}, {'CGG','GGG'}, 
         {'GTT','ATT'}, {'AAT','ATT'}, {'ACT','ATT'}, {'AGT','ATT'}, 
         {'ATT','ATA'}, {'ATC','ATT'}, {'ATG','ATT'}, {'GTC','ATC'}, 
         {'AAC','ATC'}, {'ATC','ACC'}, {'AGC','ATC'}, {'ATC','ATA'}, 
         {'ATG','ATC'}, {'GTA','ATA'}, {'AAA','ATA'}, {'ATA','ACA'}, 
         {'AGA','ATA'}, {'ATG','ATA'}, {'GTG','ATG'}, {'ATG','AAG'}, 
         {'ATG','ACG'}, {'ATG','AGG'}, {'ACT','GCT'}, {'ACT','AAT'}, 
         {'AGT','ACT'}, {'ACT','ACA'}, {'ACT','ACC'}, {'ACT','ACG'}, 
         {'GCC','ACC'}, {'AAC','ACC'}, {'AGC','ACC'}, {'ACC','ACA'}, 
         {'ACC','ACG'}, {'GCA','ACA'}, {'AAA','ACA'}, {'AGA','ACA'}, 
         {'ACG','ACA'}, {'ACG','GCG'}, {'AAG','ACG'}, {'AGG','ACG'}, 
         {'GAT','AAT'}, {'AGT','AAT'}, {'AAT','AAA'}, {'AAC','AAT'}, 
         {'AAG','AAT'}, {'AAC','GAC'}, {'AGC','AAC'}, {'AAC','AAA'}, 
         {'AAC','AAG'}, {'GAA','AAA'}, {'AGA','AAA'}, {'AAG','AAA'}, 
         {'GAG','AAG'}, {'AGG','AAG'}, {'AGT','GGT'}, {'AGT','AGA'}, 
         {'AGT','AGC'}, {'AGT','AGG'}, {'GGC','AGC'}, {'AGC','AGA'}, 
         {'AGC','AGG'}, {'GGA','AGA'}, {'AGG','AGA'}, {'AGG','GGG'}, 
         {'GAT','GTT'}, {'GTT','GCT'}, {'GTT','GGT'}, {'GTT','GTA'}, 
         {'GTC','GTT'}, {'GTG','GTT'}, {'GTC','GAC'}, {'GTC','GCC'}, 
         {'GGC','GTC'}, {'GTC','GTA'}, {'GTG','GTC'}, {'GAA','GTA'}, 
         {'GCA','GTA'}, {'GTA','GGA'}, {'GTG','GTA'}, {'GTG','GAG'}, 
         {'GTG','GCG'}, {'GTG','GGG'}, {'GAT','GCT'}, {'GGT','GCT'}, 
         {'GCA','GCT'}, {'GCC','GCT'}, {'GCT','GCG'}, {'GCC','GAC'}, 
         {'GGC','GCC'}, {'GCA','GCC'}, {'GCC','GCG'}, {'GAA','GCA'}, 
         {'GCA','GGA'}, {'GCA','GCG'}, {'GAG','GCG'}, {'GCG','GGG'}, 
         {'GAT','GGT'}, {'GAA','GAT'}, {'GAT','GAC'}, {'GAG','GAT'}, 
         {'GGC','GAC'}, {'GAA','GAC'}, {'GAG','GAC'}, {'GAA','GGA'}, 
         {'GAA','GAG'}, {'GAG','GGG'}, {'GGA','GGT'}, {'GGC','GGT'}, 
         {'GGT','GGG'}, {'GGC','GGA'}, {'GGC','GGG'}, {'GGA','GGG'}, ]

from itertools import permutations,combinations,product

class Graph():
    """
        A codon graph is constructed to search for all possible paths from 
        codon_i to codon_j.
        A depth-first search with maximum depth of current unique codon numbers
        would be able to find all possible paths.
    """
    def __init__(self,graph_dict,vertices,edges):
        """ all vertices and edges are pre-computed as above
        """
        self.graph_dict = graph_dict
        self.vertices = vertices
        self.edges = edges

    def dist_BASE(self,codon1,codon2):
        """ minmum subsitution from codon1 to codon2
        """
        d = [0 if codon1[m]==codon2[m] else 1 for m in range(3)]
        return sum(d) 

    def dist_AA(self,codon1,codon2):
        """ distance from codon1 to codon2 (base + AA substitutions)
        """
        d = 0 if codon_table[codon1]==codon_table[codon2] else 1
        d += 0 if codon1==codon2 else 1
        return d 

    def path2_dfs(self,start,end,depth,path=[]):
        """ Depth-first search with maximum depth
        """
        graph = self.graph_dict 
        path = path + [start]
        if start == end:
            return [path]
        if start not in graph:
            return []
        paths = []
        if depth>0:
            for vertex in graph[start]:
                if vertex not in path:
                    extended_paths = self.path2_dfs(vertex,end,depth-1,path)
                    for p in extended_paths: 
                        paths.append(p)
        return paths

    def dist_path(self,path):
        path_len = len(path)
        d = 0
        for i in range(path_len-1):
            d += self.dist_AA(path[i],path[i+1])
        return d

    def dist_path_BASE(self,path):
        path_len = len(path)
        d = 0
        for i in range(path_len-1):
            d += self.dist_BASE(path[i],path[i+1])
        return d

    def min_path(self,paths,codons):
        min_d = 1e10
        dist_path = {}
        for path in paths:
            if all(codon in codons for codon in path[1:-1]):
                d = 1e10
            else:
                d = self.dist_path(path)
            dist_path[tuple(path)] = d
            if d<min_d:
                min_d = d
        min_paths = []
        for path in paths:
            if dist_path[tuple(path)]  == min_d:
                min_paths += [path]
        return min_paths

    def shortest_path(self,codons,outgroup_codon):
        """ Shortest path that include all unique codons.
            The shortest path have the minimum AA and BASE substitutions.
        """
        if outgroup_codon in codons:
            codonsnew = [s for s in codons]
        else:
            codonsnew = [s for s in codons] + [outgroup_codon]
        ncodon = len(codonsnew)

        ## we only find path between closest codons and most distant codons
        ## closest codons have minimum sum of distances to all others
        ## distant codons have maximum sum of distances to all others
        max_d = 1e-10  ## find maximum sum of distances
        min_d = 1e10   ## find minimum sum of distances
        codon_dist = {} ## store sum of distances of each codon
        for i in range(ncodon):
            d = 0
            d_BASE = 0
            for j in range(ncodon):
                d += self.dist_AA(codonsnew[i],codonsnew[j])
                d_BASE += self.dist_BASE(codonsnew[i],codonsnew[j])
            codon_dist[codonsnew[i]] = d
            if d>max_d:
                max_d = d
            if d<min_d:
                min_d = d

        # raw paths contain all possible combinations of codons (in & outgroup)
        # we assume that the path starts from closest codons to distant codons
        max_dist_codons = [] ## distant codons
        min_dist_codons = [] ## closest codons
        for i in range(ncodon):
            if codon_dist[codonsnew[i]] == max_d:
                if codonsnew[i] not in max_dist_codons:
                    max_dist_codons += [codonsnew[i]] 
            if codon_dist[codonsnew[i]] == min_d:
                if codonsnew[i] not in min_dist_codons:
                    min_dist_codons += [codonsnew[i]] 
        other_codons = []
        for codon in codonsnew:
            if codon not in max_dist_codons+min_dist_codons:
                other_codons += [codon]

        # raw paths contain all possible combinations of codons (in & outgroup)
        # we assume that the path starts from closest codons to distant codons
        if max_d!=min_d:
            max_codons_permuts = list(permutations(max_dist_codons))
            min_codons_permuts = list(permutations(min_dist_codons))
            other_codons_permuts = list(permutations(other_codons))
            raw_paths  = [x+y+z for x in max_codons_permuts for y in other_codons_permuts for z in min_codons_permuts]
        else:
            max_codons_permuts = list(permutations(max_dist_codons))
            raw_paths = max_codons_permuts
        
        raw_path_dist = {}
        for raw_path in raw_paths:
            d = self.dist_path(raw_path)
            raw_path_dist[tuple(raw_path)] = d
            if d<min_d:
                min_d = d

        raw_paths_new = []
        for raw_path in raw_paths:
            if raw_path_dist[tuple(raw_path)] == min_d:
                raw_paths_new += [raw_path]

        possible_paths = []
        for raw_path in raw_paths_new:
            possible_path = []
            iters = []
            for i in range(len(raw_path)-1):
                d_BASE = self.dist_BASE(raw_path[i],raw_path[i+1])
                if d_BASE>1:
                    paths  = self.path2_dfs(raw_path[i],raw_path[i+1],d_BASE)
                    min_paths  = self.min_path(paths,codonsnew)
                    #print(raw_path,min_paths)
                    min_paths_ = [path[0:-1] for path in min_paths]
                    iters += [min_paths_]
                if d_BASE == 1:
                    #possible_path += [[raw_path[i]]]
                    iters += [[[raw_path[i]]]]
                    #print(d_BASE,raw_path[i])
            iters += [[[raw_path[-1]]]]
            possible_path  = list(product(*iters))
            npath = len(possible_path)
            #print(possible_path) 
            possible_path_ = [[] for s in range(npath)]
            for i in range(npath):
                for elements in possible_path[i]:
                    for element in elements:
                        possible_path_ [i] += [element]
            #print(possible_path_)
            for path in possible_path_:
                #print(path)
                possible_paths.append(path)
        #print(possible_paths)

        paths_contain_outgroup = []
        path_dist = {}
        for path in possible_paths:
            path_len = len(path)
            d = self.dist_path(path)
            path_dist[tuple(path)] = d

        min_d = 1e10
        for path in possible_paths:
            ### the path with minimum AA changes is selected
            d = path_dist[tuple(path)] 
            if d < min_d:
                min_d = d
                min_path = path

        min_d2 = 1e10
        for path in possible_paths:
            ### if the outgroup_codon is the start or end of the minimum path
            ### this path will be selected
            #print(path,path_dist[tuple(path)])
            if path_dist[tuple(path)] == min_d:
                #print(path,path_dist[tuple(path)])
                if path[-1] == outgroup_codon:
                    d = self.dist_path(path[:-1])
                    #print(path,path_dist[tuple(path)],d)
                    if d<min_d2:
                        min_d2 = d
                        min_path = path
                if path[0] == outgroup_codon:
                    d = self.dist_path(path[::-1][:-1])
                    #print(path,path_dist[tuple(path)],d)
                    if d<min_d2:
                        min_d2 = d
                        min_path = path[::-1]
                else:
                    pass
        #if min_path_1 and min_path_2:
        #    ## Sometimes, a path may have the same distance whether 
        #    ## the outgroup_codon is at start or end. To minimize changes,
        #    ## we will choose the path with minimum ingroup distance
        #    d1 = self.dist_path(min_path_1[:-1])
        #    d2 = self.dist_path(min_path_2[:-1])
        #    print(min_path_1,d1)
        #    print(min_path_2,d2)
        #    min_path = min_path_1 if d1<=d2 else min_path_2
        return min_path

graph = Graph(graph_dict,vertices,edges)

if __name__ == '__main__':
    #paths = graph.find_path('AAA','AAT')
    #print(paths)
    #print(paths)
    #paths = graph.path2_dfs('CCC','CAG',2)
    #print(paths)
    #mpath = graph.shortest_path(['ACG','AGG','AAT'],'AGG')
    #mpath = graph.shortest_path(['CCC','CAG'],'CAC')
    #mpath = graph.shortest_path(['AGG', 'AGC', 'AGA', 'AGT'],'AGG')
    #mpath = graph.shortest_path(['TTG', 'CTC'], 'TTC')
    #mpath = graph.shortest_path(['TTC', 'TGC'], 'AGC')
    #mpath = graph.shortest_path(['AGT'], 'TGC')
    #mpath = graph.shortest_path(['GCC'], 'ACT')
    #mpath = graph.shortest_path(['GCC'], 'ACG')
    #print(mpath)
    for i in range(100):
        mpath = graph.shortest_path(['AGG', 'AGC', 'AGA', 'AGT','CGG'],'AGG')
    #mpath = graph.shortest_path(['CGG','CGG','CAA','GAA','CGA','CAG','GAG','AGC'],'ATG')
    #mpath = graph.shortest_path(['CGG','CAA','GAA','CGA','CAG','GAG','AGC'],'ATG')
    #mpath = graph.shortest_path(['GTC'],'GCT')
    print(mpath)
