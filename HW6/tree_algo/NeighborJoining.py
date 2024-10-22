from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
import numpy as np
import pandas as pd
from dataclasses import dataclass

from .seq_distance import compute_distance_matrix

@dataclass
class NeighborJoining:
    sequences: list = None
    matrix: pd.DataFrame = None

    def __post_init__(self):
        if (self.sequences is None and self.matrix is None) or (self.sequences is not None and self.matrix is not None):
            raise ValueError('You should pass one of those: sequences or matrix')
        
        if self.sequences is not None:
            self.num_taxon = self.sequences.shape[0]
            self.matrix = compute_distance_matrix(self.sequences)
        
        elif self.matrix is not None:
            self.num_taxon = self.matrix.shape[0]
        
        self.cluster_names = list(self.matrix.columns)
        self.matrices_history = [self.matrix]
        self.Q_history = []
        self.clades = {name: Clade(name=name) for name in self.matrix.index}
    
    def find_minimum_idx(self, D):
        min_i, min_j = -1, -1
        min_dist = 1e6

        for i in range(len(D)):
            for j in range(i + 1, len(D)):
                if D.iloc[i, j] < min_dist:
                    min_dist = D.iloc[i, j]
                    min_i, min_j = i, j
    
        return min_i, min_j
    
    def compute_Q(self, D):
        Q = np.zeros_like(D)
        n = D.shape[0]

        sum_distance = D.sum(axis=1).values.repeat(D.shape[0]).reshape(D.shape)

        Q = (n-2) * D - sum_distance - sum_distance.T

        return Q

    def compute_distance(self, D, f, g):
        sum_f, sum_g = D[f].sum(), D[g].sum()
        n = D.shape[0]

        distance_fu = D.loc[f, g] / 2 + (sum_f - sum_g) / (2 * (n - 2))
        distance_gu = D.loc[f, g] / 2 + (sum_g - sum_f) / (2 * (n - 2))

        return distance_fu, distance_gu
    
    def run(self, save_history=True):
        D = self.matrix.copy()

        while len(D) > 2:
            Q = self.compute_Q(D)

            if save_history:
                self.Q_history.append(Q)

            min_i, min_j = self.find_minimum_idx(Q)
            tax_i, tax_j = D.index[min_i], D.index[min_j]

            min_dist = D.iloc[min_i, min_j]

            # add new node
            new_cluster = f'({tax_i}, {tax_j})'
            branch_len_i, branch_len_j = self.compute_distance(D, tax_i, tax_j)

            new_clade = Clade(branch_length=min_dist / 2)
            self.clades[tax_i].branch_length = branch_len_i
            self.clades[tax_j].branch_length = branch_len_j

            new_clade.clades.append(self.clades[tax_i])
            new_clade.clades.append(self.clades[tax_j])
            self.clades[new_cluster] = new_clade

            # update distance matrix
            D.loc[:, new_cluster] = np.zeros(D.shape[0])
            D.loc[new_cluster] = np.zeros(D.shape[1])

            for taxa in D.columns[:-1]:
                num_elem_i = tax_i.count(',') + 1
                num_elem_j = tax_j.count(',') + 1

                distance_i = D.loc[taxa, tax_i]
                distance_j = D.loc[taxa, tax_j]

                distance_ij = (distance_i + distance_j - min_dist) / 2

                D.loc[new_cluster, taxa] = distance_ij
                D.loc[taxa, new_cluster] = distance_ij
            
            D.drop([tax_i, tax_j], axis=0, inplace=True)
            D.drop([tax_i, tax_j], axis=1, inplace=True)

            if save_history:
                self.matrices_history.append(D.copy())

        # Return unrooted tree
        remaining_clusters = list(D.index)
        final_clade = Clade()
        tax_i, tax_j = remaining_clusters[0], remaining_clusters[1]
        branch_len_i, branch_len_j = self.compute_distance(D, tax_i, tax_j)

        self.clades[tax_i].branch_length = branch_len_i
        self.clades[tax_j].branch_length = branch_len_j

        final_clade.clades.append(self.clades[tax_i])
        final_clade.clades.append(self.clades[tax_j])

        # Create an unrooted tree (set rooted=False)
        return Tree(root=final_clade, rooted=False)

    def save_newick(self, tree, save_path='tree.nwk'):
        Phylo.write(tree, save_path, 'newick')
    
    def draw_tree(self, tree):
        # Set rooted=False to ensure the tree is visualized as unrooted
        Phylo.draw(tree, branch_labels=lambda x: round(x.branch_length, 2) if \
                    (x.branch_length is not None and not np.isnan(x.branch_length)) else '')