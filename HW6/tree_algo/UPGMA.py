from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
import numpy as np
import pandas as pd
from dataclasses import dataclass

from .seq_distance import compute_distance_matrix

@dataclass
class UPGMA:
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

    def run(self, save_history=True):
        D = self.matrix.copy()

        while len(D) > 2:
            min_i, min_j = self.find_minimum_idx(D)
            tax_i, tax_j = D.index[min_i], D.index[min_j]
            min_dist = D.iloc[min_i, min_j]

            # add new node
            new_cluster = f'({tax_i}, {tax_j})'
            branch_len_i = min_dist / 2 - (self.clades[tax_i].branch_length or 0)
            branch_len_j = min_dist / 2 - (self.clades[tax_j].branch_length or 0)

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

                distance_ij = (distance_i * num_elem_i + distance_j * num_elem_j) / (num_elem_i + num_elem_j)

                D.loc[new_cluster, taxa] = distance_ij
                D.loc[taxa, new_cluster] = distance_ij
            
            D.drop([tax_i, tax_j], axis=0, inplace=True)
            D.drop([tax_i, tax_j], axis=1, inplace=True)

            if save_history:
                self.matrices_history.append(D.copy())

        remaining_clusters = list(D.index)
        final_clade = Clade()
        tax_i, tax_j = remaining_clusters[0], remaining_clusters[1]
        branch_len_i = D.iloc[0, 1] / 2 - (self.clades[tax_i].branch_length or 0)
        branch_len_j = D.iloc[0, 1] / 2 - (self.clades[tax_j].branch_length or 0)

        self.clades[tax_i].branch_length = branch_len_i
        self.clades[tax_j].branch_length = branch_len_j

        final_clade.clades.append(self.clades[tax_i])
        final_clade.clades.append(self.clades[tax_j])
        final_clade.branch_length = 1

        return Tree(root=final_clade)

    def save_newick(self, tree, save_path='tree.nwk'):
        Phylo.write(tree, save_path, 'newick')
    
    def draw_tree(self, tree):
        Phylo.draw(tree, branch_labels=lambda c: round(c.branch_length, 2))
