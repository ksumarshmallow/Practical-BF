import numpy as np
import pandas as pd

def compute_distance_matrix(sequences):
    num_sequences = sequences.shape[0]
    D = np.zeros((num_sequences, num_sequences), dtype=int)

    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            distance = np.sum([sequences[i][pos]!=sequences[j][pos] for pos in range(len(sequences[i]))])
            D[i, j] = distance
            D[j, i] = distance

    return pd.DataFrame(D, index=sequences, columns=sequences)