import numpy as np
from pre_train import data_load


def load_data(data_name):
    gene_expression = data_load(data_name)
    Label = gene_expression
    # Img = np.reshape(Img, (Img.shape[0], 32, 32, 1))
    Img = np.reshape(gene_expression, [gene_expression.shape[0], 1, gene_expression.shape[1], 1])
    n_input = [1, gene_expression.shape[1]]

    return gene_expression, Img, Label, n_input

