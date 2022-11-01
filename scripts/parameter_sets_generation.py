import numpy as np
import pandas as pd


def generate_parameters(n_sets, n_types, p_intra, k_bound, h_bound, beta_min, beta_max):
    """
    Randomly generate n_sets of parameter sets.
    Each parameter set includes:
     - n_types uniformly distribute transmission probabilities in the interval (beta_min, beta_max)
     - n_types**2 interaction parameters for acquisition uniformly distributed in the log-scale
     - n_types**2 interaction parameters for clearance uniformly distributed in the log-scale
    :param n_sets: number of parameter sets to generate
    :param n_types: number of types
    :param p_intra: probability of having non-trivial interaction between a pair of types
    :param k_bound: the bound of strength of interaction between types
    :param h_bound: the bound of strength of interaction between types
    :param beta_min: lower bound of transmission probabilities
    :param beta_max: upper bound of transmission probabilities
    :return: np.array of dimension (n_sets, n_types + 2 * n_types**2)
    """

    n_pairs = int(n_types * (n_types - 1) / 2)
    output = np.ones((n_sets, int(n_types + 2 * n_pairs)))

    # Generate transmission probabilities
    output[:, :n_types] = np.random.uniform(beta_min, beta_max, (n_sets, n_types))

    # Generate interaction parameters
    for i in range(n_sets):
        # Generate interaction parameters in acquisition
        k_matrix = k_bound ** np.random.uniform(-1, 1, n_pairs)
        h_matrix = np.ones(n_pairs)

        # Generate interaction parameters in clearance
        if h_bound != 1:
            x_matrix = (k_bound * 2) ** np.random.uniform(-1, 1, n_pairs)
            h_matrix = k_matrix / x_matrix

        # Keep only proportion p of the k_matrix and h_matrix
        drop_pairs = np.random.choice(n_pairs, np.random.binomial(n_pairs, 1 - p_intra), replace=False)
        k_matrix[drop_pairs] = 1
        h_matrix[drop_pairs] = 1

        # Reshape the matrices in one-dimensional arrays
        output[i, n_types:(n_types + n_pairs)] = k_matrix
        output[i, (n_types + n_pairs):] = h_matrix

    return output


if __name__ == '__main__':
    n_types = 5

    p_intra = 0.25
    # p_intra = 0.5

    # k_bound = 10
    k_bound = 3
    # k_bound = 3/2

    h_bound = 1
    # h_bound = 9/2

    beta_min, beta_max = 1/2, 1
    n_sets = 1000

    beta_names = ['beta_{}'.format(i + 1) for i in range(n_types)]
    n_pairs = int(n_types * (n_types - 1) / 2)
    indices = np.triu_indices(5, 1)
    k_names = ['k_{}_{}'.format(indices[0][i] + 1, indices[1][i] + 1) for i in range(n_pairs)]
    h_names = ['h_{}_{}'.format(indices[0][i] + 1, indices[1][i] + 1) for i in range(n_pairs)]
    df = pd.DataFrame(generate_parameters(n_sets, n_types, p_intra, k_bound, h_bound, beta_min, beta_max), columns = beta_names + k_names + h_names)

    # write the generated parameter set in a csv file
    folder = 'C:/Users/mani/PycharmProjects/network/simulation/data_parameter_sets/'
    filename = 'parameter_sets_nTypes{}_k{}_h{}_pIntra{}'.format(n_types, k_bound, h_bound, p_intra).replace(".", "_") + '.csv'
    print(folder + filename)
    df.to_csv(folder + filename)
