import numpy as np
import pandas as pd
from ode_model import get_steady_state
from infection_state import computeStates


def read_parameters(n_types, k_bound, h_bound, p_intra):

    folder = 'C:/Users/mani/PycharmProjects/network/simulation/data_parameter_sets/'
    filename = 'parameter_sets_nTypes{}_k{}_h{}_pIntra{}'.format(n_types, k_bound, h_bound, p_intra).replace(".", "_") + '.csv'
    df = pd.read_csv(folder + filename, sep=',')
    df = np.array(df)[:, 1:]

    return df

def simulate_parameters(df, n_types, p_intra, c, analysis):
    # """
    # Set fixed parameters
    # :param vt, nvt: lists of vaccine types and non-vaccine types
    # :param n_types: number of types
    # :param states: list of state objects
    # :param n_states: number of states
    # :param mu: the baseline type-specific clearance rate
    # :param VEs: np.array of vaccine efficacy. VEs[0]=0 represents the pre-vaccination scenario
    # :param VE_operators: list of np.array a shape (n_types, n_types)
    #             left-multiplication with this matrix on beta gives the transmission probability with vaccine
    # :param marginal operator: matrix for left-multiplying with the state variable to get the marginal prevalence
    # :param num_simulations: number parameter combination to simulate from
    # :param simulation: simulation time per simulation
    # :param plot: if True, plot the trajectory of the marginals in time
    # :param k_matrices = list of length num_simulations
    #         each element is a symmetric np.array of shape (n_types, n_types)
    #         with value between 0 and 2
    # :param h_matrices = list of length num_simulations
    #         each element is a symmetric np.array of shape (n_types, n_types)
    #         with value between 0 and 2
    # """
    n_sets = df.shape[0]                            # number of available parameter sets
    states = computeStates(n_types)                 # infection states
    mu = np.ones(n_types)                           # baseline clearance rate
    model = 'multiplicative'                        # interaction structure
    marginal_operator = np.array([[1 if (type + 1) in s._types else 0 for s in states] for type in range(n_types)])

    # Set data frame to save output of simulated equilibrium
    n_required_success_sets = 500
    n_pairs = int(n_types * (n_types - 1) / 2)
    success_sets = np.zeros((n_required_success_sets, n_pairs * 2))

    # Simulate equilibrium for each set of parameters
    n_success_sets = 0
    for j in range(n_sets):
        # Try a new parameter set until the number of successful sets is n_required_success_sets
        if n_success_sets < n_required_success_sets:

            beta = df[j, :n_types]
            k_matrix = np.ones((n_types, n_types))
            h_matrix = np.ones((n_types, n_types))
            pair_indices = np.triu_indices(n_types, 1)
            k_matrix[pair_indices[0], pair_indices[1]] = df[j, n_types:(n_types + n_pairs)]
            k_matrix[pair_indices[1], pair_indices[0]] = df[j, n_types:(n_types + n_pairs)]
            h_matrix[pair_indices[0], pair_indices[1]] = df[j, (n_types + n_pairs):]
            h_matrix[pair_indices[1], pair_indices[0]] = df[j, (n_types + n_pairs):]

            this_eq = get_steady_state(n_types, states, k_matrix, h_matrix, c, beta, mu, marginal_operator, model)
            marginals = marginal_operator.dot(this_eq)
            this_eq[this_eq < 0] = 0
            this_eq = this_eq / sum(this_eq)

            eps = 0.001
            if min(marginal_operator.dot(this_eq)) > eps:

                success_sets[n_success_sets] = np.log(df[j, n_types:])

                # Loop through different numbers of objects
                for num_objects in [100, 1000, 10000, 100000]:
                    # Sample objects
                    realisations = np.random.multinomial(num_objects, this_eq)
                    output = np.zeros((num_objects, n_types), dtype=int)
                    for k in range(1, len(states)):
                        start = np.sum(realisations[:k])
                        end = np.sum(realisations[:(k+1)])
                        state = np.zeros(n_types)
                        state[np.array(states[k]._types) - 1] = 1
                        output[start:end] = state

                    # Write the realisations in a csv file
                    col_names = ['x_{}'.format(i + 1) for i in range(n_types)]
                    output_df = pd.DataFrame(output, columns =  col_names)
                    folder = 'C:/Users/mani/PycharmProjects/network/simulation/data_equilibria/{}/'.format(analysis)
                    filename = 'equilibrium_nTypes{}_k{}_h{}_pIntra{}_c{}_nObjects{}_set{}'.format(n_types, k_bound, h_bound, p_intra, c, num_objects, n_success_sets).replace(".", "_") + '.csv'
                    output_df.to_csv(folder + filename)

                n_success_sets += 1
                print(n_success_sets)
        else:
            print('Found ', n_success_sets, 'successful sets after ', j, 'trials.')
            break

    # Write parameter sets with successful simulation to a csv file
    indices = np.triu_indices(n_types, 1)
    logk_names = ['logk_{}_{}'.format(indices[0][i] + 1, indices[1][i] + 1) for i in range(n_pairs)]
    logh_names = ['logh_{}_{}'.format(indices[0][i] + 1, indices[1][i] + 1) for i in range(n_pairs)]
    success_sets_df = pd.DataFrame(success_sets, columns=logk_names + logh_names)
    folder = '/C:/Users/mani/PycharmProjects/network/simulation/data_parameter_sets/'
    filename = '{}_success_parameter_sets_log_nTypes{}_k{}_h{}_pIntra{}_c{}'.format(analysis, n_types, k_bound, h_bound, p_intra, c).replace(".", "_") + '.csv'
    success_sets_df.to_csv(folder + filename)


if __name__ == '__main__':
    n_types = 5
    k_bound, h_bound = 3, 1
    p_intra = 0.25
    c = 3
    analysis = 'basecase'
    df = read_parameters(n_types, k_bound, h_bound, p_intra)
    simulate_parameters(df, n_types, p_intra, c, analysis)

