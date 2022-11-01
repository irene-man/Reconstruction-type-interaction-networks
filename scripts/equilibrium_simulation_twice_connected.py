import numpy as np
import pandas as pd
from ode_model import get_steady_state
from infection_state import computeStates
import time



def read_parameters(n_types, k_bound, h_bound, p_intra):

    folder = '/home/iman/Research/Graphical/simulation/data_parameter_sets/'
    filename = 'parameter_sets_nTypes{}_k{}_h{}_pIntra{}'.format(n_types, k_bound, h_bound, p_intra).replace(".", "_") + '.csv'
    df = pd.read_csv(folder + filename, sep=',')
    df = np.array(df)[:, 1:]

    return df

def simulate_parameters(df, n_types_sub, p_intra, p_inter, c, analysis):
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
    model = 'multiplicative'                        # interaction structure
    n_types_full = 2 * n_types_sub
    states_sub = computeStates(n_types_sub)                 # infection states
    states_full = computeStates(n_types_full)
    mu_full = np.ones(n_types_full)
    mu_sub = np.ones(n_types_sub)
    marginal_operator_sub = np.array([[1 if (type + 1) in s._types else 0 for s in states_sub] for type in range(n_types_sub)])
    marginal_operator_full = np.array([[1 if (type + 1) in s._types else 0 for s in states_full] for type in range(n_types_full)])

    # Set data frame to save output of simulated equilibrium
    n_required_success_sets = 100
    n_pairs_sub = int(n_types_sub * (n_types_sub - 1) / 2)
    n_pairs_full = int(n_types_full * (n_types_full - 1) / 2)
    success_sets = np.zeros((n_required_success_sets, n_pairs_full * 2))

    # Simulate equilibrium for each set of parameters
    n_success_sets = 0
    for j in range(int(n_sets / 2)):
        # Try a new parameter set until the number of successful sets is n_required_success_sets
        print('Trial ', j)
        if n_success_sets < n_required_success_sets:
            i = int(j + n_sets / 2)
            beta_A = df[j, :n_types_sub]
            beta_B = df[i, :n_types_sub]
            k_matrix_A = np.ones((n_types_sub, n_types_sub))
            k_matrix_B = np.ones((n_types_sub, n_types_sub))
            h_matrix_A = np.ones((n_types_sub, n_types_sub))
            h_matrix_B = np.ones((n_types_sub, n_types_sub))
            pair_indices = np.triu_indices(n_types_sub, 1)
            k_matrix_A[pair_indices[0], pair_indices[1]] = df[j, n_types_sub:(n_types_sub + n_pairs_sub)]
            k_matrix_A[pair_indices[1], pair_indices[0]] = df[j, n_types_sub:(n_types_sub + n_pairs_sub)]
            k_matrix_B[pair_indices[0], pair_indices[1]] = df[i, n_types_sub:(n_types_sub + n_pairs_sub)]
            k_matrix_B[pair_indices[1], pair_indices[0]] = df[i, n_types_sub:(n_types_sub + n_pairs_sub)]
            h_matrix_A[pair_indices[0], pair_indices[1]] = df[j, (n_types_sub + n_pairs_sub):]
            h_matrix_A[pair_indices[1], pair_indices[0]] = df[j, (n_types_sub + n_pairs_sub):]
            h_matrix_B[pair_indices[0], pair_indices[1]] = df[i, (n_types_sub + n_pairs_sub):]
            h_matrix_B[pair_indices[1], pair_indices[0]] = df[i, (n_types_sub + n_pairs_sub):]

            this_eq_A = get_steady_state(n_types_sub, states_sub, k_matrix_A, h_matrix_A, c, beta_A, mu_sub, marginal_operator_sub, model)
            this_eq_B = get_steady_state(n_types_sub, states_sub, k_matrix_B, h_matrix_B, c, beta_B, mu_sub, marginal_operator_sub, model)

            this_eq_A[this_eq_A < 0] = 0
            this_eq_B[this_eq_B < 0] = 0
            this_eq_A = this_eq_A / sum(this_eq_A)
            this_eq_B = this_eq_B / sum(this_eq_B)

            eps = 0.001

            # If a subgroup does not allow coexistence, continue with a new pair of parameter sets
            if min(marginal_operator_sub.dot(this_eq_A)) < eps or min(marginal_operator_sub.dot(this_eq_B)) < eps:
                continue

            # If both subgroups allow coexistence, try connecting them
            beta_full = np.hstack((beta_A, beta_B))
            k_matrix_full = np.ones((n_types_full, n_types_full))
            h_matrix_full = np.ones((n_types_full, n_types_full))
            k_matrix_full[:n_types_sub, :n_types_sub] = k_matrix_A
            k_matrix_full[n_types_sub:, n_types_sub:] = k_matrix_B
            h_matrix_full[:n_types_sub, :n_types_sub] = h_matrix_A
            h_matrix_full[n_types_sub:, n_types_sub:] = h_matrix_B

            presence = np.random.choice([0, 1], (n_types_sub, n_types_sub), replace=True, p=[1 - p_inter, p_inter])
            k_AB = presence * k_bound ** np.random.uniform(-1, 1, (n_types_sub, n_types_sub)) + (1 - presence)
            k_matrix_full[:n_types_sub, n_types_sub:] = k_AB
            k_matrix_full[n_types_sub:, :n_types_sub] = k_AB.transpose()

            print('Started solving ode...')
            start = time.time()
            this_eq = get_steady_state(n_types_full, states_full, k_matrix_full, h_matrix_full, c, beta_full, mu_full, marginal_operator_full, model)
            print(int(time.time() - start), 'sec')

            # Continue if the entire system allows coexistence
            if min(marginal_operator_full.dot(this_eq)) > eps :

                success_sets[n_success_sets] = np.log(np.hstack((k_matrix_full[np.triu_indices(n_types_full, 1)], h_matrix_full[np.triu_indices(n_types_full, 1)])))

                # Loop through different numbers of objects
                for num_objects in [100, 1000, 10000, 100000]:
                    # Sample objects
                    realisations = np.random.multinomial(num_objects, this_eq)
                    output = np.zeros((num_objects, n_types_full), dtype=int)
                    for k in range(1, len(states_full)):
                        start = np.sum(realisations[:k])
                        end = np.sum(realisations[:(k + 1)])
                        state = np.zeros(n_types_full)
                        state[np.array(states_full[k]._types) - 1] = 1
                        output[start:end] = state

                    # Write the realisations in a csv file
                    col_names = ['x_{}'.format(i + 1) for i in range(n_types_full)]
                    output_df = pd.DataFrame(output, columns=col_names)
                    folder = '/home/iman/Research/Graphical/simulation/data_equilibria/{}/'.format(analysis)
                    filename = '{}_equilibrium_nTypes{}_k{}_h{}_pIntra{}_pInter{}_c{}_nObjects{}_set{}'.format(analysis, n_types_full, k_bound,
                                                                                                   h_bound, p_intra, p_inter, c,
                                                                                                   num_objects,
                                                                                                   n_success_sets).replace(
                        ".", "_") + '.csv'
                    output_df.to_csv(folder + filename)

                n_success_sets += 1
                print(n_success_sets)
        else:
            print('Found ', n_success_sets, 'successful sets after ', j, 'trials.')
            break

    # Write parameter sets with successful simulation to a csv file
    indices = np.triu_indices(n_types_full, 1)
    logk_names = ['logk_{}_{}'.format(indices[0][i] + 1, indices[1][i] + 1) for i in range(n_pairs_full)]
    logh_names = ['logh_{}_{}'.format(indices[0][i] + 1, indices[1][i] + 1) for i in range(n_pairs_full)]
    success_sets_df = pd.DataFrame(success_sets, columns=logk_names + logh_names)
    folder = '/home/iman/Research/Graphical/simulation/data_parameter_sets/'
    filename = '{}_success_parameter_sets_log_nTypes{}_k{}_h{}_pIntra{}_pInter{}_c{}'.format(analysis, n_types_sub * 2, k_bound, h_bound, p_intra, p_inter, c).replace(".", "_") + '.csv'
    success_sets_df.to_csv(folder + filename)


if __name__ == '__main__':
    n_types_sub = 5
    k_bound, h_bound = 3, 1
    p_intra = 0.25
    p_inter = 0.1
    c = 3
    analysis = 'sens6'
    df = read_parameters(n_types_sub, k_bound, h_bound, p_intra)
    simulate_parameters(df, n_types_sub, p_intra, p_inter, c, analysis)

