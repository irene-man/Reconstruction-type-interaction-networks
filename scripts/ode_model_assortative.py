import numpy as np
from scipy.integrate import odeint
from infection_state import getState


def get_steady_state(n_types, states, k_matrix, h_matrix, c, beta, mu, marginal_operator, multi, partition=[], plot=False):
    """
    Compute the prevalence of each state at the equilibrium
    :param n_types: number of types
    :param states: a list of infection states
    :param k_matrix: n_types**2 np.array of the interaction parameters in acquisition
    :param h_matrix: n_types**2 np.array of the interaction parameters in clearance
    :param beta: the type-specific transmission probability
    :param mu: the baseline type-specific clearance rate
    :param marginal_operator: matrix for left-multiplying with the state variable to get the marginal prevalence
    :param pre_vacc: boolean whether it is a pre- or post-vaccination simulation
    :param multi: 'typewise' or 'groupwise' multiplicative interaction structure
    :param partition: partition of groups(clusters) of types
    :return: steady_state: np.array of the prevalence of states at the end of simulation (equilibrium)
              trans_matrix: np.array of the transition rate at the end of simulation(equilibrium)
    """
    p_l, p_h = 0.8, 0.2

    # initialize state
    state_0_l = p_l * np.ones(len(states)) / len(states)  # initial state
    state_0_h = p_h * np.ones(len(states)) / len(states)  # initial state
    state_0 = np.hstack((state_0_l, state_0_h))
    len(state_0)
    # compute the time invariant part of the transition matrices
    acq_matrix = compute_acquisition_rates_time_invariant_part(n_types, states, k_matrix)
    cl_matrix = compute_clearance_rates(n_types, states, h_matrix, mu)

    parameters = (n_types, states, p_l, p_h, c, beta, acq_matrix, cl_matrix, marginal_operator)
    simulation_time = 300
    t = np.linspace(0, simulation_time, 50)

    # simulate by solving the systems of ordinary differential equations numerically
    solution = odeint(solve_ode, state_0, t, args=parameters)
    steady_state = solution[-1, :]

    # return the simulated equilibrium prevalence (and transition rates)
    return steady_state


def compute_transition_matrix(state, n_types, p_l, p_h, c, beta, acq_matrix, cl_matrix, marginal_operator):
    """
    Compute the comprehensive transition matrix including all transition rates
    :param state: state variable
    :param n_types: number of types
    :param beta: the type-specific transmission probability
    :param acq_matrix: a list of acquisition matrices per type excluding the time varying part, the foi
    :param cl_matrix: a matrix for the clearance rate
    :param marginal_operator: matrix for left-multiplying with the state variable to get the marginal prevalence
    :return: dD: np.array of length len(states) of the derivatives of each state
    """
    cv = 0.8
    c_l = c * (1 - cv * np.sqrt(p_h / p_l))
    c_h = c * (1 + cv * np.sqrt(p_l / p_h))
    eps_l = 0.5  # proportion contact to be distributed random
    eps_h = 0.5  # proportion contact to be distributed random
    r_lh = eps_l * p_h * c_h / (p_l * c_l + p_h * c_h)
    r_hl = eps_h * p_l * c_l / (p_l * c_l + p_h * c_h)
    r_ll = 1 - r_lh
    r_hh = 1 - r_hl
    state_l, state_h = np.array_split(state, 2)
    foi_l = beta * c_l * (r_ll * marginal_operator.dot(state_l) / p_l + r_lh * marginal_operator.dot(state_h) / p_h)
    foi_h = beta * c_h * (r_hl * marginal_operator.dot(state_l) / p_l + r_hh * marginal_operator.dot(state_h) / p_h)

    z = np.zeros((len(state_l), len(state_l)))
    transition_matrix = np.block([[cl_matrix, z], [z, cl_matrix]])
    for i in range(n_types):
        transition_matrix += np.block([[acq_matrix[i] * foi_l[i], z], [z, acq_matrix[i] * foi_h[i]]])
    return transition_matrix


def compute_acquisition_rates_time_invariant_part(n_types, states, k_matrix):
    """
    Compute the time invariant part of the acquisition rate
    :param n_types: number of types
    :param states: a list of infection states
    :param k_matrix: np.array of the pairwise interaction parameters in acquisition
    :return: output: list of length n_types. Each item in the list is a n_state**2 np.array.
       [j,k]-element of i-th matrix: acquisition interaction parameter from state j to state k.
    """
    n_states = len(states)
    output = []
    for t in range(1, n_types + 1):
        output_i = np.zeros((n_states, n_states))
        for end_state in states:
            if t not in end_state.types:
                continue
            start_state = getState([x for x in end_state.types if x != t], states)
            d = 1
            if len(start_state.types) > 0:
                d = np.prod([k_matrix[start_t - 1, t - 1] for start_t in start_state.types])
            output_i[start_state.state_id, end_state.state_id] = d
        output.append(output_i)
    return output


def compute_clearance_rates(n_types, states, h_matrix, mu):
    """
    Compute the clearance rate
    :param n_types: number of types
    :param states: a list of state objects
    :param h_matrix: np.array of the pairwise interaction parameters in clearance
    :param mu: the baseline type-specific clearance rate
    :return: output: a n_state**2 np.array.
        The [j,k]-element the clearance rate from state j to state k including the interaction parameters.
    """
    num_states = len(states)
    output = np.zeros((num_states, num_states))
    for t in range(1, n_types + 1):
        for start_state in states:
            if t not in start_state.types:
                continue
            end_state = getState([x for x in start_state.types if x != t], states)
            d = mu[t - 1]
            if len(end_state.types) > 0:
                d = np.prod([h_matrix[end_t - 1, t - 1] for end_t in end_state.types]) * mu[t - 1]
            output[start_state.state_id, end_state.state_id] = d
    return output


def solve_ode(state, t, n_types, states, p_l, p_h, c, beta, acq_matrix, cl_matrix, marginal_operator):
    """
    Loop for solving the Ode system in which the derivative is given by dD
    :param state: state variable
    :param t:
    :param n_types: number of types
    :param states: a list of state objects
    :param beta: the type-specific transmission probability
    :param acq_matrix: a list of acquisition matrices per type excl. the time varying part, the foi
    :param cl_matrix: a matrix for the clearance rate
    :param marginal_operator: matrix for left-multiplying with the state variable to get the marginal prevalence
    :return: dD: np.array of length len(states) of the derivatives of each state
    """
    q = compute_transition_matrix(state, n_types, p_l, p_h, c, beta, acq_matrix, cl_matrix, marginal_operator)
    n_states = len(states) * 2
    transition_matrix = np.dot(np.transpose(q), state) - np.dot(q, np.ones(n_states)) * state
    return transition_matrix

