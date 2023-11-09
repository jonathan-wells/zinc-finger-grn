#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import grn

def sample_discrete(probabilities):
    """Randomly sample an index with probability given by probs.""" 
    u = np.random.uniform()
    # Find index of item in probabilities which u is bounded by
    i = 0
    prob_sum = 0.0
    while prob_sum < u:
        prob_sum += probabilities[i]
        i += 1
    return i - 1

def znf_te_propensity(propensities, zf_grn):
    """Updates an array of propensities given a gene regulatory network.

    Args:
        propensities: an array of propensities for possible events
        zf_grn: an instance of ZincFingerGRN
    """

    for i, node in enumerate(zf_grn.nodes):
        idx = 2*i # Accounting for 2 propensities per node
    
        # Production rate calculated as product of Hill functions. This is equivalent to AND logic
        production_rate = 1.0
        
        for edge in node.input:
            production_rate *= edge.hill()
        degradation_rate = node.gamma*node.pop

        propensities[idx] = production_rate
        propensities[idx + 1] = degradation_rate

def gillespie_draw(propensities):
    """Draws an event and reaction time according to propensities.

    Args:
        propensities:  an array of propensities for possible events

    Returns:
        event: the relevant index of the event
        tau: the time interval between this event and previous
    """
    propsum = propensities.sum()
    tau = np.random.exponential(scale=1.0/propsum)
    probabilities = propensities/propsum
    event_idx = sample_discrete(probabilities)
    return event_idx, tau

def gillespie_ssa(zf_grn, time_points, replicates):
    """Run Gillespie stochastic simulation algorithm."""

    # Initialize and update propensities array
    propensities = np.zeros(2*len(zf_grn.nodes))
    znf_te_propensity(propensities, zf_grn)
    
    # Generate lookup table for possible changes in population to each node
    events = np.zeros((propensities.size, len(zf_grn.nodes)))
    j = 0
    for i in range(0, events.shape[0], 2):
        events[i, j] = 1
        events[i+1, j] = -1
        j += 1
    
    # Initialize time and population storage arrays
    time_log = np.zeros((replicates, time_points))
    pop_log = np.zeros((replicates, time_points, len(zf_grn.nodes)))
    pop_log[:, 0, :] = np.array([node.pop for node in zf_grn.nodes]) 
    
    initial_state = pop_log[0, 0, :].copy()
    for rep in range(replicates):
        if rep % 5 == 0:
            print(f'rep: {rep}')
        for i, node in enumerate(zf_grn.nodes):
            node.pop = initial_state[i]
        # Run Gillespie SSA loop
        for t in range(1, time_points):
            event_idx, tau = gillespie_draw(propensities)
            time_log[rep, t] = time_log[rep, t-1] + tau
            for j, node in enumerate(zf_grn.nodes):
                node.pop += events[event_idx, j]
                pop_log[rep, t, j] = node.pop
            znf_te_propensity(propensities, zf_grn)

    return time_log, pop_log

def plot_data(time_log, pop_log, zf_grn):
    fig, ax = plt.subplots(figsize=(10, 7))
    for i, node in enumerate(zf_grn.nodes):
        if node.label.startswith('TF'):
            color = 'orange'
        elif node.label.startswith('ZF'):
            color = 'dodgerblue'
        elif node.label.startswith('TE'):
            color = 'grey'
        else:
            continue
            # color = 'lightgrey'
        # for rep in range(pop_log.shape[0])[:5]:
        #     plt.plot(time_log[rep, :], pop_log[rep, :, i], color=color, lw=0.3)
        plt.plot(time_log.mean(axis=0), pop_log[:, :, i].mean(axis=0), color=color)
    ax.set_xlabel('time')
    ax.set_ylabel('counts')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.show()

def main():
    zf_grn = grn.ZincFingerGRN(5, 15, 15)
    zf_grn.generate_erdos_renyi(0.1)
    zf_grn.nodes[0].pop = 1
    tlog, plog = gillespie_ssa(zf_grn, 1500, 100)
    plot_data(tlog, plog, zf_grn)

if __name__ == '__main__':
    main()
