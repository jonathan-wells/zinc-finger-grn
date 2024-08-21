#!/usr/bin/env python3

from os import error
import numpy as np
import matplotlib.pyplot as plt
from zfnetwork import grn


class GillespieSSA:

    def __init__(self, zf_grn):
        self.zf_grn = zf_grn
        self.n_nodes = len(zf_grn.nodes)
        self.propensities = np.zeros(2*self.n_nodes)
        self.update_propensities()
        # Generate lookup table for possible events for each node (+1, -1)
        self.events = np.zeros((self.propensities.size, self.n_nodes))
        j = 0
        for i in range(0, self.events.shape[0], 2):
            self.events[i, j] = 1
            self.events[i+1, j] = -1
            j += 1

    def step_tf(self, pop=None, beta=None, gamma=None):
        for tf in self.zf_grn.tfs:
            if pop:
                tf.pop = pop
            if beta:
                tf.beta = beta
            if gamma:
                tf.gamma = gamma

    def sample_discrete(self, probabilities):
        """Randomly sample an index with probability given by probs.""" 
        u = np.random.uniform()
        # Find index of item in probabilities which u is bounded by
        i = 0
        prob_sum = 0.0
        while prob_sum < u:
            prob_sum += probabilities[i]
            i += 1
        return i - 1

    def update_propensities(self):
        """Updates propensities of gene regulatory network."""

        for i, node in enumerate(self.zf_grn.nodes):
            idx = 2*i # Accounting for 2 propensities per node
        
            if node.ntype == 'TF':
                production_rate = node.beta
            elif len(node.input) == 0:
                production_rate = 0.0
            else:
                # Production rate calculated as product of Hill functions. This is equivalent to AND logic
                production_rate = np.array([edge.hill() for edge in node.input]).prod()
            degradation_rate = node.gamma*node.pop

            self.propensities[idx] = production_rate
            self.propensities[idx + 1] = degradation_rate

    def gillespie_draw(self):
        """Draws an event and reaction time according to propensities.

        Returns:
            event: the relevant index of the event
            tau: the time interval between this event and previous
        """
        propsum = self.propensities.sum()
        tau = np.random.exponential(scale=1.0/propsum)
        probabilities = self.propensities/propsum
        event_idx = self.sample_discrete(probabilities)
        return event_idx, tau


    def gillespie_ssa(self, duration, user_events={}):
        """Run Gillespie stochastic simulation algorithm.
        
        Arguments:
            duration: the duration of the simulation.
            events: a dict mapping from time points to events, where an event is a function that
                will be called.

        Returns:
            time_log: array of time steps from t0 to t0+duration
            pop_log: array of population records for each node
        """

        initial_state = [node.pop for node in self.zf_grn.nodes]
        time_log = np.zeros(duration)
        pop_log = np.zeros((duration, self.n_nodes))
        pop_log[0, :] = initial_state
    
        # Run Gillespie SSA loop
        t, t_idx = 0, 1
        while t_idx < duration:
            
            # Call user events
            user_events.get(t_idx, lambda: None)()
            

            event_idx, tau = self.gillespie_draw()
            t += tau
            for j, node in enumerate(self.zf_grn.nodes):
                node.pop += self.events[event_idx, j]
            self.update_propensities()
            
            if t < t_idx:
                continue

            for _ in range(int(np.ceil(t-t_idx))):
                if t_idx >= duration:
                    break
                pop_log[t_idx] = [node.pop for node in self.zf_grn.nodes]
                time_log[t_idx] = time_log[t_idx-1] + 1
                t_idx += 1

        return time_log, pop_log

    def run(self, duration, replicates, user_events={}):
        """Run Gillespie stochastic simulation algorithm."""

        # Initialize time and population storage arrays
        initial_state = [node.pop for node in self.zf_grn.nodes]

        time_log = np.zeros((replicates, duration))
        pop_log = np.zeros((replicates, duration, self.n_nodes))
        
        for rep in range(replicates):
            if rep % 5 == 0:
                print(f'rep: {rep}', end='\r')
            
            # Reset node populations to original values
            for i, node in enumerate(self.zf_grn.nodes):
                node.pop = initial_state[i]

            tlog, plog = self.gillespie_ssa(duration, user_events)
            time_log[rep, :] = tlog
            pop_log[rep, :, :] = plog
        return time_log, pop_log

    def plot_data(time_log, pop_log, zf_grn):
        fig, ax = plt.subplots(figsize=(10, 7))
        for i, node in enumerate(zf_grn.nodes):
            if node.ntype == 'TF':
                color = 'orange'
            elif node.ntype == 'ZF':
                color = 'dodgerblue'
            elif node.ntype == 'TE':
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
    node_types = {1: 'TF', 2: 'ZF', 3: 'TE'}
    edges = [(1, 2), (1, 3), (2, 3)]
    zf_grn = grn.ZincFingerGRN()
    zf_grn.from_edge_list(edges, node_types)
    for tf in zf_grn.tfs:
        tf.pop += 10
    
    simulation = GillespieSSA(zf_grn)
    t, p = simulation.run(200, 1, {100: lambda: simulation.step_tf(pop=1, beta=0.1)})

if __name__ == '__main__':
    main()
