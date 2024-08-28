#!/usr/bin/env python3

import cProfile
import numpy as np
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
        """Manually update TF parameters."""
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
            
            # Take care here - what are the effects of this?
            if node.ntype == 'TF':
                production_rate = node.beta
            # elif len(node.input) == 0:
            #     production_rate = 0.0
            else:
                # Production rate calculated as product of Hill functions. This is equivalent to AND logic
                production_rate = np.prod([edge.hill() for edge in node.input])
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
        
        # This should only trigger in the event that the population of all nodes in the network is
        # zero and beta for all TFs is also set to zero. Alternatively, if the network is
        # mis-specified, e.g. no connected components.
        if propsum == 0.0:
            tau = np.inf
            event_idx = None
        else:
            tau = np.random.exponential(scale=1.0/propsum)
            probabilities = self.propensities/propsum
            event_idx = self.sample_discrete(probabilities)
        
        return event_idx, tau


    def gillespie_ssa(self, duration, initial_pop, user_events={}):
        """Run Gillespie stochastic simulation algorithm.
        
        Arguments:
            duration: the duration of the simulation.
            events: a dict mapping from time points to events, where an event is a function that
                will be called.

        Returns:
            time_log: array of time steps from t0 to t0+duration
            pop_log: array of population records for each node
        """

        time_log = np.zeros(duration)
        pop_log = np.zeros((duration, self.n_nodes))
        pop_log[0, :] = initial_pop
    
        # Run Gillespie SSA loop
        t, t_idx = 0, 1
        while t_idx < duration:
            
            event_idx, tau = self.gillespie_draw()
            
            # See gillespie_draw() for trigger conditions.
            if event_idx is None:
                for k in range(t_idx, duration):
                    time_log[k] = time_log[k-1] + 1
                    pop_log[k] = [node.pop for node in self.zf_grn.nodes]
                return time_log, pop_log

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
                # Call user events. Seems this is not triggering when it should...
                if t_idx in user_events:
                    user_events[t_idx]()
                t_idx += 1

        return time_log, pop_log

    def run(self, duration, replicates, user_events={}):
        """Run Gillespie stochastic simulation algorithm."""

        # Initialize time and population storage arrays
        time_log = np.zeros((replicates, duration))
        pop_log = np.zeros((replicates, duration, self.n_nodes))
        statedict = self.zf_grn.save_state()

        for rep in range(replicates):
            if rep % 5 == 0:
                print(f'rep: {rep}', end='\r')
            
            # Reset node populations to original values
            self.zf_grn.load_state(statedict)
            tlog, plog = self.gillespie_ssa(duration, statedict['nodes']['pop'], user_events)
            time_log[rep, :] = tlog
            pop_log[rep, :, :] = plog
        return time_log, pop_log


def main():
    node_types = {1: 'TF', 2: 'ZF', 3: 'TE'}
    edges = [(1, 2), (1, 3), (2, 3), (2, 2)]
    zf_grn = grn.ZincFingerGRN()
    zf_grn.from_edge_list(edges, node_types)
    
    simulation = GillespieSSA(zf_grn)
    simulation.run(600, 50)

if __name__ == '__main__':
    cProfile.run('main()')
