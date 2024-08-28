from zfnetwork import ssa, grn
import numpy as np
import unittest


class TestGillespieMethods(unittest.TestCase):

    znf_grn = grn.ZincFingerGRN()
    node_types = {'A': 'TF', 'B': 'ZF', 'C': 'TE'}
    edges = [('A', 'B'), ('A', 'C'), ('B', 'C')]
    znf_grn.from_edge_list(edges, node_types)
    
    for edge in znf_grn.edges:
        edge.k = 1.0
        edge.n = 2.0

    # TF parameters
    znf_grn['A'].pop = 1
    znf_grn['A'].beta = 10.0
    znf_grn['A'].gamma = 5.0
    
    # ZF params
    znf_grn['B'].pop = 0
    znf_grn['B'].beta = 10.0
    znf_grn['B'].gamma = 1.0
    
    # TE params
    znf_grn['C'].pop = 1
    znf_grn['C'].beta = 2.0
    znf_grn['C'].gamma = 1.0
    
    # Het params
    znf_grn['Het_1'].pop = 0
    znf_grn['Het_1'].beta = 1.0
    znf_grn['Het_1'].gamma = 1.0

    simulation = ssa.GillespieSSA(znf_grn)
   
    def test_step_tf(self):
        self.simulation.step_tf(pop=2, beta=5.0, gamma=2.5)
        self.assertEqual(self.znf_grn['A'].pop, 2)
        self.assertEqual(self.znf_grn['A'].beta, 5.0)
        self.assertEqual(self.znf_grn['A'].gamma, 2.5)

    def test_update_propensities(self):
        """Manually verify correctness of propensity values given starting params."""
        self.assertEqual(len(self.simulation.propensities), 8)

        self.assertEqual(self.simulation.propensities[0], 10.0) # TF prod 
        self.assertEqual(self.simulation.propensities[1], 5.0)  # TF deg 
        
        self.assertEqual(self.simulation.propensities[2], 5.0)  # ZF prod
        self.assertEqual(self.simulation.propensities[3], 0.0)  # ZF deg 
        
        self.assertEqual(self.simulation.propensities[4], 0.0)  # Het prod 
        self.assertEqual(self.simulation.propensities[5], 0.0)  # Het deg 
        
        self.assertEqual(self.simulation.propensities[6], 2.0)  # TE prod 
        self.assertEqual(self.simulation.propensities[7], 1.0)  # TE deg 
        
        # Update propensities and check okay after running test_step_tf
        self.simulation.update_propensities()
        self.assertEqual(len(self.simulation.propensities), 8)
        self.assertEqual(self.simulation.propensities[0], 5.0) 
        self.assertEqual(self.simulation.propensities[1], 5.0) 
    
    def test_sample_discrete(self):
        probabilities = [0.0, 1.0]
        for _ in range(1000):
            idx = self.simulation.sample_discrete(probabilities)
            if idx != 1:
                self.fail()

        probabilities = [1.0, 0.0]
        for _ in range(1000):
            idx = self.simulation.sample_discrete(probabilities)
            if idx != 0:
                self.fail()
        
        probabilities = [0.5, 0.5]
        successes = 0
        for _ in range(100000):
            successes += self.simulation.sample_discrete(probabilities)
        self.assertAlmostEqual(successes/100000, 0.5, delta=0.01)

        probabilities = [0.0, 0.0]
        with self.assertRaises(IndexError):
            self.simulation.sample_discrete(probabilities)
        
    def test_gillespie_draw(self):
        for i in range(8):
            self.simulation.propensities[i] = 0.0
        event_idx, tau = self.simulation.gillespie_draw()
        self.assertEqual(event_idx, None)
        self.assertEqual(tau, np.inf)
        
        self.simulation.update_propensities()
        event_idx, tau = self.simulation.gillespie_draw()
        self.assertNotEqual(event_idx, None)
        self.assertNotEqual(tau, np.inf)


class TestGillespieOutcomes(unittest.TestCase):
    
    def test_exponential_decay(self):
        znf_grn = grn.ZincFingerGRN()
        znf_grn.tfs.append(grn.Node(label='A', 
                                    pop=1000, 
                                    beta=0.0, 
                                    gamma=np.log(2)/100.0))
        
        znf_grn['A'].pop = 1000
        znf_grn['A'].beta = 0.0
        znf_grn['A'].gamma = np.log(2)/100.0

        simulation = ssa.GillespieSSA(znf_grn)
        tlog, plog = simulation.run(100, 100)
        
        self.assertAlmostEqual(plog.mean(axis=0)[-1], 500, delta=10)
         

    def test_steady_state(self):
        znf_grn = grn.ZincFingerGRN()
        node_types = {'A': 'TF', 'B': 'ZF'}
        edges = [('A', 'B')]
        znf_grn.from_edge_list(edges, node_types)
        
        # TF parameters
        znf_grn['A'].pop = 10
        znf_grn['A'].beta = 10.0
        znf_grn['A'].gamma = 1.0
        
        # ZF params
        znf_grn['B'].pop = 0
        znf_grn['B'].beta = 5.0
        znf_grn['B'].gamma = 2.0
        
        simulation = ssa.GillespieSSA(znf_grn)
        tlog, plog = simulation.run(25, 100)
        
        self.assertAlmostEqual(plog.mean(axis=0)[-1, 0], 10.0, delta=2)
        self.assertAlmostEqual(plog.mean(axis=0)[-1, 1], 2.5, delta=1)

if __name__ == '__main__':
    unittest.main(verbosity=3)
