from zfnetwork import grn, ssa
import unittest


class TestNode(unittest.TestCase):

    a = grn.Node('A')
    b = grn.Node('B')
    a.add_edge(b, k_xy=1.0, n=2.0)

    def test_edge_creation(self):
        self.assertEqual(len(self.a.input), 0)
        self.assertEqual(len(self.a.output), 1)
        self.assertEqual(len(self.b.input), 1)
        self.assertEqual(len(self.b.output), 0)
    
    def test_edge_parameter_assignment(self):
        self.assertEqual(self.a.output[0].k, 1.0)
        self.assertEqual(self.b.input[0].k, 1.0)
        self.a.output[0].k = 10.0
        self.assertEqual(self.a.output[0].k, 10.0)
        self.assertEqual(self.b.input[0].k, 10.0)
        self.b.input[0].k = 5.0
        self.assertEqual(self.a.output[0].k, 5.0)
        self.assertEqual(self.b.input[0].k, 5.0)

    def test_degree(self):
        self.assertEqual(self.a.degree, 1)
        self.assertEqual(self.b.degree, 1)


class TestEdge(unittest.TestCase):
    
    a = grn.Node('A', beta=1.0)
    b = grn.Node('B', beta=1.0)
    edge = grn.Edge(a, b, k_xy=10.0, n=2.0)

    def test_edge_constructor(self):
        self.assertEqual(self.a.output[0].k, 10.0)
        self.assertEqual(self.b.input[0].k, 10.0)
        self.assertEqual(self.edge.k, self.b.input[0].k)
        self.assertEqual(self.edge.k, self.a.output[0].k)

    def test_hill_activator(self):
        self.edge.k = 1.0
        self.a.mode = 'activator'
        self.a.pop = 0
        self.assertEqual(self.edge.hill(), 0.0)
        self.a.pop = 1
        self.assertEqual(self.edge.hill(), 0.5)
        self.a.pop = 10000
        self.assertAlmostEqual(self.edge.hill(), self.b.beta, delta=1e-5)

    def test_hill_repressor(self):
        self.edge.k = 1.0
        self.a.mode = 'repressor'
        self.a.pop = 0
        self.assertAlmostEqual(self.edge.hill(), self.b.beta, delta=1e-5)
        self.a.pop = 1
        self.assertEqual(self.edge.hill(), self.b.beta/2)
        self.a.pop = 10000
        self.assertAlmostEqual(self.edge.hill(), 0.0, delta=1e-5)


class TestZincFingerGRN(unittest.TestCase):

    znf_grn = grn.ZincFingerGRN()
    node_types = {1: 'TF', 2: 'ZF', 3: 'TE'}
    edges = [(1, 2), (1, 3), (2, 3), (2, 2)]
    znf_grn.from_edge_list(edges, node_types)
    
    def test_from_edge_list(self):
        
        # Test standard nodes
        self.assertEqual(len(self.znf_grn.tfs), 1)
        self.assertEqual(len(self.znf_grn.zfs), 1)
        self.assertEqual(len(self.znf_grn.tes), 1)

        # Check for creation of implicit heterochromatin node
        self.assertEqual(len(self.znf_grn.het), 2)
        self.assertEqual(len(self.znf_grn.nodes), 5) # +1 from implicit heterochromatin
        self.assertEqual(len(self.znf_grn.edges), 6) # +2 from implicit heterochromatin
        
        # Test modes
        self.assertEqual(self.znf_grn.tfs[0].mode, 'activator')
        self.assertEqual(self.znf_grn.zfs[0].mode, 'activator')
        self.assertEqual(self.znf_grn.tes[0].mode, None)
        self.assertEqual(self.znf_grn.het[0].mode, 'repressor')

    def test_nodes(self):
        self.assertEqual(len(self.znf_grn.nodes), 5)

    def test_networkx(self):
        G = self.znf_grn.to_digraph()
        self.assertEqual(len(G.nodes), 3)
        self.assertEqual(len(G.edges), 4)


if __name__ == '__main__':
    unittest.main()

