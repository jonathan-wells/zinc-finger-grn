#!/usr/bin/env python3

import grn
import unittest


class TestNode(unittest.TestCase):

    def test_edge_creation(self):
        a = grn.Node('A')
        b = grn.Node('B')
        a.add_edge(b)
        self.assertEqual(len(a.input), 0)
        self.assertEqual(len(a.output), 1)
        self.assertEqual(len(b.input), 1)
        self.assertEqual(len(b.output), 0)

    def test_degree(self):
        a = grn.Node('A')
        b = grn.Node('B')
        a.add_edge(b)
        self.assertEqual(a.degree, 1)
        self.assertEqual(b.degree, 1)

if __name__ == '__main__':
    unittest.main()

