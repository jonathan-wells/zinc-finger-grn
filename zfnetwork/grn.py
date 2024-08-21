#!/usr/bin/env python3

import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
from matplotlib.patches import ArrowStyle

class Node:
    """Implementation of node in graph.

    Used to represent a single gene, transposable element or dimensionless unit of heterochromatin.
    """
    def __init__(self, 
                 label,
                 ntype=None,
                 pop=0.0, 
                 beta=1.0, 
                 gamma=0.1, 
                 mode=None,
                 input=None, 
                 output=None):
        """Node constructor

        Args:
            label: name of noded
            ntype: type of node, from 'TF', 'TE' or 'ZF'
            pop: current population of the node
            beta: maximum production rate associated with this node
            gamma: degradation rate parameter
            mode: activator or repressor
            input: list of nodes that regulate this node
            output: list of nodes regulated by this node

        Returns:
            Node instance

        """
        self.label = label
        self.ntype = ntype
        self.pop = pop
        self.beta = beta
        self.gamma = gamma
        self.mode = mode
        if input:
            self.input = input
        else:
            self.input = []
        if output:
            self.output = output
        else:
            self.output = []

    def add_edge(self, other, k_xy=1.0, n=2.0):
        """Add directed edge between self and other node"""
        edge = Edge(self, other, k_xy=k_xy, n=n)
        self.output.append(edge)
        other.input.append(edge)

    @property
    def degree(self):
        """Total number of edges connected to this node"""
        return len(self.input) + len(self.output)
    
    def __repr__(self):
        data = [f'{self.label}',
                f'\tntype: {self.ntype}',
                f'\tpop: {self.pop}',
                f'\tbeta: {self.beta}',
                f'\tgamma: {self.gamma}',
                f'\tmode: {self.mode}',
                f'\tinput: {self.input}',
                f'\toutput: {self.output}']
        return '\n'.join(data)


class Edge:
    """Edge between two Nodes."""

    def __init__(self, node_x, node_y, k_xy=1.0, n=2.0):
        """Edge constructor

        Args:
            node_x: starting node
            node_y: ending node
            k_xy: hill function activation threshold - amount of x needed to activate y
            n: hill function cooperativity coefficient
        """
        self.x = node_x
        self.y = node_y
        self.k = k_xy
        self.n = n
    
    def hill(self):
        """Return output of Hill equation for node_x acting on node_y."""
        if self.x.mode == 'activator':
            return (self.y.beta*self.x.pop**self.n)/(self.k**self.n + self.x.pop**self.n)
        elif self.x.mode == 'repressor':
            return self.y.beta/(1.0 + (self.x.pop/self.k)**self.n)
        else:
            return None

    def __repr__(self):
        return f'({self.x.label}, {self.y.label})'


class ZincFingerGRN:
    """Representation of a gene regulatory network including TFs, ZFs and TEs."""

    def __init__(self, n_tfs=0, n_zfs=0, n_tes=0):
        """Network constructor

        Args:
            n_tfs: number of transcription factors
            n_zfs: number of zinc fingers
            n_tes: number of transposable elements

        Returns:
            ZincFingerGRN instance
        """
        self.n_tfs, self.n_zfs, self.n_tes = n_tfs, n_zfs, n_tes
        self.tfs = [Node(f'TF_{i}', mode='activator') for i in range(self.n_tfs)]
        # ZF is a repressor but "activates" heterochromatin
        self.zfs = [Node(f'ZF_{i}', mode='activator') for i in range(self.n_zfs)]
        self.tes = [Node(f'TE_{i}') for i in range(self.n_tes)]
        self.het = []
        self.edges = []

    def from_edge_list(self, edge_list, node_types):
        """Constructs network from list of edges.
        
        Args:
            edge_list: list of directed edges from A -> B, represented as tuples (A, B)
            node_types: dictionary mapping each node label to its type, from 'TF', 'TE' or 'ZF'.
        """
        node_labels = set(label for pair in edge_list for label in pair)
        nodes = {label: Node(label, node_types[label])  for label in node_labels}
        for node_i_label, node_j_label in edge_list:
            node_i, node_j = nodes[node_i_label], nodes[node_j_label]
            
            for node in node_i, node_j:
                if node.ntype == 'TF' and node.label not in [n.label for n in self.tfs]:
                    node.mode = 'activator'
                    self.tfs.append(node)
                    self.n_tfs += 1
                elif node.ntype == 'ZF' and node.label not in [n.label for n in self.zfs]:
                    node.mode = 'activator'
                    self.zfs.append(node)
                    self.n_zfs += 1
                elif node.ntype == 'TE' and node.label not in [n.label for n in self.tes]:
                    self.tes.append(node)
                    self.n_tes += 1

            if node_i.ntype == 'TF':
                self.add_tf_edge(node_i, node_j)
            elif node_i.ntype == 'ZF':
                self.add_zf_edge(node_i, node_j)
            else:
                raise ValueError(f'Node_i must be either TF or ZF')

    def add_tf_edge(self, node_i, node_j):
        """Adds an edge from a TF to something else."""
        node_i.add_edge(node_j)
        self.edges.append(node_i.output[-1])
        
    def add_zf_edge(self, node_i, node_j):
        """Adds an edge from a ZF to something else, via heterochromatin unit."""
        het_count = len(self.het) + 1
        self.het.append(Node(f'Het_{het_count}', ntype='Het', beta=0.1, gamma=0.01, mode='repressor'))
        
        node_i.add_edge(self.het[-1])
        self.edges.append(node_i.output[-1])
        
        self.het[-1].add_edge(node_j, k_xy=1)
        self.edges.append(node_j.input[-1])

    def generate_erdos_renyi(self, p):
        """Generate Erdos-Renyi-like graph.

        Edges between nodeds are generated with uniform probability, with the contstraint that edges
        must follow the usual rules for TF/ZF/TE behavior. E.g. No TEs acting as TFs. 
        """
        # Reinit graph if edges already present
        self.edges = []
        self.het = []
        for node in self.nodes:
            node.input = []
            node.output = []
        
        # Generate Erdos-Renyi graph
        for node_i in self.tfs + self.zfs:
            for node_j in self.zfs + self.tes:
                u = np.random.uniform()
                if u < p:
                    if node_i.ntype == 'TF':
                        self.add_tf_edge(node_i, node_j)
                    elif node_i.ntype == 'ZF':
                        self.add_zf_edge(node_i, node_j)
    
    @property
    def nodes(self):
        return self.tfs + self.zfs + self.het + self.tes
    
    def _extract_zf_edges(self):
        edges = []
        for zf_node in self.zfs:
            for edge in zf_node.output:
                for edge in edge.y.output:
                    edges.append((zf_node.label, edge.y.label))
        return edges
    
    def _extract_tf_edges(self):
        edges = []
        for tf_node in self.tfs:
            for edge in tf_node.output:
                edges.append((edge.x.label, edge.y.label))
        return edges
    
    def to_digraph(self):
        """Converts ZFNetwork to Networkx Digraph"""
        G = nx.DiGraph()
        G.add_nodes_from([n.label for n in self.tfs + self.zfs + self.tes])
        G.add_edges_from(self._extract_tf_edges() + self._extract_zf_edges())
        return G

    def draw(self):
        """Draw graphical representation of the GRN"""
        G = nx.DiGraph()
        G.add_nodes_from([n.label for n in self.tfs + self.zfs + self.tes])
        G.add_edges_from(self._extract_tf_edges() + self._extract_zf_edges())
        layout = nx.random_layout(G)
        tf_nodes = [node.label for node in self.tfs]
        zf_nodes = [node.label for node in self.zfs]
        te_nodes = [node.label for node in self.tes]
        nx.draw_networkx_nodes(G, nodelist=tf_nodes, pos=layout, node_color='orange',
                               edgecolors='black')
        nx.draw_networkx_nodes(G, nodelist=zf_nodes, pos=layout, node_color='dodgerblue',
                               edgecolors='black')
        nx.draw_networkx_nodes(G, nodelist=te_nodes, pos=layout, node_color='grey',
                               edgecolors='black')
        nx.draw_networkx_edges(G, edgelist=self._extract_tf_edges(), pos=layout, edge_color='grey')
        nx.draw_networkx_edges(G, edgelist=self._extract_zf_edges(), pos=layout, edge_color='red',
                               arrowstyle=ArrowStyle.BarAB(widthA=0.0, widthB=0.5))
        plt.show()

    def __repr__(self):
        node_string = ''
        for node in self.tfs + self.zfs + self.het + self.tes:
            if node.degree == 0:
                node_string += f'{node.label}\n'
        edge_string = ''
        for edge in self.edges:
            edge_string += f'{edge[0]}\t{edge[1]}\n'
        return node_string + edge_string


if __name__ == '__main__':
    node_types = {1: 'TF', 2: 'ZF', 3: 'TE'}
    edges = [(1, 2), (1, 3), (2, 3)]
    
    zf_grn = ZincFingerGRN()
    zf_grn.from_edge_list(edges, node_types)
    for tf in zf_grn.tfs:
        tf.pop += 5
    for node in zf_grn.nodes:
        print(node)


