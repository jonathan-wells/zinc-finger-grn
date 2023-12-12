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
                 k_xy=1.0, 
                 n=2.0, 
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
            n: hill function cooperativity coefficient
            mode: activator or repressor
            input: list of nodes that regulate this node
            output: list of nodes regulated by this node

        Returns:
            Node instance

        """
        self.label = label
        if ntype == None:
            self.ntype = self.label.split('_')[0]
        else:
            self.ntype = ntype
        self.pop = pop
        self.beta = beta
        self.gamma = gamma
        self.k_xy = k_xy
        self.n = n
        self.mode = mode
        if input is None:
            self.input = []
        else:
            self.input = input
        if output is None:
            self.output = []
        else:
            self.output = output

    def add_edge(self, other, k_xy=1.0):
        """Add directed edge between self and other node"""
        self.output.append(Edge(self, other, k_xy=k_xy))
        other.input.append(Edge(self, other, k_xy=k_xy))

    @property
    def degree(self):
        """Total number of edges connected to this node"""
        return len(self.input) + len(self.output)
    
    def __repr__(self):
        return self.label


class Edge:
    """Edge between two Nodes."""

    def __init__(self, node_x, node_y, k_xy=1.0):
        """Edge constructor

        Args:
            node_x: starting node
            node_y: ending node
            k_xy: hill function activation threshold - amount of x needed to activate y

        """
        self.x = node_x
        self.y = node_y
        self.k = k_xy
    
    def hill(self):
        """Return output of Hill equation for node_x acting on node_y"""
        if self.x.mode == 'activator':
            return (self.y.beta*self.x.pop**self.x.n)/(self.k**self.x.n + self.x.pop**self.x.n)
        elif self.x.mode == 'repressor':
            return self.y.beta/(1 + (self.x.pop/self.k)**self.x.n)
        else:
            return None

    def __repr__(self):
        return f'{self.x}\t{self.y}'


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
        self.zfs = [Node(f'ZF_{i}', mode='activator') for i in range(self.n_zfs)]
        self.tes = [Node(f'TE_{i}') for i in range(self.n_tes)]
        self.het = []
        self.edges = []

    def from_edge_list(edge_list, node_types):
        """Constructs network from list of edges.
        
        Args:
            edge_list: list of directed edges from A -> B, represented as tuples (A, B)
            node_types: dictionary mapping each node label to its type, from 'TF', 'TE' or 'ZF'.
        """
        for node_i_label, node_j_label in edge_list:
            nit, njt = node_types[node_i_label], node_types[node_j_label]
            node_i = Node(node_i_label, nit)
            node_j = Node(node_j_label, njt)
            if nit == 'TF':
                self.add_tf_edge(node_i, node_j)
            elif nit == 'ZF':
                self.add_zf_edge(node_i, node_j)
            else:
                raise ValueError(f'Parent node must be type TF or ZF. {node_i_label} is type {nit}')
         
    def add_tf_edge(self, node_i, node_j):
        """Adds an edge from a TF to something else."""
        node_i.add_edge(node_j)
        self.edges.append(Edge(node_i, node_j))
        
    def add_zf_edge(self, node_i, node_j):
        """Adds an edge from a ZF to something else, via heterochromatin unit."""
        het_count = len(self.het)
        self.het.append(Node(f'Het_{het_count}', beta=0.1, gamma=0.01, mode='repressor'))
        node_i.add_edge(self.het[-1])
        self.het[-1].add_edge(node_j, k_xy=1)
        self.edges.append(Edge(node_i, self.het[-1]))
        self.edges.append(Edge(self.het[-1], node_j, k_xy=1))

    def generate_erdos_renyi(self, p):
        for node_i in self.tfs + self.zfs:
            for node_j in self.zfs + self.tes:
                u = np.random.uniform()
                if u < p:
                    if node_i.ntype == 'TF':
                        self.add_tf_edge(node_i, node_j)
                    elif node_i.ntype == 'ZF':
                        self.add_zf_edge(node_i, node_j)


    def generate_random_edges(self):
        """Randomly connect edges between TFs, ZFs and tes
        
        This should be done in a biologically plausible manner, so as to realistically simulate
        different potential GRN architectures.
        """
        # Generate edges from TFs to ZFs and TEs
        for node_x in self.tfs:
            p = 0.75 # Need a better way of doing this later
            for node_y in self.zfs:
                u = np.random.uniform(0.0, 1.0)
                if u <= p:
                    node_x.add_edge(node_y)
                    self.edges.append(Edge(node_x, node_y))
            for node_z in self.tes:
                u = np.random.uniform(0.0, 1.0)
                if u <= p:
                    node_x.add_edge(node_z)
                    self.edges.append(Edge(node_x, node_z))

        het_count = 0
        for node_x in self.zfs:
            p_zf = 1/self.n_zfs
            for node_y in self.zfs:
                u = np.random.uniform(0.0, 1.0)
                if u <= p_zf:
                    self.het.append(Node(f'Het_{het_count}', beta=0.1, gamma=0.01, mode='repressor'))
                    node_x.add_edge(self.het[-1])
                    self.het[-1].add_edge(node_y, k_xy=1)
                    self.edges.append(Edge(node_x, self.het[-1]))
                    self.edges.append(Edge(self.het[-1], node_y, k_xy=1))
                    het_count += 1
            p_te = 2/self.n_tes # Need a better way of doing this later
            for node_z in self.tes:
                u = np.random.uniform(0.0, 1.0)
                if u <= p_te:
                    self.het.append(Node(f'Het_{het_count}', beta=0.1, gamma=0.01, mode='repressor'))
                    node_x.add_edge(self.het[-1])
                    self.het[-1].add_edge(node_z, k_xy=1)
                    self.edges.append(Edge(node_x, self.het[-1]))
                    self.edges.append(Edge(self.het[-1], node_z, k_xy=1))
                    het_count += 1
    
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

    def draw(self):
        """Draw graphical represation of the GRN"""
        G = nx.DiGraph()
        G.add_nodes_from([n.label for n in self.tfs + self.zfs + self.tes])
        G.add_edges_from(self._extract_tf_edges() + self._extract_zf_edges())
        layout = nx.random_layout(G)
        tf_nodes = [n for n in G.nodes if n.startswith('TF')]
        zf_nodes = [n for n in G.nodes if n.startswith('ZF')]
        te_nodes = [n for n in G.nodes if n.startswith('TE')]
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
    grn = ZincFingerGRN(5, 15, 15)
    grn.generate_erdos_renyi(0.33)
    grn.draw()


