"""This module has objects for non-directed donors (NDDs) and chains initiated by NDDs.

In contrast to many kidney-exchange formulations, we do not include in the main directed
graph a vertex for each NDD. Rather, the program maintains a separate array
of Ndd objects. Each of these Ndd objects maintains a list of outgoing edges, and
each of these edges points towards a vertex (which represents a donor-patient pair)
in the directed graph.
"""

from kidney_digraph import KidneyReadException


class Ndd:
    """A non-directed donor"""

    def __init__(self, id=None):
        self.edges = []
        self.chain_weight = 0
        self.id = id

    def add_edge(self, ndd_edge):
        """Add an edge representing compatibility with a patient who appears as a
        vertex in the directed graph."""
        self.edges.append(ndd_edge)


class NddEdge:
    """An edge pointing from an NDD to a vertex in the directed graph"""

    def __init__(self, tgt, weight, fail=False, src_id=None):
        self.tgt = tgt
        self.tgt_id = tgt.id
        self.src_id = src_id
        self.weight = weight  # edge weight
        self.fail = fail

    def __str__(self):
        return "NDD edge to V{}".format(self.tgt.id)

    def __eq__(self, other):
        return (self.tgt_id == other.tgt_id) and (self.src_id == other.src_id)

    def display(self):
        return "NDD Edge: tgt=%d, weight=%f" % (self.tgt.id, self.weight,)


class Chain(object):
    """A chain initiated by an NDD.
    
    Data members:
        ndd_index: The index of the NDD
        vtx_indices: The indices of the vertices in the chain, in order
        weight: the chain's weight
        edges: ordered list of the chain's edges
    """

    def __init__(self, ndd_index, vtx_indices, weight, chain_edges):
        self.ndd_index = ndd_index
        self.vtx_indices = vtx_indices
        self.weight = weight
        self.edges = chain_edges

    @property
    def length(self):
        return len(self.vtx_indices)

    def __repr__(self):
        return (
            "Chain NDD{} ".format(self.ndd_index)
            + " ".join(str(v) for v in self.vtx_indices)
            + " with weight "
            + str(self.weight)
        )

    def display(self):
        vtx_str = " ".join(str(v) for v in self.vtx_indices)
        return "Ndd %d, vtx = %s; weight = %f" % (self.ndd_index, vtx_str, self.weight)

    def __cmp__(self, other):
        # Compare on NDD ID, then chain length, then on weight, then
        # lexicographically on vtx indices
        if self.ndd_index < other.ndd_index:
            return -1
        elif self.ndd_index > other.ndd_index:
            return 1
        elif len(self.vtx_indices) < len(other.vtx_indices):
            return -1
        elif len(self.vtx_indices) > len(other.vtx_indices):
            return 1
        elif self.weight < other.weight:
            return -1
        elif self.weight > other.weight:
            return 1
        else:
            for i, j in zip(self.vtx_indices, other.vtx_indices):
                if i < j:
                    return -1
                elif i > j:
                    return 1
        return 0

    # return chain weight with (possibly) different edge weights
    def get_weight(self, digraph, ndds, edge_success_prob):
        # chain weight is equal to e1.weight * p + e2.weight * p**2 + ... + en.weight * p**n
        ndd = ndds[self.ndd_index]
        # find the vertex that the NDD first donates to
        tgt_id = self.vtx_indices[0]
        e1 = []
        for e in ndd.edges:
            if e.tgt.id == tgt_id:
                # get edge...
                e1 = e
                break
        if e1 == []:
            print("chain.update_weight: could not find vertex id")
        weight = e1.weight * edge_success_prob  # add weight from ndd to first pair
        for j in range(len(self.vtx_indices) - 1):
            weight += digraph.adj_mat[self.vtx_indices[j]][
                self.vtx_indices[j + 1]
            ].weight * edge_success_prob ** (j + 2)
        return weight
