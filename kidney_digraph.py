"""A Digraph class which can be used for representing donor-patient pairs
(as vertices) and their compatibilities (as weighted edges), along with
some related methods.
Modified from https://github.com/jamestrimble/kidney_solver
"""

from collections import deque

import numpy as np


class KidneyReadException(Exception):
    pass


class Vertex:
    """A vertex in a directed graph (see the Digraph class)."""

    def __init__(self, id):
        self.id = id
        self.edges = []

    def __str__(self):
        return "V{}".format(self.id)


class Cycle:
    """A cycle in a directed graph.

    Contains:
    - list of vertices, in order
    - list of edges
    - cycle weight
    """

    def __init__(self, vs):
        self.vs = vs
        self.weight = 0
        self.length = len(vs)
        self.edges = []

    def __cmp__(self, other):
        if min(self.vs) < min(other.vs):
            return -1
        elif min(self.vs) > min(other.vs):
            return 1
        else:
            return 1

    def __len__(self):
        return self.length

    def display(self):
        vtx_str = " ".join(str(v) for v in self.vs)
        return "L=%2d; %15s; weight = %f" % (self.length, vtx_str, self.weight,)

    def contains_edge(self, e):
        if e.src in self.vs:
            i = self.vs.index(e.src)
            if e.tgt == self.vs[(i + 1) % self.length]:
                return True
            else:
                return False
        return False

    def add_edges(self, es):
        # create an unordered list of edges in the cycle
        self.edges = [e for e in es if self.contains_edge(e)]


class Edge:
    """An edge in a directed graph (see the Digraph class)."""

    def __init__(self, id, weight, src, tgt, fail=False):
        self.id = id
        self.weight = weight  # edge weight
        self.fail = fail  # whether edge failed
        self.src = src  # source vertex
        self.src_id = src.id
        self.tgt = tgt  # target vertex
        self.tgt_id = tgt.id

    def __str__(self):
        return "V" + str(self.src.id) + "-V" + str(self.tgt.id)

    def __eq__(self, other):
        return (self.tgt_id == other.tgt_id) and (self.src_id == other.src_id)

    def to_dict(self):
        e_dict = {
            "type": "pair_edge",
            "id": self.id,
            "weight": self.weight,
            "src_id": self.src_id,
            "tgt_id": self.tgt_id,
        }
        return e_dict

    @classmethod
    def from_dict(cls, e_dict, digraph):
        return digraph.adj_mat[e_dict["src_id"]][e_dict["tgt_id"]]

    def display(self,):
        return "src=%d, tgt=%d, weight=%f" % (self.src.id, self.tgt.id, self.weight,)


class Digraph:
    """A directed graph, in which each edge has a numeric weight.

    Data members:
        n: the number of vertices in the digraph
        vs: an array of Vertex objects, such that vs[i].id == i
        es: an array of Edge objects, such that es[i].id = i
    """

    def __init__(self, n, vs_ids=None, adj_mat=None):
        """Create a Digraph with n vertices"""
        self.n = n

        # add vertices with specified IDs
        if vs_ids is not None:
            # check length of vs
            if len(vs_ids) != n:
                raise Warning("length of vs_ids is %d but n is %d" % (len(vs_ids), n))
            else:
                self.vs = [Vertex(id) for id in vs_ids]
        else:
            self.vs = [Vertex(i) for i in range(n)]

        # add an existing adj_mat
        if adj_mat is not None:
            # check dimensions
            if np.array(adj_mat).shape != (n, n):
                raise Warning(
                    "adj_mat is not the correct size for n=%d: %s"
                    % (n, str(np.array(adj_mat).shape))
                )
            # check diag (no self-edges)
            for i in range(n):
                if adj_mat[i][i] is not None:
                    raise Warning(
                        f"adj_mat has a diagonal entry at ({i},{i}): {str(adj_mat[i][i])}"
                    )
        else:
            self.adj_mat = [[None for x in range(n)] for x in range(n)]

        self.es = []
        self.cycles = []

    def add_edge(self, weight, source, tgt):
        """Add an edge to the digraph

        Args:
            weight: the edge's weight, as a float
            source: the source Vertex
            tgt: the edge's target Vertex
        """

        id = len(self.es)
        e = Edge(id, weight, source, tgt)
        self.es.append(e)
        source.edges.append(e)
        self.adj_mat[source.id][tgt.id] = e

    def remove_edge(self, e):
        self.es.remove(e)
        e.src.edges.remove(e)
        self.adj_mat[e.src.id][e.tgt.id] = None

    def find_cycles(self, max_length):
        """Find cycles of length up to max_length in the digraph.

        Returns:
            a list of cycles. Each cycle is represented as a list of
            vertices, with the first vertex _not_ repeated at the end.
        """
        cycle_list = [cycle for cycle in self.generate_cycles(max_length)]

        return cycle_list

    def generate_cycles(self, max_length):
        """Generate cycles of length up to max_length in the digraph.

        Each cycle yielded by this generator is represented as a list of
        vertices, with the first vertex _not_ repeated at the end.
        """

        vtx_used = [False] * len(
            self.vs
        )  # vtx_used[i]==True iff vertex i is in current path

        def cycle(current_path):
            last_vtx = current_path[-1]
            if self.edge_exists(last_vtx, current_path[0]):
                yield current_path[:]
            if len(current_path) < max_length:
                for e in last_vtx.edges:
                    v = e.tgt
                    if (
                        len(current_path) + shortest_paths_to_low_vtx[v.id]
                        <= max_length
                        and not vtx_used[v.id]
                    ):
                        current_path.append(v)
                        vtx_used[v.id] = True
                        for c in cycle(current_path):
                            yield c
                        vtx_used[v.id] = False
                        del current_path[-1]

        # Adjacency lists for transpose graph
        transp_adj_lists = [[] for v in self.vs]
        for edge in self.es:
            transp_adj_lists[edge.tgt.id].append(edge.src)

        for v in self.vs:
            shortest_paths_to_low_vtx = self.calculate_shortest_path_lengths(
                v,
                max_length - 1,
                lambda u: (w for w in transp_adj_lists[u.id] if w.id > v.id),
            )
            vtx_used[v.id] = True
            for c in cycle([v]):
                yield c
            vtx_used[v.id] = False

    def calculate_shortest_path_lengths(
        self, from_v, max_dist, adj_list_accessor=lambda v: (e.tgt for e in v.edges)
    ):
        """Calculate the length of the shortest path from vertex from_v to each
        vertex with a greater or equal index, using paths containing
        only vertices indexed greater than or equal to from_v.

        Return value: a list of distances of length equal to the number of vertices.
        If the shortest path to a vertex is greater than max_dist, the list element
        will be 999999999.

        Args:
            from_v: The starting vertex
            max_dist: The maximum distance we're interested in
            adj_list_accessor: A function taking a vertex and returning an
                iterable of out-edge targets
        """
        # Breadth-first search
        q = deque([from_v])
        distances = [999999999] * len(self.vs)
        distances[from_v.id] = 0

        while q:
            v = q.popleft()
            # Note: >= is used instead of == on next line in case max_dist<0
            if distances[v.id] >= max_dist:
                break
            for w in adj_list_accessor(v):
                if distances[w.id] == 999999999:
                    distances[w.id] = distances[v.id] + 1
                    q.append(w)

        return distances

    def edge_exists(self, v1, v2):
        """Returns true if and only if an edge exists from Vertex v1 to Vertex v2."""

        return self.adj_mat[v1.id][v2.id] is not None


def failure_aware_cycle_weight(cycle, digraph, edge_success_prob):
    """Calculate a cycle's total weight, with edge failures and no backarc recourse.

    Args:
        cycle: A list of Vertex objects in the cycle, with the first Vertex not repeated.
            UPDATE: cycle = a Cycle object
        digraph: The digraph in which this cycle appears.
        edge_success_prob: The problem that any given edge will NOT fail
    """

    return sum(
        digraph.adj_mat[cycle[i - 1].id][cycle[i].id].weight for i in range(len(cycle))
    ) * edge_success_prob ** len(cycle)


def cycle_weight(cycle, digraph):
    """Calculate the sum of a cycle's edge weights.

    Args:
        cycle: A list of Vertex objects in the cycle, with the first Vertex not repeated.
        digraph: The digraph in which this cycle appears.
    """

    return sum(
        digraph.adj_mat[cycle[i - 1].id][cycle[i].id].weight for i in range(len(cycle))
    )
