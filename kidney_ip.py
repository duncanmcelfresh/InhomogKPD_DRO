"""Solving the kidney exchange problem with PICEF, PC-TSP, and PI-TSP"""

import kidney_utils
from gurobi_functions import optimize, create_mip_model
from kidney_digraph import Cycle, failure_aware_cycle_weight, cycle_weight
from gurobipy import *

###################################################################################################
#                                                                                                 #
#                                  Code used by all formulations                                  #
#                                                                                                 #
###################################################################################################


class OptConfig(object):
    """The inputs (problem instance and parameters) for an optimisation run

    Data members:
        digraph
        ndds
        max_cycle
        max_chain
        verbose: True if and only if Gurobi output should be written to screen and log file
        timelimit
        edge_success_prob
        lp_file: The name of a .lp file to write, or None if the file should not be written
        relax: True if and only if the LP relaxation should be solved also
    """

    def __init__(
        self,
        digraph,
        ndds,
        max_cycle,
        max_chain,
        verbose=False,
        timelimit=None,
        edge_success_prob=1,
        lp_file=None,
        relax=False,
        cardinality_restriction=None,
        protection_level=0.1,
        chain_restriction=None,
        cycle_restriction=None,
        name=None,
    ):
        self.digraph = digraph
        self.ndds = ndds
        self.max_cycle = max_cycle
        self.max_chain = max_chain
        self.verbose = verbose
        self.timelimit = timelimit
        self.edge_success_prob = edge_success_prob
        self.edge_failure_prob = 1.0 - self.edge_success_prob
        self.lp_file = lp_file
        self.relax = relax
        self.cardinality_restriction = cardinality_restriction
        self.chain_restriction = chain_restriction
        self.cycle_restriction = cycle_restriction
        self.protection_level = protection_level  # for variable budget uncertainty: probability that U-set does not contain true edge weights
        self.name = name

        # are chains used?
        self.use_chains = (
            True  # (self.max_chain > 0) and (self.chain_restriction < len(self.ndds))
        )


class OptSolution(object):
    """An optimal solution for a kidney-exchange problem instance.

    Data members:
        ip_model: The Gurobi Model object
        cycles: A list of cycles in the optimal solution, each represented
            as a list of vertices
        chains: A list of chains in the optimal solution, each represented
            as a Chain object
        total_weight: The total weight of the solution
    """

    def __init__(
        self,
        ip_model,
        cycles,
        chains,
        digraph,
        edge_success_prob=1,
        infeasible=False,
        robust_weight=0,
        optimistic_weight=0,
        cycle_obj=None,
        chain_restriction=None,
        cycle_restriction=None,
        cardinality_restriction=None,
        cycle_cap=None,
        chain_cap=None,
        matching_edges=None,
        alpha_var=None,
    ):
        self.ip_model = ip_model
        self.cycles = cycles
        self.chains = chains
        self.digraph = digraph
        self.infeasible = infeasible
        if self.infeasible:
            self.total_weight = 0
        else:
            self.total_weight = sum(c.weight for c in chains) + sum(
                failure_aware_cycle_weight(c, digraph, edge_success_prob)
                for c in cycles
            )
        self.edge_success_prob = edge_success_prob
        self.cycle_obj = cycle_obj
        self.matching_edges = matching_edges
        self.robust_weight = robust_weight
        self.optimistic_weight = optimistic_weight
        self.cycle_restriction = cycle_restriction
        self.chain_restriction = chain_restriction
        self.cardinality_restriction = cardinality_restriction
        self.cycle_cap = cycle_cap
        self.chain_cap = chain_cap
        self.alpha_var = alpha_var

        if ip_model is not None:
            self.timeout = ip_model.status == GRB.TIME_LIMIT
        else:
            self.timeout = False

    def same_matching_edges(self, other):
        if len(self.matching_edges) != len(other.matching_edges):
            return False
        for self_e in self.matching_edges:
            edge_found = False
            for other_e in other.matching_edges:
                if (self_e.src_id == other_e.src_id) and (
                    self_e.tgt.id == other_e.tgt.id
                ):
                    edge_found = True
                    break
            if not edge_found:
                return False
        return True

    def add_matching_edges(self, ndds):
        """Set attribute 'matching_edges' using self.cycle_obj, self.chains, and self.digraph"""

        matching_edges = []

        for ch in self.chains:
            chain_edges = []
            tgt_id = ch.vtx_indices[0]
            for e in ndds[ch.ndd_index].edges:
                if e.tgt.id == tgt_id:
                    chain_edges.append(e)
            if len(chain_edges) == 0:
                raise Warning("NDD edge not found")
            for i in range(len(ch.vtx_indices) - 1):
                chain_edges.append(
                    self.digraph.adj_mat[ch.vtx_indices[i]][ch.vtx_indices[i + 1]]
                )
            if len(chain_edges) != (len(ch.vtx_indices)):
                raise Warning(
                    "Chain contains %d edges, but only %d edges found"
                    % (len(ch.vtx_indices), len(chain_edges))
                )
            matching_edges.extend(chain_edges)

        for cy in self.cycle_obj:
            cycle_edges = []
            for i in range(len(cy.vs) - 1):
                cycle_edges.append(self.digraph.adj_mat[cy.vs[i].id][cy.vs[i + 1].id])
            # add final edge
            cycle_edges.append(self.digraph.adj_mat[cy.vs[-1].id][cy.vs[0].id])
            if len(cycle_edges) != len(cy.vs):
                raise Warning(
                    "Cycle contains %d vertices, but only %d edges found"
                    % (len(cy.vs), len(cycle_edges))
                )
            matching_edges.extend(cycle_edges)

        self.matching_edges = matching_edges

    def display(self):
        """Print the optimal cycles and chains to standard output."""

        print("cycle_count: {}".format(len(self.cycles)))
        print("chain_count: {}".format(len(self.chains)))
        print("cycles:")
        # # cs is a list of cycles, with each cycle represented as a list of vertex IDs
        # # Sort the cycles
        if len(self.cycle_obj) > 0:
            for c in sorted(self.cycle_obj):
                print(c.display())
        else:
            cs = [[v.id for v in c] for c in self.cycles]
            # Put the lowest-indexed vertex at the start of each cycle
            for i in range(len(cs)):
                min_index_pos = cs[i].index(min(cs[i]))
                cs[i] = cs[i][min_index_pos:] + cs[i][:min_index_pos]
                print("\t".join(str(v_id) for v_id in cs[i]))
        print("chains:")
        for c in self.chains:
            print(c.display())

        print("edges:")
        for e in sorted(self.matching_edges, key=lambda x: x.weight, reverse=True):
            print(e.display())

        print("total weight:")
        print(self.total_weight)

    # get weight using a digraph with (possibly) different weights
    def get_weight(self, digraph, ndds, edge_success_prob=1.0):
        weight = sum(
            c.get_weight(digraph, ndds, edge_success_prob) for c in self.chains
        ) + sum(
            failure_aware_cycle_weight(c, digraph, edge_success_prob)
            for c in self.cycles
        )
        return weight


###################################################################################################
#                                                                                                 #
#                  Chain vars and constraints (used by HPIEF', HPIEF'' and PICEF)                 #
#                                                                                                 #
###################################################################################################


def add_chain_vars_and_constraints(
    digraph,
    ndds,
    use_chains,
    max_chain,
    m,
    vtx_to_vars,
    store_edge_positions=False,
    add_o_vars=False,
    num_o_vars=0,
    check_edge_success=False,
):
    """Add the IP variables and constraints for chains in PICEF and HPIEF'.

    Args:
        ndds: a list of NDDs in the instance
        use_chains: boolean: True if chains should be used
        max_chain: the chain cap
        m: The Gurobi model
        vtx_to_vars: A list such that for each Vertex v in the Digraph,
            vtx_to_vars[v.id] will contain the Gurobi variables representing
            edges pointing to v.
        store_edge_positions: if this is True, then an attribute grb_var_positions
            will be added to edges that have associated Gurobi variables.
            edge.grb_var_positions[i] will indicate the position of the edge respresented
            by edge.grb_vars[i]. (default: False)
        add_o_vars: (bool). specific to the heterogeneous edge success probability formulation. add o_vars and big_o_vars
            properties to each edge, and add big_o_vars_in and big_o_vars_out properties to each vertex (analagous to
            v.edges_in and v.edges_out. this requires that e.success_prob is set for each edge
        num_o_vars: (int). if set to >1, then create this number of o and big_o variables, for the DRO-SAA formulation.
            if this is set, each edge needs the property e.realizations - a binary vector of length num_o_vars
    """

    if use_chains:  # max_chain > 0:
        for v in digraph.vs:
            v.grb_vars_in = [[] for i in range(max_chain - 1)]
            v.grb_vars_out = [[] for i in range(max_chain - 1)]
            v.edges_in = [[] for i in range(max_chain - 1)]
            v.edges_out = [[] for i in range(max_chain - 1)]
            if add_o_vars:
                if num_o_vars < 2:
                    v.big_o_vars_in = [[] for i in range(max_chain - 1)]
                    v.big_o_vars_out = [[] for i in range(max_chain - 1)]
                else:
                    # v.big_o_vars[i][n] is the n^th instance of big_o_vars for the i^th edge position...
                    v.big_o_vars_in = [
                        [[] for j in range(num_o_vars)] for i in range(max_chain - 1)
                    ]
                    v.big_o_vars_out = [
                        [[] for j in range(num_o_vars)] for i in range(max_chain - 1)
                    ]

        for ndd in ndds:
            ndd_edge_vars = []
            for e in ndd.edges:
                edge_var = m.addVar(vtype=GRB.BINARY)
                if check_edge_success:
                    if not e.success:
                        m.addConstr(edge_var == 0)
                e.edge_var = edge_var
                ndd_edge_vars.append(edge_var)

                if add_o_vars:
                    if num_o_vars < 2:
                        e.o_var = m.addVar(lb=0, ub=(1 - e.success_prob))
                        e.big_o_var = m.addVar(lb=0, ub=(1 - e.success_prob))
                        m.addConstr(e.big_o_var <= e.edge_var)
                        m.addConstr(e.big_o_var <= e.o_var)
                    else:
                        e.o_vars = m.addVars(num_o_vars, lb=0, ub=1)
                        e.big_o_vars = m.addVars(num_o_vars, lb=0, ub=1)
                        for n in range(num_o_vars):
                            m.addConstr(e.big_o_vars[n] <= e.edge_var)
                            m.addConstr(e.big_o_vars[n] <= e.o_vars[n])
                            m.addConstr(e.o_vars[n] <= e.realizations[n])

                vtx_to_vars[e.tgt.id].append(edge_var)
                if max_chain > 1:
                    e.tgt.grb_vars_in[0].append(edge_var)
                    e.tgt.edges_in[0].append(e)
                    if add_o_vars:
                        if num_o_vars < 2:
                            e.tgt.big_o_vars_in[0].append(e.big_o_var)
                        else:
                            for n in range(num_o_vars):
                                e.tgt.big_o_vars_in[0][n].append(e.big_o_vars[n])

            m.update()
            m.addConstr(quicksum(ndd_edge_vars) <= 1)

        dists_from_ndd = kidney_utils.get_dist_from_nearest_ndd(digraph, ndds)

        # Add pair->pair edge variables, indexed by position in chain
        # e.grb_var are the chain variables for each edge.
        for e in digraph.es:
            e.grb_vars = []
            if add_o_vars:
                if num_o_vars < 2:
                    e.o_vars = []
                    e.big_o_vars = []
                else:
                    # e.o_vars[n][i] is the o-var for the n^th edge realization
                    e.o_vars = []
                    e.big_o_vars = []
            if store_edge_positions:
                e.grb_var_positions = []
            for i in range(max_chain - 1):
                if dists_from_ndd[e.src.id] <= i + 1:
                    edge_var = m.addVar(vtype=GRB.BINARY)
                    if check_edge_success:
                        if not e.success:
                            m.addConstr(edge_var == 0)
                    e.grb_vars.append(edge_var)
                    if add_o_vars:
                        if num_o_vars < 2:
                            o_var = m.addVar(lb=0, ub=(1 - e.success_prob))
                            big_o_var = m.addVar(lb=0, ub=(1 - e.success_prob))
                            e.o_vars.append(o_var)
                            e.big_o_vars.append(big_o_var)
                            m.addConstr(big_o_var <= edge_var)
                            m.addConstr(big_o_var <= o_var)
                        else:
                            o_vars = m.addVars(num_o_vars, lb=0, ub=1)
                            big_o_vars = m.addVars(num_o_vars, lb=0, ub=1)
                            e.o_vars.append(o_vars)
                            e.big_o_vars.append(big_o_vars)
                            for n in range(num_o_vars):
                                m.addConstr(big_o_vars[n] <= edge_var)
                                m.addConstr(big_o_vars[n] <= o_vars[n])
                                m.addConstr(o_vars[n] <= e.realizations[n])
                    if store_edge_positions:
                        e.grb_var_positions.append(i + 1)
                    vtx_to_vars[e.tgt.id].append(edge_var)
                    e.src.grb_vars_out[i].append(edge_var)
                    e.src.edges_out[i].append(e)
                    if add_o_vars:
                        if num_o_vars < 2:
                            e.src.big_o_vars_out[i].append(big_o_var)
                        else:
                            for n in range(num_o_vars):
                                e.src.big_o_vars_out[i][n].append(big_o_vars[n])

                    if i < max_chain - 2:
                        e.tgt.grb_vars_in[i + 1].append(edge_var)
                        e.tgt.edges_in[i + 1].append(e)
                        if add_o_vars:
                            if num_o_vars < 2:
                                e.tgt.big_o_vars_in[i + 1].append(big_o_var)
                            else:
                                for n in range(num_o_vars):
                                    e.tgt.big_o_vars_in[i + 1][n].append(big_o_vars[n])

            # m.update()

        # At each chain position, sum of edges into a vertex must be >= sum of edges out
        for i in range(max_chain - 1):
            for v in digraph.vs:
                m.addConstr(quicksum(v.grb_vars_in[i]) >= quicksum(v.grb_vars_out[i]))
                if add_o_vars:
                    if num_o_vars < 2:
                        assert len(v.edges_out[i]) == len(v.big_o_vars_out[i])
                        m.addConstr(
                            quicksum(v.big_o_vars_in[i])
                            >= quicksum(
                                (1.0 / e.success_prob) * var
                                for e, var in zip(v.edges_out[i], v.big_o_vars_out[i])
                            )
                        )
                    else:
                        for n in range(num_o_vars):
                            m.addConstr(
                                quicksum(v.big_o_vars_in[i][n])
                                >= quicksum(v.big_o_vars_out[i][n])
                            )

        m.update()


###################################################################################################
#                                                                                                 #
#                                              PICEF                                              #
#                                                                                                 #
###################################################################################################


def create_picef_model(cfg, check_edge_success=False, add_o_vars=False, num_o_vars=0):
    """Optimise using the PICEF formulation.

    Args:
        cfg: an OptConfig object
        check_edge_success: (bool). if True, check if each edge has e.success = False. if e.success=False, the edge cannot
            be used.
        add_o_vars/num_o_vars: passed to add_chain_vars_and_constraints

    Returns:
        an OptSolution object
    """

    cycles = cfg.digraph.find_cycles(cfg.max_cycle)

    m = create_mip_model(time_lim=cfg.timelimit, verbose=cfg.verbose)
    m.params.method = -1

    cycle_vars = [m.addVar(vtype=GRB.BINARY) for __ in cycles]

    vtx_to_vars = [[] for __ in cfg.digraph.vs]

    add_chain_vars_and_constraints(
        cfg.digraph,
        cfg.ndds,
        cfg.use_chains,
        cfg.max_chain,
        m,
        vtx_to_vars,
        store_edge_positions=True,
        add_o_vars=add_o_vars,
        num_o_vars=num_o_vars,
        check_edge_success=check_edge_success,
    )  # cfg.edge_success_prob != 1)

    for i, c in enumerate(cycles):
        for v in c:
            vtx_to_vars[v.id].append(cycle_vars[i])

    for l in vtx_to_vars:
        if len(l) > 0:
            m.addConstr(quicksum(l) <= 1)

    # add variables for each pair-pair edge indicating whether it is used in a cycle or chain
    for e in cfg.digraph.es:
        used_in_cycle = []
        for var, c in zip(cycle_vars, cycles):
            if kidney_utils.cycle_contains_edge(c, e):
                used_in_cycle.append(var)

        used_var = m.addVar(vtype=GRB.INTEGER)
        if check_edge_success:
            if not e.success:
                m.addConstr(used_var == 0)

        if cfg.use_chains:
            m.addConstr(used_var == quicksum(used_in_cycle) + quicksum(e.grb_vars))
        else:
            m.addConstr(used_var == quicksum(used_in_cycle))
        e.used_var = used_var

    # number of edges in the matching
    num_edges_var = m.addVar(vtype=GRB.INTEGER)
    pair_edge_count = [e.used_var for e in cfg.digraph.es]
    if cfg.use_chains:
        ndd_edge_count = [e.edge_var for ndd in cfg.ndds for e in ndd.edges]
        m.addConstr(num_edges_var == quicksum(pair_edge_count + ndd_edge_count))
    else:
        m.addConstr(num_edges_var == quicksum(pair_edge_count))

    # add a cardinality restriction if necessary
    if cfg.cardinality_restriction is not None:
        m.addConstr(num_edges_var <= cfg.cardinality_restriction)

    m.update()

    return m, cycles, cycle_vars, num_edges_var


def optimize_picef(cfg, check_edge_success=False):
    m, cycles, cycle_vars, _ = create_picef_model(
        cfg, check_edge_success=check_edge_success
    )

    # add cycle objects
    cycle_list = []
    for c, var in zip(cycles, cycle_vars):
        c_obj = Cycle(c)
        c_obj.add_edges(cfg.digraph.es)
        c_obj.weight = failure_aware_cycle_weight(
            c_obj.vs, cfg.digraph, cfg.edge_success_prob
        )
        c_obj.grb_var = var
        cycle_list.append(c_obj)

    if not cfg.use_chains:
        obj_expr = quicksum(
            failure_aware_cycle_weight(c, cfg.digraph, cfg.edge_success_prob) * var
            for c, var in zip(cycles, cycle_vars)
        )
    elif cfg.edge_success_prob == 1:
        obj_expr = (
            quicksum(
                cycle_weight(c, cfg.digraph) * var for c, var in zip(cycles, cycle_vars)
            )
            + quicksum(e.weight * e.edge_var for ndd in cfg.ndds for e in ndd.edges)
            + quicksum(e.weight * var for e in cfg.digraph.es for var in e.grb_vars)
        )
    else:
        obj_expr = (
            quicksum(
                failure_aware_cycle_weight(c, cfg.digraph, cfg.edge_success_prob) * var
                for c, var in zip(cycles, cycle_vars)
            )
            + quicksum(
                e.weight * cfg.edge_success_prob * e.edge_var
                for ndd in cfg.ndds
                for e in ndd.edges
            )
            + quicksum(
                e.weight * cfg.edge_success_prob ** (pos + 1) * var
                for e in cfg.digraph.es
                for var, pos in zip(e.grb_vars, e.grb_var_positions)
            )
        )
    m.setObjective(obj_expr, GRB.MAXIMIZE)

    optimize(m)

    pair_edges = [e for e in cfg.digraph.es if e.used_var.x > 0.5]

    if cfg.use_chains:
        matching_chains = kidney_utils.get_optimal_chains(
            cfg.digraph, cfg.ndds, cfg.edge_success_prob
        )
        ndd_chain_edges = [
            e for ndd in cfg.ndds for e in ndd.edges if e.edge_var.x > 0.5
        ]
    else:
        ndd_chain_edges = []
        matching_chains = []

    matching_edges = pair_edges + ndd_chain_edges

    if cfg.cardinality_restriction is not None:
        if len(matching_edges) > cfg.cardinality_restriction:
            raise Warning(
                "cardinality restriction is violated: restriction = %d edges, matching uses %d edges"
                % (cfg.cardinality_restriction, len(matching_edges))
            )

    cycles_used = [c for c, v in zip(cycles, cycle_vars) if v.x > 0.5]
    cycle_obj = [c for c in cycle_list if c.grb_var.x > 0.5]

    sol = OptSolution(
        ip_model=m,
        cycles=cycles_used,
        cycle_obj=cycle_obj,
        chains=matching_chains,
        digraph=cfg.digraph,
        edge_success_prob=cfg.edge_success_prob,
        chain_restriction=cfg.chain_restriction,
        cycle_restriction=cfg.cycle_restriction,
        cycle_cap=cfg.max_chain,
        chain_cap=cfg.max_cycle,
        cardinality_restriction=cfg.cardinality_restriction,
    )
    sol.add_matching_edges(cfg.ndds)
    kidney_utils.check_validity(
        sol, cfg.digraph, cfg.ndds, cfg.max_cycle, cfg.max_chain
    )
    return cycle_obj, matching_chains


def max_cycles(cfg):
    """
    Use PICEF to find the maximum number of cycles in a matching...
    """
    m, _, cycle_vars, _ = create_picef_model(cfg)

    num_cycles = quicksum(cycle_vars)

    m.setObjective(num_cycles, GRB.MAXIMIZE)

    optimize(m)
    if cfg.verbose:
        print("maximum number of cycles = %d" % m.objVal)
    if m.objVal != int(m.objVal):
        raise Warning("number of cycles is not integer")
    return int(m.objVal)


###################################################################################################
#                                                                                                 #
#                       PICEF with heterogeneous edge probability                                 #
#                                                                                                 #
###################################################################################################


def expected_cycle_weight(cycle):
    """
    calculate the expected cycle weight, using the cycle's vertices and edge.success_prob properties

    args:
        cycle: (kidney_digrpah.Cycle)
    """

    cycle_prob = 1
    for e in cycle.edges:
        cycle_prob = cycle_prob * e.success_prob

    return sum(e.weight for e in cycle.edges) * cycle_prob


def optimize_picef_heterogeneous_edge_prob(cfg):
    """
    solve the PICEF model with heterogeneous edge success probabilities of Ren, McElfresh, Bidkhori, Dickerson (2020)

    this requires that each edge has the property edge.success_prob
    """
    m, cycles, cycle_vars, _ = create_picef_model(cfg, add_o_vars=True)

    # add cycle objects
    cycle_list = []
    for c, var in zip(cycles, cycle_vars):
        c_obj = Cycle(c)
        c_obj.add_edges(cfg.digraph.es)
        c_obj.weight = expected_cycle_weight(c_obj)
        c_obj.grb_var = var
        cycle_list.append(c_obj)

    # add objective
    chain_weight = quicksum(
        e.weight * big_o_var
        for i in range(cfg.max_chain - 1)
        for v in cfg.digraph.vs
        for e, big_o_var in zip(v.edges_in[i], v.big_o_vars_in[i])
    )

    obj_expr = chain_weight + quicksum(c.weight * c.grb_var for c in cycle_list)

    m.setObjective(obj_expr, GRB.MAXIMIZE)

    optimize(m)

    pair_edges = [e for e in cfg.digraph.es if e.used_var.x > 0.5]

    if cfg.use_chains:
        matching_chains = kidney_utils.get_optimal_chains(
            cfg.digraph, cfg.ndds, cfg.edge_success_prob
        )
        ndd_chain_edges = [
            e for ndd in cfg.ndds for e in ndd.edges if e.edge_var.x > 0.5
        ]
    else:
        ndd_chain_edges = []
        matching_chains = []

    matching_edges = pair_edges + ndd_chain_edges

    if cfg.cardinality_restriction is not None:
        if len(matching_edges) > cfg.cardinality_restriction:
            raise Warning(
                "cardinality restriction is violated: restriction = %d edges, matching uses %d edges"
                % (cfg.cardinality_restriction, len(matching_edges))
            )

    cycles_used = [c for c, v in zip(cycles, cycle_vars) if v.x > 0.5]
    cycle_obj = [c for c in cycle_list if c.grb_var.x > 0.5]

    sol = OptSolution(
        ip_model=m,
        cycles=cycles_used,
        cycle_obj=cycle_obj,
        chains=matching_chains,
        digraph=cfg.digraph,
        edge_success_prob=cfg.edge_success_prob,
        chain_restriction=cfg.chain_restriction,
        cycle_restriction=cfg.cycle_restriction,
        cycle_cap=cfg.max_chain,
        chain_cap=cfg.max_cycle,
        cardinality_restriction=cfg.cardinality_restriction,
    )
    sol.add_matching_edges(cfg.ndds)
    kidney_utils.check_validity(
        sol, cfg.digraph, cfg.ndds, cfg.max_cycle, cfg.max_chain
    )

    return cycle_obj, matching_chains


###################################################################################################
#                                                                                                 #
#                       PICEF with DRO-SAA edge existence                                         #
#                                                                                                 #
###################################################################################################


def optimize_dro_saa_edge_existence(cfg, num_measurements, gamma, alpha):
    """
    solve the PICEF model with DRO-SAA objective for edge existence uncertainty

    this requires that each edge has the property e.realizations: a binary vector of length num_measurements
    """
    m, cycles, cycle_vars, _ = create_picef_model(
        cfg, add_o_vars=True, num_o_vars=num_measurements
    )

    # add cycle objects
    cycle_list = []
    for c, var in zip(cycles, cycle_vars):
        c_obj = Cycle(c)
        c_obj.add_edges(cfg.digraph.es)
        c_obj.weight = expected_cycle_weight(c_obj)
        c_obj.grb_var = var
        # for DRO-SAA - keep track of cycle realizations
        c_obj.realizations = [
            min(e.realizations[n] for e in c_obj.edges) for n in range(num_measurements)
        ]
        cycle_list.append(c_obj)

    # add weight variables for each realization
    w_vars = m.addVars(num_measurements, lb=-GRB.INFINITY, ub=GRB.INFINITY)

    for n in range(num_measurements):
        # objective for the n^th realization
        m.addConstr(
            w_vars[n]
            == (
                -quicksum(
                    e.weight * big_o_var
                    for i in range(cfg.max_chain - 1)
                    for v in cfg.digraph.vs
                    for e, big_o_var in zip(v.edges_in[i], v.big_o_vars_in[i][n])
                )
                - quicksum(c.weight * c.realizations[n] * c.grb_var for c in cycle_list)
            )
        )

    d_var = m.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY)

    # define pi variables
    pi_vars = m.addVars(num_measurements, lb=0, ub=GRB.INFINITY)
    for n in range(num_measurements):
        m.addConstr(pi_vars[n] >= w_vars[n] - d_var)

    # define objective
    obj_expr = (1.0 / float(num_measurements)) * quicksum(w_vars) + gamma * (
        d_var + (1.0 / (alpha * float(num_measurements))) * quicksum(pi_vars)
    )

    m.setObjective(obj_expr, GRB.MINIMIZE)

    optimize(m)

    pair_edges = [e for e in cfg.digraph.es if e.used_var.x > 0.5]

    if cfg.use_chains:
        matching_chains = kidney_utils.get_optimal_chains(
            cfg.digraph, cfg.ndds, cfg.edge_success_prob
        )
        ndd_chain_edges = [
            e for ndd in cfg.ndds for e in ndd.edges if e.edge_var.x > 0.5
        ]
    else:
        ndd_chain_edges = []
        matching_chains = []

    matching_edges = pair_edges + ndd_chain_edges

    if cfg.cardinality_restriction is not None:
        if len(matching_edges) > cfg.cardinality_restriction:
            raise Warning(
                "cardinality restriction is violated: restriction = %d edges, matching uses %d edges"
                % (cfg.cardinality_restriction, len(matching_edges))
            )

    cycles_used = [c for c, v in zip(cycles, cycle_vars) if v.x > 0.5]
    cycle_obj = [c for c in cycle_list if c.grb_var.x > 0.5]

    sol = OptSolution(
        ip_model=m,
        cycles=cycles_used,
        cycle_obj=cycle_obj,
        chains=matching_chains,
        digraph=cfg.digraph,
        edge_success_prob=cfg.edge_success_prob,
        chain_restriction=cfg.chain_restriction,
        cycle_restriction=cfg.cycle_restriction,
        cycle_cap=cfg.max_chain,
        chain_cap=cfg.max_cycle,
        cardinality_restriction=cfg.cardinality_restriction,
    )
    sol.add_matching_edges(cfg.ndds)
    kidney_utils.check_validity(
        sol, cfg.digraph, cfg.ndds, cfg.max_cycle, cfg.max_chain
    )

    return cycle_obj, matching_chains
