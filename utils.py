import os
import time


def generate_filepath(output_dir, name, extension):
    # generate filepath, of the format <name>_YYYYMMDD_HHMMDD<extension>
    timestr = time.strftime("%Y%m%d_%H%M%S")
    output_string = (name + '_%s.' + extension) % timestr
    return os.path.join(output_dir, output_string)


def score_after_edge_realization(cycles, chains):
    '''
    calculate the score of cycles/chains after edge success/failure, using property e.success

    args:
    - cycles: list(Cycle).
    - chains: list(Chain)
    '''
    # if the cycles does not have failed edges, add it
    score = 0.0
    for c in cycles:
        if all([e.success for e in c.edges]):
            score += sum([e.weight for e in c.edges])

    # calculate the truncated chain score if an edge failed
    for ch in chains:
        for e in ch.edges:
            if e.success:
                score += e.weight
            else:
                break

    return score


def initialize_edges(digraph, ndd_list, rs):
    '''initialized edge success probabilities on U[0.1, 0.9], and set all edge weights to 1'''

    low_p = 0.1
    high_p = 0.9

    for n in ndd_list:
        for e in n.edges:
            e.success_prob = rs.uniform(low=low_p, high=high_p)

    for e in digraph.es:
        e.success_prob = rs.uniform(low=low_p, high=high_p)


def realize_edge_existence(digraph, ndd_list, rs):
    '''realize each edge- set the property e.success to True or False, according to e.success_prob. the true weight will
    need to account for edge failures.'''

    for n in ndd_list:
        for e in n.edges:
            e.success = rs.uniform(low=0.0, high=1.0) < e.success_prob

    for e in digraph.es:
        e.success = rs.uniform(low=0.0, high=1.0) < e.success_prob


def set_edge_realizations(digraph, ndd_list, num_measurements, rs):
    '''set the property e.realizations to a list of num_measurements realizations (binary list)

    this requires each edge to have the property e.success_prob'''

    for n in ndd_list:
        for e in n.edges:
            e.realizations = rs.choice([0, 1], num_measurements, p=[1 - e.success_prob, e.success_prob])

    for e in digraph.es:
        e.realizations = rs.choice([0, 1], num_measurements, p=[1 - e.success_prob, e.success_prob])
