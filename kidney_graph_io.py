import glob
import os
import pandas as pd

from kidney_digraph import Digraph
from kidney_ndds import Ndd, NddEdge


def read_cmu_format(details_filename, maxcard_filename, frac_edges=None, seed=101):
    # read a "cmu" format exchange graph, using the details and maxcard files
    #
    # optional : frac_edges in (0, 1) adds only a fraction of the edges to the Digraph.

    name = os.path.basename(maxcard_filename)

    # read details.input file
    col_names = [
        "id",
        "abo_patient",
        "abo_fonor",
        "wife_patient",
        "pra",
        "in_deg",
        "out_deg",
        "is_ndd",
        "is_marginalized",
    ]
    df_details = pd.read_csv(
        details_filename, names=col_names, skiprows=1, delim_whitespace=True
    )

    pair_details = df_details.loc[df_details["is_ndd"] == 0]
    pair_id = list(pair_details["id"].unique())

    # vtx_index[id] gives the index in the digraph
    vtx_index = dict(zip(pair_id, range(len(pair_id))))

    vtx_count = len(vtx_index)
    digraph = Digraph(vtx_count)

    # label sensitized pairs
    for index, row in pair_details.iterrows():
        if row["is_marginalized"]:
            digraph.vs[vtx_index[row["id"]]].sensitized = True

    # read maxcard.inuput file (edges)
    col_names = ["src_id", "tgt_id", "weight", "c4", "c5"]
    df_edges = pd.read_csv(
        maxcard_filename, names=col_names, skiprows=1, delim_whitespace=True
    )

    # drop the last column
    df_edges.drop(df_edges.index[-1])

    # take only nonzero edges
    nonzero_edges = df_edges.loc[df_edges["weight"] > 0]

    # optional: sample from the edges
    if frac_edges is not None:
        assert (frac_edges < 1.0) and (frac_edges > 0.0)
        nonzero_edges = nonzero_edges.sample(frac=frac_edges, random_state=seed)

    # ind ndds if they exist
    ndd_details = df_details.loc[df_details["is_ndd"] == 1]
    ndd_count = len(ndd_details)

    if ndd_count > 0:
        ndd_list = [Ndd(id=i) for i in range(ndd_count)]
        ndd_id = list(ndd_details["id"].unique())

        # ndd_index[id] gives the index in the ndd list
        ndd_index = dict(zip(ndd_id, range(len(ndd_id))))
    else:
        ndd_list = []
        ndd_index = []

    use_ndds = ndd_count > 0

    # add edges to pairs and ndds
    for index, row in nonzero_edges.iterrows():
        src = row["src_id"]
        tgt_id = vtx_index[row["tgt_id"]]
        weight = row["weight"]
        if use_ndds and (src in ndd_index.keys()):  # this is an ndd edge
            src_id = ndd_index[src]
            ndd_list[src_id].add_edge(
                NddEdge(digraph.vs[tgt_id], weight, src_id=ndd_list[src_id].id)
            )
        else:  # this edge is a pair edge
            src_id = vtx_index[src]
            digraph.add_edge(weight, digraph.vs[src_id], digraph.vs[tgt_id])

    return digraph, ndd_list, name


def get_cmu_graphs(directory):
    # create a generator that produces kidney exchange graphs, given a directory containing "cmu" format exchange
    # graph files.

    # find all *maxcard.input files in the directory -- each corresponds to an exchange graph
    maxcard_files = glob.glob(os.path.join(directory, "*maxcard.input"))

    for maxcard_file in maxcard_files:

        file_base = "_".join(maxcard_file.split("_")[:-1])

        # find the details file; there can be only one
        details_files = glob.glob(
            os.path.join(directory, file_base + "_*details.input")
        )
        assert len(details_files) == 1
        details_file = details_files[0]

        digraph, ndd_list, name = read_cmu_format(details_file, maxcard_file)

        yield digraph, ndd_list, name
