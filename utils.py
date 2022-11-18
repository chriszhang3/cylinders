import networkx as nx

def find_generic_pants(cyl_diag):
    """Finds cases when n cylinders are all complete attached to the side of a
    single cylinder.
    
    Input: A cylinder diagram.
    
    Output: A list of lists repsenting homology restrictions."""
    cylinders = cyl_diag.cylinders()
    digraph_data = [[None, None] for _ in range(cyl_diag.degree())]
    for i, (bot, top) in enumerate(cyl_diag.cylinders()):
        for separatrix in bot:
            digraph_data[separatrix][0] = i
        for separatrix in top:
            digraph_data[separatrix][1] = i

    digraph = nx.DiGraph()
    digraph.add_nodes_from(range(len(cylinders)))
    for source, dest in digraph_data:
        digraph.add_edge(source, dest)

    output = []
    for n in digraph:
        successors = list(digraph.successors(n))
        if all([list(digraph.predecessors(suc)) == [n] for suc in successors]):
            successors.append(n)
            output.append(successors)

        predecessors = list(digraph.predecessors(n))
        if all([list(digraph.successors(pre)) == [n] for pre in predecessors]):
            predecessors.append(n)
            # check if this relation already exists
            if not set(successors) == set(predecessors):
                output.append(predecessors)
    return output
