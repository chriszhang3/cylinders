from sage.all import DiGraph, matrix

class CylinderGraph:
    """
    Turn a cylinder diagram to a directed graph (an instance of nx.Digraph)
    as follows:

    Every cylinder is a vertex of the graph. For each horizontal saddle
    connection `s`, there is a directed edge from the cylinder under `s` to the
    cylinder above `s`. 
    """

    def __init__(self, cd) -> None:
        """Input is a cylinder diagram.
        self.digraph is the graph described above."""
        cylinders = cd.cylinders()

        # Using `cylinders`, create a list of edges.
        # The i-th element of saddle_data is [under, above], where `under` is
        # the cylinder under saddle i and `above` is the saddle above it
        saddle_data = [[None, None] for _ in range(cd.degree())]
        for i, (bot, top) in enumerate(cylinders):
            for saddle in bot:
                saddle_data[saddle][1] = i
            for saddle in top:
                saddle_data[saddle][0] = i
        
        # Turns the list of edges into an adjacency matrix. The reason we're
        # doing this is that we want a weighted graph.
        adjacency_matrix = matrix(len(cylinders))
        for pre, suc in saddle_data:
            adjacency_matrix[pre, suc] += 1

        self.digraph = DiGraph(adjacency_matrix, loops=True, weighted = True)
