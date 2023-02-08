import networkx as nx

class CylinderGraph:
    """
    Turn a cylinder diagram to a directed graph as follows:
    Every cylinder is a vertex of the graph. For each horizontal saddle
    connection `s`, there is a directed edge from the cylinder under `s` to the
    cylinder above `s`. 
    """

    def __init__(self, cd) -> None:
        """Input is a cylinder diagram.
        self.digraph is the graph described above."""
        cylinders = cd.cylinders()

        # The i-th element of saddle_data is [under, above], where `under` is
        # the cylinder under saddle_data[i] and `above` is the saddle above it
        saddle_data = [[None, None] for _ in range(cd.degree())]
        for i, (bot, top) in enumerate(cylinders):
            for saddle in bot:
                saddle_data[saddle][1] = i
            for saddle in top:
                saddle_data[saddle][0] = i

        self.digraph = nx.DiGraph()
        self.digraph.add_nodes_from(range(len(cylinders)))
        for source, dest in saddle_data:
            self.digraph.add_edge(source, dest)
        
    # TODO: Add unittests for this
    def find_generic_pants(self):
        """Finds cases when n cylinders are all only attached to the side
        of a single cylinder. This is like a generalized version of a
        topological pair of pants allows n pant legs.
        
        Input: A cylinder diagram.
        
        Output: A list of pants. Each pants is a frozenset consisting of the 
        cylinders in the pants."""

        pants_list = []
        for n in self.digraph:
            # If every cylinder above `C` is only adjacent to `C` along its
            # bottom, this is a generic pants.
            successors = list(self.digraph.successors(n))
            if all([list(self.digraph.predecessors(suc)) == [n] 
                    for suc in successors]):
                pants = frozenset(successors + [n])
                if pants not in pants_list:
                    pants_list.append(pants)

            # If every cylinder below `C` is only adjacent to `C` along its
            # top, this is a generic pants.
            predecessors = list(self.digraph.predecessors(n))
            if all([list(self.digraph.successors(pre)) == [n] 
                    for pre in predecessors]):
                pants = frozenset(predecessors + [n])
                if pants not in pants_list:
                    pants_list.append(pants)
        return pants_list
    
    def find_leaves(self):
        """Return the leaves of self.digraph.
        
        A leaf is a vertex such that it only has one neighbor, where we 
        count a vertex as a neighbor if there is either an edge coming from it
        or an edge going to it."""
        leaves = []
        for n in self.digraph:
            neighbors = set(self.digraph.successors(n)) | \
                        set(self.digraph.predecessors(n))
            if len(neighbors) == 1:
                leaves.append(n)
        return leaves
