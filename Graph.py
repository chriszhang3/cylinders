import networkx as nx

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

        # The i-th element of saddle_data is [under, above], where `under` is
        # the cylinder under saddle i and `above` is the saddle above it
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
    def find_generalized_pants(self):
        """Finds cases when n cylinders are all only attached to the side
        of a single cylinder. This is like a generalized version of a
        topological pair of pants allows n pant legs.
        
        Input: A cylinder diagram.
        
        Output: A set of pants. Each pants is a frozenset consisting of the 
        cylinders in the pants.
        
        Note: For Python 3.7 and higher, frozensets should maintain insertion
        order. Thus, the first element of each pants is the waist curve and the
        rest are the pant legs. However, we do not use this in our code."""

        pants_set = set()
        for n in self.digraph:
            # If every cylinder above `C` is only adjacent to `C` along its
            # bottom, this is a generalized pants.
            successors = list(self.digraph.successors(n))
            if all([list(self.digraph.predecessors(suc)) == [n] 
                    for suc in successors]):
                
                pants = frozenset([n] + successors)
                pants_set.add(pants)

            # If every cylinder below `C` is only adjacent to `C` along its
            # top, this is a generalized pants.
            predecessors = list(self.digraph.predecessors(n))
            if all([list(self.digraph.successors(pre)) == [n] 
                    for pre in predecessors]):

                pants = frozenset([n] + predecessors)
                pants_set.add(pants)
        return pants_set
    
    def find_leaves(self):
        """Return the tuples (leaf, neighbor) for all leaves of self.digraph.
        
        A leaf is a vertex such that it only has one neighbor, where we 
        count a vertex as a neighbor if there is either an edge coming from it
        or an edge going to it.
        
        For a given leaf, `neighbor` is it's unique neighbor."""
        leaf_neighbers = []
        for n in self.digraph:
            neighbors = set(self.digraph.successors(n)) | \
                        set(self.digraph.predecessors(n))
            if len(neighbors) == 1:
                leaf_neighbers.append((n, next(iter(neighbors))))
        return leaf_neighbers
