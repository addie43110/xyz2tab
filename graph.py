import mod
import networkx as nx
import time

from prettify import red, warn

# each vertex has attribute "label" which stores string label
# each edge has attribute "bond" which stores bond type (-, =, etc.)

class Graph:
    def __init__(self, nxGraph=None, modGraph: mod.Graph=None, verbose=False):
        self._verbose = verbose
        self._mod_graph = modGraph if modGraph else self.nx_graph_to_mod(nxGraph)
        self._nx_graph = nxGraph if nxGraph else self.mod_to_nx_graph(modGraph)

    
    def mod_to_nx_graph(self, modGraph: mod.Graph):
        g = nx.Graph()

        for v in modGraph.vertices:
            g.add_node(v.id, label=v.stringLabel, modID=v.id)

        for e in modGraph.edges:
            g.add_edge(e.source.id, e.target.id, bond=e.bondType)
        
        return g

    def nx_graph_to_mod(self, nxGraph):
        try:
            return mod.graphGMLString(self.nx_graph_to_GML_string(nxGraph))
        except mod.libpymod.InputError: # graph is not connected probably
            if self._verbose:
                print(warn("Error converting nxGraph to mod graph. Likely graph is not connected. This will not affect rule generation."))
            return None


    def nx_graph_to_GML_string(self, nxGraph):
        out = []
        out.append("graph [")

        out.extend([f"\t\tnode [ id {nxGraph.nodes[node]['modID']} label \"{nxGraph.nodes[node]['label']}\" ]" for node in nxGraph.nodes])
        out.extend([f"\t\tedge [ source {u} target {v} label \"{nxGraph[u][v]['bond']}\" ]"
            for (u,v) in nxGraph.edges])

        out.append("]")

        return "\n".join(out)

    @property
    def nx_graph(self):
        return self._nx_graph

    @property
    def mod_graph(self):
        return self._mod_graph

    @property
    def edges(self):
        return self._nx_graph.edges

    @property
    def nodes(self):
        return self._nx_graph.nodes

    def mod_print(self):
        self._mod_graph.print()

    def __str__(self) -> str:
        return self.nx_graph_to_GML_string(self.nx_graph)

    def __hash__(self) -> int:
        return hash(self.nx_graph_to_GML_string(self.nx_graph))
