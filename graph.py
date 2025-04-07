import mod
import networkx as nx
import time
import re

from prettify import red, warn
from lookup_tables import valence_electrons

# each vertex has attribute "label" which stores string label
# each edge has attribute "bond" which stores bond order (-, =, etc.)

class Graph:
    def __init__(self, graph, verbose=False):
        self._verbose = verbose
        
        mod_lg_class = mod.libpymod.Rule.LeftGraph
        mod_rg_class = mod.libpymod.Rule.RightGraph
        mod_g_class = mod.libpymod.Graph
        update_mod = True
        if isinstance(graph, mod_lg_class) or isinstance(graph, mod_rg_class) or isinstance(graph, mod_g_class):
            self._mod_graph = graph
            self._nx_graph = self.mod_to_nx_graph(graph)
            self._gml_string = self.nx_graph_to_GML_string(self._nx_graph)
            update_mod = False
        elif isinstance(graph, nx.Graph):
            self._nx_graph = graph
            self._gml_string = self.nx_graph_to_GML_string(self._nx_graph)
        elif isinstance(graph, str):
            self._gml_string = graph
            self._nx_graph = self.GML_to_nx_graph(self._gml_string)
        else:
            print(red(f"ERROR: Graph class cannot identify graph type in constructor: {type(graph)}"))

        # up to the caller to check the number of components created
        ccs = [self._nx_graph.subgraph(c).copy() for c in nx.connected_components(self._nx_graph)]
        ccs = sorted(ccs, key=len, reverse=True)
        self._num_components = len(ccs)
        if self._num_components > 1:
            self._components = [Graph(c) for c in ccs]
        else:
            self._components = [self]

        if self._num_components == 1 and update_mod:
            self._mod_graph = self.nx_graph_to_mod(self._nx_graph)

    
    def mod_to_nx_graph(self, modGraph: mod.Graph):
        g = nx.Graph()

        for v in modGraph.vertices:
            g.add_node(int(v.id), label=str(v.stringLabel), modID=int(v.id))

        for e in modGraph.edges:
            g.add_edge(int(e.source.id), int(e.target.id), bond=str(e.bondType))
        
        return g

    def nx_graph_to_mod(self, nxGraph):
        try:
            return mod.graphGMLString(self.nx_graph_to_GML_string(nxGraph))
        except mod.libpymod.InputError: # graph is not connected probably
            if self._verbose:
                print(warn("Error converting nxGraph to mod graph. Likely graph is not connected. This will not affect rule generation."))
            return None

    def GML_to_nx_graph(self, gml_string):
        g = nx.Graph()
        lines = gml_string.split("\n")
        for line in lines:
            tokens = line.split()
            if len(tokens) < 2:
                continue
            if tokens[0] == "node":
                (_,_,_, mid, _, l, _) = tokens
                g.add_node(int(mid), label=l[1:-1], modID=int(mid))
            elif tokens[0] == "edge":
                (_, _, _, u, _, v, _, bondOrder, _ ) = tokens
                g.add_edge(int(u), int(v), bond=bondOrder[1:-1])
        return g

    def nx_graph_to_GML_string(self, nxGraph):
        out = []
        out.append("graph [")

        out.extend([f"\t\tnode [ id {nxGraph.nodes[node]['modID']} label \"{nxGraph.nodes[node]['label']}\" ]" for node in nxGraph.nodes])
        out.extend([f"\t\tedge [ source {u} target {v} label \"{nxGraph[u][v]['bond']}\" ]"
            for (u,v) in nxGraph.edges])

        out.append("]")

        return "\n".join(out)

    """ def calculate_formal_charge(self):
        if self._num_components > 1:
            for c in self._components:
                c.calculate_charge()
            return

        for v, attrs1 in self._nx_graph.nodes(data=True):
            full_label = attrs1['label']
            m = re.match(r"^([a-zA-Z]+)([-\+])?$", full_label)
            atom_name = m.group(1)

            # if the label already has a charge, assume it is correct and skip
            if m.group(2):
                continue
            
            valence_electrons[atom_name] - count_bonds(self._nx_graph, v) """


    @property
    def nx_graph(self):
        return self._nx_graph

    @property
    def mod_graph(self):
        return self._mod_graph

    @property  
    def gml_string(self):
        return self._gml_string

    @property
    def edges(self):
        return self._nx_graph.edges

    @property
    def nodes(self):
        return self._nx_graph.nodes

    @property
    def num_components(self):
        return self._num_components

    @property
    def connected_components(self):
        return self._components

    def mod_print(self):
        self._mod_graph.print()

    def __str__(self) -> str:
        return self._gml_string

    def __hash__(self) -> int:
        return hash(self._gml_string)
