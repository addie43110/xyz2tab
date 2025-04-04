import mod
import networkx as nx

from mod import BondType
from graph import Graph
from prettify import header, blue, green, red, warn, bold, underline
from typing import Dict, List

class Reaction:
    '''
    Can either call constructor with Reaction(leftGraph, rightGraph) or 
                                     Reaction(rule) (where rule is type mod.Rule)
    '''
    def __init__(self, rule=None, educts=None, products=None, isMinimal=False, name="no name", context=0):
        self._name: str = rule.name if rule else name

        if educts and products:
            self._educts = educts
            self._products = products
        else:
            self._educts: List[Graph] = Graph(rule.left).connected_components
            self._products: List[Graph] = Graph(rule.right).connected_components

        if len(self._educts) == 1:
            self._leftGraph = self._educts[0]
        else:
            g = nx.Graph()
            [g.update(s.nx_graph) for s in self._educts]
            self._leftGraph = Graph(g)

        if len(self._products) == 1:
            self._rightGraph = self._products[0]
        else:
            g = nx.Graph()
            [g.update(s.nx_graph) for s in self._products]
            self._rightGraph = Graph(g)

        self._mod_rule: mod.Rule = rule if rule else self.set_mod_rule()

        self._is_minimal: bool = isMinimal
        self._maximum_subrule: Reaction = self
        self._minimum_subrule: Reaction = self if isMinimal else self.set_minimum_subrule(add_context=int(context))

    def set_mod_rule(self):
        return mod.ruleGMLString(self.to_ruleGML_string(get_from_mod=False))

    def to_ruleGML_string(self, get_from_mod=False):

            if get_from_mod:
                return self._mod_rule.getGMLString()

            leftGraph = self._leftGraph.nx_graph
            rightGraph = self._rightGraph.nx_graph

            # NODES
            label_changing_nodes = [node for node in leftGraph.nodes if node in rightGraph.nodes and \
                                    leftGraph.nodes[node]['label']!=rightGraph.nodes[node]['label']]
            shared_nodes = [node for node in leftGraph.nodes if node in rightGraph.nodes and node not in label_changing_nodes]
            only_right_nodes = [node for node in rightGraph.nodes if node not in shared_nodes and node not in label_changing_nodes]
            only_left_nodes = [node for node in leftGraph.nodes if node not in shared_nodes and node not in label_changing_nodes]

            # EDGES
            shared_edges = []
            for (u, v, attrs1) in leftGraph.edges(data=True):
                for (p,q,attrs2) in rightGraph.edges(data=True):
                    if {u,v}=={p,q} and attrs1['bond']==attrs2['bond']:
                        shared_edges.append({u,v})
            shared_edges_list = [list(s) for s in shared_edges]

            bond_changing_edges = []
            for (u,v,attrs) in leftGraph.edges(data=True):
                for (p,q,attrs2) in rightGraph.edges(data=True):
                    if {u,v}=={p,q} and attrs['bond']!=attrs2['bond']:
                        bond_changing_edges.append({u,v})
            bond_changing_edges_list = [list(s) for s in bond_changing_edges]

            only_left_edges = []
            for (u,v) in leftGraph.edges:
                if {u,v} not in shared_edges and {u,v} not in bond_changing_edges:
                    only_left_edges.append({u,v})
            only_left_edges_list = [list(s) for s in only_left_edges]

            only_right_edges = []
            for (u,v) in rightGraph.edges:
                if {u,v} not in shared_edges and {u,v} not in bond_changing_edges:
                    only_right_edges.append({u,v})
            only_right_edges_list = [list(s) for s in only_right_edges]

            if self._name == "p1f1!!p1f1p7":
                print(f"shared_edges: {shared_edges_list}")
                print(f"only_left_edges: {only_left_edges_list}")
            """ only_left_edges = []
            rightGraph_set_edges = [{p,q} for (p,q) in rightGraph.edges]
            for (u,v) in leftGraph.edges:
                if u in only_left_nodes or v in only_left_nodes:
                    only_left_edges.append({u,v})
                elif {u,v} not in rightGraph_set_edges:
                    only_left_edges.append({u,v})
            only_left_edges_list = [list(s) for s in only_left_edges]

            only_right_edges = []
            leftGraph_set_edges = [{p,q} for (p,q) in leftGraph.edges]
            for (u,v) in rightGraph.edges:
                if u in only_right_nodes or v in only_right_nodes:
                    only_right_edges.append({u,v})
                elif {u,v} not in leftGraph_set_edges:
                    only_right_edges.append({u,v})
            only_right_edges_list = [list(s) for s in only_right_edges]
            shared_edges = [{u,v} for (u,v) in leftGraph.edges if {u,v} not in only_left_edges and {u,v} not in only_right_edges and \
                                                                {u,v} not in bond_changing_edges]
            shared_edges_list = [list(s) for s in shared_edges] """

            out = []
            out.append("rule [")
            out.append(f'\truleID "{self._name}"')

            out.append("\tleft [")
            out.extend([f"\t\tnode [ id {leftGraph.nodes[node]['modID']} label \"{leftGraph.nodes[node]['label']}\" ]" for node in only_left_nodes])
            out.extend([f"\t\tnode [ id {leftGraph.nodes[node]['modID']} label \"{leftGraph.nodes[node]['label']}\" ]" for node in label_changing_nodes])
            out.extend([f"\t\tedge [ source {leftGraph.nodes[u]['modID']} target {leftGraph.nodes[v]['modID']} label \"{leftGraph.edges[u,v]['bond']}\" ]" for [u,v] in only_left_edges_list])
            out.extend([f"\t\tedge [ source {leftGraph.nodes[u]['modID']} target {leftGraph.nodes[v]['modID']} label \"{leftGraph.edges[u,v]['bond']}\" ]" for [u,v] in bond_changing_edges_list])
            out.append("\t]")


            out.append("\tcontext [")
            out.extend([f"\t\tnode [ id {leftGraph.nodes[node]['modID']} label \"{leftGraph.nodes[node]['label']}\" ]" for node in shared_nodes])
                                                                # could take from left or right, arbitrarily chose leftGraph
            out.extend([f"\t\tedge [ source {leftGraph.nodes[u]['modID']} target {leftGraph.nodes[v]['modID']} label \"{leftGraph.edges[u,v]['bond']}\" ]" for [u,v] in shared_edges_list])
            out.append("\t]")

            out.append("\tright [")
            out.extend([f"\t\tnode [ id {rightGraph.nodes[node]['modID']} label \"{rightGraph.nodes[node]['label']}\" ]" for node in only_right_nodes])
            out.extend([f"\t\tnode [ id {rightGraph.nodes[node]['modID']} label \"{rightGraph.nodes[node]['label']}\" ]" for node in label_changing_nodes])
            out.extend([f"\t\tedge [ source {rightGraph.nodes[u]['modID']} target {rightGraph.nodes[v]['modID']} label \"{rightGraph.edges[u,v]['bond']}\" ]" for [u,v] in only_right_edges_list])
            out.extend([f"\t\tedge [ source {rightGraph.nodes[u]['modID']} target {rightGraph.nodes[v]['modID']} label \"{rightGraph.edges[u,v]['bond']}\" ]" for [u,v] in bond_changing_edges_list])
            out.append("\t]")

            out.append("]")

            return "\n".join(out)

    def set_minimum_subrule(self, add_context=0):
        
        core_left_nodes = []
        core_right_nodes = []
        left_graph = self._leftGraph.nx_graph
        left_nodes = list(left_graph.nodes)
        right_graph = self._rightGraph.nx_graph
        right_nodes = list(right_graph.nodes)

        for (u,v) in left_graph.edges:
            core_left_nodes_list = [t[0] for t in core_left_nodes]
            core_right_nodes_list = [t[0] for t in core_right_nodes]

            # bond in leftGraph, no bond in rightGraph
            try:
                right_graph[u][v]
            except KeyError:
                # since u-v is a bond in leftGraph, u and v must be part of the left core

                if u not in core_left_nodes_list:
                    core_left_nodes.append((u, {'label': left_graph.nodes[u]['label'],
                                           'modID': left_graph.nodes[u]['modID']}))
                if v not in core_left_nodes_list:
                    core_left_nodes.append((v, {'label': left_graph.nodes[v]['label'],
                                           'modID': left_graph.nodes[v]['modID']}))
                
                if u in right_nodes and u not in core_right_nodes_list:
                    core_right_nodes.append((u, {'label': right_graph.nodes[u]['label'],
                                           'modID': left_graph.nodes[u]['modID']}))
                if v in right_nodes and u not in core_right_nodes_list:
                    core_right_nodes.append((v, {'label': right_graph.nodes[v]['label'],
                                           'modID': left_graph.nodes[v]['modID']}))

                continue

            # bond type changes; both left and right must have these nodes
            if left_graph[u][v]['bond'] != right_graph[u][v]['bond']:
                if u not in core_left_nodes_list:
                    core_left_nodes.append((u, {'label': left_graph.nodes[u]['label'],
                                           'modID': left_graph.nodes[u]['modID']}))
                if v not in core_left_nodes_list:
                    core_left_nodes.append((v, {'label': left_graph.nodes[v]['label'],
                                           'modID': left_graph.nodes[v]['modID']}))
                if u not in core_right_nodes_list:
                    core_right_nodes.append((u, {'label': right_graph.nodes[u]['label'],
                                           'modID': right_graph.nodes[u]['modID']}))
                if v not in core_right_nodes_list:
                    core_right_nodes.append((v, {'label': right_graph.nodes[v]['label'],
                                           'modID': right_graph.nodes[v]['modID']}))
    
        for (u,v) in right_graph.edges:
            core_left_nodes_list = [t[0] for t in core_left_nodes]
            core_right_nodes_list = [t[0] for t in core_right_nodes]
            # bond in rightGraph, no bond in leftGraph
            try:
                left_graph[u][v]
            except KeyError:
                # since u-v is a bond in rightGraph, u and v must be part of the right core

                if u not in core_right_nodes_list:
                    core_right_nodes.append((u, {'label': right_graph.nodes[u]['label'],
                                           'modID': right_graph.nodes[u]['modID']}))
                if v not in core_right_nodes_list:
                    core_right_nodes.append((v, {'label': right_graph.nodes[v]['label'],
                                           'modID': right_graph.nodes[v]['modID']}))
                
                if u in left_nodes and u not in core_left_nodes_list:
                    core_left_nodes.append((u, {'label': left_graph.nodes[u]['label'],
                                           'modID': left_graph.nodes[u]['modID']}))
                if v in left_nodes and v not in core_left_nodes_list:
                    core_left_nodes.append((v, {'label': left_graph.nodes[v]['label'],
                                           'modID': left_graph.nodes[v]['modID']}))

        if core_left_nodes==[] and core_right_nodes==[]:
            print(red(f"Cannot set minimum subrule for reaction {self._name} as no bonds change."))
            return False


        # ADDING EXTRA NODES FOR CONTEXT
        while add_context>0:
            added=False
            int_nodes = [v for v,_ in core_left_nodes]
            for v in int_nodes:
                if left_graph.degree[v] > 1:
                    nbs = list(left_graph.neighbors(v))
                    for neighbour in nbs:
                        if neighbour not in int_nodes:
                            core_left_nodes.append((neighbour, {'label':left_graph.nodes[neighbour]['label'], 'modID':left_graph.nodes[neighbour]['modID']}))
                            core_right_nodes.append((neighbour, {'label':right_graph.nodes[neighbour]['label'], 'modID':right_graph.nodes[neighbour]['modID']}))
                            add_context-=1
                            added=True
            if added==False:
                break

        core_left_graph = nx.Graph()
        core_left_graph.add_nodes_from(core_left_nodes)
        int_nodes = core_left_graph.nodes
        core_left_edges = [(u, v, {'bond': left_graph[u][v]['bond']}) for (u,v) in left_graph.edges 
                            if u in int_nodes and v in int_nodes]
        core_left_graph.add_edges_from(core_left_edges)
        

        core_right_graph = nx.Graph()
        core_right_graph.add_nodes_from(core_right_nodes)
        int_nodes = core_right_graph.nodes
        core_right_edges = [(u, v, {'bond': right_graph[u][v]['bond']}) for (u,v) in right_graph.edges 
                            if u in int_nodes and v in int_nodes]
        core_right_graph.add_edges_from(core_right_edges)
                    

        new_left_graph = Graph(core_left_graph)
        new_right_graph = Graph(core_right_graph)

        return Reaction(educts=new_left_graph.connected_components, products=new_right_graph.connected_components, isMinimal=True, name=(str(self._name)+"_min"))

    @property
    def minimum_subrule(self):
        return self._minimum_subrule

    @property
    def maximum_subrule(self):
        return self._maximum_subrule

    @property
    def mod_rule(self):
        return self._mod_rule

    @property
    def educts(self):
        return self._educts

    @property
    def products(self):
        return self._products

    @property
    def left_graph(self):
        return self._leftGraph

    @property
    def right_graph(self):
        return self._rightGraph

    @property
    def name(self):
        return self._name

    def is_minimal(self):
        return self._is_minimal

    def __eq__(self, other: 'Reaction') -> bool:
        return hash(self) == hash(other)

    def __str__(self) -> str:
        return self._name

    def __repr__(self) -> str:
        return self._name

    def __hash__(self) -> int:
        return (hash(self.left_graph)-hash(self.right_graph))


    def _self_check_bond_counts(self):
        for nid in self._leftGraph.nx_graph.nodes:
            left_count = count_bonds(self._leftGraph.nx_graph, nid)
            right_count = count_bonds(self._rightGraph.nx_graph, nid)
            if left_count != right_count:
                print(red(f"BOND COUNT HAS CHANGED! from {left_count} to {right_count}"))

def count_bonds(graph: nx.Graph, v_id):
    count = 0
    for (u, v) in graph.edges([v_id]):
        bt = graph[u][v]['bond']
        if bt == BondType.Single:
            count+=1
        elif bt == BondType.Double:
            count+=2
        elif bt == BondType.Triple:
            count+=3
        elif bt == BondType.Aromatic:
            count+=1.5
        else:
            print(red("BOND DETECTED OTHER THAN EXPECTED"))
    return count