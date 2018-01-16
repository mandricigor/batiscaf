
# this is the class which solves orientation problem

import sys
import cplex
import math
import networkx as nx
from cplex.exceptions import CplexSolverError
from collections import Counter
from Bio.Seq import Seq
from Queue import Queue


def orient(graph):
    cpx = cplex.Cplex()
    cpx.set_results_stream("/dev/null")

    edges = list(graph.edges())
    nodes = list(graph.nodes())

    for node in nodes:
        name = "X#%s" % node
        cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])

    for u, v in edges:
        name = "Z#%s#%s" % (u, v)
        cpx.variables.add(lb=[0], ub=[1], types=["B"], names=[name])


    for u, v in edges:
        inds1 = ["X#%s" % u, "X#%s" % v, "Z#%s#%s" % (u, v)]
        vals1 = [1, 1, -1]
        c1 = cplex.SparsePair(ind=inds1, val=vals1)
        cpx.linear_constraints.add( \
            lin_expr = [c1],\
            senses = ["G"],\
            rhs = [0],\
            names = ['repeat-edge-1']\
        )
        inds2 = ["Z#%s#%s" % (u, v), "X#%s" % u]
        vals2 = [1, -1]
        c2 = cplex.SparsePair(ind=inds2, val=vals2)
        cpx.linear_constraints.add( \
            lin_expr = [c2],\
            senses = ["G"],\
            rhs = [0],\
            names = ['repeat-edge-2']\
        )
        inds3 = ["Z#%s#%s" % (u, v), "X#%s" % v]
        vals3 = [1, -1]
        c3 = cplex.SparsePair(ind=inds3, val=vals3)
        cpx.linear_constraints.add( \
            lin_expr = [c3],\
            senses = ["G"],\
            rhs = [0],\
            names = ['repeat-edge-3']\
        )

    contigs = set()
    for node in graph.nodes():
        contigs.add(node[:-2])
    for contig in contigs:
        inds = ["X#%s_1" % contig, "X#%s_2" % contig]
        vals = [1, -1]
        c = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [c],\
            senses = ["E"],\
            rhs = [0],\
            names = ['same-node']\
        )

    for node in nodes:
        neighbors = list(graph.neighbors(node))
        inds = []
        for z in neighbors:
            if (node, z) in edges:
                ind = (node, z)
            else:
                ind = (z, node)
            inds.append("Z#%s#%s" % ind)
        vals = [1 for z in range(len(neighbors))]
        c = cplex.SparsePair(ind=inds, val=vals)
        cpx.linear_constraints.add( \
            lin_expr = [c],\
            senses = ["G"],\
            rhs = [len(neighbors) - 2],\
            names = ['repeat-node-1']\
        )
        inds_prime = inds + ["X#%s" % node]
        vals_prime = vals + [-1]
        c_prime = cplex.SparsePair(ind=inds_prime, val=vals_prime)
        cpx.linear_constraints.add( \
            lin_expr = [c_prime],\
            senses = ["L"],\
            rhs = [len(neighbors) - 1],\
            names = ['repeat-node-2']\
        )


    for node in graph.nodes():
        #cpx.objective.set_linear("X#%s" % node, graph.nodes[node]["width"] * 1.0 / (0.000000001 + graph.nodes[node]["score"]))
        #cpx.objective.set_linear("X#%s" % node, 1.0 / (1 + graph.nodes[node]["width"] * graph.nodes[node]["score"]))
        cpx.objective.set_linear("X#%s" % node, graph.nodes[node]["width"])
        #cpx.objective.set_linear("X#%s" % node, 1.0 / (1 + graph.nodes[node]["score"]))

    #for ee1, ee2 in graph.edges():
    #    cpx.objective.set_linear("Z#%s#%s" % (ee1, ee2), graph[ee1][ee2]["weight"])

    cpx.objective.set_sense(cpx.objective.sense.minimize)
    cpx.set_problem_type(cpx.problem_type.MILP)
    cpx.write("program.txt", filetype="lp")
    start = cpx.get_time()
    cpx.solve()
    end = cpx.get_time()
    print end - start, "THIS MUCH TIME"

    for x, y in zip(cpx.variables.get_names(), cpx.solution.get_values()):
        print x, y  

    print int(cpx.solution.get_objective_value()), "THIS MANY REPEATS"

    toremove = []
    for x, y in zip(cpx.variables.get_names(), cpx.solution.get_values()):
        if "X" in x and y > 0.5: # this is 1 in fact - a repeat or a bad node for now
            toremove.append(x)
            print x
    print len(toremove) / 2, "I REMOVED THIS MANY REPEATS"
    return toremove
  

#gg = nx.read_graphml("../RECOMB2018/hsapiens-97-200-stranded.graphml")
#gg = nx.read_graphml("../RECOMB2018/soap2-fijiensis-97-200-stranded.graphml")
#gg = nx.read_graphml("../RECOMB2018/fijiensis-97-200-stranded.graphml")
#gg = nx.read_graphml("../RECOMB2018/rhodobacter-97-200-stranded.graphml")
gg = nx.read_graphml("soap2-rhodobacter-97-200-stranded.graphml")
#gg = nx.read_graphml("../RECOMB2018/scaffmatch-97-200-stranded.graphml")
#gg = nx.read_graphml("../RECOMB2018/graminicola-97-200-stranded.graphml")
#gg = nx.read_graphml("../RECOMB2018/soap2-graminicola-97-200-stranded.graphml")


"""
with open("scores") as f:
    a = f.readlines()
a = map(lambda x: x.strip().split(), a)
adict = {}
for x, y, z in a:
    adict[x] = float(z)

for node in gg.nodes():
    gg.nodes[node]["score"] = adict[node[:-2]]
"""

toremove = []
for x, y in gg.edges():
    if gg[x][y]["weight"] < 5:
        toremove.append((x, y))
for x, y in toremove:
    gg.remove_edge(x, y)


gg_praorig = gg.copy()


matching = {}

wws = []
for node in gg.nodes():
    nei = {}
    for node2 in gg.neighbors(node):
        nei[node2] = gg[node][node2]["weight"]
    if nei:
        maximum = max(nei.items(), key=lambda z: z[1])
        matching[(node, maximum[0])] = maximum[1]
        matching[(maximum[0], node)] = maximum[1]
        wws.append(maximum[1])

print max(wws), min(wws), "WWS"

toremove = []
for e1, e2 in gg.edges():
    if (e1, e2) not in matching:
        toremove.append((e1, e2))
for e1, e2 in toremove:
    gg.remove_edge(e1, e2)


contigs = set()
for node in gg.nodes():
    contigs.add(node[:-2])
for contig in contigs:
    gg.add_edge(contig + "_1", contig + "_2", dist=0, weight=100000000)




gg_orig = gg.copy()


# MST REMOVAL
for x, y in gg.edges():
    gg[x][y]["weight"] = -gg[x][y]["weight"]
gg = nx.minimum_spanning_tree(gg)
for x, y in gg.edges():
    gg[x][y]["weight"] = -gg[x][y]["weight"]

# try to impute copy number scores


toremove = orient(gg)
almost_repeats = toremove

gg2 = gg.copy()
for node in toremove:
    gg2.remove_node(node[2:])

scaffolds = {}
orders = []
oris = []
distances_all = []
chains = [] # this is a list of graphs-chains
trusted_nodes = set()
extremities = []
transitively_removed_edges = {}
scaf_ends = {}
ccc = 1
for cc in nx.connected_components(gg2):
    if len(cc) <= 2:
        continue
    subgraph = nx.Graph(gg2.subgraph(cc))
    if nx.cycle_basis(subgraph) != []: # it has cycles
        edg = list(subgraph.edges())
        minimal = edg[0]
        minw = gg2[minimal[0]][minimal[1]]["weight"]
        for e in edg[1:]:
            if gg2[e[0]][e[1]]["weight"] < minw:
                minw = gg2[e[0]][e[1]]["weight"]
                minimal = e
        subgraph.remove_edge(*minimal)
    start, end = [x for x in subgraph.nodes() if subgraph.degree(x) == 1]
    path = nx.shortest_path(subgraph, start, end)
    print path, "THIS IS PATH"
    chain = nx.DiGraph()
    nodes = [path[i][:-2] for i in range(len(path)) if i % 2 == 0]
    chain_no_ones = []
    for i in range(len(path) - 1):
        if path[i][:-2] == path[i + 1][:-2]:
            continue
        chain_no_ones.append((path[i], path[i + 1]))

    order = []
    orients = []
    distances = []
    for i in range(len(chain_no_ones)):
        chain2 = chain_no_ones[i]
        x, y = chain2
        ori = (x[-1], y[-1])
        if i == 0:
            if ori == ("1", "2"):
                chain.add_node(x[:-2], ori="+", length=gg.node[x[:-2] + "_1"]["width"])
                chain.add_node(y[:-2], ori="+", length=gg.node[y[:-2] + "_1"]["width"])
                orients.extend(["+", "+"])
            elif ori == ("2", "1"):
                chain.add_node(x[:-2], ori="-", length=gg.node[x[:-2] + "_1"]["width"])
                chain.add_node(y[:-2], ori="-", length=gg.node[y[:-2] + "_1"]["width"])
                orients.extend(["-", "-"])
            elif ori == ("1", "1"):
                chain.add_node(x[:-2], ori="+", length=gg.node[x[:-2] + "_1"]["width"])
                chain.add_node(y[:-2], ori="-", length=gg.node[y[:-2] + "_1"]["width"])
                orients.extend(["+", "-"])
            elif ori == ("2", "2"):
                chain.add_node(x[:-2], ori="-", length=gg.node[x[:-2] + "_1"]["width"])
                chain.add_node(y[:-2], ori="+", length=gg.node[y[:-2] + "_1"]["width"])
                orients.extend(["-", "+"])
            pair = tuple(sorted([x, y]))
            distance = max(0, int(gg[pair[0]][pair[1]]["dist"]))
            distances.append(distance)
            chain.add_edge(x[:-2], y[:-2], dist=max(0, distance))
            order.extend([x[:-2], y[:-2]])
        else:
            if ori == ("1", "2"):
                if orients[-1] == "+":
                    chain.add_node(y[:-2], orien="+", length=gg.node[y[:-2] + "_1"]["width"]) # done
                    orients.append("+")
                else:
                    chain.add_node(y[:-2], orien="+", length=gg.node[y[:-2] + "_1"]["width"]) # this is very strange!!!!!!!!!1
                    orients.append("+")
            elif ori == ("2", "1"):
                if orients[-1] == "-":
                    chain.add_node(y[:-2], orien="-", length=gg.node[y[:-2] + "_1"]["width"]) # done
                    orients.append("-")
                else:
                    chain.add_node(y[:-2], orien="-", length=gg.node[y[:-2] + "_1"]["width"]) # hz ego znaet
                    orients.append("-")
            elif ori == ("1", "1"):
                if orients[-1] == "+":
                    chain.add_node(y[:-2], orien="-", length=gg.node[y[:-2] + "_1"]["width"]) # done
                    orients.append("-")
                else:
                    chain.add_node(y[:-2], orien="-", length=gg.node[y[:-2] + "_1"]["width"]) # Igor, check
                    orients.append("-")
            elif ori == ("2", "2"):
                if orients[-1] == "-":
                    chain.add_node(y[:-2], orien="+", length=gg.node[y[:-2] + "_1"]["width"]) # done
                    orients.append("+")
                else:
                    chain.add_node(y[:-2], orien="+", length=gg.node[y[:-2] + "_1"]["width"]) # Igor, check
                    orients.append("+")
            pair = tuple(sorted([x, y]))
            distance = max(0, int(gg[pair[0]][pair[1]]["dist"]))
            distances.append(distance)
            chain.add_edge(x[:-2], y[:-2], dist=max(0, distance))
            order.append(y[:-2])
    chains.append(chain)
    orders.append(order)
    oris.append(orients)
    distances_all.append(distances)
    print order
    print orients

    adjacencies = []
    for i in range(0, len(path[1:-1]), 2):
        adjacencies.append((path[1:-1][i], path[1:-1][i + 1]))
    print adjacencies, "ADJACENCIES"

    for a1, a2 in adjacencies:
        print gg_orig[a1][a2]["weight"], a1, a2, "WEIGHT OF EDGE"

    # make an induced subgraph

    #strands = []
    #for gm, gmgm in adjacencies:
    #    strands.append(gm)
    #    strands.append(gmgm)
    strands = path[1:-1] # this is new - not to include garbage at boundaries!
    extremities.append((path[0], path[-1]))


    conts = set()
    for st in strands:
        neighs = list(gg_orig.neighbors(st))
        print neighs, "THIS IS NEIGHS"
        for nei in neighs:
            print nei, st, gg_orig[nei][st]["weight"], nei in strands
        for nn in [n for n in neighs if gg_orig[st][n]["weight"] > 5 and n not in path]:
            conts.add(nn[:-2])

    # handle the ends of scaffolds
    scaf_ends[ccc] = {}
    scaf_ends[ccc]["start"] = set()
    scaf_ends[ccc]["end"] = set()
    neighs = list(gg_praorig.neighbors(path[0]))
    for nn in [n for n in neighs if gg_praorig[path[0]][n]["weight"] > 5 and n not in path]:
        scaf_ends[ccc]["start"].add(nn[:-2])
    scaf_ends[ccc]["start"].add(path[0][:-2])
    neighs = list(gg_praorig.neighbors(path[-1]))
    for nn in [n for n in neighs if gg_praorig[path[-1]][n]["weight"] > 5 and n not in path]:
        scaf_ends[ccc]["end"].add(nn[:-2])
    scaf_ends[ccc]["end"].add(path[-1][:-2])


    strands2 = []
    for cc in conts:
        strands2.extend([cc + "_1", cc + "_2"])
    subg = nx.Graph(gg_orig.subgraph(strands2 + strands))

    """
    conts = set()
    for st in strands:
        neighs = list(gg_praorig.neighbors(st))
        print neighs, "THIS IS NEIGHS"
        #for nei in neighs:
        #    print nei, st, gg_orig[nei][st]["weight"], nei in strands
        for nn in [n for n in neighs if gg_praorig[st][n]["weight"] > 5 and n not in path]:
            conts.add(nn[:-2])
    strands2 = []
    for cc in conts:
        strands2.extend([cc + "_1", cc + "_2"])
    subg = nx.Graph(gg_praorig.subgraph(strands2 + strands))
    """





    toremove = []
    for e1, e2 in subg.edges():
        if e1 in strands2 or e2 in strands2:
            if subg[e1][e2]["weight"] < 5:
                toremove.append((e1, e2))
    for e1, e2 in toremove:
        subg.remove_edge(e1, e2)
    for node in subg.nodes():
        if node in strands:
            subg.node[node]["color"] = "red"
        else:
            subg.node[node]["color"] = "blue"
    print subg.nodes(), "SUBG NODES" 
    # make directed chains
    ordict = {}
    for u, v in zip(order, orients):
        ordict[u] = v
    dch = nx.DiGraph()
    print order, "ORDER HERE"
    for i in range(1, len(order)):
        distanta = chain[order[i - 1]][order[i]]["dist"]
        dch.add_edge("%s:%s" % (order[i - 1], orients[i - 1]), "%s:%s" % (order[i], orients[i]))
        dch["%s:%s" % (order[i - 1], orients[i - 1])]["%s:%s" % (order[i], orients[i])]["dist"] = distanta
        dch.nodes["%s:%s" % (order[i - 1], orients[i - 1])]["color"] = "red"
        dch.nodes["%s:%s" % (order[i], orients[i])]["color"] = "red"


    #nx.write_graphml(dch, "pregraph_%s.graphml" % ccc)

    print ordict, order, orients

    print strands2

    #"""""""""""""""""""""""""""""""""""
    for ee1, ee2 in subg.edges():
        print ee1, ee2, ee1 in strands2, ee2 in strands2, "SUBG EDGE"
        if ee1 in strands2 and ee2 not in strands2:
            e1 = ee1
            e2 = ee2
        elif ee1 not in strands2 and ee2 in strands2:
            e1 = ee2
            e2 = ee1
        else:
            print "CONTINUE"
            continue

        distanta = gg_praorig[ee1][ee2]["dist"]
        if ordict[e2[:-2]] == "-" and e2[-1] == "2": # ot nego
            if e1[-1] == "1":
                dch.add_edge(e2[:-2] + ":" + ordict[e2[:-2]], e1[:-2] + ":-", dist=distanta)
            elif e1[-1] == "2":
                dch.add_edge(e2[:-2] + ":" + ordict[e2[:-2]], e1[:-2] + ":+", dist=distanta)
        elif ordict[e2[:-2]] == "-" and e2[-1] == "1": # k nemu
            if e1[-1] == "1":
                dch.add_edge(e1[:-2] + ":+", e2[:-2] + ":" + ordict[e2[:-2]], dist=distanta)
            elif e1[-1] == "2":
                dch.add_edge(e1[:-2] + ":-", e2[:-2] + ":" + ordict[e2[:-2]], dist=distanta)
        elif ordict[e2[:-2]] == "+" and e2[-1] == "1": # ot nego
            if e1[-1] == "1":
                dch.add_edge(e2[:-2] + ":" + ordict[e2[:-2]], e1[:-2] + ":-", dist=distanta)
            elif e1[-1] == "2":
                dch.add_edge(e2[:-2] + ":" + ordict[e2[:-2]], e1[:-2] + ":+", dist=distanta)
        elif ordict[e2[:-2]] == "+" and e2[-1] == "2": # k nemu
            if e1[-1] == "1":
                dch.add_edge(e1[:-2] + ":+", e2[:-2] + ":" + ordict[e2[:-2]], dist=distanta)
            elif e1[-1] == "2":
                dch.add_edge(e1[:-2] + ":-", e2[:-2] + ":" + ordict[e2[:-2]], dist=distanta)
         



    print "CYCLE BASIS", list(nx.cycle_basis(subg))
    for cycle in list(nx.cycle_basis(subg)):
        dirs = []
        for c in cycle:
            dirs.append(c[-1])
        if dirs == ['1', '2'] * (len(cycle) / 2) or dirs == ['2', '1'] * (len(cycle) / 2):
            pass
        else:
            print "PROBLEMATIC CYCLE", cycle

    #nx.write_graphml(subg, "subgraph_%s.graphml" % ccc)

    # add here widths for each node
    for nod in dch.nodes():
        dch.node[nod]["width"] = gg_orig.node[nod.split(":")[0] + "_2"]["width"]


    print dch.nodes(), "NODES OF DCH"


    is_dag = nx.dag.is_directed_acyclic_graph(dch)
    if is_dag:
        dch2 = nx.transitive_reduction(dch)
        toremove = []
        for x, y in dch.edges():
            if not dch2.has_edge(x, y):
                toremove.append((x, y))
        transitively_removed_edges[id(dch)] = {}
        for x, y in toremove:
            transitively_removed_edges[id(dch)][(x, y)] = dch[x][y]["dist"]
            dch.remove_edge(x, y)
    else:
        while True:
            try:
                cycle = list(nx.find_cycle(dch))
            except nx.exception.NetworkXNoCycle:
                cycle = []
            if cycle == []:
                break
            cycle_nodes = set()
            for x, y in cycle:
                cycle_nodes.add(x)
                cycle_nodes.add(y)
            unmarked = {}
            for nod in cycle_nodes:
                if "color" not in dch.nodes[nod]:
                    unmarked[nod] = dch.node[nod]["width"]
            #print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", cycle, unmarked
            minimal_to_remove = min(unmarked.items(), key=lambda z: z[1])[0]
            affected_edges = [(x, y) for (x, y) in cycle if minimal_to_remove in (x, y)]
            print affected_edges
            for copy, ed in zip([1, 2], affected_edges):
                ed1, ed2 = ed # this is our edge to replace with a copy of a node
                if ed1 == minimal_to_remove:
                    dch.add_node("%s:%s" % (minimal_to_remove, copy), width=dch.node[minimal_to_remove]["width"])
                    dch.add_edge("%s:%s" % (minimal_to_remove, copy), ed2, dist=dch[minimal_to_remove][ed2]["dist"])
                else:
                    dch.add_node("%s:%s" % (minimal_to_remove, copy), width=dch.node[minimal_to_remove]["width"])
                    dch.add_edge(ed1, "%s:%s" % (minimal_to_remove, copy), dist=dch[ed1][minimal_to_remove]["dist"])
            dch.remove_node(minimal_to_remove)
        # transitive reduction
        dch2 = nx.transitive_reduction(dch)
        toremove = []
        for x, y in dch.edges():
            if not dch2.has_edge(x, y):
                toremove.append((x, y))
        transitively_removed_edges[id(dch)] = {}
        for x, y in toremove:
            transitively_removed_edges[id(dch)][(x, y)] = dch[x][y]["dist"]
            dch.remove_edge(x, y)
   
    articulation_points = set(nx.articulation_points(dch.to_undirected()))
    topsort = list(nx.topological_sort(dch))
    art_order = {}
    for ii, cc in enumerate(topsort):
        if cc in articulation_points:
            art_order[cc] = ii
    art_order = map(lambda zz: zz[0], sorted(art_order.items(), key=lambda z: z[1]))
    print "CUTS CUTS", art_order, order[0], orients[0], order[-1], orients[-1]
    dch_start = order[0] + ":" + orients[0]
    dch_finish = order[-1] + ":" + orients[-1]
    if dch_start not in art_order:
        art_order = [dch_start] + art_order
    if dch_finish not in art_order:
        art_order = art_order + [dch_finish]

    for nod in dch.nodes():
        dch.node[nod]["status"] = "UNMARKED"
    # first, mark articulation points with token "ART"
    for ap in art_order:
        dch.node[ap]["status"] = "ART"
    slots = {}
    art_order = ["START_TOKEN"] + art_order + ["END_TOKEN"]
    for i in range(1, len(art_order)):
        slots[(art_order[i - 1], art_order[i])] = {} # let's make it dictionary where keys are small contigs and values are distances from the first articulation point
    
    # the strategy is like this: we first see which slots those in paths belong to
    # and then see the place of the danlgling nodes


    def dfs(mygraph, start, end, acc=[]):
        if start == end:
            return acc[:-1]
        else:
            igor = set()
            children = [child for _, child in mygraph.out_edges(start)]
            if start != end and children == []:
                return []
            else:
                for child in children:
                    igor.update(dfs(mygraph, child, end, acc + [child]))
            return igor
 


    # paths stuff - a DFS-like algorithm
    for i in range(1, len(art_order) - 2):
        intermediate = dfs(dch, art_order[i], art_order[i + 1], [])
        for nod in intermediate:
            dch.node[nod]["status"] = "MARKED"
        #slots[(art_order[i], art_order[i + 1])].extend(intermediate)
        # if there was no edge at all in the graph, we have to provide the distance estimate
        # but, we have to put this in transitively_removed_edges because we will search later on for this distance
        if not dch.has_edge(art_order[i], art_order[i + 1]) and (art_order[i], art_order[i + 1]) not in transitively_removed_edges[id(dch)]:

            from itertools import islice
            def k_shortest_paths(G, source, target, k, weight=None):
                return list(islice(nx.shortest_simple_paths(G, source, target, weight=weight), k))

            shortest_paths_to_consider = k_shortest_paths(dch, art_order[i], art_order[i + 1], 10, weight="dist") # TODO: it is biased, should be another estimate
            longest_path = shortest_paths_to_consider[-1]
            estimate = 0
            for ii in range(1, len(longest_path)):
                estimate += dch[longest_path[ii - 1]][longest_path[ii]]["dist"]
            for ii in longest_path[1:-1]:
                estimate += dch.node[ii]["width"] # this is because the gap is composed of gaps and contigs!
            transitively_removed_edges[id(dch)][(art_order[i], art_order[i + 1])] = estimate
            print longest_path, "SHORTEST_PATH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        else:
            if dch.has_edge(art_order[i], art_order[i + 1]):
                estimate = dch[art_order[i]][art_order[i + 1]]["dist"]
            else:
                estimate = transitively_removed_edges[id(dch)][(art_order[i], art_order[i + 1])]
        # now try to estimate distance between each contig and let's say the first articulation point
        for inter in intermediate:
            path_to_inter = nx.shortest_path(dch, art_order[i], inter, weight="dist")
            di = 0
            for ii in range(1, len(path_to_inter)):
                di += dch[path_to_inter[ii - 1]][path_to_inter[ii]]["dist"]
            for ii in path_to_inter[1:-1]:
                di += dch.node[ii]["width"]
            slots[(art_order[i], art_order[i + 1])][inter] = di

    print slots

    # now dangling stuff 
    # iterate through slots first
    for i in range(0, len(art_order) - 1):
        if i == 0: # this is the first slot
            #continue
            print "FIRST SLOT", art_order[i], art_order[i + 1]
            # take all incoming to the art_order[i + 1] and push them to the slot - they can not go anywhere else
            # they also must be unmarked
            for desc, _ in dch.in_edges(art_order[i + 1]):
                if dch.node[desc]["status"] == "UNMARKED":
                    # the arrow should exist in the graph, so it is easy to find the distance to the second contig of the slot
                    slots[(art_order[i], art_order[i + 1])][desc] = dch[desc][art_order[i + 1]]["dist"]
                dch.node[desc]["status"] = "MARKED"
        elif i + 1 == len(art_order) - 1: #this is the second slot
            #continue
            print "LAST SLOT", art_order[i], art_order[i + 1]
            for _, asc in dch.out_edges(art_order[i]):
                if dch.node[asc]["status"] == "UNMARKED":
                    # the same, arrow should exist in the graph
                    slots[(art_order[i], art_order[i + 1])][asc] = dch[art_order[i]][asc]["dist"]
                dch.node[asc]["status"] = "MARKED"
        else:
            # first, examine the out dangling edges
            for _, asc in dch.out_edges(art_order[i]):
                if dch.node[asc]["status"] == "UNMARKED":
                    if dch.has_edge(art_order[i], art_order[i + 1]):
                        dd2 = dch[art_order[i]][art_order[i + 1]]["dist"]
                    else: # if it is not in the chain, it was transitively removed
                        dd2 = transitively_removed_edges[id(dch)][(art_order[i], art_order[i + 1])]
                    wd2 = dch.node[art_order[i + 1]]["width"]
                    dd1 = dch[art_order[i]][asc]["dist"]
                    wd1 = dch.node[asc]["width"]
                    if wd1 > dd2 + 1800:
			#dch.remove_edge(art_order[i], asc)
			dch.remove_node(asc)
                        continue
                    # the crucial stuff! in which slot to put the contig!
                    print dd1, wd1, dd2, wd2, asc, art_order[i], art_order[i + 1]
                    if dd1 + wd1 - dd2 < 300: # to change later
                        slots[(art_order[i], art_order[i + 1])][asc] = dd1 # TODO: CHECK!
                    else: # assuming it is going to the next slot!
                        slots[(art_order[i + 1], art_order[i + 2])][asc] = dd1 - dd2 # TODO: CHECK!
                    dch.node[asc]["status"] = "MARKED"
            # now examine the in dangling edges
            for desc, _ in dch.in_edges(art_order[i + 1]):
                print desc, "DESCDESC", dch.node[desc]["status"]
                if dch.node[desc]["status"] == "UNMARKED":
                    if dch.has_edge(art_order[i], art_order[i + 1]):
                        dd2 = dch[art_order[i]][art_order[i + 1]]["dist"]
                    else: # if it is not in the chain, it was transitively removed
                        dd2 = transitively_removed_edges[id(dch)][(art_order[i], art_order[i + 1])]
                    wd2 = dch.node[art_order[i + 1]]["width"]
                    dd1 = dch[desc][art_order[i + 1]]["dist"]
                    wd1 = dch.node[desc]["width"]
                    if wd1 > dd2 + 1800:
			#dch.remove_edge(desc, art_order[i + 1])
			dch.remove_node(desc)
                        continue
                    # the crucial stuff! in which slot to put the contig!
                    if dd1 + wd1 - dd2 < 300: # to change later
                        slots[(art_order[i], art_order[i + 1])][desc] = dd2 - dd1 - wd1 # TODO: CHECK!
                    else: # assuming it is going to the next slot!
                        slots[(art_order[i - 1], art_order[i])][desc] = dd1 - dd2 - wd2 # TODO: CHECK!
                    dch.node[desc]["status"] = "MARKED"
                print desc, "DESCDESC", dch.node[desc]["status"], "VTOTOE"

    print "SLOTS ARE HERE", slots

    # now sort contigs in each slots - one of the most difficult parts
    for i in range(0, len(art_order) - 1):
        if i == 0: # these are sorted by the distance to the front articulation point
            dist_dict_sorted = sorted(slots[(art_order[i], art_order[i + 1])].items(), key=lambda z: z[1], reverse=True)
            print dist_dict_sorted, "WHEN 0"
            # remove them all for now in each slot
            for newcontig, di in dist_dict_sorted:
                dch.remove_node(newcontig)
            for newcontig, di in dist_dict_sorted:
                dch.add_node(newcontig)
            for_edges = [newc[0] for newc in dist_dict_sorted] + [art_order[i + 1]]
            for ii in range(1, len(for_edges)):
                dch.add_edge(for_edges[ii - 1], for_edges[ii], dist=100) # TODO: change 100 to the correct one!
        else: # but these ones are sorted by the distance to the back articulation point
            dist_dict_sorted = sorted(slots[(art_order[i], art_order[i + 1])].items(), key=lambda z: z[1])
            print dist_dict_sorted, "WHEN", i
            if dch.has_edge(art_order[i], art_order[i + 1]):
                print "I REMOVED EDGE 1", art_order[i], art_order[i + 1]
                dch.remove_edge(art_order[i], art_order[i + 1])
            print len(dch.nodes()), "BEFORE", dch.nodes()
            print dist_dict_sorted, "WHEN AGAIN", i
            for newcontig, di in dist_dict_sorted:
                print "REMOVING", newcontig
                dch.remove_node(newcontig)
            print len(dch.nodes()), "ADRE", dch.nodes()
            for newcontig, di in dist_dict_sorted:
                dch.add_node(newcontig)
            for_edges = [art_order[i]] + [newc[0] for newc in dist_dict_sorted]
            if i + 1 != len(art_order) - 1:
                for_edges += [art_order[i + 1]]
            for ii in range(1, len(for_edges)):
                dch.add_edge(for_edges[ii - 1], for_edges[ii], dist=100) # TODO: change 100 to the correct one!

    #"""""""""""""""""""""""""""""""""""

    scaffolds[ccc] = dch
    #nx.write_graphml(dch, "subgraph_igor_%s.graphml" % ccc)
    ccc += 1

    ###########################


print "SCAF ENDS", scaf_ends


#for number, scaf in scaffolds.items():
    #start = [x for x in scaf.nodes() if scaf.in_degree(x) == 0]
    #end = [x for x in scaf.nodes() if scaf.out_degree(x) == 0]
    #assert len(start) == 1, "FFFFFFFFFF %s %s" % (len(start), number)
    #assert len(end) == 1, "KKKKKKKKKKKi %s" % (len(end))


involved = set()
for number, scaf in scaffolds.items():
    nodes = scaf.nodes()
    for node in nodes:
        involved.add(node.split(":")[0])
print len(involved), "INVOLVED"

gnodes = set()
for node in gg_praorig.nodes():
    gnodes.add(node[:-2])

#almost = set()
#for x in almost_repeats:
#    almost.add(x[2:-2])

rest = gnodes - involved
print rest, "REST"

#for r in rest:
#    for x in ["_1", "_2"]:
#        node = r + x
#        ma = {}
#        for nei in gg_orig.neighbors(node):
#            if gg_orig[node][nei]["weight"] != 100000000:
#                ma[nei] = gg_orig[node][nei]["weight"]
#        print node, ma




#gg2 = gg_orig.copy()
#for node in involved:
#    gg2.remove_node(node + "_1")
#    gg2.remove_node(node + "_2")


#for x in (involved - almost):
#    for y in (involved - almost):


"""
new_edges = []
for x in involved:
    for y in rest:
        if x != y:
            for u in ["_1", "_2"]:
                for v in ["_1", "_2"]:
                    if gg_praorig.has_edge(x + u, y + v) and gg_praorig[x + u][y + v]["weight"] >= 25:
                        already_involved = set()
                        for nn, scaf in scaffolds.items():
                            for vertex in scaf.nodes():
                                if vertex.startswith(x + ":"):
                                    already_involved.add(vertex)
                        if len(already_involved) == 1:
                            # this is a unique one
                            print x + u, y + v, gg_praorig[x + u][y + v]["weight"], "DEEEEEEEEEECI", already_involved
                            new_edges.append(((y + v, x + u), list(already_involved)[0]))


node_scaf_dict = {}
for number, scaf in scaffolds.items():
    for node in scaf.nodes():
        node_scaf_dict[node] = number

for edge, invnode in new_edges:
    invori = invnode.split(":")[1]
    e1, e2 = edge
    distanta = gg_praorig[e1][e2]["dist"]
    ddch = scaffolds[node_scaf_dict[invnode]]
    print edge, invnode, "AAAAAAAAAAA"
    if 1 == 0 and ((invori == "-" and e2[-1] == "2") or (invori == "+" and e2[-1] == "1")): # ot nego
    #if ((invori == "-" and e2[-1] == "2") or (invori == "+" and e2[-1] == "1")): # ot nego
        children = []
        for _, desc in ddch.out_edges(invnode):
            children.append(desc)
        if e1[-1] == "1":
            newori = ":-"
        elif e1[-1] == "2":
            newori = ":+"
        if not children:
            ddch.add_edge(invnode, e1[:-2] + newori, dist=100)
        else:
            # more complicated
            ddch.remove_edge(invnode, children[0])
            ddch.add_edge(invnode, e1[:-2] + newori, dist=100)
            ddch.add_edge(e1[:-2] + newori, children[0], dist=100)
    elif 1 == 0 and ((invori == "-" and e2[-1] == "1") or (invori == "+" and e2[-1] == "2")): # k nemu
    #elif ((invori == "-" and e2[-1] == "1") or (invori == "+" and e2[-1] == "2")): # k nemu
        children = []
        for asc, _ in ddch.in_edges(invnode):
            children.append(asc)
        if e1[-1] == "1":
            newori = ":+"
        elif e1[-1] == "2":
            newori = ":-"
        if not children:
            ddch.add_edge(e1[:-2] + newori, invnode, dist=100)
        else:
            # more complicated
            ddch.remove_edge(children[0], invnode)
            ddch.add_edge(e1[:-2] + newori, invnode, dist=100)
            ddch.add_edge(children[0], e1[:-2] + newori, dist=100)
"""





# try preliminary results - make fasta scaffolds

fastas = []
for number, scaf in scaffolds.items():
    print [x for x in scaf.nodes() if scaf.in_degree(x) == 0], [x for x in scaf.nodes() if scaf.out_degree(x) == 0], "STARTS"
    print len(list(nx.weakly_connected_components(scaf))), "CONNECTED COMP"
    start = [x for x in scaf.nodes() if scaf.in_degree(x) == 0][0]
    print start, "THIS IS START HERE"
    beginning = start # this is to test if we can join chains
    name_parts = start.split(":")
    contig, orient = name_parts[:2]
    seq = Seq(gg_orig.nodes[contig + "_1"]["seq"])
    if orient == "-":
        seq = str(seq.reverse_complement())
    else:
        seq = str(seq)
    fasta = seq
    while True:
        out_edges = list(scaf.out_edges(start))
        if len(out_edges) == 0:
            break
        else:
            prevstart, start = out_edges[0]
            print start, "BUT NOW THIS IS START"
            distance = scaf[prevstart][start]["dist"]
            name_parts = start.split(":")
            contig, orient = name_parts[:2]
            seq = Seq(gg_orig.nodes[contig + "_1"]["seq"])
            if orient == "-":
                seq = str(seq.reverse_complement())
            else:
                seq = str(seq)
            fasta += "N" * int(distance)
            fasta += seq
    ending = start # this is to test if we can join them
    print beginning, ending
    print len(scaf.nodes()), len(fasta)
    fastas.append(fasta)


with open("output.scaffolds.repeat-aware.fasta", "w") as f:
    for i, scaf in enumerate(fastas):
        f.write(">scaffold_%s\n" % (i + 1))
        f.write("%s\n" % scaf)
    # add the rest of them
    j = i
    for k in rest:
        f.write(">scaffold_%s\n" % j)
        f.write("%s\n" % Seq(gg_orig.nodes[k + "_1"]["seq"]))
        j += 1

print rest

print len(scaffolds)


