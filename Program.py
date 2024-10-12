import networkx as nx
import osmnx as ox
import numpy as np
import random as rd
import math
from docplex.mp.model import Model

ox.__version__

# download/model a street network for some city then visualize it
# if running code for the first time, one has to de-comment the next three lines
'''G = ox.graph_from_place("Reinickendorf, Berlin, Deutschland", network_type="drive")
H = ox.graph_from_place("Pankow, Berlin, Deutschland", network_type="drive")
F = ox.graph_from_place("Waidmannslust, Reinickendorf, Berlin, Deutschland", network_type="drive") '''

# save/load graph as a graphml file: this is the best way to save your model
# for subsequent work later
# if running code for the first time, one has to de-comment all lines in the next paragraph
filepathg = "./data/reinickendorf.graphml"
#ox.save_graphml(G, filepathg)
filepathh = "./data/pankow.graphml"
#ox.save_graphml(H, filepathh)
filepathf = "./data/waid.graphml"
#ox.save_graphml(F, filepathf)

#load graphs
G = ox.load_graphml(filepathg)
H = ox.load_graphml(filepathh)
F = ox.load_graphml(filepathf)

def completeProb(N, anzKom, capHigh, demHigh):
    # impute speed on all edges missing data
    N = ox.add_edge_speeds(N)
    # calculate travel time (seconds) for all edges
    N = ox.add_edge_travel_times(N)
    # assign random Integer from 0 to capHigh as capacity as na attribute to all edges
    rng = np.random.default_rng()
    edgelist = [e for e in N.edges]
    for e in edgelist:
        N.add_edge(e[0], e[1], e[2], capacity=rng.integers(1,capHigh))
    #choose random nodes as terminals (returns list)
    sample = rd.sample(list(N.nodes), 2*anzKom)

    #compute longest shortest path between sinks and sources with respect to travel times between 
    Graph = nx.DiGraph(N)
    def get_tt(u,  v, e):
        return Graph.edges[u, v]["travel_time"]
    maxi = 0
    for i in range(0,len(sample)//2):
        laenge = nx.shortest_path_length(Graph, source=sample[i], target=sample[i+len(sample)//2], weight=get_tt, method='dijkstra')
        if laenge > maxi:
            maxi = laenge
    print(maxi)

    # first half of terminals become sources and second half become sinks
    for i in range(0,len(sample)//2):
        d=rng.integers(1,demHigh)
        N.add_node("source"+str(i), demand=d)
        N.add_edge("source"+str(i), sample[i], 0, capacity=d+1, travel_time=0)
        N.add_node("sink"+str(i), demand=-d)
        N.add_edge(sample[i+len(sample)//2], "sink"+str(i), 0, capacity=d+1, travel_time=0)
    return N


def constructTimeExpN(N, T, Delta):
    expG = nx.MultiDiGraph()
    expG.add_nodes_from(N.nodes(data=True))
    nx.set_node_attributes(expG, 0, "time")
    i = 1
    exp = False
    while i*Delta < T:
        exp = True
        for node in list(N.nodes):
            expG.add_node(str(node)+"-"+str(i*Delta), time = i*Delta)
        i=i+1
    letzteZeit = (i-1)*Delta
    if exp == False:
        letzteZeit = 0
    for edge in list(N.edges):
        cap = N[edge[0]][edge[1]][edge[2]]["capacity"]
        ttime = N[edge[0]][edge[1]][edge[2]]["travel_time"]
        rtime = Delta * math.ceil(ttime / Delta) #round to next multiple of Delta
        i = 1 
        while i*Delta < T:
            if str(edge[0])+"-"+str(i*Delta) in expG.nodes() and str(edge[1])+"-"+str(i*Delta+rtime) in expG.nodes(): 
                expG.add_edge(str(edge[0])+"-"+str(i*Delta), str(edge[1])+"-"+str(i*Delta+rtime), edge[2], capacity = Delta*cap)
            i=i+1
        #previous while-loop does not add arcs going from and to 0-timelayer this is fixed in next four lines
        if rtime <= letzteZeit and rtime != 0:
            expG.add_edge(edge[0], str(edge[1])+"-"+str(rtime), edge[2], capacity = Delta*cap)
        if rtime <= letzteZeit and rtime == 0:
            expG.add_edge(edge[0], edge[1], edge[2], capacity = Delta*cap)
    # holdover arcs einfügen: 
    for node in list(N.nodes):
        i = 1 
        while (i+1)*Delta < T: 
            expG.add_edge(str(node)+"-"+str(i*Delta), str(node)+"-"+str((i+1)*Delta), 0, capacity=math.inf)
            i=i+1
        #previous while-loop does not add arcs coming from 0-timelayer this is fixed in next two lines
        if 1 <= letzteZeit:
            expG.add_edge(node, str(node)+"-"+str(Delta), 0, capacity=math.inf)
    # change sinks to be the last copy of one node
    node_demand_dict = nx.get_node_attributes(expG, "demand") #Dictionary of attribute values keyed by node
    for node in node_demand_dict:
        if node_demand_dict[node] < 0 and letzteZeit >= 1:
            expG.add_node(str(node)+"-"+str(letzteZeit), demand = node_demand_dict[node])
            del expG.nodes[node]['demand']
    return expG

def solveWithCplex(N): #nicht getestet
    edgelist = [e for e in N.edges]
    nodelist = [e for e in N.nodes]
    node_demand_dict = nx.get_node_attributes(N, "demand") #Dictionary of attribute values keyed by node
    AnzK = len(node_demand_dict)
    variables = []
    for i in range(0,AnzK):
        for e in edgelist:
            variables.append((e,i))
    m = Model(name='static flow')
    def get_cap(edge):
        return N.edges[edge[0][0], edge[0][1], edge[0][2]]["capacity"]
    x = m.continuous_var_dict(variables, 0.0, get_cap)
    for node in nodelist:
        for j in range(0,AnzK):
            getout = [(edge,j) for edge in N.out_edges(nbunch = node, keys = True)]
            getin =[(edge,j) for edge in N.in_edges(nbunch = node, keys = True)]
            if node in node_demand_dict: #fullfill demands for every good
                m.add_constraint(sum(x[i] for i in getout)-sum(x[i] for i in getin)==node_demand_dict[node])
            else: #Flowconservation for every good
                m.add_constraint(sum(x[i] for i in getin)==sum(x[i] for i in getout))
    return m.solve()

def putItAllTogether(G, anzKom, capHigh, demHigh):
    G = completeProb(G, anzKom, capHigh, demHigh)
    for Delta in range(1,11):
        sol = None
        T = Delta
        while sol == None: #find upper and lower bound on T
            expG = constructTimeExpN(G, T, Delta)
            sol = solveWithCplex(expG)
            T = T*2
        upper = T//2
        lower = upper//2
        T = math.ceil((lower+upper)/2)
        T = Delta * math.ceil(T / Delta) #round to next multiple of Delta
        while T != upper: #binary search
            expG = constructTimeExpN(G, T, Delta)
            sol = solveWithCplex(expG)
            if sol == None:
                lower = T
            else:
                upper = T
            T = math.ceil((lower+upper)/2)
            T = Delta * math.ceil(T / Delta) #round to next multiple of Delta
        print("Für Delta = ", Delta, "ist der optimale Zeithorizont", upper)
    return True

