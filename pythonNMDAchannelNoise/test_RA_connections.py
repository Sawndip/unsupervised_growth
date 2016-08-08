# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 15:04:57 2016

@author: jingroup
"""

"""
Script check RA to RA connections
"""
import matplotlib.pyplot as plt
import reading
import networkx as nx

fileSuper = "/home/eugene/Output/RA_RA_super_connections.bin"
fileActive = "/home/eugene/Output/RA_RA_connections.bin"
filenameRA_I = "/home/eugene//Output/RA_I_connections.bin"
filenameI = "/home/eugene/Output/I_RA_connections.bin"
filename_xy_RA = "/home/eugene/Output/RA_xy.bin"
filename_xy_I = "/home/eugene/Output/I_xy.bin"
filename_weights = "/home/eugene/Output/weights.bin"

(N_RA, RA_super_targets, RA_super_targets_G) = reading.read_connections(fileSuper)
(N_RA, RA_active_targets, RA_active_targets_G) = reading.read_connections(fileActive)


print "RA_active_targets: ", RA_active_targets
print "RA_active_targets_G: ", RA_active_targets_G

print "RA_super_targets: ", RA_super_targets
print "RA_super_targets_G: ", RA_super_targets_G


(N_RA, RA2I_edges, RA2I_weights) = reading.get_RA2I_graph(filenameRA_I)
(N_I, I2RA_edges, I2RA_weights) = reading.get_I2RA_graph(N_RA, filenameI)
#(unused, edges, weights) = reading.get_RA2RA_graph(filenameRA_RA)
        
        #print self.I2RA_edges[0][0]
        #print self.edges        
        #print self.RA2I_edges
        # read coordinates
(xx_RA, yy_RA) = reading.read_coordinates(filename_xy_RA)
(xx_I, yy_I) = reading.read_coordinates(filename_xy_I)
    
        # fill dictionary for positions of neurons
pos = {} 
       
for i in range(N_RA):
    pos[i] = (xx_RA[i], yy_RA[i])
for i in range(N_I):   
    pos[i + N_RA] = (xx_I[i], yy_I[i])
        
node_color = []
        # assign colors to nodes
for i in range(N_RA):
    node_color.append("g")
for i in range(N_I):
    node_color.append("r")
        
        #self.weights = map(lambda x: x*10, self.weights)        
        
        # create bidirectional graph
        
G = nx.DiGraph()
G.add_nodes_from(pos.keys())
       
        # biderectional graph for fixed synapses       
G_I2RA = nx.DiGraph()
G_I2RA.add_nodes_from(pos.keys())
        
G_RA2I = nx.DiGraph()
G_RA2I.add_nodes_from(pos.keys())
                
edge_color_I2RA = []        
edge_color_RA2I = []        
        
for i in range(len(I2RA_edges)):
    G_I2RA.add_edge(I2RA_edges[i][0], I2RA_edges[i][1], weight=I2RA_weights[i])
    edge_color_I2RA.append('r')
                
for i in range(len(RA2I_edges)):
    G_RA2I.add_edge(RA2I_edges[i][0], RA2I_edges[i][1], weight=RA2I_weights[i])
    edge_color_RA2I.append(1)
            
edgewidth_I2RA = [ d['weight'] for (u,v,d) in G_I2RA.edges(data=True)]
edgewidth_RA2I = [ d['weight'] for (u,v,d) in G_RA2I.edges(data=True)]

from collections import Counter

duplicates = [k for k,v in Counter(RA2I_edges).items() if v>1]


print "RA2I edges: ", RA2I_edges
print "duplicates: ", duplicates
#print "edge_color_RA2I = ", edge_color_RA2I


print len(RA2I_edges)
print len(G_RA2I.edges())

nx.draw_networkx(G_RA2I, pos=pos, node_size = 130, width = edgewidth_RA2I, 
                      node_color=node_color, with_labels=True, font_size=5)
        #axes.set_xlim((-5,SIDE+5))
        #axes.set_ylim((-5,SIDE+5))
plt.show()