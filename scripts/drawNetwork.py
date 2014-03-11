import sys
import pydot
import numpy as np

graph   = pydot.Dot(graph_type='digraph')

network = np.loadtxt(sys.argv[1])
T, E    = np.array_split(network.flatten(), 2)

pos = [
"5, 7.5!", #Bcd
"15, 7.5!", #Cad
"7.5, 10!", #hb
"12.5, 10!", #gt
"7.5, 5!", #Kr
"12.5, 5!", #Kni
"17, 9!", #hkb
"17, 6!"] #tll


node_bcd = pydot.Node("Bcd", pos = pos[0], style="filled",
						 fillcolor = "orange", width = 1.2, height = 1.2)
node_cad = pydot.Node("Cad", pos = pos[1], style="filled",
						 fillcolor = "cyan",   width = 1.2, height = 1.2)
node_hb  = pydot.Node("Hb",  pos = pos[2], style="filled",
						 fillcolor = "red", width = 1.2, height = 1.2)
node_gt  = pydot.Node("Gt",  pos = pos[3], style="filled",
						 fillcolor = "blue",   width = 1.2, height = 1.2)
node_kr  = pydot.Node("Kr",  pos = pos[4], style="filled",
						 fillcolor = "green",  width = 1.2, height = 1.2)
node_kni = pydot.Node("Kni", pos = pos[5], style="filled",
						 fillcolor = "pink",    width = 1.2, height = 1.2)
node_tll = pydot.Node("Tll", pos = pos[6], style="filled",
						 fillcolor = "purple",   width = 1.2, height = 1.2)
node_hkb = pydot.Node("Hkb", pos = pos[7], style="filled",
						 fillcolor = "gray",   width = 1.2, height = 1.2)

graph.add_node(node_bcd)
graph.add_node(node_cad)
graph.add_node(node_tll)
graph.add_node(node_hkb)
graph.add_node(node_hb)
graph.add_node(node_gt)
graph.add_node(node_kr)
graph.add_node(node_kni)

gap_genes      = [node_hb, node_kr, node_gt, node_kni]
external_genes = [node_bcd, node_cad, node_tll, node_hkb]

ngenes = 4
for i, g1 in enumerate(gap_genes):
	for j, g2 in enumerate(gap_genes):
		if g1 == node_kni: color = 'pink'
		if g1 == node_kr: color = 'green'
		if g1 == node_gt: color = 'blue'
		if g1 == node_hb: color = 'red'
		if abs(T[j + i*ngenes]) > 0.01:
			graph.add_edge(pydot.Edge(g1, g2, 
									penwidth = str(30*abs(T[j + i*ngenes])),
									arrowhead = "normal" if T[j + i*ngenes] > 0 else "tee",
									style = "solid" if T[j + i*ngenes] > 0 else "solid",
									color = color))

	for j , ex in enumerate(external_genes):
		if ex == node_bcd: color = 'orange'
		if ex == node_cad: color = 'cyan'
		if ex == node_tll: color = 'purple'
		if ex == node_hkb: color = 'gray'
		if abs(E[j + i*ngenes]) > 0.01:
			graph.add_edge(pydot.Edge(ex, g1, 
									penwidth = str(30*abs(E[j + i*ngenes])), 
									arrowhead = "normal" if E[j + i*ngenes] > 0 else "tee", 
									style = "solid" if E[j + i*ngenes] > 0 else "solid", 
									color = color))



# graph.write_png(sys.argv[1] + '.png')#, prog='neato')
graph.write_png(sys.argv[1] + '.png', prog='neato')
