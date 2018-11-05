#!/usr/bin/env python
import os, argparse, sys
import networkx as nx

def topological_features(file):
    f = open(file,"r")
    read = f.read() 
    f.close()
    a=read.split()
    G=nx.Graph()
    for j in range(0,len(a),2):
        G.add_edge(a[j],a[j+1])
    aa=nx.degree(G)
    bb=nx.closeness_centrality(G)
    cc=nx.betweenness_centrality(G)
    dd=nx.clustering(G)
    return aa, bb, cc, dd

def main(args = sys.argv):
    parser = argparse.ArgumentParser('Calculates topological features')
    parser.add_argument('-net', '--network', default='net.txt',
                        help='residue network file')
    parser.add_argument('-o', '--out', default='out.txt',
                        help='out file')
    args = parser.parse_args()
    
    degree,closeness,betweenness,clustering = topological_features(args.network)

    f = open(args.out, 'w')
    for key in degree:
        f.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(key[0],key[1],closeness[key[0]],betweenness[key[0]],clustering[key[0]]))
    f.close()
    
if __name__ == '__main__':   
    main()
