#!/usr/bin/python -tt
# -*- coding: utf-8 -*-
# =============================================================================
# Time-stamp: <Fre 2017-08-25 13:39 juergen>
# File      : fw.py 
# Creation  : 14 Jun 2017
#
# Copyright (c) 2017 Jürgen Hackl <hackl@ibi.baug.ethz.ch>
#               http://www.ibi.ethz.ch
# $Id$ 
#
# Description : A simple traffic model based on Frank-Wolfe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
# =============================================================================

import os
import subprocess
import networkx as nx
import numpy as np
import itertools
import tempfile

class TrafficAssignment(object):
    def __init__(self, graph, od_graph,od_matrix=None):
        self.graph = graph.copy()
        self.od_graph = od_graph.copy()
        self.od_matrix = np.copy(od_matrix)
        self.lost_trips = {}

        self.print_flag = True
        temp = True
        self.temp_folder = './pytras/temp/'
        #self.temp_folder = './temp/'
        if temp:
            self.temp_dir = tempfile.mkdtemp(dir = self.temp_folder)
            self.temp_name = os.path.split(self.temp_dir)[-1]
        else:
            self.temp_dir = self.temp_folder+'test/'
            self.temp_name = 'test'

        self.alpha = 0.15
        self.beta = 4
        self.toll = 0
        self.unit_factor = 1000
        self.threshold = 1e-4
        self.large = 1e+14
        # number of iterations for the corrective factors
        self.n_iter_tm = 1000
        pass

    def get_edge_attribute(self,edge,attribute):
        return self.graph[edge[0]][edge[1]][attribute]

    def set_edge_attribute(self,edge,attribute,value):
        self.graph[edge[0]][edge[1]][attribute] = value
        pass

    def set_od_matix(self, od_matrix):
        for edge in self.od_graph.edges():
            s = edge[0]
            t = edge[1]
            self.od_graph[s][t]['demand'] = od_matrix[s,t]
        pass

    def check_network_connections(self):
        graph = self.graph.copy()

        for edge in graph.edges():
            if self.get_edge_attribute(edge,'capacity') == 1:
                graph.remove_edge(edge[0],edge[1])

        for edge in self.od_graph.edges():
            s = edge[0] # source
            t = edge[1] # target
            if not nx.has_path(graph,s,t):
                self.lost_trips[(s,t)] = self.od_graph[s][t]['demand']
                self.od_graph.remove_edge(s,t)

        for node in self.od_graph.nodes():
            if self.od_graph.degree(node) <= 2:
                mapping = {node:self.od_graph.node[node]['coordinates']}
                self.od_matrix = np.delete(np.delete(self.od_matrix,node,0),node,1)
                self.od_graph.remove_node(node)
                self.graph = nx.relabel_nodes(self.graph,mapping,copy=False)
        pass

    def process_network(self):
        nodes = self.od_graph.nodes()
        nodes.sort()
        mapping = {node:i+1 for i,node in enumerate(nodes)}
        self.graph = nx.relabel_nodes(self.graph,mapping,copy=False)
        self.od_graph = nx.relabel_nodes(self.od_graph,mapping,copy=False)

        k = nodes[-1] + 1
        mapping = {node:i+k for i,node in enumerate(list(set(self.graph.nodes())-set(self.od_graph.nodes())))}
        self.graph = nx.relabel_nodes(self.graph,mapping,copy=False)

        number_of_zones = self.od_graph.number_of_nodes()
        number_of_nodes = self.graph.number_of_nodes()
        first_thru_node = 1
        number_of_links = self.graph.number_of_edges()

        net_rows = ['<NUMBER OF ZONES> ' + str(number_of_zones),
                    '<NUMBER OF NODES> ' + str(number_of_nodes),
                    '<FIRST THRU NODE> ' + str(first_thru_node),
                    '<NUMBER OF LINKS> ' + str(number_of_links),
                    '<END OF METADATA> \n',
                    '~\t Init node \t Term node \t Capacity \t Length \t Free Flow Time \t B \t Power \t Speed limit \t Toll \t Link Type \t'
        ]

        for e in self.graph.edges(data=True):
            attr = e[2]
            line = [e[0],                                                # Init node
                    e[1],                                                # Term node
                    attr['capacity'],                                    # Capacity (veh/h)
                    attr['length']/self.unit_factor,                     # Length (km)
                    attr['length']/self.unit_factor/attr['speedlimit'],  # Free Flow Time (h)
                    self.alpha,                                          # B
                    self.beta,                                           # Power
                    attr['speedlimit'],                                  # Speed limit (km/h)
                    self.toll,                                           # Toll
                    attr['oneway']                                       # Link Type
            ]
            net_rows.append('\t'+'\t'.join(map(str,line))+'\t;')

        with open(self.temp_dir+'/'+self.temp_name+'_net.tntp','w') as net_file:
            net_file.write('\n'.join(net_rows))


        total_od_flow = np.sum(self.od_matrix)

        width_1 = len(str(number_of_zones))+2
        width_2 = len(str(int(np.max(self.od_matrix))))+5

        trips_rows = ['<NUMBER OF ZONES> '+str(number_of_zones),
                      '<TOTAL OD FLOW> '+str(total_od_flow),
                      '<END OF METADATA>'
        ]
        for i,o in enumerate(self.od_matrix):
            trips_rows.append('\nOrigin '+str(i+1))
            col = 0
            line = []
            for j,d in enumerate(o):
                line.append("{0:>{2}} :{1:>{3}.{digits}f};".format(j+1,d, width_1,width_2,digits=2))
            n = 5
            lines = [line[i:i + n] for i in range(0, len(line), n)]
            [trips_rows.append(''.join(l)) for l in lines]

        with open(self.temp_dir+'/'+self.temp_name+'_trips.tntp','w') as trips_file:
            trips_file.write('\n'.join(trips_rows))
        pass

    def ta_frank_wolf(self):
        args = [self.temp_folder,self.temp_name,str(self.threshold),str(self.n_iter_tm)]
        output = subprocess.check_output('julia ./pytras/ta.jl '+' '.join(args), shell=True).decode("utf-8").strip('\n')
        #output = subprocess.check_output('julia ta.jl '+' '.join(args), shell=True).decode("utf-8").strip('\n')
        values = [float(x) for x in output.split()]
        self.flow = values[0]
        self.traveltime = values[1]
        self.car_hours = values[2]
        self.car_distances = values[3]
        pass

    def get_traveltime(self):
        return [self.traveltime]

    def get_flow(self):
        return [self.flow]


    def get_car_hours(self):
        return [self.car_hours]

    def get_car_distances(self):
        return [self.car_distances]

    def get_lost_trips(self):
        return self.lost_trips

    def run(self):

        # print(nx.info(self.od_graph))
        # assign od matrix to od graph (if matrix is given)
        if self.od_matrix is not None:
            self.set_od_matix(self.od_matrix)

        # check network if every source and target can be reached
        self.check_network_connections()

        # # create input file for the simulation
        self.process_network()

        # run frank-wolf algorithm
        self.ta_frank_wolf()

        # clean folders
        os.system('rm -rf '+ self.temp_dir)

        pass

# =============================================================================
# eof
#
# Local Variables: 
# mode: python
# mode: linum
# End: 

 
