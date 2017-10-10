#!/usr/bin/python -tt
# -*- coding: utf-8 -*-
# =============================================================================
# Time-stamp: <Die 2017-10-10 13:02 juergen>
# File      : msa.py 
# Creation  : 08 Oct 2015
#
# Copyright (c) 2015 Jürgen Hackl <hackl@ibi.baug.ethz.ch>
#               http://www.ibi.ethz.ch
# $Id$ 
#
# Description : A simple traffic model based on MSA and BPR
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
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import itertools

class TrafficAssignment(object):
    """Traffic model to estimate the flow, velocity and travel time on a road
    network.

    Calculating the flow and other parameters for a given network using the
    'Method of Successive Averages' [1]_.

    Iterative algorithms were developed, at least partially, to overcome the
    problem of allocating too much traffic to low-capacity links. In an
    iterative assignment algorithm the ‘current’ flow on a link is calculated
    as a linear combination of the current flow on the previous iteration and
    an auxiliary flow resulting from an all-or-nothing assignment in the
    present iteration. The algorithm can be described by the following steps:

    1. Select a suitable initial set of current link costs, usually free-flow
       travel times. Initialise all flows :math:`V_a = 0`; make :math:`n = 0`.
    2. Build the set of minimum cost trees with the current costs; make
       :math:`n = n + 1`.
    3. Load the whole of the matrix :math:`T` all-or-nothing to these trees
       obtaining a set of auxiliary flows :math:`F_a`.
    4. Calculate the current flows as:
    
           .. math::
              V_a^n = (1 − \\phi)V_a^{n-1} + \\phi F_a

    5. Calculate a new set of current link costs based on the flows
       :math:`V_a^n` . If the flows (or current link costs) have not changed
       significantly in two consecutive iterations, stop; otherwise proceed to
       step 2. Another, less good but quite common, criterion for stopping is
       simply to fix the maximum number of iterations; should be calculated in
       this case as well to know how close the solution is to Wardrop’s
       equilibrium.

    Iterative assignment algorithms differ in the method used to give a value
    to :math:`\\phi`. A simple rule is to make it constant, for example
    :math:`\\phi = 0.5`. A much better approach due to Smock (1962), is to
    make :math:`\\phi = 1/n`.
    
    Parameters
    ----------
    graph : NetworkX graph
       Graph which represents the road network with the actual locations of
       the source and sink nodes and the connections to the road network.

    od_graph : NetworkX graph
       Graph contains a network between source and sink notes. The demand of
       each edge should be already assigned.

    od_matrix : numpy matrix
       Matrix with origin destination demands. If the od_graph does not
       contain this information, it can be optional loded as separate matrix.

    Examples
    --------
    >>> traffic = TrafficModel(graph,od_graph) # doctest: +SKIP
    >>> traffic.run()
    >>> graph_results = traffic.get_graph()

    References
    ----------
    .. [1] Juan de Dios Ortuzar and Luis G. Willumsen (2011) Modelling Transport, Fourth Edition. John Wiley and Sons.

    """

    def __init__(self, graph, od_graph,od_matrix=None):

        self.graph = graph
        self.od_graph = od_graph
        self.od_matrix = od_matrix
        self.lost_trips = {}
        self.cut_links = []

        self.print_flag = True
        self.unit_factor = 1000
        self.alpha = 0.15#1
        self.beta = 4#2
        self.threshold = 0.001#1e-5
        self.large = 1e+14
        # number of iterations for the corrective factors
        self.n_iter_tm = 500#1000#1001

        pass

    def get_edge_attribute(self,edge,attribute):
        return self.graph[edge[0]][edge[1]][attribute]

    def set_edge_attribute(self,edge,attribute,value):
        self.graph[edge[0]][edge[1]][attribute] = value
        pass

    def calculate_initial_traveltime(self):
        for edge in self.graph.edges():
            initial_traveltime = self.get_edge_attribute(edge,'length') / self.unit_factor / self.get_edge_attribute(edge,'speedlimit')
            self.set_edge_attribute(edge,'t_0',initial_traveltime)
        pass

    def set_initial_traveltimes(self):
        for edge in self.graph.edges():
            self.set_edge_attribute(edge,'t_k',self.get_edge_attribute(edge,'t_0'))
            self.set_edge_attribute(edge,'t_h',0)
        pass

    def set_initial_flow(self):
        for edge in self.graph.edges():
            self.set_edge_attribute(edge,'flow',0)
        pass

    def set_initial_help_flow(self):
        for edge in self.graph.edges():
            self.set_edge_attribute(edge,'help_flow',0)
        pass

    def set_od_matix(self, od_matrix):
        for edge in self.od_graph.edges():
            s = edge[0]
            t = edge[1]
            self.od_graph[s][t]['demand'] = od_matrix[s,t]
        pass

    def set_help_traveltime(self):
        for edge in self.graph.edges():
            self.set_edge_attribute(edge,'t_h',self.get_edge_attribute(edge,'t_k'))
        pass

    def calculate_auxiliary_flows(self):
        self.set_initial_help_flow()

        for edge in self.od_graph.edges():
            s = edge[0] # source
            t = edge[1] # target
            # print(' new  ')
            # print(s,t)
            sp = nx.shortest_path(self.graph, source=s, target=t, weight='t_k')
            for i in range(len(sp)-1):
                u = sp[i]
                v = sp[i+1]
                self.graph[u][v]['help_flow'] += self.od_graph[s][t]['demand'] # += correct ???
        pass

    def calculate_flows(self,phi):
        for edge in self.graph.edges():
            if edge[0] != edge[1]:
                flow = (1-phi) * self.get_edge_attribute(edge,'flow') + phi * self.get_edge_attribute(edge,'help_flow')
            else: # internal flow
                if (isinstance(edge[0],int) and isinstance(edge[1],int)): # no self loops of roads
                    flow = self.od_graph[edge[0]][edge[1]]['demand']
            self.set_edge_attribute(edge,'flow',flow)
        pass

    def calculate_traveltime(self):
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge,'capacity') < self.threshold:
                traveltime = self.large
            else:
                traveltime = self.get_edge_attribute(edge,'t_0') * (1 + self.alpha * (self.get_edge_attribute(edge,'flow') / self.get_edge_attribute(edge,'capacity'))**self.beta)
            self.set_edge_attribute(edge,'t_k',traveltime)
        pass

    def stopping_criteria(self):
        t_k_list = [value for key, value in nx.get_edge_attributes(self.graph,'t_k').items()]
        t_h_list = [value for key, value in nx.get_edge_attributes(self.graph,'t_h').items()]
        t_d_list = list(np.abs(np.array(t_k_list) - np.array(t_h_list)))
        return all(i <= self.threshold for i in t_d_list)

    def check_network_connections(self):
        graph = self.graph.copy()
        for edge in graph.edges():
            if self.get_edge_attribute(edge,'capacity') == 0:
                graph.remove_edge(edge[0],edge[1])

        cut_links = []
        for edge in self.od_graph.edges():
            s = edge[0] # source
            t = edge[1] # target
            if not nx.has_path(graph,s,t):
                self.lost_trips[(s,t)] = self.od_graph[s][t]['demand']
                self.od_graph.remove_edge(s,t)

                cut_value, partition = nx.minimum_cut(self.graph,s,t)
                reachable, non_reachable = partition

                cutset = set()
                for u, nbrs in ((n, self.graph[n]) for n in reachable):
                    cutset.update((u, v) for v in nbrs if v in non_reachable)
                cut_links.extend(list(cutset))

        for edge in list(set(cut_links)):
            if self.graph[edge[0]][edge[1]]['capacity'] == 0:
                self.cut_links.append(edge)

        pass

    def run(self):

        # print(nx.info(self.graph))
        # assign od matrix to od graph (if matrix is given)
        if self.od_matrix is not None:
            self.set_od_matix(self.od_matrix)
        # calculate traveltime at t=0
        self.calculate_initial_traveltime()

        # set traveltime equal to initial traveltime
        self.set_initial_traveltimes()

        # set initial flow = 0
        self.set_initial_flow()

        #set initial help flow = 0
        self.set_initial_help_flow()

        # check network if every source and target can be reached
        self.check_network_connections()
    
        for i in range(1,self.n_iter_tm):
            phi = 1/i
            # calculate auxiliary flows
            self.calculate_auxiliary_flows()

            # calculating the flow using auxiliary flow
            self.calculate_flows(phi)

            # save old traveltimes
            self.set_help_traveltime()

            # calculate new traveltime
            self.calculate_traveltime()

            # stopping criteria
            if self.stopping_criteria():
                break
        pass


    def get_traveltime(self):
        traveltime = []
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge,'type') != 'Misc':
                if self.get_edge_attribute(edge,'capacity') != 0:
                    traveltime.append(self.get_edge_attribute(edge,'t_k'))
        return traveltime

    def get_flow(self):
        flow = []
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge,'type') != 'Misc':
                if self.get_edge_attribute(edge,'capacity') != 0:
                    flow.append(self.get_edge_attribute(edge,'flow'))
        return flow


    def get_car_hours(self):
        car_hours = []
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge,'type') != 'Misc':
                if self.get_edge_attribute(edge,'capacity') != 0:
                    car_hour = self.get_edge_attribute(edge,'t_k') * self.get_edge_attribute(edge,'flow')
                    car_hours.append(car_hour)
        return car_hours

    def get_car_distances(self):
        """vehicle-kilometre"""
        car_distances = []
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge,'type') != 'Misc':
                if self.get_edge_attribute(edge,'capacity') != 0:
                    car_distance = self.get_edge_attribute(edge,'length') / self.unit_factor * self.get_edge_attribute(edge,'flow')
                    car_distances.append(car_distance)
        return car_distances

    def get_lost_trips(self):
        return self.lost_trips

    def get_cut_links(self):
        return self.cut_links

    def print_results(self):
        for edge in self.graph.edges():
            if self.get_edge_attribute(edge,'type') != 'Misc':
                name = self.get_edge_attribute(edge,'name')
                flow = self.get_edge_attribute(edge,'flow')
                initialtraveltime = self.get_edge_attribute(edge,'t_0')
                traveltime = self.get_edge_attribute(edge,'t_k')
                print('t_0 vs t_k (flow) ' + str(name) + ': ' +str(round(initialtraveltime,2)) +', '+ str(round(traveltime,2)) +' ('+str(round(flow,0))+')')
                # print('initial traveltime ' + str(name) + ': ' +str(initialtraveltime))
                # print('traveltime ' + str(name) + ': ' +str(traveltime))
#                print('traffic flow ' + str(name) + ': ' +str(flow))

    def get_graphs(self):
        return self.graph, self.od_graph

    def get_od_matrix(self):
        return self.od_matrix

    def get_graph(self):
        return self.graph

# =============================================================================
# eof
#
# Local Variables: 
# mode: python
# mode: linum
# End: 
