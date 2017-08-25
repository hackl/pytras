#!/usr/bin/python -tt
# -*- coding: utf-8 -*-
# =============================================================================
# File      : initialize.py 
# Creation  : 26 Oct 2015
# Time-stamp: <Mit 2017-06-14 15:22 juergen>
#
# Copyright (c) 2015 Jürgen Hackl <hackl@ibi.baug.ethz.ch>
#               http://www.ibi.ethz.ch
# $Id$ 
#
# Description : Initialize Graphs and OD
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
import itertools

def read_shp(path):
    """Generates a networkx.DiGraph from shapefiles. Point geometries are
    translated into nodes, lines into edges. Coordinate tuples are used as
    keys. Attributes are preserved, line geometries are simplified into start
    and end coordinates. Accepts a single shapefile or directory of many
    shapefiles.

    "The Esri Shapefile or simply a shapefile is a popular geospatial vector
    data format for geographic information systems software [1]_."

    Parameters
    ----------
    path : file or string
       File, directory, or filename to read.

    Returns
    -------
    G : NetworkX graph

    Examples
    --------
    >>> G=nx.read_shp('test.shp') # doctest: +SKIP

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Shapefile
    """
    try:
        from osgeo import ogr
    except ImportError:
        raise ImportError("read_shp requires OGR: http://www.gdal.org/")

    if not isinstance(path, str):
        return

    net = nx.DiGraph()
    shp = ogr.Open(path)
    for lyr in shp:
        fields = [x.GetName() for x in lyr.schema]
        for f in lyr:
            flddata = [f.GetField(f.GetFieldIndex(x)) for x in fields]
            g = f.geometry()
            attributes = dict(zip(fields, flddata))
            attributes["ShpName"] = lyr.GetName()
            if g.GetGeometryType() == 1:  # point
                residual = 2
                xy = (round(g.GetPoint_2D(0)[0],residual),round(g.GetPoint_2D(0)[1],residual))
                name = attributes['name']
                attributes['coordinates'] = (xy)
                net.add_node(name, attributes)
            if g.GetGeometryType() == 2:  # linestring
                attributes["Wkb"] = g.ExportToWkb()
                attributes["Wkt"] = g.ExportToWkt()
                attributes["Json"] = g.ExportToJson()
                last = g.GetPointCount() - 1
                residual = 2
                # calculate new true length of a road segment
                # attributes['length'] = g.Length()
                # create notes and round cooridinates
                u_node = (round(g.GetPoint_2D(0)[0],residual),round(g.GetPoint_2D(0)[1],residual))
                v_node = (round(g.GetPoint_2D(last)[0],residual),round(g.GetPoint_2D(last)[1],residual))
                net.add_edge(u_node, v_node, attributes)

    return net
def remove_oneways(graph):
    for edge in graph.edges():
        u = edge[0]
        v = edge[1]
        if graph[u][v]['oneway'] == 1:
            graph.remove_edge(u,v)
    return graph

def complete_graph_from_list(L, G=None, selfloops=False, create_using=None):
    if G == None:
        G = nx.empty_graph(len(L),create_using)
    if len(L)>1:
        if G.is_directed():
            edges = itertools.permutations(L,2)
        else:
            edges = itertools.combinations(L,2)
        G.add_edges_from(edges,demand=0)
        if selfloops:
            edges = [(i,i)for i in L]
            G.add_edges_from(edges,demand=0)
    return G

def convert_matrix_to_graph_dict(matrix):
    d = dict()
    for i in range(np.shape(matrix)[0]):
        for j in range(np.shape(matrix)[1]):
            d[(i,j)] = matrix[i,j]
    return d


def create_od_graph(filename,od_matrix=None):
    # load shp file with node locations
    graph = read_shp(filename)

    # create a complete graph for the OD nodes (with self loops)
    complete_graph_from_list(graph.nodes(),G=graph,selfloops=True)

    if od_matrix != None:
        # Transform matrix to dict
        demand = convert_matrix_to_graph_dict(od_matrix)

        # add demand values to graph
        nx.set_edge_attributes(graph, 'demand', demand)
    return graph

def create_network_graph(road_graph,od_graph,con_edges):
    # add road graph structure to graph
    graph = road_graph.copy()

    # add nodes of the od to the graph
    graph.add_nodes_from(od_graph)

    # add connections from od nodes to road nodes
    for edge in con_edges.edges():
        u = edge[0]
        v = edge[1]
        w = con_edges[u][v]['id']
        graph.add_edge(w,v,name='a_'+str(w),type='Misc',capacity=10000, length=0,speedlimit=500,oneway=9,damage=0)
        graph.add_edge(v,w,name='b_'+str(w),type='Misc',capacity=10000, length=0,speedlimit=500,oneway=0,damage=0)

    # add self loops for internal traveltime
    for node in od_graph.nodes():
        graph.add_edge(node,node,name='i_'+str(node),type='Internal',capacity=1000, length=3750,speedlimit=45,oneway=8,damage=0)
    return graph


# =============================================================================
# eof
#
# Local Variables: 
# mode: python
# End: 

 
