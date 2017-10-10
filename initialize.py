#!/usr/bin/python -tt
# -*- coding: utf-8 -*-
# =============================================================================
# File      : initialize.py 
# Creation  : 26 Oct 2015
# Time-stamp: <Die 2017-10-10 11:20 juergen>
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


def write_shp(G, outdir, filename):
    """Writes a networkx.DiGraph to two shapefiles, edges and nodes.
    Nodes and edges are expected to have a Well Known Binary (Wkb) or
    Well Known Text (Wkt) key in order to generate geometries. Also
    acceptable are nodes with a numeric tuple key (x,y).

    "The Esri Shapefile or simply a shapefile is a popular geospatial vector
    data format for geographic information systems software [1]_."

    Parameters
    ----------
    outdir : directory path
       Output directory for the two shapefiles.

    Returns
    -------
    None

    Examples
    --------
    nx.write_shp(digraph, '/shapefiles') # doctest +SKIP

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Shapefile
    """
    try:
        from osgeo import ogr
    except ImportError:
        raise ImportError("write_shp requires OGR: http://www.gdal.org/")
    # easier to debug in python if ogr throws exceptions
    ogr.UseExceptions()

    def netgeometry(key, data):
        if 'Wkb' in data:
            geom = ogr.CreateGeometryFromWkb(data['Wkb'])
        elif 'Wkt' in data:
            geom = ogr.CreateGeometryFromWkt(data['Wkt'])
        elif type(key[0]).__name__ == 'tuple':  # edge keys are packed tuples
            geom = ogr.Geometry(ogr.wkbLineString)
            _from, _to = key[0], key[1]
            try:
                geom.SetPoint(0, *_from)
                geom.SetPoint(1, *_to)
            except TypeError:
                # assume user used tuple of int and choked ogr
                _ffrom = [float(x) for x in _from]
                _fto = [float(x) for x in _to]
                geom.SetPoint(0, *_ffrom)
                geom.SetPoint(1, *_fto)
        else:
            geom = ogr.Geometry(ogr.wkbPoint)
            try:
                geom.SetPoint(0, *key)
            except TypeError:
                # assume user used tuple of int and choked ogr
                fkey = [float(x) for x in key]
                geom.SetPoint(0, *fkey)

        return geom

    # Create_feature with new optional attributes arg (should be dict type)
    def create_feature(geometry, lyr, attributes=None):
        feature = ogr.Feature(lyr.GetLayerDefn())
        feature.SetGeometry(g)
        if attributes != None:
            # Loop through attributes, assigning data to each field
            for field, data in attributes.items():
                feature.SetField(field, data)
        lyr.CreateFeature(feature)
        feature.Destroy()

    drv = ogr.GetDriverByName("ESRI Shapefile")
    shpdir = drv.CreateDataSource(outdir)
    # delete pre-existing output first otherwise ogr chokes
    try:
        shpdir.DeleteLayer(filename + '_n')
    except:
        pass
    nodes = shpdir.CreateLayer(filename + '_n', None, ogr.wkbPoint)
    for n in G:
        data = G.node[n]
        g = netgeometry(n, data)
        create_feature(g, nodes)
    try:
        shpdir.DeleteLayer(filename + '_e')
    except:
        pass
    edges = shpdir.CreateLayer(filename + '_e', None, ogr.wkbLineString)

    # New edge attribute write support merged into edge loop
    fields = {}      # storage for field names and their data types
    attributes = {}  # storage for attribute data (indexed by field names)

    # Conversion dict between python and ogr types
    OGRTypes = {int: ogr.OFTInteger, str: ogr.OFTString, float: ogr.OFTReal}

    # Edge loop
    for e in G.edges(data=True):
        data = G.get_edge_data(*e)
        g = netgeometry(e, data)
        # Loop through attribute data in edges
        for key, data in e[2].items():
            # Reject spatial data not required for attribute table
            if (key != 'Json' and key != 'Wkt' and key != 'Wkb'
                and key != 'ShpName'):
                  # For all edges check/add field and data type to fields dict
                    if key not in fields:
                  # Field not in previous edges so add to dict
                        if type(data) in OGRTypes:
                            fields[key] = OGRTypes[type(data)]
                        else:
                            # Data type not supported, default to string (char 80)
                            fields[key] = ogr.OFTString
                        # Create the new field
                        newfield = ogr.FieldDefn(key, fields[key])
                        edges.CreateField(newfield)
                        # Store the data from new field to dict for CreateLayer()
                        attributes[key] = data
                    else:
                     # Field already exists, add data to dict for CreateLayer()
                        attributes[key] = data
        # Create the feature with, passing new attribute data
        create_feature(g, edges, attributes)

    nodes, edges = None, None


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

 
