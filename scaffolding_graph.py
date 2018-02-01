#!/usr/bin/env python

'''
Created on Jan 31, 2018

@author: igorm
'''

import argparse
import sys
import os
import mmap
from copy import copy
from random import random
from collections import deque, Counter
import networkx as nx
import fasta.io as io
from math import sqrt, log
import numpy
import operator
import os.path
from Bio import SeqIO

class Link(object):
    def __init__(self, **args):
        for x, y in args.iteritems():
            self.__dict__[x] = y

    def __str__(self):
        return "distance: %s" % self.dist



class GraphConstructor(object):
    
    _settings = None
    
    def __init__(self):
        self._hits = {}
        self._sizes = {}
        self._dist = {}
        self._IGORgraph = nx.Graph() # the main scaffolding graph
        self._contigs= {}
        
    def set_settings(self, settings):
        self._settings = settings

    
    def scaffolding_graph(self):
        """Read mapping files chunk by chunk and feed the lines
        corresponding to valid read pairs to graph construction"""

        # read the contigs and determine their widths
        for record in SeqIO.parse(self._settings.get("contigs_fasta"), "fasta"):
            self._contigs[record.id] = {"seq": str(record.seq), "width": len(record.seq)}
            self._IGORgraph.add_node(record.id + "_1")
            self._IGORgraph.add_node(record.id + "_2")

        libraries = self._settings.get("libraries")

        for lib_id in libraries.keys(): # stub for multiple libraries
            lib = libraries[lib_id]
            sam1, sam2 = lib["sam1"], lib["sam2"]
            with open(sam1) as xx, open(sam2) as yy:
                with open(sam1, "r+b") as f:
                    kmap = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
                    fileiter1 = iter(kmap.readline, "")
                                       
                with open(sam2, "r+b") as f:
                    kmap = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
                    fileiter2 = iter(kmap.readline, "")
              
                fileiter1 = deque([xxx for xxx in fileiter1 if xxx[0] != "@"])
                fileiter2 = deque([xxx for xxx in fileiter2 if xxx[0] != "@"])

                while fileiter1 and fileiter2:
                    if len(fileiter1) <= 1 or len(fileiter2) <= 1:
                        break
                    reads11 = []
                    curread, curcounter1, curline = "", 0, ""
                    bbb = []
                    while True:
                        try:
                            read = fileiter1.popleft()
                        except Exception:
                            break
                        r = read.split()[0]
                        if r == curread:
                            curcounter1 += 1
                        else:
                            if curcounter1 == 0:
                                reads11.append((curread, curcounter1))
                                bbb.append(read)
                                curread, curcounter1 = r, 1
                                curline = read
                            else:
                                reads11.append((curread, curcounter1))
                                break
                    fileiter1.appendleft(read)
                    line1 = curline

                    reads22 = []
                    curread, curcounter2, curline = "", 0, ""
                    bbb = []
                    while True:
                        try:
                            read = fileiter2.popleft()
                        except Exception:
                            break
                        r = read.split()[0]
                        if r == curread:
                            curcounter2 += 1
                        else:
                            if curcounter2 == 0:
                                reads22.append((curread, curcounter2))
                                bbb.append(read)
                                curread, curcounter2 = r, 1
                                curline = read
                            else:
                                reads22.append((curread, curcounter2))
                                break
                    fileiter2.appendleft(read)
                    line2 = curline

                    if curcounter1 > 1 or curcounter2 > 1:
                        continue
                    else:
                        line1, line2 = line1.split(), line2.split()

                        if line1[2] == "*" or line2[2] == "*": # we use this with bowtie2
                            pass
                        else:
                            if line1[2] != line2[2]: # discard those who has mapped to the same contig
				mismatches1, mismatches2 = line1[13].split(":")[-1], line2[13].split(":")[-1]
				if int(mismatches1) > 1000 and int(mismatches2) > 1000:
				    continue
				else:
                                    self._paired_read_to_graph(line1, line2, lib_id)
                            else:
                                pass
        self._build_graph()
        # add the information to the scaffolding graph

        for node, data in self._contigs.items():
            self._IGORgraph.node[node + "_1"]["width"] = data["width"]
            self._IGORgraph.node[node + "_2"]["width"] = data["width"]
            self._IGORgraph.node[node + "_1"]["seq"] = data["seq"]
            self._IGORgraph.node[node + "_2"]["seq"] = data["seq"]

        #for node in self._IGORgraph.nodes():
        #    print node, node[:-2], "AAAAAAAAAA"
        #    self._IGORgraph.node[node]["width"] = self._contigs[node[:-2]]["width"]
        #    self._IGORgraph.node[node]["seq"] = self._contigs[node[:-2]]["seq"]
        return self._IGORgraph


    def _paired_read_to_graph(self, line1, line2, lib_id):
        """
        @param line1: line corresponding to first fragment
        @param line2: line corresponding to second fragment
        @param lib_id: id of the library from which pair comes
        """
        # insert size and standard deviation
        library = self._settings.get("libraries").get(lib_id)
        ins_size, std_dev, pair_mode = library["ins"], library["std"], library["pm"]
 
        # first leg of the read
        oflag1, rname1, lpos1 = line1[1], line1[2], int(line1[3])

        width1 = self._contigs[rname1]["width"]
        rpos1 = lpos1 + len(line1[9])
        
        # second leg of the read
        oflag2, rname2, lpos2 = line2[1], line2[2], int(line2[3])
        width2 = self._contigs[rname2]["width"]
        rpos2 = lpos2 + len(line2[9])

        op, oq = self._get_orientation(oflag1, oflag2, pair_mode)
        orients = (op, oq)
               
        if orients == (0, 0):
            distance = ins_size - (width1 - lpos1) - rpos2
            edge = (rname1 + "_1", rname2 + "_2")
        elif orients == (1, 1):
            distance = ins_size - (width2 - lpos2) - rpos1
            edge = (rname1 + "_2", rname2 + "_1")
        elif orients == (0, 1):
            distance = ins_size - (width1 - lpos1) - (width2 - lpos2)
            edge = (rname1 + "_1", rname2 + "_1")
        elif orients == (1, 0):
            distance = ins_size - rpos1 - rpos2
            edge = (rname1 + "_2", rname2 + "_2")

        if -std_dev * 3 <= distance <= ins_size + 3 * std_dev:
            pair = tuple(sorted(edge))
            if rname1 < rname2:
                pos1, pos2 = lpos1, lpos2
            else:
                pos1, pos2 = lpos2, lpos1
                line1, line2 = line2, line1 # used to calculate entropy
   
            para = self._dist.get(pair, [])
            link_dict = {"pos1": pos1, "pos2": pos2, "dist": distance,
                            "ins": ins_size, "std": std_dev}
            link = Link(**link_dict)
            para.append(link)
	    self._dist[pair] = para
                
    def _get_orientation(self, oflag1, oflag2, pair_mode):
        """Gets orientation from SAM object for pairs"""
        if oflag1 == "0":
            o1 = 0
        else:
            o1 = 1
        if oflag2 == "0":
            o2 = 0
        else:
            o2 = 1
        if pair_mode == "fr":
            o2 = 1 - o2
        if pair_mode == "rf":
            o1 = 1 - o1
        return o1, o2


    def _build_graph(self):
        """build the graph on the data obtained
        from parsing the libraries"""
        """
        @note: perform bundling step as described in the paper
        "The Greedy Path-Merging Algorithm for Contig Scaffolding"
        by Huson, Reinert and Myers, 2002
        """
        nodedict = {}
        for node1, node2 in self._dist:
            ar = self._dist.get((node1, node2), [])
            ar1 = nodedict.get(node1, {})
            ar2 = nodedict.get(node2, {})
          
            for link in ar:
                pos1 = link.pos1
                aaa = ar1.get(pos1, [])
                aaa.append(node2)
                ar1[pos1] = aaa

                pos2 = link.pos2
                bbb = ar2.get(pos2, [])
                bbb.append(node1)
                ar2[pos2] = bbb

            nodedict[node1] = ar1
            nodedict[node2] = ar2
                 
        for node1, node2 in self._dist:

            links = self._dist[(node1, node2)]
            # ------------------------
            linkdists = {}
            for link in links:
                linkdists[link] = link.dist

            distslinks = {} # identify the link by dist
            for link in links:
                distslinks[link.dist] = link
            sorted_x = sorted(linkdists.items(), key=operator.itemgetter(1))
            median = len(sorted_x) / 2
            p = 0
            q = 0
            newmean = 0
            newstd = 0
            size = 0
            
            lowerbound = sorted_x[median][1] - 3 * distslinks[sorted_x[median][1]].std  #std_dev
            upperbound = sorted_x[median][1] + 3 * distslinks[sorted_x[median][1]].std
            
            iterator = 0
            while iterator < len(linkdists) and sorted_x[iterator][1] < upperbound:
                std_dev = distslinks[sorted_x[iterator][1]].std
                if lowerbound < sorted_x[iterator][1] < upperbound:
                    size += 1
                    p += sorted_x[iterator][1] * 1.0 / pow(std_dev, 2)
                    q += 1.0 / pow(std_dev, 2)
                    iterator += 1
                else:
                    iterator += 1

            newmean = p / q
            self._IGORgraph.add_edge(node1, node2, weight=size, dist=newmean)


class Settings(object):
    _instance = None
    
    _settings = {}
    _defaults = copy(_settings)
    
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Settings, cls).__new__(
                                cls, *args, **kwargs)
        return cls._instance
    
    def update(self, new_settings):
        self._settings.update(new_settings)
        
    def get(self, key):
        return self._settings.get(key)
    
    def set(self, key, value):
        self._settings[key] = value
           
    def __str__(self):
        return "\n".join(["%s: %s" % x for x in self._settings.iteritems()])

    def iteritems(self):
        it = []
        for s, v in self._settings.iteritems():
            if s not in self._defaults and s != "logger":
                it.append((s, v))
        return it


def parse_args():
    main_p = argparse.ArgumentParser(description='BATISCAF scaffolding graph construction helper script. Produces a .graphml file')
    main_p.add_argument('-o', dest='output_graphml', help='output graphml file')
    main_p.add_argument('-c', dest='contigs_fasta', help='fasta file with contigs')
    main_p.add_argument('-m1', dest='mappings1', required=True, help='comma separated list of .sam files (first read in the read pair)')
    main_p.add_argument('-m2', dest='mappings2', required=True, help='comma separated list of .sam files (second read in the read pair)')  
    main_p.add_argument('-i', dest='ins_size', required=True, help='insert sizes (comma separated values)')
    main_p.add_argument('-p', dest='pair_mode', required=True, help='pair modes (fr - innie style -> <-, rf - outtie style <- ->) (comma separated values)')
    main_p.add_argument('-s', dest='std_dev', required=True, help='libraries standard deviations (comma separated values)')
    return vars(main_p.parse_args())  


def prepare_libraries(args_dict):
    """If there are more then one single
    library, pair mappings as necessary"""
    mappings1 = args_dict["mappings1"]
    args_dict.pop("mappings1", None)
    mappings2 = args_dict["mappings2"]
    args_dict.pop("mappings2", None)
    insert_sizes = args_dict["ins_size"]
    std_devs = args_dict["std_dev"]
    pair_modes = args_dict["pair_mode"]
    try:
        mappings1 = mappings1.split(",")
        mappings2 = mappings2.split(",")
        insert_sizes = insert_sizes.split(",")
        std_devs = std_devs.split(",")
        pair_modes = pair_modes.split(",")
    except:
        print "You have issues in defining mapping files"
        sys.exit()
    if not (len(mappings1) == len(mappings2) == len(insert_sizes) == len(std_devs) == len(pair_modes)):
        print "You missed some arguments, please check them carefully"
        sys.exit()
    libraries = {}
    libcount = 1
    for sam1, sam2, ins, std, pm in \
                    zip(mappings1, mappings2, insert_sizes, std_devs, pair_modes):
        libraries[libcount] = {"ins": int(ins), "std": int(std), "sam1": sam1, "sam2": sam2, "pm": pm}
        libcount += 1
    args_dict["libraries"] = libraries


if __name__ == '__main__':
    # Parse arguments
    args = parse_args()
    prepare_libraries(args)
    settings = Settings()
    settings.update(args)
    gcon = GraphConstructor()
    gcon.set_settings(settings)
    scaffolding_graph = gcon.scaffolding_graph()
    nx.write_graphml(scaffolding_graph, settings.get("output_graphml"))

