#!/usr/bin/python
# coding=utf-8

from pylab import *

from copy import *

import commands
import subprocess

import time

import os, sys

import bz2

import muller

# Muller's ratchet in finite population, genes have "good" and "bad" states (1 or 0),
# fitness is (1 + fb) ** E, where E is number of good genes

# N - number of organisms
# G - number of genes
# f - fitness of gene, F - fitness of organism (additive or multiplicative)
# M - mutation frequency
# B - probability of beneficial mutation (all other mutations are deleterious)
# fb - fitness benefit from beneficial mutation
# T - frequency of horizontal gene transfer (transformation)
# C - cost of horizontal gene transfer
# Tmut - mutation of T
# Ttransform - whether T and M gene will be transfered
# 
# Mmut - mutation of M
# 
# 

default_params = {"steps" : 200, "N" : 100, "G" : 100, "M" : 0.015, "B" : 0.1, "fb" : 0.05, "T" : 0., "Tmut" : 0., "Mmut" : 0., "Ttransform" : 1., "C" : 0., "X" : 1, "even": True, "constantX": True, "Binitial" : -1., "interval" : 1, "seed" : -1, "verbose" : False}

# extra_param_names = ["steps", "interval", "verbose"] # params that are not params of the model so should not be passed into model initialization

swigworld_param_names = ["N", "G", "B", "fb", "M", "Mmut", "T", "Tmut", "Ttransform", "C", "X", "even", "constantX", "Binitial", "seed"] # parameters for C++/swig World object initialization - unfortunately swig does not mention parameter names

exe_param_names = ["N", "G", "B", "fb", "M", "Mmut", "T", "Tmut", "Ttransform", "C", "Binitial", "interval", "seed"] # parameters needed for running c++ executable

def STATVARS(*names):
    o = []
    for name in names:
        o.extend([name + "avg", name + "std", name + "min", name + "max"])
    return o

stat_names = ["time"] + STATVARS("E", "EE", "X", "F", "M", "T", "EG") + ["Tplus"]

# stat_names = ["time", "Eavg", "Estd", "Emin", "Emax", "Favg", "Fstd", "Fmin", "Fmax", "Mavg", "Mstd", "Mmin", "Mmax", "Tavg", "Tstd", "Tmin", "Tmax", "Tplus", "EGavg", "EGstd", "EGmin", "EGmax"]


class population_exe:
    
    def __init__(self, **kwargs):
        
        # filling attributes at once instead of mentioning each one
        # if parameter is not present in arguments then it will be taken from default_params
        self.params = copy(default_params)
        self.params.update(kwargs)
        for k, v in self.params.iteritems():
            setattr(self, k, v)
        
        # self.paramline = " " + " ".join(map(str, (self.N, self.G, self.B, self.fb, self.M, self.Mmut, self.T, self.Tmut, self.Ttransform, self.C, self.Binitial, self.interval, self.seed))) # parameters needed for running c++ executable
        self.paramline = ""
        for p in exe_param_names:
            self.paramline += (" " + str(self.params[p]))
        
        # if seed < 0:
            # pass
    
    def readstat(self, text):
        # text is content of statistics file or output of modeling program
        # TODO: generally move to csv, scipy.io or another module (also use bzip2)
        #reads = [[] for i in xrange(len(self.stat))]
        s = text.splitlines()
        for i in xrange(len(s)):
            if s[i] == "Statistics begin":
                begin = i
                break
        
        # parsing header representing statistic variable names
        self.stat = {}
        self.stat_names = []
        for statname in s[begin+1].split():
            self.stat[statname] = []
            self.stat_names.append(statname)
            
        for line in s[begin+2:]:
            s = line.split()
            for i in xrange(len(s)):
                self.stat[self.stat_names[i]].append(float(s[i]))
        
        for k, v in self.stat.iteritems():
            self.stat[k] = array(v) # for memory economy!
        
            
    def run(self, steps = None):
        if steps:
            self.steps = steps
        # commandline = 
        # print self.paramline
        command = "nice -n 19 ./muller "
        if sys.platform == "win32":
            command = "muller.exe "
            self.output = subprocess.check_output(command + str(self.steps) + " " + self.paramline)
        else:
            command = "nice -n 19 ./muller "
            self.output = commands.getoutput(command + str(self.steps) + " " + self.paramline)
        if self.verbose:
            print self.output
        self.readstat(self.output)
    
    def save(self, name = "trajectories"):
        if not hasattr(self, "output"):
            print "There are no results, run calculation first!"
            return
        if os.path.isfile(name):
            print "file", name, "exists! delete it or save in another place."
            return
        if not os.path.exists(name):
            os.mkdir(name)
        f = open(name + "/" + time.strftime("%Y-%m-%d_%H-%M-%S") + self.paramline + ".txt.bz2", "w")
        f.write(bz2.compress(self.output))
        f.close()

class population_swig(population_exe):
    
    def __init__(self, **kwargs):
        # super(type(self), self).__init__(**kwargs)
        population_exe.__init__(self, **kwargs)        
        
        if self.Binitial < 0:
            self.params["Binitial"] = self.B
            self.Binitial = self.B
        
        if self.seed < 0: # then random seed is generated by numpy random
            seed = randint(2**32)
            self.params["seed"] = seed
            self.seed = seed
        
        self.model_params = []
        for p in swigworld_param_names:
            self.model_params.append(self.params[p])
        self.model = muller.World(*self.model_params)
        
        # statistics initialization
        self.stat = {s: [] for s in stat_names}
        self.stat_names = stat_names

        
    def writestat(self):
        "returns string representing statistics table"
        # self.params["steps"] = len(self.stat["time"])
        s = "Model of Muller's ratchet https://github.com/dining-philosopher/muller.git\n"
        s += str(self.params) + "\n"
        s += "Statistics begin\n"
        s += " ".join(self.stat_names) + "\n"
        for i in xrange(len(self.stat["time"])):
            s += " ".join([str(self.stat[n][i]) for n in self.stat_names]) + "\n"
            # s +=  + "\n"
        return s
    
    def append_stat(self):
        self.model.calc_stat()
        for k, v in self.stat.iteritems():
            # self.stat[k].append(getattr(self.model, k))
            v.append(getattr(self.model, k))
        
    
    def run(self, steps = None):
        if steps:
            self.params["steps"] = steps + self.model.time
            self.steps = steps
        for i in xrange(self.steps - self.model.time):
            self.model.step()
            if self.model.time % self.interval == 0:
                self.append_stat()
     


def readstat(text):
    # text is content of statistics file or output of modeling program
    # reads parameters as dictionary!
    # returns params, stats
    # TODO: generally move to csv, scipy.io or another module (also use bzip2)
    s = text.splitlines()
    params = {}
    for i in xrange(len(s)):
        if s[i][0] == "{":
            params = eval(s[i], {}, {})
        if s[i] == "Statistics begin":
            begin = i
            break
    
    # parsing header representing statistic variable names
    stat = {}
    stat_names = []
    for statname in s[begin+1].split():
        stat[statname] = []
        stat_names.append(statname)
        
    for line in s[begin+2:]:
        s = line.split()
        for i in xrange(len(s)):
            stat[stat_names[i]].append(float(s[i]))
    
    for k, v in stat.iteritems():
        stat[k] = array(v) # for memory economy!
    return params, stat


class Cache:
    
    def __init__(self, datadir = "trajectories", indexname = "dataindex.txt"):
        self.dataindex = []
        self.datadir = datadir
        self.indexname = indexname
        if os.path.isfile(indexname):
            self.readindex()
        else:
            f = open(indexname, "w")
            f.close()
    
    def readindex(self):
        self.dataindex = [] # array of dictionaries {"params" : {..}, "stat_names": {..}, "file": "smth"}
        f = open(self.indexname)
        s = f.read().splitlines()
        for i in xrange(len(s)):
            if s[i][0] == "{":
                self.dataindex.append(eval(s[i], {}, {}))
        f.close()
        
    def saveindex(self):
        f = open(self.indexname, "w")
        for d in self.dataindex:
            f.write(str(d) + "\n")
        f.close()
    
    def save(self, pop):
        "saves trajectory file and registers it in index file"
        if os.path.isfile(self.datadir):
            print "file", self.datadir, "exists! delete it or save in another place."
            return
        if not os.path.exists(self.datadir):
            os.mkdir(self.datadir)
        
        paramline = "_".join(map(lambda n: str(pop.params[n]), swigworld_param_names))
        filename = time.strftime("%Y-%m-%d_%H-%M-%S_") + paramline + ".txt.bz2"
        f = open(self.datadir + "/" + filename , "w")
        f.write(bz2.compress(pop.writestat()))
        f.close()
        
        self.readindex()
        self.dataindex.append({"params" : pop.params, "stat_names": pop.stat_names, "file": filename})
        self.saveindex()
        
    def find(self, **kwargs):
        params = copy(default_params)
        params.update(kwargs)
        steps = params.pop("steps")
        if kwargs.has_key("steps"):
            kwargs.pop("steps")
        # print "Seed:", params["seed"]
        if params["Binitial"] < 0:
            params["Binitial"] = params["B"]
            if kwargs.has_key("Binitial"):
                kwargs.pop("Binitial")
        if params["seed"] < 0: # if seed < 0 then we can choose trajectory with any seed
            seed = params.pop("seed")
            if kwargs.has_key("seed"):
                kwargs.pop("seed")
        found = self.dataindex
        for k, v in params.iteritems():
            # print k, v, len(found)
            found = filter(lambda a: a["params"][k] == v, found)
        if len(found) == 0:
            return None
        maxsteps = max(map(lambda a: a["params"]["steps"], found))
        if maxsteps < steps:
            return None
        found = filter(lambda a: a["params"]["steps"] == maxsteps, found)
        return found[randint(len(found))] # return some trajectory from those who have maximum length
    
    def load(self, **kwargs):
        found = self.find(**kwargs)
        if found == None:
            return None
        filename = found["file"]
        f = open(self.datadir + "/" + filename)
        params, stat = readstat(bz2.decompress(f.read()))
        return params, stat

cache = Cache()

class population_cached(population_swig):
    
    def __init__(self, **kwargs):
        data = cache.load(**kwargs)
        if data:
            self._params, self._stat = data
            self.params, self.stat = self._params, self._stat
            for k, v in self.params.iteritems():
                setattr(self, k, v)
            print "Loading data from file"
        else:
            population_swig.__init__(self, **kwargs)
            calc_units = self.N * self.G * self.X * self.steps
            print "Simulation will be done: ", calc_units, "calculation units or about ", 5e-8 * calc_units, " seconds"
        
    
    def run(self, steps = None):
        if hasattr(self, "model"): # if no cache available
            population_swig.run(self, steps)
            cache.save(self) # now this case is present in cache!
            print "Data saved in file"
            return
        if steps:
            if steps > self.params["steps"]:
                print "WARNING!! Requested length of simulation exceeds length of statistics in file!"
                self.params["steps"] = steps
                population_swig.__init__(self, **self.params)
                population_swig.run(self, steps)
                cache.save(self) # now this case is present in cache!
                print "Data saved in file"
            else:
                self.params, self.stat = self._params, self._stat
                for k, v in self.stat.iteritems():
                    v = v[:(steps % self.interval)]
        else:
            self.params, self.stat = self._params, self._stat


# population = population_swig
population = population_cached
