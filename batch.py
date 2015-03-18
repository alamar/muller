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
        self.statnames = []
        for statname in s[begin+1].split():
            self.stat[statname] = []
            self.statnames.append(statname)
            
        for line in s[begin+2:]:
            s = line.split()
            for i in xrange(len(s)):
                self.stat[self.statnames[i]].append(float(s[i]))
        
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
            print "File", name, "exists! Delete it or save in another place."
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
        # returns string representing statistics table
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
            self.steps = steps
        for i in xrange(self.steps):
            self.model.step()
            if self.model.time % self.interval == 0:
                self.append_stat()
     
population = population_swig


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
    statnames = []
    for statname in s[begin+1].split():
        stat[statname] = []
        statnames.append(statname)
        
    for line in s[begin+2:]:
        s = line.split()
        for i in xrange(len(s)):
            stat[statnames[i]].append(float(s[i]))
    
    for k, v in stat.iteritems():
        stat[k] = array(v) # for memory economy!
    return params, stat


class batch:
    
    def __init__(self, params = {}, constants = {}, verbose = True):
        """Generate cartesian product of parameters and initialize simulation systems.
        Params is dictionary of lists of variable parameters.
        Constants is dictionary of constant parameters. Example:

        b0 = batch.batch({"G": [10, 100, 1000], "N": [10, 100, 1000]}, {"fb": 0.04, "B": 0.1, "T": 0., "seed": 265})
        b.run(2000)"""
        self.params = params
        self.constants = constants
        
        cases = [constants]
        for k, v in params.items():
            cc = []
            for vv in v:
                for c in cases:
                    a = deepcopy(c)
                    a[k] = vv
                    cc.append(a)
            cases = cc
        self.models = [population_swig(**c) for c in cases]
        for i in xrange(len(self.models)):
            self.models[i].params = cases[i]
        if verbose:
            print len(cases), "cases: ", cases, len(cases)
        
    def run(self, steps = 20000, verbose = True):
        if verbose:
            print "Running", len(self.models), "tasks, here will be shown some average values over last half of simulation."
        start_t = time.time()
        i = 0
        for m in self.models:
            if verbose:
                print i, m.params
            i += 1
            
            start = time.time()
            m.run(steps)
            end = time.time()
            t = end - start
            
            if verbose:
                print t, "seconds, fitness:", average(m.stat["Favg"][(steps/2):]), ", good genes:", average(m.stat["Eavg"][(steps/2):]), ", transformation:",  average(m.stat["Tavg"][(steps/2):]), ", Tmax:", average(m.stat["Tmax"][(steps/2):])
        end_t = time.time()
        if verbose:
            print "Total time:", end_t - start_t, "seconds"
    
    def save(self, name = "trajectories"):
        for m in self.models:
            m.save(name)
    
    def grid(self, x, y, z, p = 0.5):
        """Construct 2d array representing averaged statistic vs two parameters. x, y - names of parameters, z - name of statistic, p - part of statistic that will be averaged"""
        xvalues = {}
        for m in self.models:
            xvalues[m.params[x]] = True
        xvalues = sorted(xvalues.keys())
        yvalues = {}
        for m in self.models:
            yvalues[m.params[y]] = True
        yvalues = sorted(yvalues.keys())
        
        zvalues = zeros((len(yvalues), len(xvalues)))
        for m in self.models:
            zvalues[yvalues.index(m.params[y]), xvalues.index(m.params[x])] = average(m.stat[z][int(len(m.stat[z]) * p):])
        return xvalues, yvalues, zvalues
        
    def index_grid(self, x, y):
        """Construct 2d array representing indices of models in batch vs two parameters. x, y - names of parameters, z - name of statistic, p - part of statistic that will be averaged"""
        xvalues = {}
        for m in self.models:
            xvalues[m.params[x]] = True
        xvalues = sorted(xvalues.keys())
        yvalues = {}
        for m in self.models:
            yvalues[m.params[y]] = True
        yvalues = sorted(yvalues.keys())
        
        ivalues = [[0 for a in xrange(len(xvalues))] for b in xrange(len(yvalues))]
        for i in xrange(len(self.models)):
            m = self.models[i]
            ivalues[yvalues.index(m.params[y])][xvalues.index(m.params[x])] = i
        return xvalues, yvalues, ivalues
    

def difference(y, di):
    o = []
    for i in xrange(di, len(y) - di):
        o.append((y[i + di] - y[i - di]) / 2.)
    return array(o)

def window_avg(y, di):
    o = []
    for i in xrange(di, len(y) - di):
        o.append(average(y[i - di : i + di]))
    return array(o)

# functions from fraction of genes
# def fitness_function(fb, G, E):
#     return pow(1 - fb, (1 - E) * G)

# def reverse_fitness_function(fb, G, F):
#     # calculates E from given F
#     return 1 - (log(F) / log(1 - fb)) / G

# functions from quantity of genes
def fitness_function(fb, G, E):
    return pow(1 - fb, G - E)

def log_fitness_function(fb, G, E):
    # logarithmic fitness (rate of replication) as it occurs in fisher's theorem
    return (G - E) * log(1 - fb)

def reverse_fitness_function(fb, G, F):
    # calculates E from given F
    return G - (log(F) / log(1 - fb))

def reverse_log_fitness_function(fb, G, logF):
    # calculates E from given F
    return G - (logF / log(1 - fb))

def fisher_plot(model, dt = 10):
    di = max(dt / model.interval, 1)
    dt = model.interval * di
    
    E0 = array(model.stat["EEavg"]) * model.G
    F0 = log_fitness_function(model.fb, model.G, E0) # logarithmic fitness
    Fvar0 = (array(model.stat["EEstd"]) * model.G * log(1 - model.fb)) ** 2 # variance of the above
    
    E = window_avg(E0, di)
    F = window_avg(F0, di)
    Fvar = window_avg(Fvar0, di)
    dFdt = difference(F0, di) / dt
    
    Fnext = F + Fvar # fitness after one selection round according to Fisher's law
    Enext = reverse_log_fitness_function(model.fb, model.G, Fnext) # E after selection
    
    dEmut = (model.B * model.G - Enext) * model.M # change of E due to mutations for one generation
    dFmut = - dEmut * log(1 - model.fb)
    E1 = Enext + dEmut # E after selection and mutation
    F1 = log_fitness_function(model.fb, model.G, E1) # fitness after selection and mutation
    
    figure()
    title("Actual and theoretical fitness change")
    xlabel("fitness variance")
    ylabel("fitness change per one step")
    plot(Fvar, dFdt, ".", label="actual dF/dt")   # total actual dF/dt
    plot(Fvar, F1 - F, label="predicted dF/dt") # predicted total dF/dt
    plot(Fvar, dFmut, ".", markersize=2, label="mutational dF/dt")  # predicted dF/dt of mutation
    plot(Fvar, Fvar, ".", markersize=2, label="selection dF/dt")   # variance of logarithmic fitness = predicted dF/dt of selection
    legend()



def fisher_test(world, steps):
    """Test for selection and mutation processes. Mutation ok (as of 3.03.2015), selection is 15 times stronger than theoretical :( """
    for i in xrange(steps):
        fvec = log([o.F for o in world]) # logarithmic fitness as it occurs in fisher's theorem
        F = average(fvec)
        # Fvar = var(fvec)
        
        world.select()
        # world.swap_pop()
        
        fvec0 = log([o.F for o in world])
        F0 = average(fvec0)
        
        print ""
        print "Step ", i
        print ""
        print "Selection:"
        print "Fitness before: ", F
        print "Fitness variance (= predicted fitness change): ", var(fvec)
        print "                        Actual fitness change: ", F0 - F
        print "Fitness after: ", F0
        
        evec0 = [o.EE for o in world]
        
        # for o in world: o.mutate()
        world.mutate()
        
        fvec1 = log([o.F for o in world])
        evec1 = [o.EE for o in world]
        
        dEpred = (world.B - average(evec0) / world.G) * world.M * world.G
        F0calc = log_fitness_function(world.fb, world.G, average(evec0)) # this is wrong cos fitness is nonlinear function of E, so average fitness is not fitness of average E
        F1calc = log_fitness_function(world.fb, world.G, average(evec0) + dEpred) # wrong ^^
        print ""
        print "Mutation:"
        print "Good genes (EE) before: ", average(evec0)
        print "               Predicted EE change: ", dEpred
        print "                  Actual EE change: ", (average(evec1) - average(evec0))
        print "           Fitness before: ", average(fvec0)
        print "Calculated fitness before: ", F0calc
        print "         Predicted fitness change: ", F1calc - F0calc
        print "            Actual fitness change: ", average(fvec1) - F0

def bools_to_int(a):
    o = 0
    for i in a:
        o = o * 2 + i
    return o

def int_to_bools(i, length = 0):
    o = []
    while i > 0:
        o = [bool(i % 2)] + o
        i = i / 2
    while length > len(o):
        o = [False] + o
    return o

def replication_test(world, steps, even = False, constantX = False):
    """test for uneven division (amitosis) with constant or variable number of chromosomes"""
    a = world[0]
    o = world.offsprings[0]
    X = a.X
    G = a.G
    a.even = even
    a.constantX = constantX
    print "Initial number of chromosomes X = ", X
    for i in xrange(steps):
        for x in xrange(a.X):
            xb = int_to_bools(x, G)
            for j in xrange(G):
                a.chromosomes[x][j] = xb[j]
        o.replicate_from(a)
        print o.X, " chromosomes: ",
        for x in xrange(o.X):
            xl = [o.chromosomes[x][j] for j in xrange(G)]
            print bools_to_int(xl),
        t = a
        a = o
        o = t
        print ""

def replication_alloc_test(world, steps, even = False, constantX = False):
    a = world[0]
    o = world.offsprings[0]
    X = a.X
    G = a.G
    a.even = even
    a.constantX = constantX
    for i in xrange(steps):
        o.replicate_from(a)
    

def selection_test(world, steps):
    """Test for selection process. Returns fitnesses, numbers of children, theoretical numbers of children. Example usage:

        p0 = batch.population(N = 1000, G = 100)
        f, c, p = batch.selection_test(w, 10000); figure(); plot(f, c, "."); plot(f, p); show()
        
        Passed ok (as of 17.03.2015)"""
    f = [o.F for o in world] # fitnesses
    children = [0 for o in world]
    for i in xrange(len(world)):
        world[i].E = i # dummy number to identify ancestors
    for i in xrange(steps):
        world.select()
        world.swap_pop() # offsprings will be always in offsprings array, pop won't change
        for i in xrange(len(world)):
            children[world.offsprings[i].E] += 1
    
    k = steps * 1. / average(f)
    children_pred = map(lambda f: f * k, f) # theoretical number of children
    return f, children, children_pred


def mutation_test(world, steps, B = None):
    """Test for mutational process. Passed ok (as of 17.03.2015)"""
    a = world[0]
    o = world.offsprings[0]
    B0 = a.B
    if B:
        a.B = B
    else:
        B = a.B
    u = 0 # number of increments of genes (beneficial mutations)
    d = 0 # number of decrements of genes (deleterious mutations)
    M = a.M
    E0 = a.E
    EE = a.EE
    G = a.G
    X = a.X
    e0 = 1. * E0 / G
    for i in xrange(steps):
        o.copy_from(a)
        o.mutate()
        for x in xrange(o.X):
            for j in xrange(o.G):
                if o.chromosomes[x][j] and not a.chromosomes[x][j]: u += 1
                if a.chromosomes[x][j] and not o.chromosomes[x][j]: d += 1
    print "Statistics of ", steps, " mutations of ", G, " genes (total ", G * steps, "mutations):"
    print "E0: ", E0
    print "e0: ", e0
    print "G: ", G
    print "X: ", X
    print "EE: ", EE
    print "M: ", M
    print "B: ", B
    print "Beneficial mutations : ", u
    print "Deleterious mutations: ", d
    print "Predicted beneficial mutations  (G * X - E0) * M * B * steps: ", (G * X - E0) * M * B * steps
    print "Predicted deleterious mutations E0 * M * (1 - B) * steps: ", E0 * M * (1 - B) * steps
    print ""
    
    a.B = B0
    o.B = B0


def draw_stats(batches, xname = None, yname = None, additional_stats = [], log = True, html = None):
    "Draw many interesting graphs about results. xname, yname - names of parameters"
    #fig_n = 0
    
    def logaxes():
        if log:
            xscale("log")
            yscale("log")
    
    def fig():
        if html:
            gcf().number += 1
            clf()
        else:
            figure()
    
    def save(ext = ".svg"):
        "save current fig if saving enabled"
        if html:
            global fig_n
            #fig_n += 1
            name = str(gcf().number) + ext
            savefig(dirname + "/" + name)
            htmlfile.write(' <img src = "{}" /> '.format(name))
        
    def htmlwrite(s):
        if html:
            htmlfile.write(s)
    
    def htmlclose():
        if html:
            htmlfile.close()
    
    
    if xname == None:
        xname = batches[0].params.keys()[0]
    if yname == None:
        yname = batches[0].params.keys()[1]
        
    if html:
        dirname = xname + " " + yname + " " + time.strftime("%Y-%m-%d_%H-%M-%S")
        #dirname = str(time.time())
        os.mkdir(dirname)
        htmlname = dirname + "/" + "results.html"
        htmlfile = open(htmlname, "w")
        htmlfile.write("""<!DOCTYPE html>
<html>
   <head>
      <title>{} {}</title>
      <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
   </head>
   <body>""".format(xname, yname))

    fitness = []
    good_genes = []
    transformation = []
    #stats = [[] for s in additional_stats]
    for b in batches:
        
        htmlwrite("<p><b>Постоянные параметры:  " + str(b.constants) + "</b></p>")
        htmlwrite("<p><b>Варьируемые параметры: " + str(b.params) + "</b></p>")
        if not hasattr(b.models[0], "Ttransform"):
            b.models[0].Ttransform = 1. # for compatibility!
        htmlwrite("Все параметры (на примере первого запуска): N = {}, G = {}, B = {}, fb = {}, M = {}, Mmut = {}, T = {}, Tmut = {}, Ttransform = {} ".format(b.models[0].N, b.models[0].G, b.models[0].B, b.models[0].fb, b.models[0].M, b.models[0].Mmut, b.models[0].T, b.models[0].Tmut, b.models[0].Ttransform))
        
        f = b.grid(xname, yname, "Favg")
        e = b.grid(xname, yname, "Eavg")
        t = b.grid(xname, yname, "Tavg")
        #f = grid(b, xname, yname, "Favg") # needed for objects lacking grid method
        #e = grid(b, xname, yname, "Eavg") # define global grid if needed
        #t = grid(b, xname, yname, "Tavg")
        fitness.append(f)
        good_genes.append(e)
        transformation.append(t)
        
        htmlwrite("<!--")
        htmlwrite("<p>Команда: b = batch.batch(params = {}, constants = {}); b.run()</p>" .format(b.params, b.constants))
        htmlwrite("<p>Результаты расчёта (рисуются с помощью contourf(*e), например):</p>")
        htmlwrite("<p>Приспособленность f = {}</p>".format(f))
        htmlwrite("<p>Доля хороших генов e = {}</p>".format(e))
        #htmlwrite("<p>Вероятность трансформации t = {}</p>".format(t))
        htmlwrite("-->")
        
        
        #for i in xrange(len(additional_stats)):
            #stats[i].append(b.grid(xname, yname, additional_stats[i]));
        
        htmlwrite("<p>")
        fig(); contourf(*f); title("F ({}, {})".format(xname, yname)); xlabel(xname); ylabel(yname); colorbar(); logaxes();save()
        fig(); contourf(*e); title("E ({}, {})".format(xname, yname)); xlabel(xname); ylabel(yname); colorbar(); logaxes();save()
        fig(); contourf(*t); title("T ({}, {})".format(xname, yname)); xlabel(xname); ylabel(yname); colorbar(); logaxes();save()
        htmlwrite("</p>")
        
        #for i in xrange(len(additional_stats)):
            #s = stats[i][-1]
            #fig(); contourf(*e); title("{} ({}, {})".format(additional_stats[i], xname, yname)); xlabel(xname); ylabel(yname); colorbar(); xscale("log"); yscale("log")
        
    if len(batches) > 1:
        f = fitness[0] # needed for enumerating xname and yname coords
        for i in xrange(1, len(batches)):
            fratio = fitness[i][2] / fitness[i-1][2]
            eratio = good_genes[i][2] / good_genes[i-1][2]
            tratio = transformation[i][2] / transformation[i-1][2]
            
            htmlwrite("<p></p>")
            htmlwrite("<p><b>Сравнение запусков" + str(batches[i].constants) + " и  " + str(batches[i-1].constants) + "</b></p>")
            htmlwrite("<!--")
            htmlwrite("<p>Соотношение приспособленностей fratio = {}</p>".format(fratio))
            htmlwrite("<p>Соотношение долей хороших генов eratio = {}</p>".format(eratio))
            #htmlwrite("<p>Соотношение вероятностей трансформации tratio = {}</p>".format(tratio))
            htmlwrite("-->")
            htmlwrite("<p>")
            fig(); contourf(f[0], f[1], log10(fratio)); title("log10 F ({}) /F ({}) ({}, {})".format(batches[i].constants, batches[i-1].constants, xname, yname)); xlabel(xname); ylabel(yname); logaxes();colorbar(); save()
            #fig(); contourf(f[0], f[1], log10(fratio)); title("log10 F ({}) /F ({}) ({}, {})".format(batches[i].constants, batches[i-1].constants, xname, yname)); xlabel(xname); ylabel(yname); colorbar(); save()
            fig(); contourf(f[0], f[1], eratio); title("E ({}) /E ({}) ({}, {})".format(batches[i].constants, batches[i-1].constants, xname, yname)); xlabel(xname); ylabel(yname); colorbar(); logaxes();save()
            fig(); contourf(f[0], f[1], tratio); title("T ({}) /T ({}) ({}, {})".format(batches[i].constants, batches[i-1].constants, xname, yname)); xlabel(xname); ylabel(yname); colorbar(); logaxes();save()
            htmlwrite("</p>")
            
    htmlwrite("""</body>
</html>""")
    htmlclose()
            
    #for b in batches:
        #i = index_grid(b, xname, yname)
        #time = b.models[0].stat["time"]
        ##print i
        #for y in xrange(len(i[0])):
            #fig(); xlabel("time"); ylabel("Eavg")
            #leg = []
            #for x in xrange(len(i[1])):
                #plot(time, b.models[i[2][y][x]].stat["Eavg"]);
                #leg.append(xname + " = " + str(i[1][x]))
            #title(yname + " = " + str(i[0][y]))
            #legend(leg)
            #fig(); xlabel("time"); ylabel("Favg")
            #leg = []
            #for x in xrange(len(i[1])):
                #semilogy(time, b.models[i[2][y][x]].stat["Favg"])
                #leg.append(xname + " = " + str(i[1][x]))
            #title(yname + " = " + str(i[0][y]))
            #legend(leg)
            #fig()
            #leg = []
            #for x in xrange(len(i[1])):
                #plot(time, b.models[i[2][y][x]].stat["Tavg"])
                #leg.append(xname + " = " + str(i[1][x]))
            #title(yname + " = " + str(i[0][y]))
            #legend(leg)
        #for x in xrange(len(i[1])):
            #fig()
            #for y in xrange(len(i[0])):
                #plot(time, b.models[i[2][y][x]].stat["Eavg"])
            #fig()
            #for y in xrange(len(i[0])):
                #semilogy(time, b.models[i[2][y][x]].stat["Favg"])
            #fig()
            #for y in xrange(len(i[0])):
                #plot(time, b.models[i[2][y][x]].stat["Tavg"])
                









