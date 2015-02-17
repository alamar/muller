#!/usr/bin/python
# coding=utf-8

from pylab import *

from copy import *

import commands
import subprocess

import time

import os, sys

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


class population:
    
    def __init__(self, steps = 200, N = 100, G = 100, M = 0.015, B = 0.1, fb = 0.05, T = 0., Tmut = 0., Mmut = 0., Ttransform = 1., C = 0., Binitial = -1., interval = 1, verbose = False):
        
        # C and seed is not implemented yet!
        
        self.steps = steps
        self.interval = interval
        self.verbose = verbose
        
        self.N = N
        self.G = G
        self.M = M
        self.B = B
        self.fb = fb
        self.T = T
        self.Tmut = Tmut
        self.Mmut = Mmut
        self.Ttransform = Ttransform
        self.C = C
        
        self.paramline = " {} {} {} {} {} {} {} {} {} {} {} ".format(N, G, B, fb, M, Mmut, T, Tmut, Ttransform, C, Binitial)
        #self.paramline = " " + N +  " " + G + " " + B + " " + fb + " " + M + " " + Mmut + " " + T + " " + Tmut + " " + Ttransform + " "
        #commandline = "./cpp_gpg " + paramline
        #self.run()
    
    def readstat(self):
        self.stat = {
            "time": [],
            "Eavg": [], "Estd": [], "Emax": [], "Emin": [],
            "Favg": [], "Fstd": [], "Fmax": [], "Fmin": [],
            "Mavg": [], "Mstd": [], "Mmax": [], "Mmin": [],
            "Tavg": [], "Tstd": [], "Tmax": [], "Tmin": [],
            "Tplus": [],
            "EGavg": [], "EGstd": [], "EGmax": [], "EGmin": [],
            }
        #reads = [[] for i in xrange(len(self.stat))]
        s = self.output.splitlines()
        for i in xrange(len(s)):
            if s[i] == "Statistics begin":
                begin = i
                break
        for line in s[begin+2:]:
            s = line.split()
            self.stat["time"].append(float(s[0]))
            
            self.stat["Eavg"].append(float(s[1]))
            self.stat["Estd"].append(float(s[2]))
            self.stat["Emin"].append(float(s[3]))
            self.stat["Emax"].append(float(s[4]))
            
            self.stat["Favg"].append(float(s[5]))
            self.stat["Fstd"].append(float(s[6]))
            self.stat["Fmin"].append(float(s[7]))
            self.stat["Fmax"].append(float(s[8]))
            
            self.stat["Mavg"].append(float(s[9]))
            self.stat["Mstd"].append(float(s[10]))
            self.stat["Mmin"].append(float(s[11]))
            self.stat["Mmax"].append(float(s[12]))
            
            self.stat["Tavg"].append(float(s[13]))
            self.stat["Tstd"].append(float(s[14]))
            self.stat["Tmin"].append(float(s[15]))
            self.stat["Tmax"].append(float(s[16]))
            
            self.stat["Tplus"].append(float(s[17]))
    
            self.stat["EGavg"].append(float(s[18]))
            self.stat["EGstd"].append(float(s[19]))
            self.stat["EGmin"].append(float(s[20]))
            self.stat["EGmax"].append(float(s[21]))
            
    def run(self, steps = None):
        if steps:
            self.steps = steps
        #commandline = 
        print self.paramline
        command = "nice -n 19 ./muller "
        if sys.platform == "win32":
            command = "muller.exe "
            self.output = subprocess.check_output(command + str(self.steps) + self.paramline + str(self.interval))
        else:
            command = "nice -n 19 ./muller "
            self.output = commands.getoutput(command + str(self.steps) + self.paramline + str(self.interval))
        if self.verbose:
            print self.output
        self.readstat()
        self.output = ""
    
    
    
class batch:
    
    def __init__(self, params = {}, constants = {}, verbose = True):
        """Generate cartesian product of parameters and initialize simulation systems.
        Params is dictionary of lists of variable parameters.
        Constants is dictionary of constant parameters. Example:
        
        b = muller_finite.batch({"T": [0., 0.5, 1.], "G": [10, 20], "fb": [0.04]})
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
        self.models = [population(**c) for c in cases]
        for i in xrange(len(self.models)):
            self.models[i].params = cases[i]
        if verbose:
            print len(cases), "cases: ", cases, len(cases)
        
    def run(self, steps = 20000, resume = False, verbose = True):
        if verbose:
            print "Running", len(self.models), "tasks, here will be shown some average values over last half of simulation."
        start_t = time.time()
        for m in self.models:
            if verbose:
                print m.params
            
            start = time.time()
            m.run(steps)
            end = time.time()
            t = end - start
            
            if verbose:
                print t, "seconds, fitness:", average(m.stat["Favg"][(steps/2):]), ", good genes:", average(m.stat["Eavg"][(steps/2):]), ", transformation:",  average(m.stat["Tavg"][(steps/2):]), ", Tmax:", average(m.stat["Tmax"][(steps/2):])
        end_t = time.time()
        if verbose:
            print "Total time:", end_t - start_t, "seconds"
    
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
    
    
    if xname == None:
        xname = batches[0].params.keys()[0]
    if yname == None:
        yname = batches[0].params.keys()[1]
        
    if html:
        dirname = xname + " " + yname + " " + time.ctime(time.time())
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
                









