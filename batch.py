#!/usr/bin/python
# coding=utf-8

from population import *

class batch:
    
    def __init__(self, params = {}, constants = {}, verbose = True):
        """Generate cartesian product of parameters and initialize simulation systems.
        Params is dictionary of lists of variable parameters.
        Constants is dictionary of constant parameters. Example:

        b0 = batch.batch({"G": [10, 100, 1000], "N": [10, 100, 1000]}, {"fb": 0.04, "B": 0.1, "T": 0., "seed": 265})
        b0.run()"""
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
        # for i in xrange(len(self.models)):
        #     self.models[i].params = cases[i]
        if verbose:
            print len(cases), "cases: ", cases, len(cases)
        
    def run(self, steps = None, verbose = True):
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
                begin = m.params["steps"] / 2
                print t, "seconds, fitness:", average(m.stat["Favg"][begin:]), ", good genes:", average(m.stat["Eavg"][begin:]), ", transformation:",  average(m.stat["Tavg"][begin:]), ", Tmax:", average(m.stat["Tmax"][begin:])
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


# old shit

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
        fig(); contourf(*f); title("F ({}) ({}, {})".format(b.constants, xname, yname)); xlabel(xname); ylabel(yname); colorbar(); logaxes();save()
        fig(); contourf(*e); title("E ({}) ({}, {})".format(b.constants, xname, yname)); xlabel(xname); ylabel(yname); colorbar(); logaxes();save()
        fig(); contourf(*t); title("T ({}) ({}, {})".format(b.constants, xname, yname)); xlabel(xname); ylabel(yname); colorbar(); logaxes();save()
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
                









