#!/usr/bin/python
# coding=utf-8

from batch import *

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
        t = a
        a = o
        o = t
    

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


