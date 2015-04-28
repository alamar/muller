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
    """test for uneven division (amitosis) with constant or variable number of chromosomes via replicate_from() function"""
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

def division_test(world, steps = 1):
    """test for binary division via divide_to() function, there should be exactly two chromosomes with the same number in genomes of two offsprings"""
    a = world[0]
    o = world.offsprings[0]
    X = a.X
    G = a.G
    print "Initial number of chromosomes X = ", X
    for i in xrange(steps):
        for x in xrange(a.X):
            xb = int_to_bools(x, G)
            for j in xrange(G):
                a.chromosomes[x][j] = xb[j]
        a.divide_to(o)
        print "First offspring has", a.X, "chromosomes: ",
        for x in xrange(a.X):
            xl = [a.chromosomes[x][j] for j in xrange(G)]
            print bools_to_int(xl),
        print ""
        print "Second offspring has", o.X, " chromosomes: ",
        for x in xrange(o.X):
            xl = [o.chromosomes[x][j] for j in xrange(G)]
            print bools_to_int(xl),
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
    

def chromosome_to_list(o, x):
    """Converts chromosome to pythonic list. Direct use of chromosome (like for g in o.chromosomes[x]) leads to crash due to lack of range checking."""
    # TODO: migrate all chromosome retrievals to this
    return [int(o.chromosomes[x][g]) for g in xrange(o.G)]

def compare_chromosomes(a, o):
    """checks if all offspring's chromosomes are copies of ancestor's"""
    w = True
    for i in xrange(o.X):
        x = [o.chromosomes[i][g] for g in xrange(o.G)]
        e = False
        for j in xrange(a.X):
            y = [a.chromosomes[j][g] for g in xrange(a.G)]
            e = e or (x == y)
        w = w and e
    return w

def uneven_selection_test(world):
    for i in xrange(len(world)):
        world[i].C = i # dummy number to identify ancestors
    world.select()
    for o in world:
        i = int(o.C)
        if i < 0 or i >= world.N:
            print "C is not in range of organism indices!", o.C
        else:
            a = world.offsprings[int(o.C)]
            if a.C != o.C:
                print "C is not index of ancestor!", o.C, " ", a.C
            if not compare_chromosomes(a, o):
                print "Bastard detected! ", o.C

def uneven_test(world, steps = 1):
    G = world.G
    for i in xrange(len(world)):
        ib = int_to_bools(i, G)
        a = world[i]
        a.C = i # dummy number to identify ancestors
        for x in xrange(a.X):
            # xb = int_to_bools(x, G)
            for j in xrange(G):
                a.chromosomes[x][j] = ib[j]
    for t in xrange(steps):
        perm = permutation(world.N)
        for i in xrange(len(world)):
            world.offsprings[i].replicate_from(world[perm[i]])
        for i in xrange(len(world)):
            o = world.offsprings[i]
            a = world[perm[i]]
            if a.C != o.C:
                print "Step ", t, ", organism ", i, ": C is not equal to ", a.C, " but ", o.C
            for x in xrange(a.X):
                xo = bools_to_int([o.chromosomes[x][j] for j in xrange(G)])
                xa = bools_to_int([a.chromosomes[x][j] for j in xrange(G)])
                if xo != xa:
                    print "Step ", t, ", organism ", i, ", chromosome ", x, "is not ", xa, " but ", xo
        world.swap_pop()

    


def selection_test_old(world, steps):
    """Test for selection process. Returns fitnesses, numbers of children, theoretical numbers of children.
    This version is for binary = False! Otherwise use selection_test().
    Example usage:

        p0 = batch.population_raw(N = 1000, G = 100)
        f, c, p = batch.selection_test(p0.model, 10000); figure(); plot(f, c, "."); plot(f, p); show()
        
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

def selection_test(world, steps, bias = 0):
    """Test for selection process. Returns fitnesses, numbers of children, theoretical numbers of children.
    This version is for binary = True! Otherwise use selection_test_old().
    Example usage:

        p0 = batch.population_raw(N = 1000, G = 100)
        f, c, p = batch.selection_test(p0.model, 10000); figure(); plot(f, c, "."); plot(f, p); show()
        
        """
    children = [0 for o in world] + [0 for i in xrange(len(world))]
    fitnesses = sorted([bias + rand() for i in xrange(2 * len(world))])
    for i in xrange(steps):
        for o in world:
            o.F = rand()
        for i in xrange(len(world)):
            world.offsprings[i].F = rand()
        for i in xrange(len(world)):
            world[i].E = i # dummy number to identify ancestors
            world.offsprings[i].E = i + len(world) # dummy number to identify ancestors
            world[i].F = fitnesses[i]
            world.offsprings[i].F = fitnesses[i + len(world)]
        world.select_half()
        for i in xrange(len(world)):
            children[world[i].E] += 1
    
    k = steps * 0.5 / average(fitnesses)
    children_pred = map(lambda f: f * k, fitnesses) # theoretical number of children
    return fitnesses, children, children_pred


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


