#include "world.h"

void World::swap_pop(){
    swap(pop, offsprings);
};

void World::select(){
    vector <real> roulette(N);
    real total = 0;
    for(int i = 0; i < N; i++){
        total += pop[i]->F;
        roulette[i] = total;
//     cout << total << "\n";
    }
    std::uniform_real_distribution<real> random_roll(0, total);
    vector <real>::iterator pos;
    for(int i = 0; i < N; i++){
        // pos = lower_bound(roulette.begin(), roulette.end(), frand() * total);
        pos = lower_bound(roulette.begin(), roulette.end(), random_roll(generator));
//         cout << pos - roulette.begin() << "\n";
        offsprings[i]->copy_from(pop[pos - roulette.begin()]);
    }
}

void World::step(){
    select();
    std::uniform_int_distribution<int> random_donor(0, N - 1);
    for(int i = 0; i < N; i++){
        offsprings[i]->mutate();
        //offsprings[i]->transform(pop[randint(N)]);
        offsprings[i]->transform(pop[random_donor(generator)]);
    }
    // swap(pop, offsprings);
    swap_pop();
    time += 1;
}

void World::run(int steps, int interval){
    write_header();
    write_stat();
    for(int i = 0; i < steps; i++){
        step();
        if (time % interval == 0) write_stat();
//         cout << Favg << " ";
    }
}

void World::calc_stat(){
    Eavg = 0; // part of good genes
    Estd = 0;
    Emin = pop[0]->E;
    Emax = pop[0]->E;
    
    Favg = 0; // fitness
    Fstd = 0;
    Fmin = pop[0]->F;
    Fmax = pop[0]->F;
    
    Mavg = 0; // mutation probability
    Mstd = 0;
    Mmin = pop[0]->M;
    Mmax = pop[0]->M;
    
    Tavg = 0; // transformation probability
    Tstd = 0;
    Tmin = pop[0]->T;
    Tmax = pop[0]->T;
    
    Tplus = 0.;
    
    EGavg = 0; // part of good genes by number of gene
    EGstd = 0;
    EGmin = N;
    EGmax = 0.;
    
    // int * EG = (int *) alloca(G * sizeof(int)); // int EG[G]
    int * EG = new int[G]; // int EG[G]
//     static int EG[MAX_GENES];
    for(int j = 0; j < G; j++) EG[j] = 0;
    
    for(int i = 0; i < N; i++){
        real t;
        
        t = pop[i]->E;
        Eavg += t;
        Estd += t * t;
        if (t < Emin) Emin = t;
        if (t > Emax) Emax = t;
        
        t = pop[i]->F;
        Favg += t;
        Fstd += t * t;
        if (t < Fmin) Fmin = t;
        if (t > Fmax) Fmax = t;
        
        t = pop[i]->M;
//         cout << t << " ";
        Mavg += t;
        Mstd += t * t;
        if (t < Mmin) Mmin = t;
        if (t > Mmax) Mmax = t;
        
        t = pop[i]->T;
        Tavg += t;
        Tstd += t * t;
        if (t < Tmin) Tmin = t;
        if (t > Tmax) Tmax = t;
        
        Tplus += (pop[i]->T > 0.);
        
        for(int j = 0; j < G; j++) EG[j] += pop[i]->genes[j];
    }
    // https://en.wikipedia.org/wiki/Algebraic_formula_for_the_variance
    Eavg = Eavg / G / N;
    Estd = sqrt(abs((Estd / G / G - N * Eavg * Eavg) / (N - 1)));
    Emin = Emin / G;
    Emax = Emax / G;
    
    Favg = Favg / N;
    Fstd = sqrt(abs((Fstd - N * Favg * Favg) / (N - 1)));
    
    Mavg = Mavg / N;
    Mstd = sqrt(abs((Mstd - N * Mavg * Mavg) / (N - 1)));
    
    Tavg = Tavg / N;
    Tstd = sqrt(abs((Tstd - N * Tavg * Tavg) / (N - 1)));
    
    Tplus = Tplus / N;
    
    for(int j = 0; j < G; j++) {
        real t = EG[j];
        EGavg += t;
        EGstd += t * t;
        if (t < EGmin) EGmin = t;
        if (t > EGmax) EGmax = t;
    }
    EGavg = EGavg / N / G;
    EGstd = sqrt(abs((EGstd / N / N - G * EGavg * EGavg) / (G - 1)));
    EGmin = EGmin / N;
    EGmax = EGmax / N;
    
    delete EG;
}

void World::write_stat(){
    calc_stat();
    cout << time << " ";
    cout << Eavg << " " << Estd << " " << Emin << " " << Emax << " ";
    cout << Favg << " " << Fstd << " " << Fmin << " " << Fmax << " ";
    cout << Mavg << " " << Mstd << " " << Mmin << " " << Mmax << " ";
    cout << Tavg << " " << Tstd << " " << Tmin << " " << Tmax << " ";
    cout << Tplus << " ";
    cout <<EGavg << " " <<EGstd << " " <<EGmin << " " <<EGmax << " ";
    cout << endl;
}

void World::write_header(){
    cout << "\n"
    "Statistics begin\n"
    "time "
    "Eavg  Estd  Emin  Emax  "
    "Favg     Fstd     Fmin     Fmax     "
    "Mavg     Mstd     Mmin     Mmax     "
    "Tavg     Tstd     Tmin     Tmax     "
    "Tplus "
    "EGavg  EGstd  EGmin  EGmax  "
    "\n";
}

World::World(int _N, int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, real _Binitial, long long int _seed){
    
    N = _N;
    G = _G;
    B = _B;
    fb = _fb;
    M = _M;
    Mmut = _Mmut;
    T = _T;
    Tmut = _Tmut;
    Ttransform = _Ttransform;
    C = _C;
    Binitial = _Binitial;
    
    time = 0;
    
    seed = _seed;
    std::seed_seq seed2 = {seed}; 
    generator.seed(seed2);
    
    pop = new Organism * [N];
    offsprings = new Organism * [N];
    
    for(int i = 0; i < N; i++){
        pop[i] = new Organism(_G, _B, _fb, _M, _Mmut, _T, _Tmut, _Ttransform, _C, _Binitial, &generator);
        offsprings[i] = new Organism(_G, _B, _fb, _M, _Mmut, _T, _Tmut, _Ttransform, _C, _Binitial, &generator);
    }
}

