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
    }
    std::uniform_real_distribution<real> random_roll(0, total);
    vector <real>::iterator pos;
    for(int i = 0; i < N; i++){
        pos = lower_bound(roulette.begin(), roulette.end(), random_roll(generator));
        offsprings[i]->replicate_from(pop[pos - roulette.begin()]);
    }
    swap_pop();
}


void World::select_half() { 
    // selection of half of population by roulette algorithm
    // each organism may be selected only once
    // population is assumed to be in pop+offsprings, chosen ones will be in pop
    bool * chosen = new bool [N * 2];
    for(int i = 0; i < N * 2; i++) chosen [i] = false;
    // select_n_from (N, N * 2, chosen, generator);
    
    real total = 0;
    std::vector <real> roulette(N * 2);
    for(int i = 0; i < N; i++){
        total += pop[i]->F;
        roulette[i] = total;
    }
    for(int i = 0; i < N; i++){
        total += offsprings[i]->F;
        roulette[i + N] = total;
    }
    // cout << "Fitness calculated" << endl;
    // roulette is now consists of sums of fitnesses of pop concatenated with offsprings
    // so we choose half of total population
    int Nnew = 0;
    vector <real>::iterator pos;
    std::uniform_real_distribution<real> random_roll(0, total);
    while (Nnew < N) {
        pos = lower_bound(roulette.begin(), roulette.end(), random_roll(generator));
        if (!(chosen [pos - roulette.begin()])) {
            chosen [pos - roulette.begin()] = true;
            Nnew += 1;
        }
    }
    // cout << "Numbers selected " << Nnew << endl;
    
    // and gather all chosen ones into pop
    // int swap_n = 0; // for debug purposes
    int j = 0; // index of member of offsprings to be interchanged
    for(int i = 0; i < N; i++) if (!(chosen [i])) {
        while (!(chosen [j + N])) j += 1;
        // cout << i << " " << j << endl;
        swap (pop[i], offsprings[j]);
        // Organism * t = pop[i];
        // pop[i] = offsprings[j];
        // offsprings[j] = t;
        j += 1;
        // swap_n += 1;
    }
    // cout << "Swaps: " << swap_n << endl;
    delete chosen;
}


void World::divide() { // binary division, offsprings will be both in pop and offsprings
    for (int i = 0; i < N; i++) {
        pop[i] -> divide_to (offsprings[i]);
    }
}


void World::mutate(){ // mutate pop[]
    for (int i = 0; i < N; i++){
        pop[i]->mutate();
    }
}

void World::transform(){ // transform pop[]
    for (int i = 0; i < N; i++){
        offsprings[i]->copy_from(pop[i]);
    }
    
    // shuffling for better randomness
    shuffle(&offsprings[0], &offsprings[N], generator);
    
    // reciprocal transformation in pairs
    // if population size is odd then the last organism is alone!
    for (int i = 0; i < N-1; i += 2){
        offsprings[i]->transform(pop[i+1]);
        offsprings[i+1]->transform(pop[i]);
    }
    swap_pop();
}

void World::step(){
    divide();
    select_half();
    // select();
    mutate();
    transform();
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
    
#define STATZERO(name, inf, sup) name##avg = 0; name##std = 0; name##min = sup; name##max = inf;
    
    STATZERO(E, 0, MAX_CHROMOSOMES*G);
    STATZERO(EE, 0, G);
    STATZERO(X, 0, MAX_CHROMOSOMES);
    STATZERO(F, 0, 1);
    STATZERO(M, 0, 1);
    STATZERO(T, 0, 1);
    STATZERO(EG, 0, MAX_CHROMOSOMES*N);
    
    Tplus = 0.;
    
    // int * EG = (int *) alloca(G * sizeof(int)); // int EG[G]
    int * EG = new int[G]; // int EG[G]
    for(int j = 0; j < G; j++) EG[j] = 0;
    
    for(int i = 0; i < N; i++){
        real t;
        
#define STATADD(name) \
        t = pop[i]->name; \
        name##avg += t; \
        name##std += t * t; \
        if (t < name##min) name##min = t; \
        if (t > name##max) name##max = t;
        
        STATADD(X);
        STATADD(E);
        STATADD(EE);
        STATADD(F);
        STATADD(M);
        STATADD(T);
         
        Tplus += (pop[i]->T > 0.);
        
        for(int j = 0; j < G; j++) 
            EG[j] += pop[i]->good_genes[j];
    }

//    name##std = sqrt(abs((name##std - name##avg * name##avg / N) / (N - 1))) / (norm); \
    // https://en.wikipedia.org/wiki/Algebraic_formula_for_the_variance
#define STATNORMALIZE(name, norm) \
    name##avg = name##avg / (norm) / N; \
    name##std = sqrt(abs((name##std / (norm * norm) - name##avg * name##avg * N) / (N - 1))); \
    name##min = name##min / (norm); \
    name##max = name##max / (norm);
    
    STATNORMALIZE(X, 1);
    STATNORMALIZE(E, G * Xavg);
    STATNORMALIZE(EE, G);
    STATNORMALIZE(F, 1);
    STATNORMALIZE(M, 1);
    STATNORMALIZE(T, 1);
    
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
    EGmin = EGmin / N / X;
    EGmax = EGmax / N / X;
    EGavg = EGavg / X;
    EGstd = EGstd / X;
    
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

World::World(int _N, int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, int _X, bool _even, bool _constantX, real _Binitial, long long int _seed){
    
    N = _N;
    G = _G;
    X = _X;
    B = _B;
    fb = _fb;
    M = _M;
    Mmut = _Mmut;
    T = _T;
    Tmut = _Tmut;
    Ttransform = _Ttransform;
    C = _C;
    even = _even;
    constantX = _constantX;
    Binitial = _Binitial;
    
    time = 0;
    
    seed = _seed;
    std::seed_seq seed2 = {seed}; 
    generator.seed(seed2);
    
    pop = new Organism * [N];
    offsprings = new Organism * [N];
    
    for(int i = 0; i < N; i++){
        pop[i] = new Organism(_G, _B, _fb, _M, _Mmut, _T, _Tmut, _Ttransform, _C, _X, _even, _constantX, _Binitial, &generator);
        offsprings[i] = new Organism(_G, _B, _fb, _M, _Mmut, _T, _Tmut, _Ttransform, _C, _X, _even, _constantX, _Binitial, &generator);
    }
}

World::~World() {
    for(int i = 0; i < N; i++) {
        delete pop[i];
        delete offsprings[i];
    }
    delete pop;
    delete offsprings;
}
