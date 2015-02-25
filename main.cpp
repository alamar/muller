#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <string>
#include <algorithm>

#define BINOMIAL

// #ifdef BINOMIAL

#include <random>

std::default_random_engine generator;
// std::ranlux24 generator;
// std::mt19937 generator;
// std::ranlux24_base generator;


// #endif

// const int MAX_GENES = 100000;

typedef double real;
// typedef float real;

using namespace std;

/*
real frand(){
    return static_cast<real>(rand()) / RAND_MAX;
};

bool randbool(real probability){
    // returns true with given probability
    return rand() < (probability * RAND_MAX);
}

int randint(int n){
    return rand() % n;
}
*/

class Organism {
    
public:
    
    bool * genes; // genome
    
    int G; // number of genes
    real B; // probability of beneficial mutation
    real fb; // benefit from one mutation
    real M; // mutation probability
    real Mmut; // mutation probability mutation
    real T; // transformation probability
    real Tmut; // transformation probability mutation
    bool Ttransform; // transformation of transformation/mutation genes
    real C; // cost of transformation
    
    int E; // calculated number of good genes
    real F; // calculated fitness
    
    void mutate();
    void transform(Organism * donor);
    void copy_from(Organism * donor);
    void calc_fitness();
    
    Organism(int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, real _Binitial);
};

void Organism::mutate(){
    
    std::bernoulli_distribution random_mutation(B);
    std::bernoulli_distribution whether_to_mutate(M);
    if (M > 0) {
        
        
#ifdef BINOMIAL
        if ((M < 0.11) & (G > 40)) {
            // std::binomial_distribution<int> distribution(G, M);
            int mutation_number = std::binomial_distribution<int>(G, M)(generator);
	    std::uniform_int_distribution<int> random_gene(0, G - 1);
            
//             bool mutation_log[G]; // incompatible with msvc?
//             vector<bool> mutation_log(G); // too slow
            // bool * mutation_log = (bool *) alloca(G * sizeof(bool)); // incompatible with mingw!
            bool * mutation_log = new bool[G]; // slower than alloca but cross platform
//             static bool mutation_log[MAX_GENES]; // fast and compatible, but not thread safe!

            for(int i = 0; i < G; i++) mutation_log[i] = false;
            
            for(int i = 0; i < mutation_number; i++){
                // int pos = randint(G);
                int pos = random_gene(generator);
                while(mutation_log[pos] == true) {
                    // pos = randint(G);
                    pos = random_gene(generator);
                }
                // genes[i] = randbool(B);
                genes[pos] = random_mutation(generator);
                mutation_log[pos] = true;
            }
            delete mutation_log;
        } else 
#endif
        {
            for(int i = 0; i < G; i++){
                // if (randbool(M)){
                if (whether_to_mutate(generator)){
                    // genes[i] = randbool(B);
                    genes[i] = random_mutation(generator);
                }
            };
        }
        
        // if ((Tmut > 0) && randbool(M)) T = max(0., T + frand() * Tmut - 0.5 * Tmut);
        if ((Tmut > 0) && whether_to_mutate(generator)) 
            T = max(0., T + std::uniform_real_distribution<real> (- 0.5 * Tmut, 0.5 * Tmut)(generator));
        // if ((Mmut > 0) && randbool(M)) M = max(0., M + frand()* Mmut - 0.5 * Mmut);
        if ((Mmut > 0) && whether_to_mutate(generator))
            M = max(0., M + std::uniform_real_distribution<real> (- 0.5 * Mmut, 0.5 * Mmut)(generator));
        calc_fitness();
    };
};

void Organism::transform(Organism * donor){
    std::bernoulli_distribution whether_to_transform(T);
    if (T > 0){
        
#ifdef BINOMIAL
        if ((T < 0.11) & (G > 40)) {
            std::binomial_distribution<int> distribution(G, T);
	    std::uniform_int_distribution<int> random_gene(0, G - 1);
            
            int transformation_number = distribution(generator);
            // bool * transformation_log = (bool *) alloca(G * sizeof(bool));
            bool * transformation_log = new bool[G];
//             static bool transformation_log[MAX_GENES];
            for(int i = 0; i < G; i++) transformation_log[i] = false;
            
            for(int i = 0; i < transformation_number; i++){
                // int pos = randint(G);
                int pos = random_gene(generator);
                while(transformation_log[pos] == true) {
                    // pos = randint(G);
                    pos = random_gene(generator);
                }
                genes[pos] = donor->genes[pos];
                transformation_log[pos] = true;
            }
            delete transformation_log;
        } else 
#endif
        {
            for(int i = 0; i < G; i++){
                // if (randbool(T)){ 
                if (whether_to_transform(generator)){ 
                    genes[i] = donor->genes[i];
                };
            };
        }
        
        if (Ttransform){
            // if (randbool(T)) T = donor->T;
            if (whether_to_transform(generator)) T = donor->T;
            // if (randbool(T)) M = donor->M;
            if (whether_to_transform(generator)) M = donor->M;
        }
        calc_fitness();    
    }
};

void Organism::copy_from(Organism * donor){
    
    // for(int i = 0; i < G; i++){
        // genes[i] = donor->genes[i];
    // };
	copy(donor->genes, donor->genes + G, genes);
    
    G = donor->G;
    B = donor->B;
    fb = donor->fb;
    M = donor->M;
    Mmut = donor->Mmut;
    T = donor->T;
    Tmut = donor->Tmut;
    Ttransform = donor->Ttransform;
    C = donor->C;
    
    E = donor->E;
    F = donor->F;
};

void Organism::calc_fitness(){
    
    E = 0;
    for(int i = 0; i < G; i++){
//         cout  << genes[i] << " ";
        E += (int)genes[i];
    };
    F = pow((1 - fb), (G - E)); // multiplicative fitness, maximum fitness is always 1
    if (T > 0.) F = F * (1 - C);
};

Organism::Organism(int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, real _Binitial){
    
    G = _G;
    B = _B;
    fb = _fb;
    M = _M;
    Mmut = _Mmut;
    T = _T;
    Tmut = _Tmut;
    Ttransform = _Ttransform;
    C = _C;
    
    genes = new bool[G];
    M = 1.;
    B = _Binitial;
//     cout << "B: " << B << " ";
    mutate();
    M = _M;
    B = _B;
}

class World {
    
public:
    
    int N; // number of individuals
    int G; // number of genes
    
    real B; // probability of beneficial mutation
    real fb; // benefit from one mutation
    real M; // initial mutation probability
    real Mmut; // initial mutation probability mutation
    real T; // initial transformation probability
    real Tmut; // initial transformation probability mutation
    bool Ttransform; // transformation of transformation/mutation genes
    real C; // cost of transformation
    real Binitial; // initial content of beneficial genes
    
    Organism ** pop; // population
    
    Organism ** offsprings; // temporary population
    
    void select();
    void step();
    void run(int steps, int interval);
    
    // void init(int _N, int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut);
    World(int _N, int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, real _Binitial);
    
    // statistics
    
    int time;
    
    real Eavg; // part of good genes
    real Estd; // by organism
    real Emin;
    real Emax;
    
    real Favg; // fitness
    real Fstd;
    real Fmin;
    real Fmax;
    
    real Mavg; // mutation probability
    real Mstd;
    real Mmin;
    real Mmax;
    
    real Tavg; // transformation probability
    real Tstd;
    real Tmin;
    real Tmax;
    
    real Tplus;
    
    real EGavg; // part of good genes
    real EGstd; // by number of gene
    real EGmin;
    real EGmax;
    
    void calc_stat();
    void write_header();
    void write_stat();
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
    swap(pop, offsprings);
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

World::World(int _N, int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, real _Binitial){
    
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
    
    pop = new Organism * [N];
    offsprings = new Organism * [N];
    
    for(int i = 0; i < N; i++){
        pop[i] = new Organism(_G, _B, _fb, _M, _Mmut, _T, _Tmut, _Ttransform, _C, _Binitial);
        offsprings[i] = new Organism(_G, _B, _fb, _M, _Mmut, _T, _Tmut, _Ttransform, _C, _Binitial);
    }
}

int main(int argc, char* argv[]){
//     cout << sizeof(bool) << endl; // sizeof(bool) == 1 in g++
    if (argc >= 10){
        int steps = atoi(argv[1]); 
        int N = atoi(argv[2]);
        int G = atoi(argv[3]); 
        real B = atof(argv[4]);
        real fb = atof(argv[5]);
        real M = atof(argv[6]);
        real Mmut = atof(argv[7]);
        real T = atof(argv[8]);
        real Tmut = atof(argv[9]);
        
        int Ttransform = (argc >= 11) ? atoi(argv[10]): 1;
        real C = (argc >= 12) ? atof(argv[11]): 1;
        real Binitial = ((argc >= 13) & (atof(argv[12]) >= 0.)) ? atof(argv[12]): B;
        int interval = (argc >= 14) ? atof(argv[13]) : 1;
        int seed = (argc >= 15) ? atof(argv[14]) : 0;
        
        cout << "calculation units to calculate             " << (real)steps * N * G << "\n";
        cout << "\n";
        cout << "steps to go                                " << steps << "\n";
        cout << "N - number of individuals                  " << N << "\n";
        cout << "G - number of genes                        " << G << "\n";
        cout << "\n";
        cout << "B - probability of beneficial mutation     " << B << "\n";
        cout << "fb - benefit from one mutation             " << fb << "\n";
        cout << "M - mutation probability                   " << M << "\n";
        cout << "Mmut - mutation probability mutation       " << Mmut << "\n";
        cout << "T - transformation probability             " << T << "\n";
        cout << "Tmut - transformation probability mutation " << Tmut << "\n";
        cout << "\n";
        cout << "Ttransform - transformation of T/M genes   " << Ttransform << "\n";
        cout << "C - cost of transformation                 " << C << "\n";
        cout << "Binitial - initial beneficial genes        " << Binitial << "\n";
        cout << "interval between statistics outputs        " << interval << "\n";
        cout << "seed - random seed                         " << seed << "\n";
        
        World w(N, G, B, fb, M, Mmut, T, Tmut, Ttransform, C, Binitial);
        
//         w.calc_stat();
//         cout << w.Favg << "\n";
        
        w.run(steps, interval);
        
//         w.calc_stat();
//         cout << w.Favg << "\n";
        
    } else {
        cout << "Too few arguments!\n"
        "\n"
        "Parameters should include:\n"
        "\n"
        "steps to go\n"
        "N - number of individuals\n"
        "G - number of genes\n"
        "B - probability of beneficial mutation\n"
        "fb - benefit from one mutation\n"
        "M - mutation probability\n"
        "Mmut - mutation probability mutation\n"
        "T - transformation probability\n"
        "Tmut - transformation probability mutation\n"
        "\n"
        "Optional parameters:"
        "Ttransform - transformation of T/M genes\n"
        "C - cost of transformation\n"
        "Binitial - initial content of beneficial genes (if -1 then Binitial = B)\n"
        "interval - interval between statistics outputs\n"
        "seed - random seed\n"
        "\n"
        "Example:\n"
        "\n"
        "./muller 2000 1000 100 0.1 0.04 0.01 0. 0. 0.05 1 0. -1. 1000\n"
        ;
    }
//     w.step();
    
//     cout << "Hello world!!!\n";
}
