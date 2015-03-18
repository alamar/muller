// #define BINOMIAL

// #ifdef BINOMIAL


// std::default_random_engine generator;
// std::ranlux24 generator;
// std::mt19937 generator;
// std::ranlux24_base generator;


// #endif

// const int MAX_GENES = 100000;
const int MAX_CHROMOSOMES = 100;

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
    
    std::default_random_engine * generator;
    
    int * good_genes; // numbers of good genes at each locus
    bool ** chromosomes; // genome
    
    int G; // number of genes
    int X; // number of chromosomes
    bool even; // even division (mitosis) or uneven (amitosis)
    bool constantX; // keep number of chromosomes constant in uneven division (in mitosis X is always constant)

    real B; // probability of beneficial mutation
    real fb; // benefit from one mutation
    real M; // mutation probability
    real Mmut; // mutation probability mutation
    real T; // transformation probability
    real Tmut; // transformation probability mutation
    bool Ttransform; // transformation of transformation/mutation genes
    real C; // cost of transformation
    
    int E; // calculated number of good genes
    int EE; // effective number of good genes
    real F; // calculated fitness
    
    void mutate();
    void transform(Organism * donor);
    void calc_fitness();
    
    void copy_from(Organism * donor); // simple copy
    void uneven_division_from(Organism * parent);
    int  change_ploidy(int Xnew);
    void replicate_from(Organism * donor); // replication according to genetic program
    
    Organism(int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, int _X, bool _even, bool _constantX, real _Binitial, std::default_random_engine * _generator);
};


class World {
    
public:
    
    std::default_random_engine generator;
    // std::ranlux24 generator;
    // std::mt19937 generator;
    // std::ranlux24_base generator;
    
    long long int seed;
    
    int N; // number of individuals
    int G; // number of genes
    int X; // number of chromosomes
    bool even; // even division (mitosis) or uneven (amitosis)
    bool constantX; // keep number of chromosomes constant in uneven division (in mitosis X is always constant)
    
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
    
    void swap_pop(); // swap pop and offsprings
    void select();
    void mutate();
    void transform();
    
    void step();
    void run(int steps, int interval);
    
    // void init(int _N, int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut);
    World(int N, int G, real B, real fb, real M, real Mmut, real T, real Tmut, bool Ttransform, real C, int _X, bool _even, bool _constantX, real _Binitial, long long int seed);
    
    // statistics
    
    int time;
    
// #define STATVARS(v) double v##avg
#define STATVARS(name) real name##avg; real name##std; real name##min; real name##max;
    
    // real Eavg; real Estd; real Emin; real Emax;
    STATVARS(E); // part of good genes by organism
    STATVARS(EE); // part of effective good genes
    STATVARS(X); // number of chromosomes
    STATVARS(F); // fitness
    STATVARS(M); // mutation probability
    STATVARS(T); // transformation probability
    STATVARS(EG); // part of good genes by number of gene
    
    /*
    real Eavg; // part of good genes
    real Estd; // by organism
    real Emin;
    real Emax;
    
    real EEavg; // part of effective good genes
    real EEstd; // by organism
    real EEmin;
    real EEmax;

    real Xavg; // number of chromosomes 
    real Xstd;
    real Xmin;
    real Xmax;
    
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
    
    
    real EGavg; // part of good genes
    real EGstd; // by number of gene
    real EGmin;
    real EGmax;
    */
    
    real Tplus;
    
    void calc_stat();
    void write_header();
    void write_stat();
};


