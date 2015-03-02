#include "world.h"


void Organism::mutate(){
    
    std::bernoulli_distribution random_mutation(B);
    std::bernoulli_distribution whether_to_mutate(M);
    if (M > 0) {
        
        
#ifdef BINOMIAL
        if ((M < 0.11) & (G > 40)) {
            // std::binomial_distribution<int> distribution(G, M);
            int mutation_number = std::binomial_distribution<int>(G, M)(*generator);
	    std::uniform_int_distribution<int> random_gene(0, G - 1);
            
//             bool mutation_log[G]; // incompatible with msvc?
//             vector<bool> mutation_log(G); // too slow
            // bool * mutation_log = (bool *) alloca(G * sizeof(bool)); // incompatible with mingw!
            bool * mutation_log = new bool[G]; // slower than alloca but cross platform
//             static bool mutation_log[MAX_GENES]; // fast and compatible, but not thread safe!

            for(int i = 0; i < G; i++) mutation_log[i] = false;
            
            for(int i = 0; i < mutation_number; i++){
                // int pos = randint(G);
                int pos = random_gene(*generator);
                while(mutation_log[pos] == true) {
                    // pos = randint(G);
                    pos = random_gene(*generator);
                }
                // genes[i] = randbool(B);
                genes[pos] = random_mutation(*generator);
                mutation_log[pos] = true;
            }
            delete mutation_log;
        } else 
#endif
        {
            for(int i = 0; i < G; i++){
                // if (randbool(M)){
                if (whether_to_mutate(*generator)){
                    // genes[i] = randbool(B);
                    genes[i] = random_mutation(*generator);
                }
            };
        }
        
        // if ((Tmut > 0) && randbool(M)) T = max(0., T + frand() * Tmut - 0.5 * Tmut);
        if ((Tmut > 0) && whether_to_mutate(*generator)) 
            T = max(0., T + std::uniform_real_distribution<real> (- 0.5 * Tmut, 0.5 * Tmut)(*generator));
        // if ((Mmut > 0) && randbool(M)) M = max(0., M + frand()* Mmut - 0.5 * Mmut);
        if ((Mmut > 0) && whether_to_mutate(*generator))
            M = max(0., M + std::uniform_real_distribution<real> (- 0.5 * Mmut, 0.5 * Mmut)(*generator));
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
            
            int transformation_number = distribution(*generator);
            // bool * transformation_log = (bool *) alloca(G * sizeof(bool));
            bool * transformation_log = new bool[G];
//             static bool transformation_log[MAX_GENES];
            for(int i = 0; i < G; i++) transformation_log[i] = false;
            
            for(int i = 0; i < transformation_number; i++){
                // int pos = randint(G);
                int pos = random_gene(*generator);
                while(transformation_log[pos] == true) {
                    // pos = randint(G);
                    pos = random_gene(*generator);
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
                if (whether_to_transform(*generator)){ 
                    genes[i] = donor->genes[i];
                };
            };
        }
        
        if (Ttransform){
            // if (randbool(T)) T = donor->T;
            if (whether_to_transform(*generator)) T = donor->T;
            // if (randbool(T)) M = donor->M;
            if (whether_to_transform(*generator)) M = donor->M;
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

Organism::Organism(int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, real _Binitial, std::default_random_engine * _generator){
    
    G = _G;
    B = _B;
    fb = _fb;
    M = _M;
    Mmut = _Mmut;
    T = _T;
    Tmut = _Tmut;
    Ttransform = _Ttransform;
    C = _C;

    generator = _generator;
    
    genes = new bool[G];
    M = 1.;
    B = _Binitial;
//     cout << "B: " << B << " ";
    mutate();
    M = _M;
    B = _B;
}


