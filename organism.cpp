#include "world.h"

// uniformly select n items from N items, out[i] is whether i-th item is selected
void select_n_from(int n, int N, bool * out, std::default_random_engine * generator) {
    for(int i = 0; i < N; i++) out[i] = false;
    std::uniform_int_distribution<int> random_position(0, N - 1);
    for(int i = 0; i < n; i++) {
        int pos = random_position(*generator);
        while(out[pos] == true) {
            // pos = randint(G);
            pos = random_position(*generator);
        }
        out[pos] = true;
    }
}

// uniformly select items from N items with probability p
// out[i] is whether i-th item is selected
int select_binomial_from(real p, int N, bool * out, std::default_random_engine * generator) {
    int selected_number = std::binomial_distribution<int>(N, p)(*generator);
        // TODO: implement dumb selection (alternative to (M < 0.11) & (G > 40) or ifdef BINOMIAL)
    select_n_from(selected_number, N, out, generator);
    return selected_number;
}


void Organism::mutate(){
    
    std::bernoulli_distribution random_mutation(B);
    std::bernoulli_distribution whether_to_mutate(M);
    if (M > 0) {
        
        // TODO: migrate to select_binomial_from()
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
            for(int x = 0; x < X; x++) for(int i = 0; i < G; i++) {
                // if (randbool(M)){
                if (whether_to_mutate(*generator)){
                    // genes[i] = randbool(B);
                    chromosomes[x][i] = random_mutation(*generator);
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
    std::uniform_int_distribution<int> random_chromosome(0, donor->X - 1);

    if (T > 0){
        
        // TODO: migrate to select_binomial_from()
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
            for(int x = 0; x < X; x++) for(int i = 0; i < G; i++){
                // if (randbool(T)){ 
                if (whether_to_transform(*generator)){ 
                    chromosomes[x][i] = donor->chromosomes[random_chromosome(*generator)][i]; // may be slow!
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

void Organism::change_ploidy(int Xnew) {
    if (Xnew != X) {
        if (Xnew > X) {
            if (Xnew > MAX_CHROMOSOMES) Xnew = MAX_CHROMOSOMES;
            for(int x = X; x < Xnew; x++)
                chromosomes[x] = new bool [G];
        } else {
            if (Xnew < 1) Xnew = 1;
            for(int x = Xnew; x < X; x++)
                delete chromosomes[x];
        }
        X = Xnew;
    }
}

void Organism::uneven_division_from(Organism * parent) {
    int Xnew = 0;
    bool * chromosomes_replicated = new bool[X * 2]; // firstly genome is replicated, then approximately half of chromosomes are chosen
    if (constantX)
        select_n_from(X, X * 2, chromosomes_replicated, generator);
    else
        Xnew = select_binomial_from(0.5, X * 2, chromosomes_replicated, generator);
    int to = 0;
    for (int from = 0; from < X * 2; from++) {
        if (chromosomes_replicated[from]) {
            copy(parent->chromosomes[from / 2], parent->chromosomes[from / 2] + G, chromosomes[to]);
            to += 1;
        }
    }
}


void Organism::copy_from(Organism * donor){
    
    // for(int i = 0; i < G; i++){
        // genes[i] = donor->genes[i];
    // };
    // copy(donor->genes, donor->genes + G, genes);

    change_ploidy(donor -> X);
    for(int x = 0; x < X; x++)
        copy(donor->chromosomes[x], donor->chromosomes[x] + G, chromosomes[x]);
    E = donor->E;
    EE = donor->EE;
    F = donor->F;

    // X = donor->X; // already set
    even = donor->even;
    G = donor->G;
    B = donor->B;
    fb = donor->fb;
    M = donor->M;
    Mmut = donor->Mmut;
    T = donor->T;
    Tmut = donor->Tmut;
    Ttransform = donor->Ttransform;
    C = donor->C;
    
};


void Organism::replicate_from(Organism * donor){
    
    if (donor->even) {
        copy_from(donor);
    } else {
        uneven_division_from(donor);
        calc_fitness();
    }
    
    even = donor->even;
    G = donor->G;
    B = donor->B;
    fb = donor->fb;
    M = donor->M;
    Mmut = donor->Mmut;
    T = donor->T;
    Tmut = donor->Tmut;
    Ttransform = donor->Ttransform;
    C = donor->C;
    
};

void Organism::calc_fitness(){
    
    E = 0;
    EE = 0;
    for(int i = 0; i < G; i++) {
        good_genes[i] = 0;
        int g = 0; // good_genes[i]
        for(int x = 0; x < X; x++)
            g += (int) chromosomes[x][i];
//         cout  << genes[i] << " ";
        EE += (g > 0);
        E += g;
        good_genes[i] = g;
    }
    F = pow((1 - fb), (G - EE)); // multiplicative fitness, maximum fitness is always 1
    if (T > 0.) F = F * (1 - C);
};

Organism::Organism(int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, int _X, bool _even, bool _constantX, real _Binitial, std::default_random_engine * _generator){
    
    G = _G;
    X = _X;
    even = _even;
    constantX = _constantX;
    
    B = _B;
    fb = _fb;
    M = _M;
    Mmut = _Mmut;
    T = _T;
    Tmut = _Tmut;
    Ttransform = _Ttransform;
    C = _C;

    generator = _generator;
    
    good_genes = new int[G];
    chromosomes = new bool * [G];
    for(int x = 0; x < X; x++) {
        chromosomes[x] = new bool [G];
    }
    
    M = 1.;
    B = _Binitial;
//     cout << "B: " << B << " ";
    mutate();
    M = _M;
    B = _B;
}


