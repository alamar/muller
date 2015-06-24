#include "world.h"

// uniformly select n items from N items, out[i] is whether i-th item is selected
// slow if part of selected items is more than ~0.11 cos needs to choose another position if selects already chosen
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
int select_some_from(real p, int N, bool * out, std::default_random_engine * generator) {
    int selected_number = 0;
    if ((p < 0.11) & (N > 40)) { // cos different algorithms are faster in different cases
        selected_number = std::binomial_distribution<int>(N, p)(*generator);
        select_n_from(selected_number, N, out, generator);
    } else {
        std::bernoulli_distribution whether_to_choose(p);
        for(int i = 0; i < N; i++) {
            bool chosen = whether_to_choose(*generator);
            out[i] = chosen;
            selected_number += (int) chosen;
        }
    }
    return selected_number;
}


void Organism::mutate(){
    
    if (M > 0) {
        
        std::bernoulli_distribution whether_to_mutate(M);
        std::bernoulli_distribution random_mutation(B);
        
        bool * mutation_array = new bool[G];
        for(int x = 0; x < X; x++) {
            select_some_from(M, G, mutation_array, generator);
            for(int i = 0; i < G; i++)
                if (mutation_array[i])
                    chromosomes[x][i] = random_mutation(*generator);
        }
        delete mutation_array;
        
        if ((Tmut > 0) && whether_to_mutate(*generator)) 
            T = max(0., T + std::uniform_real_distribution<real> (- 0.5 * Tmut, 0.5 * Tmut)(*generator));
        if ((Mmut > 0) && whether_to_mutate(*generator))
            M = max(0., M + std::uniform_real_distribution<real> (- 0.5 * Mmut, 0.5 * Mmut)(*generator));
        calc_fitness();
    };
};

void Organism::transform(Organism * donor){
    std::bernoulli_distribution whether_to_transform(T);
    std::uniform_int_distribution<int> random_chromosome(0, donor->X - 1);

    if (T > 0){
        
        bool * transformation_array = new bool[G];
        for(int x = 0; x < X; x++) {
            select_some_from(T, G, transformation_array, generator);
            for(int i = 0; i < G; i++)
                if (transformation_array[i])
                    chromosomes[x][i] = donor->chromosomes[random_chromosome(*generator)][i]; // may be slow!
        }
        delete transformation_array;
        
        if (Ttransform){
            if (whether_to_transform(*generator)) T = donor->T;
            if (whether_to_transform(*generator)) M = donor->M;
        }
        calc_fitness();
    }
};

int Organism::change_ploidy(int Xnew) {
    if (Xnew != X) {
        if (Xnew > X) {
            if (Xnew > MAX_CHROMOSOMES) Xnew = MAX_CHROMOSOMES;
            for(int x = X; x < Xnew; x++){
                chromosomes[x] = new bool [G];
                if (X > 1)
                    copy(chromosomes[0], chromosomes[0] + G, chromosomes[x]);
            }
        } else {
            //if (Xnew < 1) Xnew = 1;
            for(int x = Xnew; x < X; x++)
                delete chromosomes[x];
        }
        X = Xnew;
    }
    return X;
}

void Organism::uneven_division_from(Organism * parent) {
    int Xold = parent -> X;
    int Xnew = parent -> X;
    bool * chromosomes_replicated = new bool[Xold * 2]; // firstly genome is replicated, then approximately half of chromosomes are chosen
    //if ( ! (parent -> constantX)) {
    if ( std::bernoulli_distribution (1 - parent -> constantX) (*generator) ) {
        Xnew = std::binomial_distribution<int>(Xold * 2, 0.5)(*generator);
    }
    Xnew = change_ploidy(Xnew);
    select_n_from(Xnew, Xold * 2, chromosomes_replicated, generator);
    /*
    if (Xnew == 0) { // then we replicate only one random chromosome
        select_n_from(1, Xold * 2, chromosomes_replicated, generator);
    }
    */
    int to = 0;
    for (int from = 0; from < Xold * 2; from++) {
        if (chromosomes_replicated[from]) {
            copy(parent->chromosomes[from / 2], parent->chromosomes[from / 2] + G, chromosomes[to]);
            to += 1;
        }  
    }
    if (to != Xnew) cout << to << " " << Xnew << endl;
    delete chromosomes_replicated;
}


void Organism::copy_from(Organism * donor) {
    
    // for(int i = 0; i < G; i++){
        // genes[i] = donor->genes[i];
    // };
    copy(donor->good_genes, donor->good_genes + G, good_genes);

    change_ploidy(donor -> X);
    for(int x = 0; x < X; x++)
        copy(donor->chromosomes[x], donor->chromosomes[x] + G, chromosomes[x]);
    E = donor->E;
    EE = donor->EE;
    F = donor->F;

    // X = donor->X; // already set
    even = donor->even;
    constantX = donor->constantX;
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


void Organism::replicate_from(Organism * donor) {
    
    if (donor->even) {
        copy_from(donor);
        return;
    } else {
        uneven_division_from(donor);
        calc_fitness();
    }
    
    even = donor->even;
    constantX = donor->constantX;
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


void Organism::divide_to(Organism * offspring) {
    // binary division - offsprings will be in this and offspring
    offspring -> copy_from (this); // even division (or preparation to amitosis)
    
    if (!even) { // uneven division
        int X2 = X * 2;
        bool ** chr2 = new bool * [X2]; // double set of chromosomes
        // we have one copy of chromosomes in this and one in offspring,
        // so we put they together in chr2 in order to distribute unevenly
        for (int x = 0; x < X; x++) {
            chr2 [x * 2] = chromosomes [x];
            chr2 [x * 2 + 1] = offspring -> chromosomes [x];
        }
        
        int Xnew = X;
        // if (!constantX) {
        if ( std::bernoulli_distribution (1 - constantX) (*generator) ) {
            Xnew = std::binomial_distribution<int>(X2, 0.5)(*generator);
            if (Xnew < 1) Xnew = 1;
            if (Xnew >= X2) Xnew = X2 - 1;
        }
        bool * chromosomes_replicated = new bool[X2];
        select_n_from(Xnew, X2, chromosomes_replicated, generator);
        
        X = 0;
        offspring -> X = 0;
        for (int from = 0; from < X2; from++) {
            if (chromosomes_replicated[from]) {
                offspring -> chromosomes[offspring -> X] = chr2 [from];
                offspring -> X += 1;
            } else {
                chromosomes[X] = chr2 [from];
                X += 1;
            }
        }
        delete chr2;
        delete chromosomes_replicated;
        calc_fitness();
        offspring -> calc_fitness();
    }
}


void Organism::calc_fitness() {
    
    E = 0;
    EE = 0;
    if (X <= 0) {
        F = 0;
        return;
    }
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

Organism::Organism(int _G, real _B, real _fb, real _M, real _Mmut, real _T, real _Tmut, bool _Ttransform, real _C, int _X, bool _even, real _constantX, real _Binitial, std::default_random_engine * _generator){
    
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
    chromosomes = new bool * [MAX_CHROMOSOMES];
    for(int x = 0; x < X; x++) {
        chromosomes[x] = new bool [G];
    }
    
    M = 1.;
    B = _Binitial;
//     cout << "B: " << B << " ";
    mutate();
    M = _M;
    B = _B;
    // calc_fitness();
}

Organism::~Organism() {
    delete good_genes;
    for(int x = 0; x < X; x++) {
        delete chromosomes[x];
    }
    delete chromosomes;
}
