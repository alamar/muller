#include "world.h"


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
        long long int seed = (argc >= 15) ? atof(argv[14]) : -1;
        
        if (seed < 0){ // setting random seed from current time if seed < 0
            using namespace std::chrono;
            system_clock::time_point tp = system_clock::now();
            system_clock::duration dtn = tp.time_since_epoch();
            // seed = dtn.count(); // current time
            seed = dtn.count() * system_clock::period::num % system_clock::period::den; // subsecond part of current time
        }
        
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
        
        World w(N, G, B, fb, M, Mmut, T, Tmut, Ttransform, C, Binitial, seed);
        
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
        "seed - random seed (if -1 then seed will be random)\n"
        "\n"
        "Example:\n"
        "\n"
        "./muller 2000 1000 100 0.1 0.04 0.01 0. 0. 0.05 1 0. -1. 1000\n"
        ;
    }
//     w.step();
    
//     cout << "Hello world!!!\n";
}
