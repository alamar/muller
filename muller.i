%module muller 
%{
#include "world.h"
// typedef Spam * Spam_p;


static int myErr = 0; // flag to save error state

%}

%include "carrays.i"
%include "exception.i"
%include "macros.i"

// typedefs

// %array_class(float, floatArray);
%array_class(bool, boolArray);
%array_class(int, intArray);

%inline %{
    typedef boolArray * boolp;
%}
%array_class(boolp, boolArrayArray);

// class member wraps

ARRAYMEMBER(Organism, good_genes, intArray);
// ARRAYCLASS(Organism, genes, G, bool);
ARRAYMEMBER(Organism, chromosomes, boolArrayArray); // NO RANGE CHECKING YET!!1

%inline %{
    typedef Organism * pOrganism;
%}
%array_class(pOrganism, organismArray);

ARRAYMEMBER(World, pop, organismArray);
ARRAYMEMBER(World, offsprings, organismArray);
ARRAYCLASS(World, pop, N, Organism*);



%include "world_export.h"

