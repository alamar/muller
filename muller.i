%module muller 
%{
#include "world.h"
// typedef Spam * Spam_p;


static int myErr = 0; // flag to save error state

%}

%include "carrays.i"
%include "exception.i"
%include "macros.i"


// %array_class(float, floatArray);
%array_class(bool, boolArray);

ARRAYMEMBER(Organism, genes, boolArray);
ARRAYCLASS(Organism, genes, G, bool);

%inline %{
    typedef Organism * pOrganism;
%}
%array_class(pOrganism, organismArray);

ARRAYMEMBER(World, pop, organismArray);
ARRAYMEMBER(World, offsprings, organismArray);
ARRAYCLASS(World, pop, N, Organism*);



%include "world_export.h"

