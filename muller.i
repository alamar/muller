%module muller 
%{
#include "world.h"
// typedef Spam * Spam_p;


static int myErr = 0; // flag to save error state

%}

%include "carrays.i"
%include "exception.i"

%include "world_export.h"

