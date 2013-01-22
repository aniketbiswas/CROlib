#pragma once
#include <stdlib.h>

double randDouble()
{
    return (double)rand()/((double)RAND_MAX + 1.0);
}

int randInt(int lim)
{
    return (int)(randDouble() * (double)lim);
}
