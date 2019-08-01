#ifndef PP_H
#define PP_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


class PP
{
public:
    PP()
    {
    }
    void performCollisions(int number, double kappa);
};

#endif // PP_H
