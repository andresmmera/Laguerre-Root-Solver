/***************************************************************************
                                LaguerreMethod.h
                             -------------------
    begin                : Jan 19 2017
    copyright            : (C) 2017 by Andres Martinez-Mera
    email                : andresmartinezmera@gmail.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <complex>
#include <deque>

using namespace std;

#define MAX_IT 1000
#define MT 10
#define EPSS 1e-7
#define EPS 2e-6


/*
  This class implements the Laguerre method for root solving.

  References:
  [1] Numerical Recipes in C. The art of scientific computing. Second Edition. William H. Press et al. Cambridge University Press. 1992
*/
class LaguerreMethod
{
public:
  LaguerreMethod(deque<complex<double>>);
  deque<complex<double>> solve_roots();

private:
  complex<double> laguerre_core(deque<complex<double>>, complex<double>);
  deque<complex<double>> PolyCoeffs;// Polynomial coefficients
};
