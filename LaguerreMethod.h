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
#ifndef LAGUERREMETHOD_H
#define LAGUERREMETHOD_H

#include <complex>
#include <deque>

using namespace std;

#define MAX_IT 1000
#define MT 10
#define EPSS 1e-7
#define EPS 2e-6


/*
  This class implements the Laguerre method for polyinomial root solving:
  
  Given P(x) = (p_N-1)·x^N + (p_N-2)·x^N-2 + ... + (p_0)
  the problem is to find x_k, k = {1,N} such as P(x_k) = 0

  The Laguerre method guarantees convergence to some root regardless the initial guess. In this sense, this code is based on calling the
  Lagurre method several times until the N roots are found.

  Since P(x) can be expressed in terms of its roots, P(x) = (x-x_N)·(x-x_N-1)· ... · (x-x_1), we split the products into sums by taking
  the natural logarithm of P(x), so

  ln(P(x)) = ln(|x - x_N|) + ln(|x - x_N-1|) + ... + ln(|x - x_1|)

  Now, G(x) and H(x) are the first and the second derivatives of ln(P(x)):

  G(x) = d(ln(P(x)))/dx = 1/(x-x_N) + 1/(x-x_N-1) + ... 1/(x-x_1)
  H(x) = d^2(ln(P(x)))/dx^2 = 1/(x-x_N)^2 + 1/(x-x_N-1)^2 + ... 1/(x-x_1)^2

  At this point, the Laguerre's method assumes that each root, x_k, is isolated from the others, which are bunched. The distance between
  x_k and the other roots is denoted as dx and this can be used for deriving the following formula for finding the root. This process is
  repeated until dx -> 0 which means that the algorithm converged to the root.
 
  The price of the algorithm simplicity is that Laguerre's method has more computational load than, e.g., Jenkins-Traub method which is
  the standard method for root solving in most commercial software packages.

  Notice that the code structure as well as the way the formulas were implemented was inspired by [2]

  References:
  [1] Numerical Methods that Work. Forman S. Acton. The mathematical association of America. 1990
  [2] Numerical Recipes in C. The art of scientific computing. Second Edition. William H. Press et al. Cambridge University Press. 1992
*/
class LaguerreMethod
{
public:
  LaguerreMethod(deque<complex<double>>);//Class constructor
  deque<complex<double>> solve_roots();

private:
  complex<double> laguerre_core(deque<complex<double>>, complex<double>);//Core function for root calculation
  deque<complex<double>> PolyCoeffs;// Polynomial coefficients
};

#endif
