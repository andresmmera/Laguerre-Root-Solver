/***************************************************************************
                                LaguerreMethod.cpp
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


#include "LaguerreMethod.h"

LaguerreMethod::LaguerreMethod(deque<complex<double>> P)
{
 PolyCoeffs = P;
}


// The Laguerre method guarantees convergence to some root regardless the initial guess.
complex<double> LaguerreMethod::laguerre_core(deque<complex<double>> P, complex<double> x)
{
  deque<double> frac={0,0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.};//Fractional steps for moving the root guess
  complex<double> m (P.size(),0);
  complex<double> dx, x1, g, gp, gm, h, b, d, f;
  double err;
  for (unsigned int it = 1; it <= MAX_IT; it++)
  {
    b = P[m.real()];
    err = abs(b);
    d = complex<double>(0,0);
    f = complex<double>(0,0);
    for ( int j=m.real()-1; j>=0; j--)
    {
      f = x*f + d;//Second derivative
      d = x*d + b;//First derivative
      b = x*b + P[j];//Polynom evaluation
      err = abs(b) + abs(x)*err;//Error term
    }
    err *= EPSS;
    if (abs(b) < err) return x;//x is already a root
     
    //Laguerre method
    g = d/b;// G(x_k) = P'(x_k) / P(x_k)
    h = g*g - 2.*f/b;// H(x_k) = G(x_k)^2 - (P''(x_k)/P(x_k))
    gp = g + sqrt((m-1.)*(m*h-g*g));
    gm = g - sqrt((m-1.)*(m*h-g*g));
    if (abs(gp) < abs(gm)) gp = gm;
    if (max(abs(gp), abs(gm)) > 0.) dx = m/gp;
    else dx = (1.+abs(x))*complex<double>(cos(1.*it), sin(1.*it));
    x1 = x - dx;
    //This procedure is reapeated until dx -> 0 or the maximum number of iterations is exceeded
    if ((x.real() == x1.real()) && (x.imag() == x1.imag())) return x;
    if (it % MT) x = x1;
    else x = x - complex<double>(frac[it/MT],0)*dx;
  }
  return x;
}


// Driver routine for the Laguerre method
deque<complex<double>> LaguerreMethod::solve_roots()
{
 deque<complex<double>> temp_poly, roots;
 complex<double> x(0,0), b, c;

 temp_poly = PolyCoeffs;

 for (int j = PolyCoeffs.size()-1; j >= 1; j--)
 {
   x=complex<double>(0,0);
   x = laguerre_core(temp_poly, x);
   if (abs(x.imag()) < 2.0*EPS*abs(x.real())) x=complex<double>(x.real(),0);//Pure real
   roots.push_back(x);//Every time laguerre() is called, it retrieves a root, no matter where the initial guess is

   //Before finding the next root, it is necessary to remove the current root from the polynomial
   b = temp_poly[j];
   for (int jj = j-1; jj>=0; jj--)
   {
     c = temp_poly[jj];
     temp_poly[jj] = b;
     b = x*b + c;
   }
   temp_poly.erase(temp_poly.begin());
 }
  //Root polish. The Laguerre's method is used again using every root as starting point. This should refine the previous result
  for (int j=1; j <= PolyCoeffs.size()-1; j++) roots[j] = laguerre_core(PolyCoeffs, roots[j]);
  return roots;
}

















