/***************************************************************************
                                main.cpp
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

/*
  This function calls LaguerreMethod class for solving the root of a polynom.
  The polynom coefficients are defined at P such as {1,2,3} would represent 3x^2 + 2x + 1 = 0
*/

#include <deque>
#include <complex>
#include <iostream>
#include "LaguerreMethod.h"

using namespace std;

main()
{
  deque<complex<double>> P = {1,2,3}, R;
  LaguerreMethod L(P);
  R = L.solve_roots();

  // Display the equation to solve
  for (int i = P.size()-1; i >=1 ; i--) cout << P[i] <<"*x^" << i << " + ";
  cout << P[0] << "= 0"<< endl;

  //Display roots
  cout << "--------- ROOTS ---------" << endl;
  for (int i = 0; i < R.size(); i++) cout << R[i] << endl;
  
  //In order to gather the quality of the roots found, the lines below evaluate the polynom at every
  //root, R_i, and takes the maximum and the mean deviation to zero 
  cout << endl << endl;
  complex<double> P_x(0.,0.);
  double mean_err = 0, max_err = -1e12, abs_err;
  cout << "--------- Root analysis ---------" << endl;
  for (int j = 0; j < R.size(); j++)
  { 
    P_x = complex<double>(0,0);
    for (int i = 0; i < P.size(); i++) P_x += P[i]*pow(R[j], i);
    abs_err = abs(P_x);
    if (abs(P_x) > max_err) max_err = abs_err;
    mean_err += abs_err;
    cout << "Error{Root[" << j << "]}= " << abs_err << endl;
  }
  cout << endl;
  cout << "Mean error = sum(|P(R_i)|, i={1, N}/N = " << mean_err/R.size() << endl;
  cout << "Maximum error = max(|P(R_i)|), i={1, N}/N = " << max_err << endl;
}
