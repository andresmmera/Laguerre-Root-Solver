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
  This function calls LaguerreMethod class for solving the root of a polynom
*/

#include <deque>
#include <complex>
#include <iostream>
#include "LaguerreMethod.h"

using namespace std;

main()
{
  deque<complex<double>> P = {1,1,1}, R;
  LaguerreMethod L(P);
  R = L.solve_roots();

  for (int i = P.size()-1; i >=1 ; i--) cout << P[i] <<"*x^" << i << " + ";
  cout << P[0] << "= 0"<< endl;
  cout << "--------- ROOTS ---------" << endl;
  for (int i = 0; i < R.size(); i++) cout << R[i] << endl;
  
  complex<double> P_x(0.,0.);
  double err = 0;
  for (int j = 0; j < R.size(); j++)
  {
    for (int i = 0; i < P.size(); i++) P_x += P[i]*pow(R[j], i);
    err += abs(P_x);
  }
  cout << endl;
  cout << "Error = sum(|P(R_i)|, i={1, N}/N = " << err/R.size() << endl;
}
