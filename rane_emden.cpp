#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <time.h>
#include <assert.h>
#include <limits>

using namespace std;

void diff_eq(float, float[], float[]);
void rk4(float[], float[], int, float, float, float[], 
     void (*)(float, float[], float[]));

#define NSIZE 2
float y[NSIZE], dydx[NSIZE], yout[NSIZE];

int main(void)
{
  float x, error;
  const float h=1.0e-02;
  ofstream ofs( "output.dat");
  x   = h;
  float n = 1.5;
  y[0]= 1.0-h*h/6.0+n*h*h*h*h/120.0;
  y[1]= -h*h*h/3.0+n*h*h*h*h*h/30.0;
  while (x <= 2.*M_PI) {
    error = cos(x)-y[0];
    ofs << x <<" "<< y[0] <<" "<< y[1] <<" "<< error << endl;
    diff_eq(x, y, dydx);
    rk4(y, dydx, NSIZE, x, h, yout, diff_eq); 
    y[0] = yout[0];
    y[1] = yout[1];
    x += h;
  }
  return 0;
} 


void diff_eq(float x, float y[], float dydx[])
{
	dydx[0] = y[1]/(x*x);
   dydx[1] = -x*x*pow(y[0],1.5);
}



void rk4(float y[], float dydx[], int n, float x, float h, float yout[],
    void (*derivs)(float, float [], float []))
/* Given values for the variables y[1..n] and their derivatives dydx[1..n]
known at x, use the fourth-order Runge-Kutta method to advance the solution
over an interval h and return the incremented variables as yout[1..n], which
need not be a distinct array from y. The user supplies the routine
derivs(x,y,dydx), which returns derivatives dydx at x. */ 
{
  int i;
  float xh, hh, h6, dym[NSIZE], dyt[NSIZE], yt[NSIZE];
  hh = h*0.5;
  h6 = h/6.0;
  xh = x + hh;

  for (i=0; i<n; i++) {
    yt[i] = y[i] + hh*dydx[i];			// First step.
  }
  (*derivs)(xh,yt,dyt);				// Second step.
  for (i=0; i<n; i++) {
    yt[i] = y[i] + hh*dyt[i];
  }
  (*derivs)(xh,yt,dym);				// Third step.
  for (i=0; i<n; i++) {
    yt[i] = y[i] + h*dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(x+h,yt,dyt);			// Fourth step.

  for (i=0; i<n; i++) {
    yout[i] = y[i] + h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  }
}