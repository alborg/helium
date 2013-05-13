#include "minimise.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include "slaterdeterminant.h"
#include "correlation.h"
#include "hamiltonian.h"

using namespace  std;
using namespace arma;


minimise::minimise(int nP):
    nParticles(nP)
{
}

//This file is based on the model example from FYS4411 webpage.


double minimise::steepestDescent(vec(*func)(double alpha, double beta, int mpi_steps,
                                 slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr,
                                 double *allEnergies), double alpha, double beta, int min_steps,
                                 slaterDeterminant *slater, Hamiltonian *hamiltonian, correlation *corr) {
    double gtol = 1e-5;
    double step = 0.1;
    double* allEnergies = new double[min_steps+1];
    double gradient = 0;
    double alpha_new = alpha;
    double alpha_old;
    int max_iter = 20;
    int iter = 0;

    while(iter < max_iter && abs(alpha_new - alpha) > gtol) {

        alpha_old = alpha_new;

        vec energies = (*func)(alpha_old, beta, min_steps, slater, hamiltonian, corr, allEnergies);
        vec energiesPsi = (*func)(alpha_old, beta, min_steps, slater, hamiltonian, corr, allEnergies);

        gradient = (2*energies(2)*energies(0) - energiesPsi(3))/(min_steps * nParticles);

        alpha_new = alpha_old - step*gradient;
        cout << alpha_new<<endl;

        ++iter;
    }

    cout <<"Alpha, iter: "<<alpha_new<<" "<<iter<<endl;

    return alpha_new;
}



static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))



#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0


void minimise::dfpmin(int var,const mat &r,double var2, slaterDeterminant *slater,correlation *corr,
                      vec &p, int n, double gtol, int *iter, double *fret, Hamiltonian *hamiltonian)
{

    double alpha, beta;

    if(var == 1) {alpha = p(0); beta = var2;}
    else {beta = p(0); alpha = var2;}

  int check,i,its,j;
  double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
  vec dg = zeros<vec>(n);
  vec g = zeros<vec>(n);
  vec hdg = zeros<vec>(n);
  vec pnew = zeros<vec>(n);
  vec xi = zeros<vec>(n);
  mat hessian = zeros<vec>(n,n);


  //fp=;
  //fp = (*func)(p);
  cout<<"alpha, beta, fp: "<<alpha<<" "<<beta<<" "<<fp<<endl;
  //g(0) = ;
  //(*dfunc)(p,g);
  cout<<"g(0): "<<g(0)<<endl;
  for (i = 0;i < n;i++) {
    for (j = 0; j< n;j++) hessian(i,j)=0.0;
    hessian(i,i)=1.0;
    xi(i) = -g(i);
    sum += p(i)*p(i);
  }
  stpmax=STPMX*FMAX(sqrt(sum),(double)n);
  for (its=1;its<=ITMAX;its++) {
    *iter=its;
    lnsrch(var,r,var2, slater,corr,1,p,fp,g,xi,pnew,fret,stpmax,&check,hamiltonian);
    fp = *fret;
    for (i = 0; i< n;i++) {
      xi(i)=pnew(i)-p(i);
      p(i)=pnew(i);
    }
    test=0.0;
    for (i = 0;i< n;i++) {
      temp=fabs(xi(i))/FMAX(fabs(p(i)),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) {
      return;
    }
    for (i=0;i<n;i++) dg(i)=g(i);
    //g(0) = hamiltonian->dlocalEnergy(var,r,alpha,beta,slater,corr);
    //(*dfunc)(p,g);
    test=0.0;
    den=FMAX(*fret,1.0);
    for (i=0;i<n;i++) {
      temp=fabs(g(i))*FMAX(fabs(p(i)),1.0)/den;
      if (temp > test) test=temp;
    }
    if (test < gtol) {
      return;
    }
    for (i=0;i<n;i++) dg(i)=g(i)-dg(i);
    for (i=0;i<n;i++) {
      hdg(i)=0.0;
      for (j=0;j<n;j++) hdg(i) += hessian(i,j)*dg(j);
    }
    fac=fae=sumdg=sumxi=0.0;
    for (i=0;i<n;i++) {
      fac += dg(i)*xi(i);
      fae += dg(i)*hdg(i);
      sumdg += SQR(dg(i));
      sumxi += SQR(xi(i));
    }
    if (fac*fac > EPS*sumdg*sumxi) {
      fac=1.0/fac;
      fad=1.0/fae;
      for (i=0;i<n;i++) dg(i)=fac*xi(i)-fad*hdg(i);
      for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      hessian(i,j) += fac*xi(i)*xi(j)
        -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j);
    }
      }
    }
    for (i=0;i<n;i++) {
      xi(i)=0.0;
      for (j=0;j<n;j++) xi(i) -= hessian(i,j)*g(j);
    }
  }
  cout << "too many iterations in dfpmin" << endl;
}
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

#define ALF 1.0e-4
#define TOLX 1.0e-7

void minimise::lnsrch(int var,const mat &r,double var2, slaterDeterminant *slater,correlation *corr,
                      int n, vec &xold, double fold, vec &g, vec &p, vec &x,
                      double *f, double stpmax, int *check, Hamiltonian *hamiltonian)
{

    double alpha, beta;
  int i;
  double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
    test,tmplam;

  *check=0;
  for (sum=0.0,i=0;i<n;i++) sum += p(i)*p(i);
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=0;i<n;i++) p(i) *= stpmax/sum;
  for (slope=0.0,i=0;i<n;i++)
    slope += g(i)*p(i);
  test=0.0;
  for (i=0;i<n;i++) {
    temp=fabs(p(i))/FMAX(fabs(xold(i)),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) {
    for (i=0;i<n;i++) x(i)=xold(i)+alam*p(i);
    if(var == 1) {alpha = x(i); beta = var2;}
    else {beta = x(i); alpha = var2;}
    //*f=(*func)(x);
    //*f=hamiltonian->localEnergy(r,alpha,beta,slater,corr);
    //*f=hamiltonian->analyticEnergyHe(r,alpha,beta);
    if (alam < alamin) {
      for (i=0;i<n;i++) x(i)=xold(i);
      *check=1;
      return;
    } else if (*f <= fold+ALF*alam*slope) return;
    else {
      if (alam == 1.0)
    tmplam = -slope/(2.0*(*f-fold-slope));
      else {
    rhs1 = *f-fold-alam*slope;
    rhs2=f2-fold2-alam2*slope;
    a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
    b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
    if (a == 0.0) tmplam = -slope/(2.0*b);
    else {
      disc=b*b-3.0*a*slope;
      if (disc<0.0) cout << "Roundoff problem in lnsrch." << endl;
      else tmplam=(-b+sqrt(disc))/(3.0*a);
    }
    if (tmplam>0.5*alam)
      tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=FMAX(tmplam,0.1*alam);
  }
}
#undef ALF
#undef TOLX

////  this function defines the Energy function
//double E_function(Vector  &x)
//{
//  double value = x(0)*x(0)*0.5+1.0/(8*x(0)*x(0));
//  return value;
//} // end of function to evaluate

////  this function defines the derivative of the energy
//void dE_function(Vector &x, Vector &g)
//{

//  g(0) = x(0)-1.0/(4*x(0)*x(0)*x(0));

//} // end of function to evaluate


//   Main function begins here
//int main()
//{
//     int n, iter;
//     double gtol, fret;
//     double alpha;
//     n = 1;
////   reserve space in memory for vectors containing the variational
////   parameters
//     Vector g(n), p(n);
//     cout << "Read in guess for alpha" << endl;
//     cin >> alpha;
//     gtol = 1.0e-5;
////   now call dfmin and compute the minimum
//     p(0) = alpha;
//     dfpmin(p, n, gtol, &iter, &fret, E_function, dE_function);
//     cout << "Value of energy minimum = " << fret << endl;
//     cout << "Number of iterations = " << iter << endl;
//     cout << "Value of alpha at minimum = " << p(0) << endl;
//      return 0;
//}  // end of main program



