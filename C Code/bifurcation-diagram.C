/* Solves the model equations of Qu et al. (2018) Siam J Appl Math, pp. 826--852
  that accounts for sexes and age class of mosquitoes, adapted to
  include two diffusively coupled habitats (A and B).  This program computes the 
  steady-state proportion of Wolbachia-infected mosquitoes in each habitat for 
  various choices of the migration parameter m.  The output can then be used to
  generate a bifurcation diagram showing infected populations versus m. */

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "parameters.h"

using namespace std;


ofstream ofsA ("infA-vs-m.txt");    // steady-state proportion of Wolbachia carriers in habitat A versus parameter m
ofstream ofsB ("infB-vs-m.txt");    // steady-state proportion of Wolbachia carriers in Habitat B versus parameter m


int main() {


  double t = 0.0;
  double dt = 0.001;
  int counter = 0;
  double m = 0.0;               //migration paramter

  double Aua, Auaold, Aub, Aubold, Awa, Awaold, Awb, Awbold;          //aquatic stage populations
  double Fua, Fuaold, Fub, Fubold, Fwa, Fwaold, Fwb, Fwbold;          //non-pregnant female populations
  double Fpua, Fpuaold, Fpub, Fpubold, Fpwa, Fpwaold, Fpwb, Fpwbold;  //pregnant female populations
  double Mua, Muaold, Mub, Mubold, Mwa, Mwaold, Mwb, Mwbold;          //male populations
  double Fpsa, Fpsaold, Fpsb, Fpsbold;                                //sterile female populations

  double infectedA, uninfectedA, infectedB, uninfectedB;              //for computing overall percent of infected ADULTS in each habitat


  /*  forward Euler ODE solver */


while (m <= 0.005) {

  /*  initial conditions */

  Auaold = 200000.0, Awaold = 300000.0, Fuaold = 300000.0, Fwaold = 500000.0, Fpuaold = 100000.0;
  Fpwaold = 100000.0, Muaold = 400000.0, Mwaold = 0.0, Fpsaold = 0.0;

  Aubold = 200000.0, Awbold = 0.0, Fubold = 800000.0, Fwbold = 0.0, Fpubold = 100000.0;
  Fpwbold = 0.0, Mubold = 400000.0, Mwbold = 0.0, Fpsbold = 0.0;

  t = 0.0;
  counter = 0;

  while (t < 5000.0) {

    t = t + dt;
    counter++;
    Aua = Auaold + dt*((phiu*Fpuaold + vu*phiw*Fpwaold)*(1.0 - (Auaold + Awaold)/Ka) - (mmua + psi)*Auaold);
    Awa = Awaold + dt*(vw*phiw*Fpwaold*(1.0 - (Auaold + Awaold)/Ka) - (mmua + psi)*Awaold);
    Fua = Fuaold + dt*(bf*psi*Auaold - (sigma + mufu)*Fuaold + m*(Fubold - Fuaold));
    Fwa = Fwaold + dt*(bf*psi*Awaold - (sigma + mufw)*Fwaold + m*(Fwbold - Fwaold));
    Fpua = Fpuaold + dt*(sigma*Fuaold*Muaold/(Muaold + Mwaold) - mufu*Fpuaold + m*(Fpubold - Fpuaold));
    Fpwa = Fpwaold + dt*(sigma*Fwaold - mufw*Fpwaold + m*(Fpwbold - Fpwaold));
    Mua = Muaold + dt*(bm*psi*Auaold - mumu*Muaold + m*(Mubold - Muaold));
    Mwa = Mwaold + dt*(bm*psi*Awaold - mumw*Mwaold + m*(Mwbold - Mwaold));
    Fpsa = Fpsaold + dt*(sigma*Fuaold*Mwaold/(Muaold + Mwaold) - mufu*Fpsaold + m*(Fpsbold - Fpsaold));

    Aub = Aubold + dt*((phiu*Fpubold + vu*phiw*Fpwbold)*(1.0 - (Aubold + Awbold)/Ka) - (mmua + psi)*Aubold);
    Awb = Awbold + dt*(vw*phiw*Fpwbold*(1.0 - (Aubold + Awbold)/Ka) - (mmua + psi)*Awbold);
    Fub = Fubold + dt*(bf*psi*Aubold - (sigma + mufu)*Fubold + m*(Fuaold - Fubold));
    Fwb = Fwbold + dt*(bf*psi*Awbold - (sigma + mufw)*Fwbold + m*(Fwaold - Fwbold));
    Fpub = Fpubold + dt*(sigma*Fubold*Mubold/(Mubold + Mwbold) - mufu*Fpubold + m*(Fpuaold - Fpubold));
    Fpwb = Fpwbold + dt*(sigma*Fwbold - mufw*Fpwbold + m*(Fpwaold - Fpwbold));
    Mub = Mubold + dt*(bm*psi*Aubold - mumu*Mubold + m*(Muaold - Mubold));
    Mwb = Mwbold + dt*(bm*psi*Awbold - mumw*Mwbold + m*(Mwaold - Mwbold));
    Fpsb = Fpsbold + dt*(sigma*Fubold*Mwbold/(Mubold + Mwbold) - mufu*Fpsbold + m*(Fpsaold - Fpsbold));

    infectedA   = Awa + Fwa + Fpwa + Mwa;        // total infected with Wolbachia in habitat A
    uninfectedA = Aua + Fua + Fpua + Mua + Fpsa; // total uninfected with Wolbachia in habitat A
    infectedB   = Awb + Fwb + Fpwb + Mwb;        // total infected with Wolbachia in habitat B
    uninfectedB = Aub + Fub + Fpub + Mub + Fpsb; // total uninfected with Wolbachia in habitat B


    Auaold = Aua,  Awaold = Awa,  Fuaold = Fua,  Fwaold = Fwa,  Fpuaold = Fpua,  Fpwaold = Fpwa,  Muaold = Mua,  Mwaold = Mwa,  Fpsaold = Fpsa;
    Aubold = Aub,  Awbold = Awb,  Fubold = Fub,  Fwbold = Fwb,  Fpubold = Fpub,  Fpwbold = Fpwb,  Mubold = Mub,  Mwbold = Mwb,  Fpsbold = Fpsb;

  }

  ofsA << m << "\t" << (infectedA/(infectedA + uninfectedA)) << endl; 
  ofsB << m << "\t" << (infectedB/(infectedB + uninfectedB)) << endl; 
  m = m + 0.00001;
  cout << "m = " << m << endl;

}


  /* Comment out the next line if you prefer to use software other than gnuplot to
     generate figures from the .txt files generated by this program.  If you do have
     gnuplot installed, you may need to change the path for calling gnuplot.  */

  system ("/Applications/Gnuplot.app/Contents/Resources/bin/gnuplot bifplot");

  return 0;

}
