/* Solve the model equations of Hu et al., Bulletin of Mathematical Biology 
   (2021) 83:58, adapted to include two diffusively-coupled habitats (A and B).  */

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>


#define bI      0.285000
#define bU      0.300000
#define deltaI  0.079000
#define deltaU  0.071000
#define CCA	400.0000
#define CCB	500.0000
#define m	0.001000

using namespace std;


ofstream ofsxA1 ("xA1.txt");    // first set of initial conditions
ofstream ofsyA1 ("yA1.txt");
ofstream ofsxB1 ("xB1.txt");
ofstream ofsyB1 ("yB1.txt");

ofstream ofsxA2 ("xA2.txt");    // second set of initial conditions
ofstream ofsyA2 ("yA2.txt");
ofstream ofsxB2 ("xB2.txt");
ofstream ofsyB2 ("yB2.txt");

ofstream ofsxA3 ("xA3.txt");    // third set of initial conditions
ofstream ofsyA3 ("yA3.txt");
ofstream ofsxB3 ("xB3.txt");
ofstream ofsyB3 ("yB3.txt");


int main() {


  double t = 0.0;
  double dt = 0.01;
  int counter = 0;
  double dA = ((bU - deltaU)*CCA - m*(CCA - CCB))/(CCA*CCA);
  double dB = ((bU - deltaU)*CCB + m*(CCA - CCB))/(CCB*CCB);
  cout << "dA = " << dA << endl;
  cout << "dB = " << dB << endl;
  double ratio = (bI-deltaI+deltaU)/(bU-bI+deltaI-deltaU);
  cout << "ratio = " << ratio << endl;

  double xA, xAold, yA, yAold;          // Habitat A:  x = infected with Wolbachia, y = uninfected
  double xB, xBold, yB, yBold;          // Habitat B:  x = infected with Wolbachia, y = uninfected


  /* first set of initial conditions  */

  t = 0.0;
  counter = 0.0;
  xAold = 28.0,  yAold = 400.0,  xBold = 36.0,  yBold = 500.0;

  /*  forward Euler ODE solver */

  while (t < 250.0) {

    t = t + dt;
    counter++;
    xA = xAold + dt*((bI - deltaI)*xAold - dA*xAold*(xAold + yAold) + m*(xBold - xAold));
    yA = yAold + dt*(bU*yAold*yAold/(xAold + yAold) - deltaU*yAold - dA*yAold*(xAold + yAold) + m*(yBold - yAold));
    xB = xBold + dt*((bI - deltaI)*xBold - dB*xBold*(xBold + yBold) + m*(xAold - xBold));
    yB = yBold + dt*(bU*yBold*yBold/(xBold + yBold) - deltaU*yBold - dB*yBold*(xBold + yBold) + m*(yAold - yBold));

    if ((counter%10) == 0) {
      ofsxA1 << t << "\t" << xA << endl;
      ofsxB1 << t << "\t" << xB << endl;
      ofsyA1 << t << "\t" << yA << endl;
      ofsyB1 << t << "\t" << yB << endl;
    }

    xAold = xA;
    yAold = yA;
    xBold = xB;
    yBold = yB;

  }

  /* second set of initial conditions  */

  t = 0.0;
  counter = 0.0;
  xAold = 33.213,  yAold = 400.0,  xBold = 41.516,  yBold = 500.0;

  /*  forward Euler ODE solver */

  while (t < 250.0) {

    t = t + dt;
    counter++;
    xA = xAold + dt*((bI - deltaI)*xAold - dA*xAold*(xAold + yAold) + m*(xBold - xAold));
    yA = yAold + dt*(bU*yAold*yAold/(xAold + yAold) - deltaU*yAold - dA*yAold*(xAold + yAold) + m*(yBold - yAold));
    xB = xBold + dt*((bI - deltaI)*xBold - dB*xBold*(xBold + yBold) + m*(xAold - xBold));
    yB = yBold + dt*(bU*yBold*yBold/(xBold + yBold) - deltaU*yBold - dB*yBold*(xBold + yBold) + m*(yAold - yBold));

    if ((counter%10) == 0) {
      ofsxA2 << t << "\t" << xA << endl;
      ofsxB2 << t << "\t" << xB << endl;
      ofsyA2 << t << "\t" << yA << endl;
      ofsyB2 << t << "\t" << yB << endl;
    }

    xAold = xA;
    yAold = yA;
    xBold = xB;
    yBold = yB;

  }

  /* third set of initial conditions  */

  t = 0.0;
  counter = 0.0;
  xAold = 34.0,  yAold = 400.0,  xBold = 42.0,  yBold = 500.0;

  /*  forward Euler ODE solver */

  while (t < 280.0) {

    t = t + dt;
    counter++;
    xA = xAold + dt*((bI - deltaI)*xAold - dA*xAold*(xAold + yAold) + m*(xBold - xAold));
    yA = yAold + dt*(bU*yAold*yAold/(xAold + yAold) - deltaU*yAold - dA*yAold*(xAold + yAold) + m*(yBold - yAold));
    xB = xBold + dt*((bI - deltaI)*xBold - dB*xBold*(xBold + yBold) + m*(xAold - xBold));
    yB = yBold + dt*(bU*yBold*yBold/(xBold + yBold) - deltaU*yBold - dB*yBold*(xBold + yBold) + m*(yAold - yBold));

    if ((counter%10) == 0) {
      ofsxA3 << t << "\t" << xA << endl;
      ofsxB3 << t << "\t" << xB << endl;
      ofsyA3 << t << "\t" << yA << endl;
      ofsyB3 << t << "\t" << yB << endl;
    }

    xAold = xA;
    yAold = yA;
    xBold = xB;
    yBold = yB;

  }

  /* Comment out the next line if you prefer to use software other than gnuplot to
     generate figures from the .txt files generated by this program.  If you do have
     gnuplot installed, you may need to change the path for calling gnuplot.  */
     

  system ("/Applications/Gnuplot.app/Contents/Resources/bin/gnuplot huplot");


  return 0;

}



