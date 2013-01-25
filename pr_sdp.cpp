/**
 * @mainpage 
 * This is an implementation of the dual only, potential reduction interior point method
 * for optimizing the second order density matrix using the P, Q, G and T_1 N-representability conditions.
 * Compiling can be done with the options PQ, PQG, PQGT1, PQGT2 and PQGT (for all conditions active) with logical consequences for the program.
 * @author Brecht Verstichel, Ward Poelmans
 * @date 22-02-2010
 */

#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * In the main the actual program is run.\n 
 * We start from the unity density matrix normed on the particle number and minimize the 
 * ojective function:\n\n
 * Tr (Gamma H) - t * ln(det P(Gamma)) \n\n
 * Once the minimum is found the parameter t is reduced and a new search is initiated,
 * this goes on until convergence is reached.\n
 * The potential is minimized using the Newton-Raphson method and the resulting linear system
 * is solved via the linear conjugate gradient method.
 */
int main(void) {

   srand(time(NULL));

   cout.precision(10);

   const int L = 20;//dim sp hilbert space
   const int N = 20;//nr of particles

   Tools::init(L,N);

   TPM::init();
   PHM::init();
   DPM::init();

   TPTPM::init();

   Gradient::init();

   DPM dpm;
   dpm.fill_Random();

   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   double **dparray = new double * [2*L];

   for(int B = 0;B < L;++B)//S = 1/2
      dparray[B] = new double [4*L4];

   for(int B = L;B < 2*L;++B)//S = 3/2
      dparray[B] = new double [L4];
/*
   dpm.convert(dparray);

   ofstream out("first.out");
   out.precision(10);

   for(int K = 0;K < L;++K)
      for(int a = 0;a < L;++a)
         for(int b = 0;b < L;++b)
            for(int c = 0;c < L;++c)
               for(int d = 0;d < L;++d)
                  out << K << "\t" << a << "\t" << b << "\t" << c << "\t" << d << "\t" << dparray[K + L][a + b*L + c*L2 + d*L3] << endl;
*/

   dpm.convert_fast(dparray);
/*
   for(int K = 0;K < L;++K)
      for(int a = 0;a < L;++a)
         for(int b = 0;b < L;++b)
            for(int c = 0;c < L;++c)
               for(int d = 0;d < L;++d)
                  cout << K << "\t" << a << "\t" << b << "\t" << c << "\t" << d << "\t" << dparray[K + L][a + b*L + c*L2 + d*L3] << endl;
*/
   //remove the array
   for(int B = 0;B < 2*L;++B)
      delete [] dparray[B];

   delete [] dparray;

/*
   Newton newton;

   //hamiltoniaan
   TPM ham;
   ham.hubbard(1.0);

   TPM rdm;
   rdm.unit();

   TPM backup_rdm(rdm);

   double t = 1.0;
   double tolerance = 1.0e-5;

   int tot_iter = 0;

   //outer iteration: scaling of the potential barrier
   //while(t > 1.0e-12){

      cout << t << "\t" << rdm.trace() << "\t" << rdm.ddot(ham) << "\t";

      int nr_newton_iter = 0;

      double convergence = 1.0;

      //inner iteration: 
      //Newton's method for finding the minimum of the current potential
      //while(convergence > tolerance){

         ++nr_newton_iter;

         SUP P;

         P.fill(rdm);

         P.invert();

         //fill the Newton object with the correct information, and solve for Delta
         newton.construct(t,ham,P);

         //dit wordt de stap:
         TPM delta;
         delta.convert(newton.gGradient());

         //line search
         double a = delta.line_search(t,P,ham);

         //rdm += a*delta;
         rdm.daxpy(a,delta);

         convergence = a*a*delta.ddot(delta);

     // }

      cout << nr_newton_iter << endl;

      t /= 2.0;

      //what is the tolerance for the newton method?
      tolerance = 1.0e-5*t;

      if(tolerance < 1.0e-12)
         tolerance = 1.0e-12;

      //extrapolatie:
      TPM extrapol(rdm);

      extrapol -= backup_rdm;

      //overzetten voor volgende stap
      backup_rdm = rdm;

      double b = extrapol.line_search(t,rdm,ham);

      rdm.daxpy(b,extrapol);

      tot_iter += nr_newton_iter;

//   }

   cout << endl;

   cout << "Final Energy:\t" << ham.ddot(rdm) << endl;
   cout << endl;
   cout << "Final Spin:\t" << rdm.S_2() << endl;

   cout << endl;
   cout << "Total nr of Newton steps = " << tot_iter << endl;
*/
   Gradient::clear();

   TPTPM::clear();

   DPM::clear();
   PHM::clear();
   TPM::clear();

   Tools::clear();

   return 0;

}
