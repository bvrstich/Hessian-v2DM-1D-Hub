#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

/**
 * standard constructor:
 */
TPSPM::TPSPM() : RecMat(TPTPM::gn(),Tools::gL()) { }

/**
 * copy constructor: constructs RecMat object
 * @param tpspm_c object that will be copied into this.
 */
TPSPM::TPSPM(const TPSPM &tpspm_c) : RecMat(tpspm_c){ }

/**
 * destructor
 */
TPSPM::~TPSPM(){ }

/**
 * construct a TPSPM by tracing the direct product of two TPM's
 */
void TPSPM::dpt(double scale,const TPM &Q){

  int B,I,J;

  int S;

  int K,a,b,c,d;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);
      K = TPM::gblock_char(B,1);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int k = 0;k < Tools::gL();++k){

         int l = (K - k + Tools::gL())%Tools::gL();

         if(k == l)
            (*this)(i,k) = 4.0 * scale * Q(S,a,b,k,l) * Q(S,c,d,k,l);
         else
            (*this)(i,k) = 2.0 * scale * Q(S,a,b,k,l) * Q(S,c,d,k,l);

      }
   }

}

/**
 * construct the singly-traced antisymmetrized direct product of two PHM matrices
 */
void TPSPM::dpt(double scale,const PHM &phm){

   int L = Tools::gL();

   double **pharray = new double * [2*L];

   for(int B = 0;B < 2*L;++B)
      pharray[B] = new double [L*L];

   phm.convert(pharray);

   int B,I,J;

   int S;

   int sign;

   int a,b,c,d;
   int c_,d_;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      sign = 1 - 2*S;

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      c_ = Tools::par(c);
      d_ = Tools::par(d);

      for(int k = 0;k < Tools::gL();++k){

         (*this)(i,k) = 0.0;

         for(int Z = 0;Z < 2;++Z){

            //(a,d,c,b)
            int P = (a + d_)%L;

            double ward = pharray[P + Z*L][a + k*L] * pharray[P + Z*L][c + k*L];

            //(b,d,c,a)
            P = (b + d_)%L;

            ward += sign * pharray[P + Z*L][b + k*L] * pharray[P + Z*L][c + k*L];

            //(a,c,d,b)
            P = (a + c_)%L;

            ward += sign * pharray[P + Z*L][a + k*L] * pharray[P + Z*L][d + k*L];

            //(b,c,d,a)
            P = (b + c_)%L;

            ward += pharray[P + Z*L][b + k*L] * pharray[P + Z*L][d + k*L];

            (*this)(i,k) += 2.0 * (2.0*Z + 1.0) * Tools::g6j(0,0,Z,S) * ward;

         }

         (*this)(i,k) *= scale;

      }

   }

   for(int B = 0;B < 2*L;++B)
      delete [] pharray[B];

   delete [] pharray;

}
