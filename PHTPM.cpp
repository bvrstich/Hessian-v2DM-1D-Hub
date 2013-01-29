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
PHTPM::PHTPM() : RecMat(PHPHM::gn(),TPTPM::gn()) { }

/**
 * copy constructor: constructs RecMat object
 * @param tpspm_c object that will be copied into this.
 */
PHTPM::PHTPM(const PHTPM &tpspm_c) : RecMat(tpspm_c){ }

/**
 * destructor
 */
PHTPM::~PHTPM(){ }

/**
 * construct a TPPHM by once tracing and once skew-tracing the direct product of two PPHM matrices, already translated to 'array' for for faster access
 */
void PHTPM::dptw(double **ppharray){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int a_,b_;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   int sign;

   for(int i = 0;i < gm();++i){

      B = PHPHM::gphmm2ph(i,0);

      S = PHM::gblock_char(B,0);

      I_i = PHPHM::gphmm2ph(i,1);
      J_i = PHPHM::gphmm2ph(i,2);

      a = PHM::gph2s(B_,I_i,0);
      b = PHM::gph2s(B_,I_i,1);
      c = PHM::gph2s(B_,J_i,0);
      d = PHM::gph2s(B_,J_i,1);

      for(int j = 0;j < gn();++j){

         B_ = TPTPM::gtpmm2t(j,0);

         S_ = TPM::gblock_char(B_,0);

         K_i = TPTPM::gtpmm2t(j,1);
         L_i = TPTPM::gtpmm2t(j,2);

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         (*this)(i,j) = 0.0;

         for(int k = 0;k < L;++k){

            int K_pph = (k + a + b)%L;

            //first S'' = 1/2
            (*this)(i,j) += 2.0 * (ppharray[K_pph][k + a*L + e*L2 + z*L3 + S*L4 + 2*S_*L4] * ppharray[K_pph][k + d*L + t*L2 + h*L3 + S*L4 + 2*S_*L4]

               + ppharray[K_pph][k + a*L + t*L2 + h*L3 + S*L4 + 2*S_*L4] * ppharray[K_pph][k + d*L + e*L2 + z*L3 + S*L4 + 2*S_*L4] );

            //first S'' = 3/2
            (*this)(i,j) += 4.0 * (ppharray[K_pph + L][k + a*L + e*L2 + z*L3] * ppharray[K_pph + L][k + d*L + t*L2 + h*L3]

                  + ppharray[K_pph + L][k + a*L + t*L2 + h*L3] * ppharray[K_pph + L][k + d*L + e*L2 + z*L3] );

         }

      }
   }

}
