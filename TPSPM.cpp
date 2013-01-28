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

/**
 * construct a TPSPM object by triple-tracing the direct product of two DPM objects
 */
void TPSPM::dpt3(double scale,const DPM &dpm){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   double **dparray = new double * [2*L];

   for(int B = 0;B < L;++B)//S = 1/2
      dparray[B] = new double [4*L4];

   for(int B = L;B < 2*L;++B)//S = 3/2
      dparray[B] = new double [L4];

   dpm.convert(dparray);

   int B,I,J;

   int S;

   int a,b,c,d;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int e = 0;e < L;++e){

         double ward = 0.0;

         //first S = 1/2
         for(int k = 0;k < L;++k){//abk and cdk

            int K = (a + b + k)%L;

            //first S_ep = 0
            for(int p = 0;p < e;++p)
               ward += dparray[K][a + b*L + e*L2 + p*L3 + S*L4] * dparray[K][c + d*L + e*L2 + p*L3 + S*L4];

            ward += 2.0 * dparray[K][a + b*L + e*L2 + e*L3 + S*L4] * dparray[K][c + d*L + e*L2 + e*L3 + S*L4];

            for(int p = e + 1;p < L;++p)
               ward += dparray[K][a + b*L + e*L2 + p*L3 + S*L4] * dparray[K][c + d*L + e*L2 + p*L3 + S*L4];


            //then S_ep = 1
            for(int p = 0;p < e;++p)
               ward += dparray[K][a + b*L + e*L2 + p*L3 + S*L4 + 2*L4] * dparray[K][c + d*L + e*L2 + p*L3 + S*L4 + 2*L4];

            for(int p = e + 1;p < L;++p)
               ward += dparray[K][a + b*L + e*L2 + p*L3 + S*L4 + 2*L4] * dparray[K][c + d*L + e*L2 + p*L3 + S*L4 + 2*L4];

         }

         (*this)(i,e) = 2.0/(2*S + 1.0) * ward;

         //then S = 3/2, only when
         if(S == 1){

            ward = 0.0;

            for(int k = 0;k < L;++k){//abk and cdk

               int K = (a + b + k)%L;

               //only S_ep = 1 term possible
               for(int p = 0;p < e;++p)
                  ward += dparray[K + L][a + b*L + e*L2 + p*L3] * dparray[K + L][c + d*L + e*L2 + p*L3];

               for(int p = e + 1;p < L;++p)
                  ward += dparray[K + L][a + b*L + e*L2 + p*L3] * dparray[K + L][c + d*L + e*L2 + p*L3];

            }

            (*this)(i,e) += 4.0/3.0 * ward;

         }

         (*this)(i,e) *= scale;

      }

   }

   //remove the array
   for(int B = 0;B < 2*L;++B)
      delete [] dparray[B];

   delete [] dparray;


}

/**
 * construct a TPSPM object by triple-tracing the direct product of two DPM objects, here the DPM has already been transformed to a double **array
 */
void TPSPM::dpt3(double scale,double **dparray){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   int B,I,J;

   int S;

   int a,b,c,d;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int e = 0;e < L;++e){

         double ward = 0.0;

         //first S = 1/2
         for(int k = 0;k < L;++k){//abk and cdk

            int K = (a + b + k)%L;

            //first S_ep = 0
            for(int p = 0;p < e;++p)
               ward += dparray[K][a + b*L + e*L2 + p*L3 + S*L4] * dparray[K][c + d*L + e*L2 + p*L3 + S*L4];

            ward += 2.0 * dparray[K][a + b*L + e*L2 + e*L3 + S*L4] * dparray[K][c + d*L + e*L2 + e*L3 + S*L4];

            for(int p = e + 1;p < L;++p)
               ward += dparray[K][a + b*L + e*L2 + p*L3 + S*L4] * dparray[K][c + d*L + e*L2 + p*L3 + S*L4];


            //then S_ep = 1
            for(int p = 0;p < e;++p)
               ward += dparray[K][a + b*L + e*L2 + p*L3 + S*L4 + 2*L4] * dparray[K][c + d*L + e*L2 + p*L3 + S*L4 + 2*L4];

            for(int p = e + 1;p < L;++p)
               ward += dparray[K][a + b*L + e*L2 + p*L3 + S*L4 + 2*L4] * dparray[K][c + d*L + e*L2 + p*L3 + S*L4 + 2*L4];

         }

         (*this)(i,e) = 2.0/(2*S + 1.0) * ward;

         //then S = 3/2, only when
         if(S == 1){

            ward = 0.0;

            for(int k = 0;k < L;++k){//abk and cdk

               int K = (a + b + k)%L;

               //only S_ep = 1 term possible
               for(int p = 0;p < e;++p)
                  ward += dparray[K + L][a + b*L + e*L2 + p*L3] * dparray[K + L][c + d*L + e*L2 + p*L3];

               for(int p = e + 1;p < L;++p)
                  ward += dparray[K + L][a + b*L + e*L2 + p*L3] * dparray[K + L][c + d*L + e*L2 + p*L3];

            }

            (*this)(i,e) += 4.0/3.0 * ward;

         }

         (*this)(i,e) *= scale;

      }

   }

}

/**
 * construct a TPSPM object by triple skew-tracing the direct product of two PPHM objects, here the PPHM has already been transformed to a double **array
 */
void TPSPM::dpw3(double scale,double **ppharray){

}

/**
 * construct a TPSPM object by once tracing, and twice skew-tracing the direct product of two PPHM objects, here the PPHM has already been transformed to a double **array
 */
void TPSPM::dptw2(double scale,double **ppharray){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   int B,I,J;

   int S;

   int a,b,c,d;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int e = 0;e < L;++e){

         double ward = 0.0;

         //first S = 1/2
         for(int S_kl = 0;S_kl < 2;++S_kl)
            for(int k = 0;k < L;++k)
               for(int l = k + S_kl;l < k;++l){

                  int K_pph = (k + l + e)%L;

                  ward += ppharray[K_pph][a + b*L + k*L2 + l*L3 + S*L4 + 2*S_kl*L4] * ppharray[K_pph][c + d*L + k*L2 + l*L3 + S*L4 + 2*S_kl*L4];

            }

         (*this)(i,e) = 2.0/(2.0*S + 1.0) * ward;

         if(S == 1){

            ward = 0.0;

            for(int k = 0;k < L;++k)
               for(int l = k + 1;l < k;++l){

                  int K_pph = (k + l + e)%L;

                  //first S_kl = 0
                  ward += ppharray[K_pph + L][a + b*L + k*L2 + l*L3] * ppharray[K_pph + L][c + d*L + k*L2 + l*L3];

               }

            (*this)(i,e) += 4.0/9.0 * ward;

         }

         (*this)(i,e) *= scale;

      }
   }

}
