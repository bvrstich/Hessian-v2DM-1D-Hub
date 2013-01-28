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
Hessian::Hessian() : Matrix(TPTPM::gn() + 1) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix hess_c
 * @param hess_c object that will be copied into this.
 */
Hessian::Hessian(const Hessian &hess_c) : Matrix(hess_c){ }

/**
 * destructor
 */
Hessian::~Hessian(){ }

ostream &operator<<(ostream &output,const Hessian &hess_p){

   int B,I,J,B_,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);
      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);

      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int j = i;j < TPTPM::gn();++j){

         B_ = TPTPM::gtpmm2t(j,0); 
         K = TPTPM::gtpmm2t(j,1);
         L = TPTPM::gtpmm2t(j,2);

         e = TPM::gt2s(B_,K,0);
         z = TPM::gt2s(B_,K,1);

         t = TPM::gt2s(B_,L,0);
         h = TPM::gt2s(B_,L,1);

         if(fabs(hess_p(i,j)) < 1.0e-14){

            output << i << "\t" << j << "\t|\t(" << B << ")\t" << I << "\t" << J << "\t(" << B_ << ")\t" << K << "\t" << L << "\t|\t" << 

               "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << 0 << endl;

         }
         else{
            output << i << "\t" << j << "\t|\t(" << B << ")\t" << I << "\t" << J << "\t(" << B_ << ")\t" << K << "\t" << L << "\t|\t" << 

               "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << hess_p(i,j) << endl;

         }

      }

   }

   return output;

}

/**
 * construct the I part of the hessian matrix
 */
void Hessian::I(const TPM &tpm){

   int B,I,J,B_,K,L;

   int S;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);
      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      for(int j = i;j < TPTPM::gn();++j){

         B_ = TPTPM::gtpmm2t(j,0);
         K = TPTPM::gtpmm2t(j,1);
         L = TPTPM::gtpmm2t(j,2);

         if(B == B_)
            (*this)(i,j) = 2.0 * (2.0*S + 1.0) * Gradient::gnorm(i) * Gradient::gnorm(j) * ( tpm(B,I,K) * tpm(B,J,L) + tpm(B,I,L) * tpm(B,J,K) );
         else
            (*this)(i,j) = 0.0;

      }

   }

}

/**
 * construct the lagrange multiplier part of the Hessian
 */
void Hessian::lagr(){

   int B,I,J;
   int S;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);
      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      if(I == J)
         (*this)(i,TPTPM::gn()) = std::sqrt(2.0*S + 1.0);
      else
         (*this)(i,TPTPM::gn()) = 0.0;

   }

}

/**
 * construct the Q part of the hessian
 */
void Hessian::Q(const TPM &Q){

   int N = Tools::gN();

   TPM Q2;
   Q2.squaresym(Q);

   SPM Q2bar;
   Q2bar.bar(8.0/(N*(N - 1.0)*(N - 1.0)),Q2);

   double Q2trace = 16 * Q2.trace()/ (N*N*(N - 1.0)*(N - 1.0));

   TPSPM dpt;
   dpt.dpt(1.0/(N - 1.0),Q);

   SPSPM dpt2;
   dpt2.dpt2(1.0/((N - 1.0)*(N - 1.0)),Q);

   int B,I,J,B_,K,L;

   int S,S_;

   //first store everything in ward, then multiply with norms and add to (*this)!
   double ward;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      S = TPM::gblock_char(B,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int j = i;j < TPTPM::gn();++j){

         B_ = TPTPM::gtpmm2t(j,0);

         S_ = TPM::gblock_char(B_,0);

         K = TPTPM::gtpmm2t(j,1);
         L = TPTPM::gtpmm2t(j,2);

         e = TPM::gt2s(B_,K,0);
         z = TPM::gt2s(B_,K,1);
         t = TPM::gt2s(B_,L,0);
         h = TPM::gt2s(B_,L,1);

         ward = 0.0;

         if(B == B_)
            ward += 2.0 / ( 2.0*S + 1.0 ) * ( Q(B,I,K) * Q(B,J,L) +  Q(B,I,L) * Q(B,J,K) );

         if(I == J){

            if(K == L){

               ward += Q2trace - Q2bar[a] - Q2bar[b] - Q2bar[e] - Q2bar[z];

               ward += dpt2(a,e) + dpt2(a,z) + dpt2(b,e) + dpt2(b,z);

            }

            ward -= dpt(j,a) + dpt(j,b);

            ward += 8.0/ ( Tools::gN() * (Tools::gN() - 1.0) ) * Q2(S_,e,z,t,h);

         }

         if(K == L){

            ward -= dpt(i,e) + dpt(i,z);

            ward += 8.0/ ( Tools::gN() * (Tools::gN() - 1.0) ) * Q2(S,a,b,c,d);

         }

         //finally
         (*this)(i,j) += ward * Gradient::gnorm(i) * Gradient::gnorm(j) * (2.0*S + 1.0) * (2.0*S_ + 1.0);

      } 
   }

}

/**
 * construct the G part of the Hessian
 */
void Hessian::G(const PHM &G){

   TPTPM dp;
   dp.dp(G);

   TPSPM dpt;
   dpt.dpt(1.0/(Tools::gN() - 1.0),G);

   SPSPM dpt2;
   dpt2.dpt2(1.0/((Tools::gN() - 1.0)*(Tools::gN() - 1.0)),G);

   int B,I,J,B_,K,L;

   int S,S_;

   //first store everything in ward, then multiply with norms and add to (*this)!
   double ward;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      S = TPM::gblock_char(B,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int j = i;j < TPTPM::gn();++j){

         B_ = TPTPM::gtpmm2t(j,0);

         S_ = TPM::gblock_char(B_,0);

         K = TPTPM::gtpmm2t(j,1);
         L = TPTPM::gtpmm2t(j,2);

         e = TPM::gt2s(B_,K,0);
         z = TPM::gt2s(B_,K,1);
         t = TPM::gt2s(B_,L,0);
         h = TPM::gt2s(B_,L,1);

         ward = 2.0 * TPM::gnorm(a,b) * TPM::gnorm(c,d) * TPM::gnorm(e,z) * TPM::gnorm(t,h) * dp(i,j);

         if(I == J){

            if(K == L)
               ward += dpt2(a,e) + dpt2(a,z) + dpt2(b,e) + dpt2(b,z);

            ward -= TPM::gnorm(e,z) * TPM::gnorm(t,h) * ( dpt(j,a) +  dpt(j,b) );

         }

         if(K == L)
            ward -= TPM::gnorm(a,b) * TPM::gnorm(c,d) * ( dpt(i,e) +  dpt(i,z) );

         //the norms
         (*this)(i,j) += ward * Gradient::gnorm(i) * Gradient::gnorm(j) * (2.0*S + 1.0) * (2.0*S_ + 1.0);

      }
   }

}

/**
 * construct the T1 part of the Hessian matrix
 */
void Hessian::T(const DPM &T){

   int N = Tools::gN();

   DPM T2;
   T2.squaresym(T);

   TPM T2bar;
   T2bar.bar(8.0/(N*(N - 1.0)),T2);

   SPM T2barbar;
   T2barbar.bar(1.0/(2*(N - 1.0)),T2bar);

   double T2trace = 16 * T2.trace()/ (N*N*(N - 1.0)*(N - 1.0));

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   double **dparray = new double * [2*L];

   for(int B = 0;B < L;++B)//S = 1/2
      dparray[B] = new double [4*L4];

   for(int B = L;B < 2*L;++B)//S = 3/2
      dparray[B] = new double [L4];

   T.convert_fast(dparray);

   TPTPM dpt2;
   dpt2.dpt2(dparray);

   TPSPM dpt3;
   dpt3.dpt3(1.0/(N - 1.0),dparray);

   SPSPM dpt4;
   dpt4.dpt4(0.5/( (N - 1.0) * (N - 1.0) ),dparray);

   int B,I_i,J_i,B_,K_i,L_i;

   int S,S_;

   //first store everything in ward, then multiply with norms and add to (*this)!
   double ward;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      S = TPM::gblock_char(B,0);

      I_i = TPTPM::gtpmm2t(i,1);
      J_i = TPTPM::gtpmm2t(i,2);

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      for(int j = i;j < TPTPM::gn();++j){

         B_ = TPTPM::gtpmm2t(j,0);

         S_ = TPM::gblock_char(B_,0);

         K_i = TPTPM::gtpmm2t(j,1);
         L_i = TPTPM::gtpmm2t(j,2);

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         ward = 2.0 * dpt2(i,j);

         if(I_i == J_i){

            if(K_i == L_i){

               ward += T2trace; 
               
               ward -= T2barbar[a] + T2barbar[b] + T2barbar[e] + T2barbar[z];

               ward += dpt4(a,e) + dpt4(b,e) + dpt4(a,z) + dpt4(b,z);

            }

            ward += T2bar(S_,e,z,t,h) - dpt3(j,a) - dpt3(j,b);

         }

         if(K_i == L_i)
            ward += T2bar(S,a,b,c,d) - dpt3(i,e) - dpt3(i,z);

         //the norms
         (*this)(i,j) += ward * Gradient::gnorm(i) * Gradient::gnorm(j) * (2.0*S + 1.0) * (2.0*S_ + 1.0);

      }
   }

   //remove the array
   for(int B = 0;B < 2*L;++B)
      delete [] dparray[B];

   delete [] dparray;

}
