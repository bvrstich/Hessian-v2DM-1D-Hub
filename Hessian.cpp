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
