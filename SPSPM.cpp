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
SPSPM::SPSPM() : Matrix(Tools::gL()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpmm_c
 * @param spmm_c object that will be copied into this.
 */
SPSPM::SPSPM(const SPSPM &spmm_c) : Matrix(spmm_c){ }

/**
 * destructor
 */
SPSPM::~SPSPM(){ }

ostream &operator<<(ostream &output,const SPSPM &spmm_p){

   for(int a = 0;a < Tools::gL();++a)
      for(int e = a;e < Tools::gL();++e)
         output << a << "\t" << e << spmm_p(a,e) << endl;

   return output;

}

/**
 * construct a SPSPM by doubly-tracing out the direct product of two TPM's
 */
void SPSPM::dpt2(double scale,const TPM &Q){

   for(int a = 0;a < Tools::gL();++a)
      for(int e = a;e < Tools::gL();++e){

         (*this)(a,e) = 0.0;

         //first S = 0
         for(int k = 0;k < Tools::gL();++k){

            int l = (a + k - e + Tools::gL())%Tools::gL();

            (*this)(a,e) += (Q(0,a,k,e,l) * Q(0,a,k,e,l) )/ ( TPM::gnorm(a,k) * TPM::gnorm(a,k) * TPM::gnorm(e,l) * TPM::gnorm(e,l) );

         }

         double ward = 0.0;

         //then S = 1
         for(int k = 0;k < Tools::gL();++k){

            int l = (a + k - e + Tools::gL())%Tools::gL();

            ward += Q(1,a,k,e,l) * Q(1,a,k,e,l);

         }

         (*this)(a,e) += 3.0 * ward;

         (*this)(a,e) *= scale;


      }

   this->symmetrize();

}

/**
 * construct the doubly-traced direct product of two PHM matrices
 */
void SPSPM::dpt2(double scale,const PHM &phm){

   for(int a = 0;a < Tools::gL();++a)
      for(int e = a;e < Tools::gL();++e){

         (*this)(a,e) = 0.0;

         //first S = 0
         for(int k = 0;k < Tools::gL();++k){

            int l = (a + k - e + Tools::gL())%Tools::gL();

            (*this)(a,e) += phm(0,a,k,e,l) * phm(0,a,k,e,l);

         }

         double ward = 0.0;

         //then S = 1
         for(int k = 0;k < Tools::gL();++k){

            int l = (a + k - e + Tools::gL())%Tools::gL();

            ward += phm(1,a,k,e,l) * phm(1,a,k,e,l);

         }

         (*this)(a,e) += 3.0 * ward;

         (*this)(a,e) *= scale;


      }

   this->symmetrize();

}

/**
 * construct a SPSPM by quadruple tracing the direct product of two DPM's
 */
void SPSPM::dpt4(double scale,const DPM &dpm){

   int L = Tools::gL();

   for(int a = 0;a < L;++a)
      for(int e = a;e < L;++e){

         (*this)(a,e) = 0.0;

         double ward = 0.0;

         //first S = 1/2
         for(int l = 0;l < L;++l)
            for(int k = 0;k < L;++k){

               int K = (a + l + k)%L;

               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_en = 0;S_en < 2;++S_en)
                     for(int n = 0;n < L;++n){

                        int p = (K - e - n + 2*L)%L;

                        ward += dpm(0,S_al,a,l,k,S_en,e,n,p) * dpm(0,S_al,a,l,k,S_en,e,n,p) / ( TPM::gnorm(a,l) * TPM::gnorm(a,l) * TPM::gnorm(e,n) * TPM::gnorm(e,n) );

                     }

            }

         (*this)(a,e) = 2.0 * ward;

         ward = 0.0;

         //then S = 3/2
         for(int l = 0;l < L;++l)
            for(int k = 0;k < L;++k){

               int K = (a + l + k)%L;

               for(int n = 0;n < L;++n){

                  int p = (K - e - n + 2*L)%L;

                  ward += dpm(1,1,a,l,k,1,e,n,p) * dpm(1,1,a,l,k,1,e,n,p);

               }

            }

         (*this)(a,e) += 4.0 * ward;

         //scale
         (*this)(a,e) *= 0.5 * scale;

      }

   this->symmetrize();

}

/**
 * construct a SPSPM by quadruple tracing the direct product of two DPM's, input the array
 */
void SPSPM::dpt4(double scale,double **dparray){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   for(int a = 0;a < L;++a)
      for(int e = a;e < L;++e){

         (*this)(a,e) = 0.0;

         double ward = 0.0;

         //first S = 1/2
         for(int l = 0;l < L;++l)
            for(int k = 0;k < L;++k){

               int K = (a + l + k)%L;

               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_en = 0;S_en < 2;++S_en)
                     for(int n = 0;n < L;++n){

                        ward += dparray[K][a + l*L + e*L2 + n*L3 + S_al*L4 + 2*S_en*L4] * dparray[K][a + l*L + e*L2 + n*L3 + S_al*L4 + 2*S_en*L4]

                           / ( TPM::gnorm(a,l) * TPM::gnorm(a,l) * TPM::gnorm(e,n) * TPM::gnorm(e,n) );


                     }

            }

         (*this)(a,e) = 2.0 * ward;

         ward = 0.0;

         //then S = 3/2
         for(int l = 0;l < L;++l)
            for(int k = 0;k < L;++k){

               int K = (a + l + k)%L;

               for(int n = 0;n < L;++n)
                  ward += dparray[K + L][a + l*L + e*L2 + n*L3] * dparray[K + L][a + l*L + e*L2 + n*L3];

            }

         (*this)(a,e) += 4.0 * ward;

         //scale
         (*this)(a,e) *= 0.5 * scale;

      }

   this->symmetrize();

}
