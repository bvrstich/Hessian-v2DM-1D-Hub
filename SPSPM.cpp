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

            (*this)(a,e) += Q(0,a,k,e,l) * Q(0,a,k,e,l) / ( TPM::gnorm(a,k) * TPM::gnorm(a,k) * TPM::gnorm(e,l) * TPM::gnorm(e,l) );

         }

         double ward = 0.0;

         //then S = 1
         for(int k = 0;k < Tools::gL();++k){

            int l = (a + k - e + Tools::gL())%Tools::gL();

            (*this)(a,e) += Q(1,a,k,e,l) * Q(1,a,k,e,l);

         }

         (*this)(a,e) += 3.0 * ward;

         (*this)(a,e) *= scale;


      }

}
