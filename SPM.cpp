#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * constructor, makes vector of dimension L, the SPM is completely diagonal in momentum space.
 */
SPM::SPM() {

   spm = new double [Tools::gL()];

}

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) {

   spm = new double [Tools::gL()];

   for(int k = 0;k < Tools::gL();++k)
      spm[k] = spm_copy[k];

}

/**
 * destructor
 */
SPM::~SPM(){

   delete [] spm;

}

/**
 * read access to the number
 */
double SPM::operator[](int k) const{

   return spm[k];

}

/**
 * read and write access to the number
 */
double &SPM::operator[](int k){

   return spm[k];

}

ostream &operator<<(ostream &output,const SPM &spm_p){

   for(int k = 0;k < Tools::gL();++k)
      output << k << "\t" << spm_p[k] << endl;

   return output;

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,const TPM &tpm){

   for(int k = 0;k < Tools::gL();++k){

      spm[k] = 0.0;

      //first S = 0

      //p < k
      for(int p = 0;p < k;++p)
         spm[k] += tpm(0,k,p,k,p);

      //p == k: factor 2 for norm basis
      spm[k] += 2.0 * tpm(0,k,k,k,k);

      for(int p = k + 1;p < Tools::gL();++p)
         spm[k] += tpm(0,k,p,k,p);

      //then S = 1

      //p < k
      for(int p = 0;p < k;++p)
         spm[k] += 3.0 * tpm(1,k,p,k,p);

      for(int p = k + 1;p < Tools::gL();++p)
         spm[k] += 3.0 * tpm(1,k,p,k,p);

      //spin norm and scale!
      spm[k] *= 0.5 * scale;

   }

}

/**
 * Trace out a set of indices to create the "bar" matrix of a PHM, slight difference from the bar(TPM) function (normalization of the tp basisset).
 * @param scale the factor u want the SPM to be scaled with
 * @param phm the PHM out of which the SPM will be filled
 */
void SPM::bar(double scale,const PHM &phm){

   for(int k = 0;k < Tools::gL();++k){

      (*this)[k] = 0.0;

      //S = 0
      for(int p = 0;p < Tools::gL();++p)
         (*this)[k] += phm(0,k,p,k,p);

      //S = 1
      for(int p = 0;p < Tools::gL();++p)
         (*this)[k] += 3.0 * phm(1,k,p,k,p);

      (*this)[k] *= 0.5 * scale;

   }

}

/** 
 * This bar function maps a PPHM object directly onto a SPM object, scaling it with a factor scale
 * @param scale the scalefactor
 * @param pphm Input PPHM object
 */
void SPM::bar(double scale,const PPHM &pphm){

   for(int k = 0;k < Tools::gL();++k){

      (*this)[k] = 0.0;

      //first S = 1/2 part
      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < Tools::gL();++a){

            for(int b = 0;b < a;++b)//b < a
               (*this)[k] += pphm(0,S_ab,a,b,k,S_ab,a,b,k);

            //a == b norm correction
            (*this)[k] += 2.0 * pphm(0,S_ab,a,a,k,S_ab,a,a,k);

            for(int b = a + 1;b < Tools::gL();++b)//b > a
               (*this)[k] += pphm(0,S_ab,a,b,k,S_ab,a,b,k);

         }
      }

      //then S = 3/2 part:
      for(int a = 0;a < Tools::gL();++a)
         for(int b = 0;b < Tools::gL();++b)
            (*this)[k] += 2.0 * pphm(1,1,a,b,k,1,a,b,k);

      //scaling
      (*this)[k] *= scale;

   }

}

/* vim: set ts=3 sw=3 expandtab :*/
