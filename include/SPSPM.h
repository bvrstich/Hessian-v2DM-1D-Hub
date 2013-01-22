#ifndef SPSPM_H
#define SPSPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class PHM;

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class SPSPM is a class written for matrices of two particle matrices, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */

class SPSPM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param spmm_p the SPSPM you want to print
    */
   friend ostream &operator<<(ostream &output,const SPSPM &spmm_p);

   public:
      
      //constructor
      SPSPM();

      //copy constructor
      SPSPM(const SPSPM &);

      //destructor
      virtual ~SPSPM();

      using Matrix::operator=;

      using Matrix::operator();

      void dpt2(double,const TPM &);

      void dpt2(double,const PHM &);

   private:

};

#endif
