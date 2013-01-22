#ifndef HESSIAN_H
#define HESSIAN_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class Hessian is a class written for matrices of two particle matrices, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */

class Hessian : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param hess_p the Hessian you want to print
    */
   friend ostream &operator<<(ostream &output,const Hessian &hess_p);

   public:
      
      //constructor
      Hessian();

      //copy constructor
      Hessian(const Hessian &);

      //destructor
      virtual ~Hessian();

      using Matrix::operator=;

      using Matrix::operator();

      void I(const TPM &);

      void lagr();

      void Q(const TPM &);

      static int gn();
      
      static void init();

      static void clear();

   private:

};

#endif
