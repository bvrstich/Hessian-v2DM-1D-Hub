#ifndef PHPHM_H
#define PHPHM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class PHM;
class DPM;

/**
 * @author Brecht Verstichel
 * @date 23-11-2012\n\n
 * This class PHPHM is a class written for matrices of two particle matrices, it inherits alle the function from its mother 
 * Matrix, some special member functions and two lists that give the relationship between the sp and the tp 
 * basis.
 */

class PHPHM : public Matrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param phmm_p the PHPHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PHPHM &phmm_p);

   public:
      
      //constructor
      PHPHM();

      //copy constructor
      PHPHM(const PHPHM &);

      //destructor
      virtual ~PHPHM();

      using Matrix::operator=;

      using Matrix::operator();

      static int gn();

      static int gphmm2ph(int,int);

      static int gph2phmm(int,int,int);

      static void init();

      static void clear();

   private:

      //!list relating the particle-hole space to the PHPHM basis
      static vector< vector<int> > phmm2ph;

      //!list relating the particle-hole space to the PHPHM basis
      static int ***ph2phmm;

};

#endif
