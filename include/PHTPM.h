#ifndef PHTPM_H
#define PHTPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "RecMat.h"

class PHM;
class DPM;

/**
 * @author Brecht Verstichel
 * @date 14-01-2013\n\n
 * This class PHTPM is a class written for the singly-traced TPTPM matrix object,
 * being a rectangular matrix on TP for the rows, and SP for the columns, it inherits alle the function from its mother 
 * RecMat, a rectangular matrix, some special member functions and two lists 
 */
class PHTPM : public RecMat {

   public:
      
      //constructor
      PHTPM();

      //copy constructor
      PHTPM(const PHTPM &);

      //destructor
      virtual ~PHTPM();

      using RecMat::operator=;

      using RecMat::operator();

      void dptw(double **);

   private:

};

#endif
