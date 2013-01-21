#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 11-05-2010\n\n
 * This class SPM was written for single particle matrices in a spinsymmetrical and translationally invariant system. 
 */
class SPM {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << spm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << spm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de SPM you want to print
    */
   friend ostream &operator<<(ostream &output,const SPM &spm_p);

   public:
      
      //constructor
      SPM();

      //copy constructor
      SPM(const SPM &);

      //destructor
      virtual ~SPM();

      double operator[](int) const;

      double &operator[](int);

      void bar(double,const TPM &);

   private:

      double *spm;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
