#ifndef TPM_H
#define TPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::ifstream;
using std::vector;

#include "BlockMatrix.h"

class SUP;
class Gradient;

/**
 * @author Brecht Verstichel
 * @date 10-05-2010\n\n
 * This class TPM is a class written for two particle matrices with spinsymmetry and translational symemtry included, it inherits alle the function from its mother 
 * BlockMatrix, some special member functions and two lists that give the relationship between the sp and the tp basis.
 */
class TPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << tpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << tpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param tpm_p the TPM you want to print
    */
   friend ostream &operator<<(ostream &output,const TPM &tpm_p);

   public:

      //constructor
      TPM();

      //copy constructor
      TPM(const TPM &);

      //destructor
      virtual ~TPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode and blockindex
      double operator()(int S,int a,int b,int c,int d) const;

      void hubbard(double U);

      void unit();

      double S_2() const;

      double line_search(double t,SUP &P,const TPM &ham);

      double line_search(double t,const TPM &rdm,const TPM &ham);

      void convert(const Gradient &);

      static int gdim(int);

      static int gt2s(int,int,int);

      static int gs2t(int,int,int);

      static int gblock_char(int,int);

      static double gnorm(int,int);

      static void init();

      static void clear();
      
   private:

      //!static list of dimension [nr][dim[i]][2] that takes in a tp index i and a blockindex B, and returns two momentum indices: k_a = t2s[B][i][0] and k_b = t2s[B][i][1]
      static vector< vector<int> > *t2s;

      //!static list of dimension [nr][L][L] that takes two sp momentum indices k_a,k_b and a blockindex B, and returns a tp index i: i = s2t[B][k_a][k_b]
      static int ***s2t;

      //!static list that takes a blockindex B and returns the tp spin S and the tp momentum K. S = block_char[B][0] , K = block_char[B][1]
      static int **block_char;

      //!list containing the norms arising because of the symmetry between the sp's in the S = 0 block
      static double **norm;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
