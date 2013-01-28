#ifndef PPHM_H
#define PPHM_H

#include <iostream>
#include <vector>

using std::ostream;
using std::vector;

#include "BlockMatrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 03-06-2010\n\n
 * This class, PPHM, is a class written for spinsymmetrical, translationally invaraint two-particle-one-hole matrices. It is written specially for the T_2 condition. 
 * It inherits all the functions from its mother class BlockMatrix, some special member functions and two lists that give
 * the relationship between the pph (two-particle one hole) and the sp basis. This matrix has M blocks, M/2 for S = 1/2 block with degeneracy 2
 * and M/2 for S = 3/2 block with degeneracy 4.
 */
class PPHM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << pphm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << pphm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param pphm_p the PPHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PPHM &pphm_p);

   public:
      
      //constructor
      PPHM();

      //copy constructor
      PPHM(const PPHM &);

      //destructor
      virtual ~PPHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int S,int S_ab,int k_a,int k_b,int k_c,int S_de,int k_d,int k_e,int k_z) const;

      int get_inco(int B,int S_ab,int k_a,int k_b,int k_c,int &i) const;

      //maak een PPHM van een TPM via de T2 conditie
      void T(const TPM &);

      void convert(double **) const;

      static int gblock_char(int,int);

      static int gpph2s(int,int,int);

      static int gs2pph(int,int,int,int,int);

      static void init();

      static void clear();
      
   private:

      //!static list of dimension [M][dim[B]][4] that takes in a pph index i and a blockindex B for spin and momentum, and returns three momentum sp indices: k_a = pph2s[B][i][1], k_b = pph2s[B][i][2] and k_c = pph2s[B][i][3] and an intermediate spin S_ab = pph2s[B][i][0]
      static vector< vector<int> > *pph2s;

      //!static list of dimension [M][2][M/2][M/2][M/2] that takes three momentum sp indices k_a, k_b and k_c, a blockindex B for total spin and momentum, and an intermediate spinindex S_ab, and returns a pph index i: i = s2pph[B][S_ab][k_a][k_b][k_c]
      static int *****s2pph;

      //!list of block characteristics: block_char[B][0] = S, block_char[B][1] = K
      static int **block_char;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
