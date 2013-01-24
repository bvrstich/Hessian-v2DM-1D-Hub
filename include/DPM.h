#ifndef DPM_H
#define DPM_H

#include <iostream>
#include <vector>

using std::ostream;
using std::vector;

#include "BlockMatrix.h"
#include "TPM.h"

/**
 * @author Brecht Verstichel
 * @date 31-05-2010\n\n
 * This class, DPM, is a class written for spinsymmetrical, translationally invariant three-particle matrices 
 * (name comes from drie-particle matrix). It is written specially for the T_1 condition. 
 * It inherits all the functions from its mother class BlockMatrix, some special member functions and two lists that
 * give the relationship between the dp (three-particle) and the sp basis. This matrix falls apart in M blocks: S = 1/2 with degeneracy 2 and S = 3/2 with degeneracy 4,
 * both have a blockstructure of (k_a + k_b + k_b) % M == K. The basis is determined by the dp-momentum K, the sp-momenuntum qn's k_a,k_b,k_c,
 * an intermediate spincoupling quantumnumber S_ab = 0 or 1, and the total spin S.
 */
class DPM : public BlockMatrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << dpm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << dpm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param dpm_p the DPM you want to print
    */
   friend ostream &operator<<(ostream &output,const DPM &dpm_p);

   public:
      
      //constructor
      DPM();

      //copy constructor
      DPM(const DPM &);

      //destructor
      virtual ~DPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //easy to access the numbers, in sp mode
      double operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int f) const;

      //generalized T1 map
      void T(double,double,double,const TPM &);

      //maak een DPM van een TPM via de T1 conditie
      void T(const TPM &);

      //maak een DPM van een TPM via de hat functie
      void hat(const TPM &);

      static int get_inco(int S,int S_ab,int a,int b,int c,int *i,double *coef);

      static int gblock_char(int,int);
      
      static int gdp2s(int,int,int);

      static int gs2dp(int,int,int,int,int);
      
      static void init();

      static void clear();

   private:

      //!list that takes in a dp index i for block B and returns an intermediate spin: S_ab and three sp-momentum indices: a, b and c = dp2s[B][i][2]
      static vector< vector<int> > *dp2s;

      //!list that takes a block index B, an intermediate spin-index S_ab and three sp-momentum indices a, b and c "and" returns a dp index i for block B: i
      static int *****s2dp;

      //!static list that takes a blockindex B and returns the dp spin S and the dp momentum K. S = block_char[B][0] , K = block_char[B][1]
      static int **block_char;

};

#endif

/* vim: set ts=3 sw=3 expandtab :*/
