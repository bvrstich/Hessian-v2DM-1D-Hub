#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using std::ostream;
using std::ofstream;
using std::vector;
using std::endl;

#include "include.h"

vector< vector<int> > *DPM::dp2s;
int *****DPM::s2dp;

int **DPM::block_char;

/**
 * initialize the static lists and variables
 */
void DPM::init(){

   int L = Tools::gL();

   dp2s = new vector< vector<int> > [2*L];

   s2dp = new int **** [2*L];//nr of blocks

   for(int B = 0;B < L;++B){//first loop over the S = 1/2 blocks

      s2dp[B] = new int *** [2];//for the S = 1/2, we have that S_ab can be 0 or 1

      for(int S_ab = 0;S_ab < 2;++S_ab){//loop and allocate

         s2dp[B][S_ab] = new int ** [L];

         for(int a = 0;a < L;++a){

            s2dp[B][S_ab][a] = new int * [L];

            for(int b = 0;b < L;++b)
               s2dp[B][S_ab][a][b] = new int [L];

         }

      }

   }

   for(int B = L;B < 2*L;++B){//then loop over the S = 3/2 blocks

      s2dp[B] = new int *** [1];//for the S = 3/2, we have that S_ab can be only 1

      s2dp[B][0] = new int ** [L];

      for(int a = 0;a < L;++a){//loop and allocate

         s2dp[B][0][a] = new int * [L];

         for(int b = 0;b < L;++b)
            s2dp[B][0][a][b] = new int [L];

      }

   }

   //allocate the block_char list
   block_char = new int * [2*L];

   for(int B = 0;B < 2*L;++B)
      block_char[B] = new int [2];

   //initialize:
   for(int K = 0;K < L;++K){

      block_char[K][0] = 0;//S (by 0 I mean 1/2)
      block_char[K][1] = K;//K

      block_char[K + L][0] = 1;//S (by 1 I mean 3/2)
      block_char[K + L][1] = K;//K

   }

   //initialize the lists
   int dp;

   vector<int> v(4);

   //first S == 1/2
   for(int K = 0;K < L;++K){

      dp = 0;

      for(int a = 0;a < L;++a){//S_ab == 0 and a == b != c

         for(int c = 0;c < a;++c){

            if( (2*a + c)%L == K){

               v[0] = 0;//S_ab
               v[1] = a;
               v[2] = a;
               v[3] = c;

               dp2s[K].push_back(v);

               s2dp[K][0][a][a][c] = dp;

               ++dp;

            }

         }

         for(int c = a + 1;c < L;++c){

            if( (2*a + c)%L == K){

               v[0] = 0;//S_ab

               v[1] = a;
               v[2] = a;
               v[3] = c;

               dp2s[K].push_back(v);

               s2dp[K][0][a][a][c] = dp;

               ++dp;

            }

         }

      }

      //S=1/2 S_ab=0 a != b != c
      for(int a = 0;a < L;++a)
         for(int b = a + 1;b < L;++b)
            for(int c = b + 1;c < L;++c){

               if( (a + b + c)%L == K ){

                  v[0] = 0;//S_ab

                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  dp2s[K].push_back(v);

                  s2dp[K][0][a][b][c] = dp;

                  ++dp;

               }

            }

      //S == 0, S_ab == 1, a != b != c
      for(int a = 0;a < L;++a)
         for(int b = a + 1;b < L;++b)
            for(int c = b + 1;c < L;++c){

               if( (a + b + c)%L == K ){

                  v[0] = 1;//S_ab

                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  dp2s[K].push_back(v);

                  s2dp[K][1][a][b][c] = dp;

                  ++dp;

               }

            }

   }

   for(int K = 0;K < L;++K){//now for the S=3/2 blocks:

      dp = 0;

      //S == 3/2, S_ab == 1, a != b != c
      for(int a = 0;a < L;++a)
         for(int b = a + 1;b < L;++b)
            for(int c = b + 1;c < L;++c){

               if( (a + b + c)%L == K ){

                  v[0] = 1;//S_ab

                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  dp2s[K + L].push_back(v);

                  s2dp[K + L][0][a][b][c] = dp;//watch out! 0 S_ab index means S_ab == 1 here!

                  ++dp;

               }

            }

   }

}

/**
 * deallocate the static lists
 */
void DPM::clear(){

   //first delete S = 1/2 blocks
   for(int B = 0;B < Tools::gL();++B){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < Tools::gL();++a){

            for(int b = 0;b < Tools::gL();++b)
               delete [] s2dp[B][S_ab][a][b];

            delete [] s2dp[B][S_ab][a];

         }

         delete [] s2dp[B][S_ab];

      }

      delete [] s2dp[B];

   }

   //then the S = 3/2 part
   for(int B = Tools::gL();B < Tools::gM();++B){

      for(int a = 0;a < Tools::gL();++a){

         for(int b = 0;b < Tools::gL();++b)
            delete [] s2dp[B][0][a][b];

         delete [] s2dp[B][0][a];

      }

      delete [] s2dp[B][0];

      delete [] s2dp[B];

   }

   delete [] s2dp;

   delete [] dp2s;

   for(int B = 0;B < Tools::gM();++B)
      delete [] block_char[B];

   delete [] block_char;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 * Tools::gL() blocks, Tools::gL() for S = 1/2 and Tools::gL() 3/2.
 */
DPM::DPM() : BlockMatrix(Tools::gM()) {

   //set the dimension and the degeneracies of the blocks
   for(int K = 0;K < Tools::gL();++K){

      this->setMatrixDim(K,dp2s[K].size(),2);
      this->setMatrixDim(K + Tools::gL(),dp2s[K + Tools::gL()].size(),4);

   }

}

/**
 * copy constructor: constructs BlockMatrix object with 2 * Tools::gL() blocks, on Tools::gL() for S=1/2 and Tools::gL() for S=3/2, and copies the content of the dpm_c blocks into it,
 * @param dpm_c DPM to be copied into (*this)
 */
DPM::DPM(const DPM &dpm_c) : BlockMatrix(dpm_c) { }

/**
 * destructor
 */
DPM::~DPM(){ }

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * DPM(B,S_ab,a,b,c,S_de,d,e,z) = sum_S_ac (some terms dependent on spin) DPM(B,S_ac,a,c,b,S_de,d,e,z) etc...
 * @param S The block index
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param a first sp index that forms the dp row index i of block B together with b, c and S_ab
 * @param b second sp index that forms the dp row index i of block B together with a, c and S_ab
 * @param c third sp index that forms the dp row index i of block B together with a, b and S_ab
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param d first sp index that forms the dp column index j of block B together with e, z and S_de
 * @param e second sp index that forms the dp column index j of block B together with d, z and S_de
 * @param z third sp index that forms the dp column index j of block B together with d, e and S_de
 * @return the number on place DPM(B,i,j) with the right phase and forefactor.
 */
double DPM::operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   int K = (a + b + c)%Tools::gL();

   if( (d + e + z)%Tools::gL() != K )
      return 0.0;

   //blockindex:
   int B = S*Tools::gL() + K;

   int *i = new int [2];
   double *coef_i = new double [2];

   int dim_i = get_inco(B,S_ab,a,b,c,i,coef_i);

   if(dim_i == 0){

      delete [] i;
      delete [] coef_i;

      return 0.0;

   }

   int *j = new int [2];
   double *coef_j = new double [2];

   int dim_j = get_inco(B,S_de,d,e,z,j,coef_j);

   if(dim_j == 0){

      delete [] i;
      delete [] j;

      delete [] coef_i;
      delete [] coef_j;

      return 0.0;

   }

   double ward = 0.0;

   for(int I = 0;I < dim_i;++I)
      for(int J = 0;J < dim_j;++J)
         ward += coef_i[I] * coef_j[J] * (*this)(B,i[I],j[J]);

   delete [] i;
   delete [] j;

   delete [] coef_i;
   delete [] coef_j;

   return ward;

}

/** 
 * Static member function that gets the dp-indices and their coefficients of the (s and t )-p indices S_ab,a,b,c.
 * @param B block index of the state
 * @param S_ab intermediate spincoupling of a and b.
 * @param a first sp-momentum orbital
 * @param b second sp-momentum orbital
 * @param c third sp-momentum orbital
 * @param i pointer of dim 1 or 2 containing the indices occuring in the expansion of this particular dp state in the normal basis (a==b,c a < b < c).
 * @param coef pointer of dim 1 or 2 containing the coefficients occuring in the expansion.
 * @return the number of terms in the expansion (1 or 2), also the dim of pointers i and coef. When zero is returned this is not a valid element.
 */
int DPM::get_inco(int B,int S_ab,int a,int b,int c,int *i,double *coef){

   //they cannot all be equal
   if(a == b && b == c)
      return 0;

   if(B < Tools::gL()){//spin 1/2 block:

      //if normal basis:
      if(a == b){

         if(S_ab == 1)//spin has to be zero for a == b
            return 0;

         i[0] = s2dp[B][0][a][b][c];
         coef[0] = 1;

         return 1;

      }
      else if (a < b && b < c){

         i[0] = s2dp[B][S_ab][a][b][c];
         coef[0] = 1;

         return 1;

      }
      else{//anomal basis:

         int min,max,phase;

         //first order a and b for code saving reasons
         if(a < b){

            min = a;
            max = b;

            phase = 1;

         }
         else{

            min = b;
            max = a;

            phase = 1 - 2*S_ab;

            if(c > max){//we still have one simple dim = 1 term left: b < a < c

               i[0] = s2dp[B][S_ab][b][a][c];
               coef[0] = phase;

               return 1;

            }

         }

         //now we have four possibilities left:
         //don't forget to multiply every result by phase to get the right a and b for min and max!
         // 1) c < min < max
         // 2) c == min < max
         // 3) min < c < max
         // 4) min < max == c
         if(c < min){//c < min < max

            //the S_ca == 0 part:
            i[0] = s2dp[B][0][c][min][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            //the S_ca == 1 part:
            i[1] = s2dp[B][1][c][min][max];
            coef[1] = phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * std::sqrt(3.0) * Tools::g6j(0,0,1,S_ab);

            return 2;

         }
         else if(c == min){//c == min < max: this will also be a 1 dim list, because S_ac can only be 0 if a == c.

            i[0] = s2dp[B][0][c][min][max];
            coef[0] = std::sqrt(2.0) * phase * std::sqrt(2.0*S_ab + 1) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            return 1;

         }
         else if(c < max){//min < c < max

            //S_ac == 0 part:
            i[0] = s2dp[B][0][min][c][max];
            coef[0] = phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * Tools::g6j(0,0,0,S_ab);

            //S_ac == 1 part:
            i[1] = s2dp[B][1][min][c][max];
            coef[1] = - phase * std::sqrt(2.0*S_ab + 1.0) * (1 - 2*S_ab) * std::sqrt(3.0) * Tools::g6j(0,0,1,S_ab);

            return 2;

         }
         else{// min < c == max: also a 1 dim list, S_bc can only be 0 if b == c

            i[0] = s2dp[B][0][max][c][min];
            coef[0] = phase * std::sqrt(2.0) * std::sqrt(2.0*S_ab + 1.0) * Tools::g6j(0,0,0,S_ab);

            return 1;

         }

      }

   }
   else{//spin 3/2 block, totally antisymmetrical in the spatial sp orbs.

      //only S_ab == 1 can couple to 3/2's.
      if(S_ab == 0)
         return 0;

      //if any of the sp orbs are equal, antisymmetry leads to zero:
      if(a == b || b == c || c == a)
         return 0;

      if(a < b){

         if(b < c){//a < b < c

            i[0] = s2dp[B][0][a][b][c];
            coef[0] = 1;

         }
         else if(c < a){//c < a < b

            i[0] = s2dp[B][0][c][a][b];
            coef[0] = 1;

         }
         else{//a < c < b

            i[0] = s2dp[B][0][a][c][b];
            coef[0] = -1;

         }

      }
      else{//b < a

         if(a < c){//b < a < c

            i[0] = s2dp[B][0][b][a][c];
            coef[0] = -1;

         }
         else if(c < b){//c < b < a

            i[0] = s2dp[B][0][c][b][a];
            coef[0] = -1;

         }
         else{//b < c < a

            i[0] = s2dp[B][0][b][c][a];
            coef[0] = 1;

         }

      }

      return 1;

   }

}

/**
 * The spincoupled, translationally invariant T1-like (generalized T1) map: maps a TPM object (tpm) on a DPM object (*this)
 * @param A term before the tp part of the map
 * @param B term before the np part of the map
 * @param C term before the sp part of the map
 * @param tpm input TPM
 */
void DPM::T(double A,double B,double C,const TPM &tpm){

   //make sp matrix out of tpm
   SPM spm;
   spm.bar(C,tpm);

   double ward = 2.0*B*tpm.trace();

   int a,b,c,d,e,z;
   int S_ab,S_de;

   int sign_ab,sign_de;

   double norm_ab,norm_de;

   double hard;

   //start with the S = 1/2 blocks, these are the most difficult:
   for(int B = 0;B < Tools::gL();++B){

      for(int i = 0;i < gdim(B);++i){

         S_ab = dp2s[B][i][0];

         a = dp2s[B][i][1];
         b = dp2s[B][i][2];
         c = dp2s[B][i][3];

         sign_ab = 1 - 2*S_ab;

         norm_ab = 1.0;

         if(a == b)
            norm_ab /= std::sqrt(2.0);

         for(int j = i;j < gdim(B);++j){

            S_de = dp2s[B][j][0];

            d = dp2s[B][j][1];
            e = dp2s[B][j][2];
            z = dp2s[B][j][3];

            sign_de = 1 - 2*S_de;

            norm_de = 1.0;

            if(d == e)
               norm_de /= std::sqrt(2.0);

            hard = std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * Tools::g6j(0,0,S_ab,S_de);

            //init
            (*this)(B,i,j) = 0.0;

            //the np + sp part
            if(i == j)
               (*this)(B,i,j) = ward - spm[a] - spm[b] - spm[c];

            //other parts are a bit more difficult.

            //tp(1)
            if(c == z)
               if(S_ab == S_de)
                  (*this)(B,i,j) += A * tpm(S_ab,a,b,d,e);

            //tp(2)
            if(b == z){

               if(a == c)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);
               else
                  (*this)(B,i,j) += A * norm_ab * sign_ab * sign_de * hard * tpm(S_de,a,c,d,e);

            }

            //tp(3)
            if(a == z){

               if(b == c)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);
               else
                  (*this)(B,i,j) += A * norm_ab * sign_de * hard * tpm(S_de,b,c,d,e);

            }

            //tp(4)
            if(c == e){

               if(d == z)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);
               else
                  (*this)(B,i,j) += A * norm_de * sign_ab * sign_de * hard * tpm(S_ab,a,b,d,z);

            }

            //tp(5)
            if(b == e){

               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,a,c,d,z);

               //correct for norms of the tpm
               if(a == c)
                  hulp *= std::sqrt(2.0);

               if(d == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * norm_ab * norm_de * sign_ab * sign_de * std::sqrt( (2*S_ab + 1.0) * (2*S_de + 1.0) ) * hulp;

            }

            //tp(6)
            if(a == e){

               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,b,c,d,z);

               if(b == c)
                  hulp *= std::sqrt(2.0);

               if(d == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * sign_de * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

            }

            //tp(7)
            if(c == d){

               if(e == z)
                  (*this)(B,i,j) += std::sqrt(2.0) * A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);
               else
                  (*this)(B,i,j) += A * norm_de * sign_ab * hard * tpm(S_ab,a,b,e,z);

            }

            //tp(8)
            if(b == d){

               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,a,c,e,z);

               if(a == c)
                  hulp *= std::sqrt(2.0);

               if(e == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

            }

            //tp(9)
            if(a == d){

               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_de) * tpm(Z,b,c,e,z);

               if(b == c)
                  hulp *= std::sqrt(2.0);

               if(e == z)
                  hulp *= std::sqrt(2.0);

               (*this)(B,i,j) += A * std::sqrt( (2*S_ab + 1) * (2*S_de + 1.0) ) * norm_ab * norm_de * hulp;

            }

         }
      }

   }

   //then the S = 3/2 blocks, this should be easy, totally antisymmetrical 
   for(int B = Tools::gL();B < Tools::gM();++B){

      for(int i = 0;i < gdim(B);++i){

         a = dp2s[B][i][1];
         b = dp2s[B][i][2];
         c = dp2s[B][i][3];

         for(int j = i;j < gdim(B);++j){

            d = dp2s[B][j][1];
            e = dp2s[B][j][2];
            z = dp2s[B][j][3];

            (*this)(B,i,j) = 0.0;

            //np + sp part:
            if(i == j)
               (*this)(B,i,j) = ward - spm[a] - spm[b] - spm[c];

            //tp(1)
            if(c == z)
               (*this)(B,i,j) += A * tpm(1,a,b,d,e);

            //tp(2)
            if(b == z)
               (*this)(B,i,j) -= A * tpm(1,a,c,d,e);

            //tp(4)
            if(c == e)
               (*this)(B,i,j) -= A * tpm(1,a,b,d,z);

            //tp(5)
            if(b == e)
               (*this)(B,i,j) += A * tpm(1,a,c,d,z);

            //tp(7)
            if(c == d)
               (*this)(B,i,j) += A * tpm(1,a,b,e,z);

            //tp(8)
            if(b == d)
               (*this)(B,i,j) -= A * tpm(1,a,c,e,z);

            //tp(9)
            if(a == d)
               (*this)(B,i,j) += A * tpm(1,b,c,e,z);

         }
      }

   }

   this->symmetrize();

}

/**
 * The T1-map: maps a TPM object (tpm) on a DPM object (*this). 
 * @param tpm input TPM
 */
void DPM::T(const TPM &tpm){

   double a = 1.0;
   double b = 1.0/(Tools::gN()*(Tools::gN() - 1.0));
   double c = 1.0/(Tools::gN() - 1.0);

   this->T(a,b,c,tpm);

}

/** 
 * The hat function maps a TPM object tpm to a DPM object (*this) so that bar(this) = tpm,
 * The inverse of the TPM::bar function. It is a T1-like map.
 * @param tpm input TPM
 */
void DPM::hat(const TPM &tpm){

   double a = 1.0/(Tools::gM() - 4.0);
   double b = 1.0/((Tools::gM() - 4.0)*(Tools::gM() - 3.0)*(Tools::gM() - 2.0));
   double c = 1.0/((Tools::gM() - 4.0)*(Tools::gM() - 3.0));

   this->T(a,b,c,tpm);

}

ostream &operator<<(ostream &output,const DPM &dpm_p){

   for(int B = 0;B < dpm_p.gnr();++B){

      output << B << "\t(" << dpm_p.gblock_char(B,0) << "," << dpm_p.gblock_char(B,1) << ")\t" << dpm_p.gdim(B) << "\t" << dpm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < dpm_p.gdim(B);++i)
         for(int j = 0;j < dpm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << 

               dpm_p.dp2s[B][i][0] << "\t" << dpm_p.dp2s[B][i][1] << "\t" << dpm_p.dp2s[B][i][2] << "\t" << dpm_p.dp2s[B][i][3] << 

               "\t" << dpm_p.dp2s[B][j][0] << "\t" << dpm_p.dp2s[B][j][1] << "\t" << dpm_p.dp2s[B][j][2] << "\t" << dpm_p.dp2s[B][j][3] << "\t" << dpm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * access to the lists from outside of the class
 */
int DPM::gs2dp(int B,int S_ab,int a,int b,int c){

   return s2dp[B][S_ab][a][b][c];

}

/**
 * access to the lists from outside of the class
 */
int DPM::gdp2s(int B,int i,int option){

   return dp2s[B][i][option];

}
/**
 * access to the lists from outside of the class
 */
int DPM::gblock_char(int B,int option){

   return block_char[B][option];

}

/**
 * convert a PHM matrix to a double array for faster access to the number
 */
void DPM::convert(double **array) const {

   int K;

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   for(int B = 0;B < L;++B){//first S = 1/2

      K = block_char[B][1];

      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int a = 0;a < L;++a)
            for(int b = 0;b < L;++b)
               for(int S_de = 0;S_de < 2;++S_de)
                  for(int d = 0;d < L;++d)
                     for(int e = 0;e < L;++e)
                        array[B][a + b*L + d*L2 + e*L3 + S_ab*L4 + 2*S_de*L4] = (*this)(0,S_ab,a,b,(K - a - b + 2*L)%L,S_de,d,e,(K - d - e + 2*L)%L);

   }

   for(int B = L;B < 2*L;++B){//then S = 3/2

      K = block_char[B][1];

      for(int a = 0;a < L;++a)
         for(int b = 0;b < L;++b)
            for(int d = 0;d < L;++d)
               for(int e = 0;e < L;++e)
                  array[B][a + b*L + d*L2 + e*L3] = (*this)(1,1,a,b,(K - a - b + 2*L)%L,1,d,e,(K - d - e + 2*L)%L);

   }

}

/**
 * convert a DPM matrix to a double array for faster access to the number, fast conversion
 */
void DPM::convert_fast(double **array) const {

   int K;

   int i,j;
   int i1,i2,j1,j2;
   double coef_i1,coef_i2,coef_j1,coef_j2;

   double SQ_2 = std::sqrt(2.0);
   double SQ_3 = std::sqrt(3.0);

   int sign_ab,sign_de;
   double sqab,sqde;

   int c,z;

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   //first S = 1/2
   for(int K = 0;K < L;++K){

      //S_ab can be 0 or 1
      for(int a = 0;a < L;++a){

         // (1) a == b only when S_ab == 0
         c = (K - 2*a + 2*L)%L;

         if(c == a){//everything zero!

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_de = 0;S_de < 2;++S_de)
                  for(int d = 0;d < L;++d)
                     for(int e = 0;e < L;++e)
                        array[K][a*(1 + L) + d*L2 + e*L3 + S_ab*L4 + 2*S_de*L4] = 0.0;

         }
         else{

            //when S_ab == 1: all zero if a == b
            for(int S_de = 0;S_de < 2;++S_de)
               for(int d = 0;d < L;++d)
                  for(int e = 0;e < L;++e)
                     array[K][a*(1 + L) + d*L2 + e*L3 + L4 + 2*S_de*L4] = 0.0;

            //index for S_ab == 0
            i = s2dp[K][0][a][a][c];

            //only S_ab == 0 terms remain to be filled
            for(int d = 0;d < L;++d){

               //(1) d == e only when S_de == 0
               z = (K - 2*d + 2*L)%L;

               if(z == d){//everything zero

                  for(int S_de = 0;S_de < 2;++S_de)
                     array[K][a*(1 + L) + d*(L2 + L3) + 2*S_de*L4] = 0.0;

               }
               else{

                  //when S_de == 1: zero!
                  array[K][a*(1 + L) + d*(L2 + L3) + 2*L4] = 0.0;

                  //index for S_de == 0
                  j = s2dp[K][0][d][d][z];

                  array[K][a*(1 + L) + d*(L2 + L3)] = (*this)(K,i,j);

               }

               //(2) d < e, both S_de = 0 and 1
               for(int S_de = 0;S_de < 2;++S_de){

                  sqde = std::sqrt(2*S_de + 1.0);
                  sign_de = 1 - 2*S_de;

                  for(int e = d + 1;e < L;++e){

                     z = (K - d - e + 2*L)%L;

                     if(z < d){//z < d < e

                        j1 = s2dp[K][0][z][d][e];
                        j2 = s2dp[K][1][z][d][e];

                        coef_j1 = sqde * sign_de * Tools::g6j(0,0,S_de,0);
                        coef_j2 = sqde * sign_de * SQ_3 * Tools::g6j(0,0,S_de,1);

                        array[K][a*(1 + L) + d*L2 + e*L3 + 2*S_de*L4] = coef_j1 * (*this)(K,i,j1) + coef_j2 * (*this)(K,i,j2);
                        array[K][a*(1 + L) + e*L2 + d*L3 + 2*S_de*L4] = sign_de * array[K][a*(1 + L) + d*L2 + e*L3 + 2*S_de*L4];

                     }
                     else if(z == d){//z == d < e

                        j1 = s2dp[K][0][z][d][e];

                        coef_j1 = SQ_2 * sqde * sign_de * Tools::g6j(0,0,S_de,0);

                        array[K][a*(1 + L) + d*L2 + e*L3 + 2*S_de*L4] = coef_j1 * (*this)(K,i,j1);
                        array[K][a*(1 + L) + e*L2 + d*L3 + 2*S_de*L4] = sign_de * array[K][a*(1 + L) + d*L2 + e*L3 + 2*S_de*L4];

                     }
                     else if(z < e){//d < z < e

                        j1 = s2dp[K][0][d][z][e];
                        j2 = s2dp[K][1][d][z][e];

                        coef_i1 = sqde * sign_de * Tools::g6j(0,0,S_de,0);
                        coef_i2 = - sqde * sign_de * SQ_3 * Tools::g6j(0,0,S_de,1);

                        array[K][a*(1 + L) + d*L2 + e*L3 + 2*S_de*L4] = coef_j1 * (*this)(K,i,j1) + coef_j2 * (*this)(K,i,j2);
                        array[K][a*(1 + L) + e*L2 + d*L3 + 2*S_de*L4] = sign_de * array[K][a*(1 + L) + d*L2 + e*L3 + 2*S_de*L4];

                     }
                     else if(z == e){//d < z == e

                        j1 = s2dp[K][0][z][e][d];

                        coef_j1 = SQ_2 * sqde * Tools::g6j(0,0,S_de,0);

                        array[K][a*(1 + L) + d*L2 + e*L3 + 2*S_de*L4] = coef_j1 * (*this)(K,i,j1);
                        array[K][a*(1 + L) + e*L2 + d*L3 + 2*S_de*L4] = sign_de * array[K][a*(1 + L) + d*L2 + e*L3 + 2*S_de*L4];

                     }
                     else{//d < e < z

                        j = s2dp[K][S_de][d][e][z];

                        array[K][a*(1 + L) + d*L2 + e*L3 + 2*S_de*L4] = (*this)(K,i,j);
                        array[K][a*(1 + L) + e*L2 + d*L3 + 2*S_de*L4] = sign_de * (*this)(K,i,j);

                     }

                  }

               }

            }

         }//end else of if a == b == c

         //(2) a < b, both S_ab = 0 and 1
         for(int S_ab = 0;S_ab < 2;++S_ab){

            sqab = std::sqrt(2*S_ab + 1.0);
            sign_ab = 1 - 2*S_ab;

            for(int b = a + 1;b < L;++b){

               c = (K - a - b + 2*L)%L;

               if(c < a){//c < a < b

                  i1 = s2dp[K][0][c][a][b];
                  i2 = s2dp[K][1][c][a][b];

                  coef_i1 = sqab * sign_ab * Tools::g6j(0,0,S_ab,0);
                  coef_i2 = sqab * sign_ab * SQ_3 * Tools::g6j(0,0,S_ab,1);

               }
               else if(c == a){//a a < b

                  i1 = s2dp[K][0][c][a][b];

                  coef_i1 = SQ_2 * sqab * sign_ab * Tools::g6j(0,0,S_ab,0);

               }
               else if(c < b){//a < c < b

                  i1 = s2dp[K][0][a][c][b];
                  i2 = s2dp[K][1][a][c][b];

                  coef_i1 = sqab * sign_ab * Tools::g6j(0,0,S_ab,0);
                  coef_i2 = - sqab * sign_ab * SQ_3 * Tools::g6j(0,0,S_ab,1);

               }
               else if(c == b){//a < b  b

                  i1 = s2dp[K][0][c][b][a];

                  coef_i1 = SQ_2 * sqab * Tools::g6j(0,0,S_ab,0);


               }
               else{//a < b < c

                  i = s2dp[K][S_ab][a][b][c];

               }

            }

         }

      }

   }

   //then S = 3/2
   for(int B = L;B < 2*L;++B){

      K = block_char[B][1];

      for(int a = 0;a < L;++a){

         // (1) first a == b
         for(int d = 0;d < L;++d)
            for(int e = 0;e < L;++e)
               array[B][a*(L + 1) + d*L2 + e*L3] = 0.0;

         // (2) a < b
         for(int b = a + 1;b < L;++b){

            c = (K - a - b + 2*L)%L;

            if(c < a){

               i = s2dp[B][0][c][a][b];

               for(int d = 0;d < L;++d){

                  //(1) first d == e
                  array[B][a + b*L + d*(L2 + L3)] = array[B][b + a*L + d*(L2 + L3)] = 0.0;

                  //(2) d < e
                  for(int e = d + 1;e < L;++e){

                     z = (K - d - e + 2*L)%L;

                     if(z < d){//z < d < e

                        j = s2dp[B][0][z][d][e];

                        array[B][a + b*L + d*L2 + e*L3] = (*this)(B,i,j);
                        array[B][b + a*L + d*L2 + e*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][a + b*L + e*L2 + d*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][b + a*L + e*L2 + d*L3] = array[B][a + b*L + d*L2 + e*L3];

                     }
                     else if(z == d){

                        array[B][a + b*L + d*L2 + e*L3] = 0.0;
                        array[B][b + a*L + d*L2 + e*L3] = 0.0;
                        array[B][a + b*L + e*L2 + d*L3] = 0.0;
                        array[B][b + a*L + e*L2 + d*L3] = 0.0;

                     }
                     else if(z < e){//d < z < e

                        j = s2dp[B][0][d][z][e];

                        array[B][a + b*L + d*L2 + e*L3] = -(*this)(B,i,j);
                        array[B][b + a*L + d*L2 + e*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][a + b*L + e*L2 + d*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][b + a*L + e*L2 + d*L3] = array[B][a + b*L + d*L2 + e*L3];

                     }
                     else if(z == e){

                        array[B][a + b*L + d*L2 + e*L3] = 0.0;
                        array[B][b + a*L + d*L2 + e*L3] = 0.0;
                        array[B][a + b*L + e*L2 + d*L3] = 0.0;
                        array[B][b + a*L + e*L2 + d*L3] = 0.0;

                     }
                     else{//d < e < z

                        j = s2dp[B][0][d][e][z];

                        array[B][a + b*L + d*L2 + e*L3] = (*this)(B,i,j);
                        array[B][b + a*L + d*L2 + e*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][a + b*L + e*L2 + d*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][b + a*L + e*L2 + d*L3] = array[B][a + b*L + d*L2 + e*L3];

                     }

                  }

               }

            }
            else if(c == a){

               for(int d = 0;d < L;++d)
                  for(int e = 0;e < L;++e)
                     array[B][a + b*L + d*L2 + e*L3] = array[B][b + a*L + d*L2 + e*L3] = 0.0; 

            }
            else if(c < b){

               //add minus!
               i = s2dp[B][0][a][c][b];

               for(int d = 0;d < L;++d){

                  //(1) first d == e
                  array[B][a + b*L + d*(L2 + L3)] = array[B][b + a*L + d*(L2 + L3)] = 0.0;

                  //(2) d < e
                  for(int e = d + 1;e < L;++e){

                     z = (K - d - e + 2*L)%L;

                     if(z < d){//z < d < e

                        j = s2dp[B][0][z][d][e];

                        array[B][a + b*L + d*L2 + e*L3] = -(*this)(B,i,j);
                        array[B][b + a*L + d*L2 + e*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][a + b*L + e*L2 + d*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][b + a*L + e*L2 + d*L3] = array[B][a + b*L + d*L2 + e*L3];

                     }
                     else if(z == d){

                        array[B][a + b*L + d*L2 + e*L3] = 0.0;
                        array[B][b + a*L + d*L2 + e*L3] = 0.0;
                        array[B][a + b*L + e*L2 + d*L3] = 0.0;
                        array[B][b + a*L + e*L2 + d*L3] = 0.0;

                     }
                     else if(z < e){//d < z < e

                        j = s2dp[B][0][d][z][e];

                        array[B][a + b*L + d*L2 + e*L3] = (*this)(B,i,j);
                        array[B][b + a*L + d*L2 + e*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][a + b*L + e*L2 + d*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][b + a*L + e*L2 + d*L3] = array[B][a + b*L + d*L2 + e*L3];

                     }
                     else if(z == e){

                        array[B][a + b*L + d*L2 + e*L3] = 0.0;
                        array[B][b + a*L + d*L2 + e*L3] = 0.0;
                        array[B][a + b*L + e*L2 + d*L3] = 0.0;
                        array[B][b + a*L + e*L2 + d*L3] = 0.0;

                     }
                     else{//d < e < z

                        j = s2dp[B][0][d][e][z];

                        array[B][a + b*L + d*L2 + e*L3] = -(*this)(B,i,j);
                        array[B][b + a*L + d*L2 + e*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][a + b*L + e*L2 + d*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][b + a*L + e*L2 + d*L3] = array[B][a + b*L + d*L2 + e*L3];

                     }

                  }

               }

            }
            else if(c == b){

               for(int d = 0;d < L;++d)
                  for(int e = 0;e < L;++e)
                     array[B][a + b*L + d*L2 + e*L3] = array[B][b + a*L + d*L2 + e*L3] = 0.0;

            }
            else{//a < b < c

               i = s2dp[B][0][a][b][c];

               for(int d = 0;d < L;++d){

                  //(1) first d == e
                  array[B][a + b*L + d*(L2 + L3)] = array[B][b + a*L + d*(L2 + L3)] = 0.0;

                  //(2) d < e
                  for(int e = d + 1;e < L;++e){

                     z = (K - d - e + 2*L)%L;

                     if(z < d){//z < d < e

                        j = s2dp[B][0][z][d][e];

                        array[B][a + b*L + d*L2 + e*L3] = (*this)(B,i,j);
                        array[B][b + a*L + d*L2 + e*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][a + b*L + e*L2 + d*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][b + a*L + e*L2 + d*L3] = array[B][a + b*L + d*L2 + e*L3];

                     }
                     else if(z == d){

                        array[B][a + b*L + d*L2 + e*L3] = 0.0;
                        array[B][b + a*L + d*L2 + e*L3] = 0.0;
                        array[B][a + b*L + e*L2 + d*L3] = 0.0;
                        array[B][b + a*L + e*L2 + d*L3] = 0.0;

                     }
                     else if(z < e){//d < z < e

                        j = s2dp[B][0][d][z][e];

                        array[B][a + b*L + d*L2 + e*L3] = -(*this)(B,i,j);
                        array[B][b + a*L + d*L2 + e*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][a + b*L + e*L2 + d*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][b + a*L + e*L2 + d*L3] = array[B][a + b*L + d*L2 + e*L3];

                     }
                     else if(z == e){

                        array[B][a + b*L + d*L2 + e*L3] = 0.0;
                        array[B][b + a*L + d*L2 + e*L3] = 0.0;
                        array[B][a + b*L + e*L2 + d*L3] = 0.0;
                        array[B][b + a*L + e*L2 + d*L3] = 0.0;

                     }
                     else{//d < e < z

                        j = s2dp[B][0][d][e][z];

                        array[B][a + b*L + d*L2 + e*L3] = (*this)(B,i,j);
                        array[B][b + a*L + d*L2 + e*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][a + b*L + e*L2 + d*L3] = -array[B][a + b*L + d*L2 + e*L3];
                        array[B][b + a*L + e*L2 + d*L3] = array[B][a + b*L + d*L2 + e*L3];

                     }

                  }

               }

            }

         }

      }

   }

}

/* vim: set ts=3 sw=3 expandtab :*/
