#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>

using std::ostream;
using std::ofstream;
using std::cout;
using std::endl;
using std::vector;

#include "include.h"

vector< vector<int> > *PPHM::pph2s;
int *****PPHM::s2pph;

int **PPHM::block_char;

/**
 * initialize the static lists and variables
 */
void PPHM::init(){

   int L = Tools::gL();

  //allocate block_char
   block_char = new int * [2*L];

   for(int B = 0;B < 2*L;++B)
      block_char[B] = new int [2];

   //initiate
   for(int K = 0;K < L;++K){

      block_char[K][0] = 0;//S = 1/2
      block_char[K + L][0] = 1;//S = 3/2

      block_char[K][1] = K;
      block_char[K + L][1] = K;

   }

   //first allocation
   pph2s = new vector< vector<int> > [2*L];//nr of blocks

   s2pph = new int **** [2*L];//nr of blocks

   for(int B = 0;B < L;++B){//loop over the S=1/2 blocks:

      s2pph[B] = new int *** [2];//for the S = 1/2 blocks, we have that S_ab can be 0 or 1

      for(int S_ab = 0;S_ab < 2;++S_ab){//loop and allocate

         s2pph[B][S_ab] = new int ** [L];

         for(int a = 0;a < L;++a){

            s2pph[B][S_ab][a] = new int * [L];

            for(int b = 0;b < L;++b)
               s2pph[B][S_ab][a][b] = new int [L];

         }

      }

   }

   for(int B = L;B < 2*L;++B){//loop over the S=3/2 blocks:

      s2pph[B] = new int *** [1];//for the S = 3/2, we have that S_ab can be only 1

      s2pph[B][0] = new int ** [L];

      for(int a = 0;a < L;++a){//loop and allocate

         s2pph[B][0][a] = new int * [L];

         for(int b = 0;b < L;++b)
            s2pph[B][0][a][b] = new int [L];

      }

   }

   //initialize the lists
   int pph;

   vector<int> v(4);

   for(int K = 0;K < L;++K){//first loop over the S = 1/2 blocks

      //re-init teller for blocks
      pph = 0;

      //S = 1/2 S_ab = 0: a <= b, c
      for(int a = 0;a < L;++a)
         for(int b = a;b < L;++b)
            for(int c = 0;c < L;++c){

               if( (a + b + c)%L == K ){

                  v[0] = 0;//S_ab

                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  pph2s[K].push_back(v);

                  s2pph[K][0][a][b][c] = pph;

                  ++pph;

               }

            }

      //S = 1/2, S_ab = 1, a < b ,c
      for(int a = 0;a < L;++a)
         for(int b = a + 1;b < L;++b)
            for(int c = 0;c < L;++c){

               if( (a + b + c)%L == K ){

                  v[0] = 1;//S_ab

                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  pph2s[K].push_back(v);

                  s2pph[K][1][a][b][c] = pph;

                  ++pph;

               }

            }

   }

   for(int B = L;B < 2*L;++B){//then loop over the S = 3/2 blocks

      //re-init teller for blocks
      pph = 0;

      //S = 3/2, S_ab = 1, a < b, c
      for(int a = 0;a < L;++a)
         for(int b = a + 1;b < L;++b)
            for(int c = 0;c < L;++c){

               if( (a + b + c)%L == block_char[B][1] ){

                  v[0] = 1;//S_ab: watch it, onlyh one S_ab here!

                  v[1] = a;
                  v[2] = b;
                  v[3] = c;

                  pph2s[B].push_back(v);

                  s2pph[B][0][a][b][c] = pph;

                  ++pph;

               }

            }

   }

}
/**
 * deallocate the static lists
 */
void PPHM::clear(){

   int L = Tools::gL();

   //delete block_char:
   for(int B = 0;B < 2*L;++B)
      delete [] block_char[B];

   delete [] block_char;

   //first delete S = 1/2 parts
   for(int B = 0;B < L;++B){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < L;++a){

            for(int b = 0;b < L;++b)
               delete [] s2pph[B][S_ab][a][b];

            delete [] s2pph[B][S_ab][a];

         }

         delete [] s2pph[B][S_ab];

      }

      delete [] s2pph[B];

   }

   //then the S = 3/2 part
   for(int B = L;B < 2*L;++B){

      for(int a = 0;a < L;++a){

         for(int b = 0;b < L;++b)
            delete [] s2pph[B][0][a][b];

         delete [] s2pph[B][0][a];

      }

      delete [] s2pph[B][0];
      delete [] s2pph[B];

   }

   delete [] s2pph;

   delete [] pph2s;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 1/2 and 3/2.
 */
PPHM::PPHM() : BlockMatrix(Tools::gM()) {

   int L = Tools::gL();

   for(int K = 0;K < L;++K){

      //set the dimension and the degeneracies of the blocks
      this->setMatrixDim(K,pph2s[K].size(),2);//S=1/2 blocks
      this->setMatrixDim(K + Tools::gL(),pph2s[K + Tools::gL()].size(),4);//S=3/2 blocks

   }

}

/**
 * copy constructor: constructs BlockMatrix object with M blocks, M/2 for S=1/2 and M/2 for S=3/2, and copies the content of the pphm_c blocks into it,
 * @param pphm_c PPHM object to be copied into (*this)
 */
PPHM::PPHM(const PPHM &pphm_c) : BlockMatrix(pphm_c) { }

/**
 * Destructor
 */
PPHM::~PPHM(){ }

/**
 * access the elements of the matrix in sp mode, special symmetry and antisymmetry relations are automatically accounted for:\n\n
 * @param S The pphm-spin index, when == 0 then access the block S = 1/2, for spinindex == 1 we access the S = 3/2.
 * @param S_ab The intermediate spinquantumnumber of a and b.
 * @param a first sp index that forms the pph row index i together with b, c and S_ab in block B
 * @param b second sp index that forms the pph row index i together with a, c and S_ab in block B
 * @param c third sp index that forms the pph row index i together with a, b and S_ab in block B
 * @param S_de The intermediate spinquantumnumber of d and e.
 * @param d first sp index that forms the pph column index j together with e, z and S_de in block B
 * @param e second sp index that forms the pph column index j together with d, z and S_de in block B
 * @param z third sp index that forms the pph column index j together with d, e and S_de in block B
 * @return the number on place PPHM(B,i,j) with the right phase.
 */
double PPHM::operator()(int S,int S_ab,int a,int b,int c,int S_de,int d,int e,int z) const {

   int K = (a + b + c)%Tools::gL();

   //check the momentum
   if( (d + e + z)%Tools::gL() != K)
      return 0;

   //associated blockindex
   int B = Tools::gL()*S + K;

   int i,j;

   int phase_i = get_inco(B,S_ab,a,b,c,i);

   if(phase_i == 0)
      return 0;

   int phase_j = get_inco(B,S_de,d,e,z,j);

   if(phase_j == 0)
      return 0;

   return phase_i*phase_j* (*this)(B,i,j);

}

/** 
 * Member function that gets the pph-index and phase corresponding to the sp indices S, K, S_ab, a, b, c.
 * @param S spin-index of the state, 0 -> S = 1/2, 1 -> S = 3/2
 * @param K momentum-index of the state, 0 <= K < M/2
 * @param S_ab intermediate spincoupling of a and b. = 0 or 1
 * @param a first sp orbital
 * @param b second sp orbital
 * @param c third sp orbital
 * @param i the corresponding pph index will be stored in this int after calling the function
 * @return the phase needed to get to a normal ordering of indices that corresponds to a pph index i
 */
int PPHM::get_inco(int B,int S_ab,int a,int b,int c,int &i) const{

   if(B < Tools::gL()){//S = 1/2

      if(S_ab == 0){//symmetric in spatial sp's

         if(a <= b)
            i = s2pph[B][0][a][b][c];
         else
            i = s2pph[B][0][b][a][c];

         return 1;

      }
      else{//antisymmetric in spatial sp's

         if(a == b)
            return 0;

         if(a < b){

            i = s2pph[B][1][a][b][c];

            return 1;

         }
         else{

            i = s2pph[B][1][b][a][c];

            return -1;

         }

      }

   }
   else{//S = 3/2

      if(S_ab == 0)//no possibile for S = 3/2
         return 0;

      if(a == b)//no possibile for S = 3/2
         return 0;

      if(a < b){

         i = s2pph[B][0][a][b][c];

         return 1;

      }
      else{

         i = s2pph[B][0][b][a][c];

         return -1;

      }

   }

}

/**
 * The spincoupled, translationally invariant T2 map, maps a TPM onto a PPHM object. See notes for more info
 * be aware that the c and z in the T2 notation become -c and -z in TPM space (remember the G-map)
 * @param tpm input TPM matrix
 */
void PPHM::T(const TPM &tpm){

   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   int a,b,c,d,e,z;
   int S_ab,S_de;

   double norm_ab,norm_de;
   int sign_ab,sign_de;

   //first the S = 1/2 blocks, these should be the most difficult ones.
   for(int B = 0;B < Tools::gL();++B){

      for(int i = 0;i < gdim(B);++i){

         S_ab = pph2s[B][i][0];

         a = pph2s[B][i][1];
         b = pph2s[B][i][2];

         //change to tp-notation:
         c = Tools::par(pph2s[B][i][3]);

         sign_ab = 1 - 2*S_ab;

         norm_ab = 1.0;

         if(a == b)
            norm_ab /= std::sqrt(2.0);

         for(int j = i;j < gdim(B);++j){

            S_de = pph2s[B][j][0];

            d = pph2s[B][j][1];
            e = pph2s[B][j][2];

            //change to tp-notation:
            z = Tools::par(pph2s[B][j][3]);

            sign_de = 1 - 2*S_de;

            norm_de = 1.0;

            if(d == e)
               norm_de /= std::sqrt(2.0);

            //start the map: init
            (*this)(B,i,j) = 0.0;

            //sp term becomes diagonal here:
            if(i == j)
               (*this)(B,i,j) += spm[c];

            //tp(1)
            if(c == z)
               if(S_ab == S_de)
                  (*this)(B,i,j) += tpm(S_ab,a,b,d,e);

            //tp(2)
            if(a == d){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,e,z,b);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == e)
                  ward *= std::sqrt(2.0);

               if(z == b)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= ward;

            }

            //tp(3)
            if(b == d){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,e,z,a);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == e)
                  ward *= std::sqrt(2.0);

               if(z == a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= sign_ab * ward;

            }

            //tp(4)
            if(a == e){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,d,z,b);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == d)
                  ward *= std::sqrt(2.0);

               if(z == b)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= sign_de * ward;

            }

            //tp(5)
            if(b == e){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_de) * tpm(Z,c,d,z,a);

               ward *= norm_ab * norm_de * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) );

               if(c == d)
                  ward *= std::sqrt(2.0);

               if(z == a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= sign_ab * sign_de * ward;

            }

         }

      }

   }

   //the easier S = 3/2 part:
   for(int B = Tools::gL();B < Tools::gM();++B){

      for(int i = 0;i < gdim(B);++i){

         a = pph2s[B][i][1];
         b = pph2s[B][i][2];

         //change to correct sp-momentum
         c = Tools::par(pph2s[B][i][3]);

         for(int j = i;j < gdim(B);++j){

            d = pph2s[B][j][1];
            e = pph2s[B][j][2];

            //change to correct sp-momentum
            z = Tools::par(pph2s[B][j][3]);

            //init
            (*this)(B,i,j) = 0.0;

            //sp part is diagonal
            if(i == j)
               (*this)(B,i,j) += spm[c];

            //tp(1)
            if(c == z)
               (*this)(B,i,j) += tpm(1,a,b,d,e);

            //tp(2)
            if(a == d){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g6j(0,0,1,Z) * tpm(Z,c,e,z,b);

               if(c == e)
                  ward *= std::sqrt(2.0);

               if(z == b)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= ward;

            }

            //tp(3)
            if(b == d){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g6j(0,0,1,Z) * tpm(Z,c,e,z,a);

               if(c == e)
                  ward *= std::sqrt(2.0);

               if(z == a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) += ward;

            }

            //tp(5)
            if(b == e){

               double ward = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  ward += (2*Z + 1.0) * Tools::g6j(0,0,1,Z) * tpm(Z,c,d,z,a);

               if(c == d)
                  ward *= std::sqrt(2.0);

               if(z == a)
                  ward *= std::sqrt(2.0);

               (*this)(B,i,j) -= ward;

            }

         }

      }

   }

   this->symmetrize();

}

ostream &operator<<(ostream &output,const PPHM &pphm_p){

   for(int B = 0;B < pphm_p.gnr();++B){

      output << pphm_p.gblock_char(B,0) << "\t" << pphm_p.gblock_char(B,1) << "\t" << pphm_p.gdim(B) << "\t" << pphm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < pphm_p.gdim(B);++i)
         for(int j = 0;j < pphm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << 

               pphm_p.pph2s[B][i][0] << "\t" << pphm_p.pph2s[B][i][1] << "\t" << pphm_p.pph2s[B][i][2] << "\t" << pphm_p.pph2s[B][i][3] << 

               "\t" << pphm_p.pph2s[B][j][0] << "\t" << pphm_p.pph2s[B][j][1] << "\t" << pphm_p.pph2s[B][j][2] << "\t" << pphm_p.pph2s[B][j][3] 

               << "\t" << pphm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * access to the lists from outside of the class
 */
int PPHM::gpph2s(int B,int i,int option){

   return pph2s[B][i][option];

}

int PPHM::gs2pph(int B,int S_ab,int a,int b,int c){

   return s2pph[B][S_ab][a][b][c];

}

int PPHM::gblock_char(int B,int option){

   return block_char[B][option];

}

/**
 * convert a PPHM matrix to a double array for faster access to the number, fast conversion
 */
void PPHM::convert(double **array) const {

   int L = Tools::gL();
   int L2 = L * L;
   int L3 = L * L2;
   int L4 = L * L3;

   int i,j;
   int c,z;

   int sign_ab,sign_de;

   //first S = 1/2
   for(int K = 0;K < L;++K){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         sign_ab = 1 - 2*S_ab;

         for(int a = 0;a < L;++a)
            for(int b = a + S_ab;b < L;++b){

               c = (K - a - b + 2*L)%L;
               i = s2pph[K][S_ab][a][b][c];

               for(int S_de = 0;S_de < 2;++S_de){

                  sign_de = 1 - 2*S_de;

                  for(int d = 0;d < L;++d)
                     for(int e = d + S_de;e < L;++e){

                        z = (K - d - e + 2*L)%L;
                        j = s2pph[K][S_de][d][e][z];

                        array[K][a + b*L + d*L2 + e*L3 + S_ab*L4 + 2*S_de*L4] = (*this)(K,i,j);
                        array[K][b + a*L + d*L2 + e*L3 + S_ab*L4 + 2*S_de*L4] = sign_ab * array[K][a + b*L + d*L2 + e*L3 + S_ab*L4 + 2*S_de*L4];
                        array[K][a + b*L + e*L2 + d*L3 + S_ab*L4 + 2*S_de*L4] = sign_de * array[K][a + b*L + d*L2 + e*L3 + S_ab*L4 + 2*S_de*L4];
                        array[K][b + a*L + e*L2 + d*L3 + S_ab*L4 + 2*S_de*L4] = sign_ab * sign_de * array[K][a + b*L + d*L2 + e*L3 + S_ab*L4 + 2*S_de*L4];

                     }

               }

            }

      }

   }//end of S=1/2 block loop

   //S = 3/2 is easier
   for(int B = L;B < 2*L;++B){

      int K = block_char[B][1];

      for(int a = 0;a < L;++a)
         for(int b = a + 1;b < L;++b){

            c = (K - a - b + 2*L)%L;
            i = s2pph[B][0][a][b][c];

            for(int d = 0;d < L;++d)
               for(int e = d + 1;e < L;++e){

                  z = (K - d - e + 2*L)%L;
                  j = s2pph[B][0][d][e][z];

                  array[B][a + b*L + d*L2 + e*L3] = (*this)(K,i,j);
                  array[B][b + a*L + d*L2 + e*L3] =  -array[K][a + b*L + d*L2 + e*L3];
                  array[B][a + b*L + e*L2 + d*L3] =  -array[K][a + b*L + d*L2 + e*L3];
                  array[B][b + a*L + e*L2 + d*L3] =  array[K][a + b*L + d*L2 + e*L3];

               }

         }

   }

}

/* vim: set ts=3 sw=3 expandtab :*/
