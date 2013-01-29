#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::endl;

#include "include.h"

vector< vector<int> > *PHM::ph2s;
int ***PHM::s2ph;

int **PHM::block_char;

/** 
 * initialize the statics
 */
void PHM::init(){

   int L = Tools::gL();

   //allocation of block_char
   block_char = new int * [2*L];

   for(int B = 0;B < 2*L;++B)
      block_char[B] = new int [2];

   //initialization of block_char
   for(int K = 0;K < L;++K){

      block_char[K][0] = 0;
      block_char[K][1] = K;

      block_char[K + L][0] = 1;
      block_char[K + L][1] = K;

   }

   //allocation of s2ph
   s2ph = new int ** [2*L];

   for(int B = 0;B < 2*L;++B){

      s2ph[B] = new int * [L];

      for(int i = 0;i < L;++i)
         s2ph[B][i] = new int [L];

   }

   //allocation of ph2s
   ph2s = new vector< vector<int> > [2*L];

   vector<int> v(2);

   //initialization of the lists
   int ph;

   int K;

   for(int B = 0;B < 2*L;++B){

      ph = 0;

      K = block_char[B][1];

      for(int a = 0;a < L;++a)
         for(int b = 0;b < L;++b){

            //if momentumconservation is satisfied, add elements to list
            if( (a + b)%L == K){

               v[0] = a;
               v[1] = b;

               ph2s[B].push_back(v);

               s2ph[B][a][b] = ph;

               ++ph;

            }

         }

   }

}

/**
 * deallocate lists
 */
void PHM::clear(){

   int L = Tools::gL();

   for(int B = 0;B < 2*L;++B){

      for(int k = 0;k < L;++k)
         delete [] s2ph[B][k];

      delete [] s2ph[B];

      delete [] block_char[B];

   }

   delete [] s2ph;
   delete [] ph2s;

   delete [] block_char;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 0 and 1 of dimension M*M/4.
 */
PHM::PHM() : BlockMatrix(2*Tools::gL()) {

   //set the dimension of the blocks
   for(int K = 0;K < Tools::gL();++K){

      this->setMatrixDim(K,ph2s[K].size(),1);
      this->setMatrixDim(K + Tools::gL(),ph2s[K + Tools::gL()].size(),3);

   }

}

/**
 * copy constructor
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(const PHM &phm_c) : BlockMatrix(phm_c){ }

/**
 * destructor
 */
PHM::~PHM(){ }

/**
 * access the elements of the matrix in sp mode, 
 * @param S The blockindex of the block you want to access
 * @param a first sp momentum index that forms the ph row index i in block B together with b
 * @param b second sp momentum index that forms the ph row index i in block B together with a
 * @param c first sp momentum index that forms the ph column index j in block B together with d
 * @param d second sp momentum index that forms the ph column index j in block B together with c
 * @return the number on place PHM(i,j)
 */
 
double PHM::operator()(int S,int a,int b,int c,int d) const{

   int K = (a + b)%Tools::gL();

   //momentum conservation:
   if( (c + d)%Tools::gL() != K)
      return 0;

   int B = K + S*Tools::gL();

   int i = s2ph[B][a][b];
   int j = s2ph[B][c][d];

   return (*this)(B,i,j);

}

ostream &operator<<(ostream &output,const PHM &phm_p){

   int S,K;

   for(int B = 0;B < phm_p.gnr();++B){

      S = phm_p.block_char[B][0];
      K = phm_p.block_char[B][1];

      output << "S =\t" << S << "\tK =\t" << K << "\tdimension =\t" << phm_p.gdim(B) << "\tdegeneracy =\t" << phm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < phm_p.gdim(B);++i)
         for(int j = 0;j < phm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << phm_p.ph2s[B][i][0] << "\t" << phm_p.ph2s[B][i][1]

               << "\t" << phm_p.ph2s[B][j][0] << "\t" << phm_p.ph2s[B][j][1] << "\t" << phm_p(B,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * The G map, maps a TPM object on a PHM object.
 * @param tpm input TPM
 */
void PHM::G(const TPM &tpm){

   int L = Tools::gL();

   //construct the SPM corresponding to the TPM
   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   int a,b,c,d;

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         a = ph2s[B][i][0];

         //transform b to tpm sp-momentum:
         b = (-ph2s[B][i][1] + L)%L;

         for(int j = i;j < gdim(B);++j){

            c = ph2s[B][j][0];

            //transform d to tpm sp-momentum:
            d = (-ph2s[B][j][1] + L)%L;

            (*this)(B,i,j) = -Tools::g6j(0,0,0,S) * tpm(0,a,d,c,b) / ( TPM::gnorm(a,d) * TPM::gnorm(c,b) ) - 3.0 * Tools::g6j(0,0,1,S) * tpm(1,a,d,c,b);

         }

         (*this)(B,i,i) += spm[a];

      }

   }

   this->symmetrize();

}

/**
 * convert a PHM matrix to a double array for faster access to the number
 */
void PHM::convert(double **array) const {

   int K,S;

   int L = Tools::gL();

   for(int B = 0;B < 2*L;++B){

      S = block_char[B][0];
      K = block_char[B][1];

      for(int a = 0;a < L;++a)
         for(int c = 0;c < L;++c)
            array[B][a + c*L] = (*this)(S,a,(K - a + L)%L,c,(K - c + L)%L);

   }

}

/**
 * The bar function that maps a PPHM object onto a PHM object by tracing away the first pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void PHM::bar(double scale,const PPHM &pphm){

   int a,b,c,d;

   double ward,hard;

   int S;

   for(int B = 0;B < gnr();++B){//loop over the blocks PHM

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         a = ph2s[B][i][0];
         b = ph2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = ph2s[B][j][0];
            d = ph2s[B][j][1];

            //init
            (*this)(B,i,j) = 0.0;

            //first the S = 1/2 block of the PPHM matrix
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_de = 0;S_de < 2;++S_de){

                  ward = 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_de + 1.0) ) * Tools::g6j(0,0,S,S_ab) * Tools::g6j(0,0,S,S_de);

                  for(int p = 0;p < Tools::gL();++p){

                     hard = ward * pphm(0,S_ab,p,a,b,S_de,p,c,d);

                     //norms
                     if(p == a)
                        hard *= std::sqrt(2.0);

                     if(p == c)
                        hard *= std::sqrt(2.0);

                     (*this)(B,i,j) += hard;

                  }

               }

            //then the S = 3/2 block
            if(S == 1)
               for(int p = 0;p < Tools::gL();++p)
                  (*this)(B,i,j) += 4.0/3.0 * pphm(1,1,p,a,b,1,p,c,d);

            (*this)(B,i,j) *= scale;

         }
      }

   }

   this->symmetrize();

}

/**
 * @return the dimension of the block corresponding to index B
 */
int PHM::gdim(int B) {

   return ph2s[B].size();

}

/**
 * access to the lists from outside of the class
 */
int PHM::gph2s(int B,int i,int option){

   return ph2s[B][i][option];

}

/**
 * access to the lists from outside of the class
 */
int PHM::gs2ph(int B,int a,int b){

   return s2ph[B][a][b];

}

/**
 * access to the lists from outside of the class
 */
int PHM::gblock_char(int block,int option){

   return block_char[block][option];

}

/* vim: set ts=3 sw=3 expandtab :*/
