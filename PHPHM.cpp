#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

int ***PHPHM::ph2phmm;
vector< vector<int> > PHPHM::phmm2ph;

/**
 * initialize the static lists
 */
void PHPHM::init(){

   int L = Tools::gL();

   ph2phmm = new int ** [2*L];

   for(int B = 0;B < 2*L;++B){

      ph2phmm[B] = new int * [PHM::gdim(B)];

      for(int i = 0;i < PHM::gdim(B);++i)
      ph2phmm[B][i] = new int [PHM::gdim(B)];

   }

   vector<int> v(3);

   int phmm = 0;

   for(int B = 0;B < 2*L;++B){

      for(int i = 0;i < PHM::gdim(B);++i)
         for(int j = i;j < PHM::gdim(B);++j){

            v[0] = B;
            v[1] = i;
            v[2] = j;

            phmm2ph.push_back(v);

            ph2phmm[B][i][j] = phmm;
            ph2phmm[B][j][i] = phmm;

            ++phmm;

         }

   }

}

/**
 * deallocate the static lists
 */
void PHPHM::clear(){

   for(int B = 0;B < 2*Tools::gL();++B){

      for(int i = 0;i < PHM::gdim(B);++i)
         delete [] ph2phmm[B][i];

      delete [] ph2phmm[B];

   }

   delete [] ph2phmm;

}

/**
 * standard constructor:
 */
PHPHM::PHPHM() : Matrix(phmm2ph.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix phmm_c
 * @param phmm_c object that will be copied into this.
 */
PHPHM::PHPHM(const PHPHM &phmm_c) : Matrix(phmm_c){ }

/**
 * destructor
 */
PHPHM::~PHPHM(){ }

ostream &operator<<(ostream &output,const PHPHM &phmm_p){

   int B,I,J,B_,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < PHPHM::gn();++i){

      B = phmm_p.phmm2ph[i][0];
      I = phmm_p.phmm2ph[i][1];
      J = phmm_p.phmm2ph[i][2];

      a = PHM::gph2s(B,I,0);
      b = PHM::gph2s(B,I,1);

      c = PHM::gph2s(B,J,0);
      d = PHM::gph2s(B,J,1);

      for(int j = i;j < PHPHM::gn();++j){

         B_ = phmm_p.phmm2ph[j][0]; 
         K = phmm_p.phmm2ph[j][1];
         L = phmm_p.phmm2ph[j][2];

         e = PHM::gph2s(B_,K,0);
         z = PHM::gph2s(B_,K,1);

         t = PHM::gph2s(B_,L,0);
         h = PHM::gph2s(B_,L,1);

         output << i << "\t" << j << "\t|\t(" << B << ")\t" << I << "\t" << J << "\t(" << B_ << ")\t" << K << "\t" << L << "\t|\t" << 

            "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << phmm_p(i,j) << endl;

      }

   }

   return output;

}

/**
 * @return the dimension of a PHPHM matrix
 */
int PHPHM::gn(){

   return phmm2ph.size();

}

/**
 * access to the lists from outside the class
 */
int PHPHM::gph2phmm(int B,int I,int J){

   return ph2phmm[B][I][J];

}

/**
 * access to the lists from outside the class
 * @param option == 0 return B, == 1 return I, == 2 return J
 */
int PHPHM::gphmm2ph(int i,int option){

   return phmm2ph[i][option];

}
