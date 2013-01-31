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

int ***TPTPM::t2tpmm;
vector< vector<int> > TPTPM::tpmm2t;

/**
 * initialize the static lists
 */
void TPTPM::init(){

   int L = Tools::gL();

   t2tpmm = new int ** [2*L];

   for(int B = 0;B < 2*L;++B){

      t2tpmm[B] = new int * [TPM::gdim(B)];

      for(int i = 0;i < TPM::gdim(B);++i)
      t2tpmm[B][i] = new int [TPM::gdim(B)];

   }

   vector<int> v(3);

   int tpmm = 0;

   for(int B = 0;B < 2*L;++B){

      for(int i = 0;i < TPM::gdim(B);++i)
         for(int j = i;j < TPM::gdim(B);++j){

            v[0] = B;
            v[1] = i;
            v[2] = j;

            tpmm2t.push_back(v);

            t2tpmm[B][i][j] = tpmm;
            t2tpmm[B][j][i] = tpmm;

            ++tpmm;

         }

   }

}

/**
 * deallocate the static lists
 */
void TPTPM::clear(){

   for(int B = 0;B < 2*Tools::gL();++B){

      for(int i = 0;i < TPM::gdim(B);++i)
         delete [] t2tpmm[B][i];

      delete [] t2tpmm[B];

   }

   delete [] t2tpmm;

}

/**
 * standard constructor:
 */
TPTPM::TPTPM() : Matrix(tpmm2t.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpmm_c
 * @param tpmm_c object that will be copied into this.
 */
TPTPM::TPTPM(const TPTPM &tpmm_c) : Matrix(tpmm_c){ }

/**
 * destructor
 */
TPTPM::~TPTPM(){ }

/**
 * access the elements of the matrix in tp mode
 * @param B block index of the first two indices
 * @param I first tp index that forms the tpmm row index i together with J
 * @param J second tp index that forms the tpmm row index i together with I
 * @param B_ block index of the second two indices
 * @param K first tp index that forms the tpmm column index j together with L
 * @param L second tp index that forms the tpmm column index j together with K
 * @return the number on place TPTPM(i,j)
 */
double TPTPM::operator()(int B,int I,int J,int B_,int K,int L) const{

   int i = t2tpmm[B][I][J];
   int j = t2tpmm[B_][K][L];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,const TPTPM &tpmm_p){

   int B,I,J,B_,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      B = tpmm_p.tpmm2t[i][0];
      I = tpmm_p.tpmm2t[i][1];
      J = tpmm_p.tpmm2t[i][2];

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);

      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int j = i;j < TPTPM::gn();++j){

         B_ = tpmm_p.tpmm2t[j][0]; 
         K = tpmm_p.tpmm2t[j][1];
         L = tpmm_p.tpmm2t[j][2];

         e = TPM::gt2s(B_,K,0);
         z = TPM::gt2s(B_,K,1);

         t = TPM::gt2s(B_,L,0);
         h = TPM::gt2s(B_,L,1);

         output << i << "\t" << j << "\t|\t(" << B << ")\t" << I << "\t" << J << "\t(" << B_ << ")\t" << K << "\t" << L << "\t|\t" << 

            "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << tpmm_p(i,j) << endl;

      }

   }

   return output;

}

/**
 * @return the dimension of a TPTPM matrix
 */
int TPTPM::gn(){

   return tpmm2t.size();

}

/**
 * access to the lists from outside the class
 */
int TPTPM::gt2tpmm(int B,int I,int J){

   return t2tpmm[B][I][J];

}

/**
 * access to the lists from outside the class
 * @param option == 0 return B, == 1 return a, == 2 return b
 */
int TPTPM::gtpmm2t(int i,int option){

   return tpmm2t[i][option];

}

/**
 * construct the antisymmetrized direct product of two PHM matrices
 */
void TPTPM::dp(const PHM &phm){

   int L = Tools::gL();

   double **pharray = new double * [2*L];

   for(int B = 0;B < 2*L;++B)
      pharray[B] = new double [L*L];

   phm.convert(pharray);

   int B,B_;

   int sign,sign_;

   int a,b,c,d;
   int e,z,t,h;

   int d_;
   int t_,h_;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      sign = 1 - 2*S;

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      d_ = Tools::par(d); 

      for(int j = i;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         sign_ = 1 - 2*S_;

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         t_ = Tools::par(t); 
         h_ = Tools::par(h); 

         (*this)(i,j) = 0.0;

         int P = (a + d_)%L;
         int P_ = Tools::par(P);

         //(a,d,c,b)_(e,h,t,z) and (b,c,d,a)_(z,t,h,e)
         if(P == (e + h_)%L){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L][a + e*L] * pharray[P + Z*L][c + t*L]

                     + pharray[P + Z*L][a + t*L]* pharray[P + Z*L][c + e*L] + pharray[P_ + Z*L][b + z*L] * pharray[P_ + Z*L][d + h*L]

                     + pharray[P_ + Z*L][b + h*L] * pharray[P_ + Z*L][d + z*L] ); 

            }

         }

         //(a,d,c,b)_(z,h,t,e) and (b,c,d,a)_(e,t,h,z)
         if(P == (z + h_)%L){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign_ * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L][a + z*L] * pharray[P + Z*L][c + t*L]

                     + pharray[P + Z*L][a + t*L]* pharray[P + Z*L][c + z*L] + pharray[P_ + Z*L][b + e*L] * pharray[P_ + Z*L][d + h*L]

                     + pharray[P_ + Z*L][b + h*L] * pharray[P_ + Z*L][d + e*L] ); 

            }

         }

         //(a,d,c,b)_(e,t,h,z) and (b,c,d,a)_(z,h,t,e)
         if(P == (e + t_)%L){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign_ * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L][a + e*L] * pharray[P + Z*L][c + h*L]

                     + pharray[P + Z*L][a + h*L]* pharray[P + Z*L][c + e*L]  + pharray[P_ + Z*L][b + z*L] * pharray[P_ + Z*L][d + t*L]

                     + pharray[P_ + Z*L][b + t*L] * pharray[P_ + Z*L][d + z*L] ); 

            }

         }

         //(a,d,c,b)_(z,t,h,e) and (b,c,d,a)_(e,h,t,z)
         if(P == (z + t_)%L){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L][a + z*L] * pharray[P + Z*L][c + h*L]

                     + pharray[P + Z*L][a + h*L]* pharray[P + Z*L][c + z*L] + pharray[P_ + Z*L][b + e*L] * pharray[P_ + Z*L][d + t*L]

                     + pharray[P_ + Z*L][b + t*L] * pharray[P_ + Z*L][d + e*L] ); 

            }

         }

         P = (b + d_)%L;
         P_ = Tools::par(P);

         //(b,d,c,a)_(e,h,t,z) and (a,c,d,b)_(z,t,h,e)
         if(P == (e + h_)%L){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L][b + e*L] * pharray[P + Z*L][c + t*L]

                     + pharray[P + Z*L][b + t*L]* pharray[P + Z*L][c + e*L] + pharray[P_ + Z*L][a + z*L] * pharray[P_ + Z*L][d + h*L]

                     + pharray[P_ + Z*L][a + h*L] * pharray[P_ + Z*L][d + z*L] ); 

            }

         }

         //(b,d,c,a)_(z,h,t,e) and (a,c,d,b)_(e,t,h,z)
         if(P == (z + h_)%L){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign * sign_ * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L][b + z*L] * pharray[P + Z*L][c + t*L]

                     + pharray[P + Z*L][b + t*L]* pharray[P + Z*L][c + z*L] + pharray[P_ + Z*L][a + e*L] * pharray[P_ + Z*L][d + h*L]

                     + pharray[P_ + Z*L][a + h*L] * pharray[P_ + Z*L][d + e*L] ); 

            }

         }

         //(b,d,c,a)_(e,t,h,z) and (a,c,d,b)_(z,h,t,e)
         if(P == (e + t_)%L){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign * sign_ * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L][b + e*L] * pharray[P + Z*L][c + h*L]

                     + pharray[P + Z*L][b + h*L]* pharray[P + Z*L][c + e*L] + pharray[P_ + Z*L][a + z*L] * pharray[P_ + Z*L][d + t*L]

                     + pharray[P_ + Z*L][a + t*L] * pharray[P_ + Z*L][d + z*L] ); 

            }

         }

         //(b,d,c,a)_(z,t,h,e) and (a,c,d,b)_(e,h,t,z)
         if(P == (z + t_)%L){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L][b + z*L] * pharray[P + Z*L][c + h*L]

                     + pharray[P + Z*L][b + h*L]* pharray[P + Z*L][c + z*L] + pharray[P_ + Z*L][a + e*L] * pharray[P_ + Z*L][d + t*L]

                     + pharray[P_ + Z*L][a + t*L] * pharray[P_ + Z*L][d + e*L] ); 

            }

         }

      }
   }

   for(int B = 0;B < 2*L;++B)
      delete [] pharray[B];

   delete [] pharray;

}

/**
 * construct a TPTPM object by double tracing the direct product of two DPM's
 */
void TPTPM::dpt2(const DPM &dpm){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   double **dparray = new double * [2*L];

   for(int B = 0;B < L;++B)//S = 1/2
      dparray[B] = new double [4*L4];

   for(int B = L;B < 2*L;++B)//S = 3/2
      dparray[B] = new double [L4];

   dpm.convert(dparray);

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      for(int j = i;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         double ward = 0.0;

         //first S = 1/2 part
         for(int k = 0;k < L;++k){

            int K_dp = (a + b + k)%L;

            ward += dparray[K_dp][a + b*L + e*L2 + z*L3 + S*L4 + 2*S_*L4] * dparray[K_dp][c + d*L + t*L2 + h*L3 + S*L4 + 2*S_*L4]

               + dparray[K_dp][a + b*L + t*L2 + h*L3 + S*L4 + 2*S_*L4] * dparray[K_dp][c + d*L + e*L2 + z*L3 + S*L4 + 2*S_*L4];

         }

         (*this)(i,j) = 2.0 / ( (2*S + 1.0)*(2*S_ + 1.0) ) * ward;

         if(S == 1 && S_ == 1){//only then contribution from S = 3/2 part

            ward = 0.0;

            for(int k = 0;k < L;++k){

               int K_dp = (a + b + k)%L;

               ward += dparray[K_dp + L][a + b*L + e*L2 + z*L3] * dparray[K_dp + L][c + d*L + t*L2 + h*L3]

                  + dparray[K_dp + L][a + b*L + t*L2 + h*L3] * dparray[K_dp + L][c + d*L + e*L2 + z*L3];

            }

            (*this)(i,j) += 4.0 / 9.0  * ward;

         }

      }
   }

   //remove the array
   for(int B = 0;B < 2*L;++B)
      delete [] dparray[B];

   delete [] dparray;

}

/**
 * construct a TPTPM object by double tracing the direct product of two DPM's, already transformed to a double **array
 */
void TPTPM::dpt2(double **dparray){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      for(int j = i;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         double ward = 0.0;

         //first S = 1/2 part
         for(int k = 0;k < L;++k){

            int K_dp = (a + b + k)%L;

            ward += dparray[K_dp][a + b*L + e*L2 + z*L3 + S*L4 + 2*S_*L4] * dparray[K_dp][c + d*L + t*L2 + h*L3 + S*L4 + 2*S_*L4]

               + dparray[K_dp][a + b*L + t*L2 + h*L3 + S*L4 + 2*S_*L4] * dparray[K_dp][c + d*L + e*L2 + z*L3 + S*L4 + 2*S_*L4];

         }

         (*this)(i,j) = 2.0 / ( (2*S + 1.0)*(2*S_ + 1.0) ) * ward;

         if(S == 1 && S_ == 1){//only then contribution from S = 3/2 part

            ward = 0.0;

            for(int k = 0;k < L;++k){

               int K_dp = (a + b + k)%L;

               ward += dparray[K_dp + L][a + b*L + e*L2 + z*L3] * dparray[K_dp + L][c + d*L + t*L2 + h*L3]

                  + dparray[K_dp + L][a + b*L + t*L2 + h*L3] * dparray[K_dp + L][c + d*L + e*L2 + z*L3];

            }

            (*this)(i,j) += 4.0 / 9.0  * ward;

         }

      }
   }

 }

/**
 * construct a TPTPM by double tracing the direct product of two PPHM matrices, already translated to 'array' for for faster access
 */
void TPTPM::dpt2_pph(double **ppharray){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      for(int j = i;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         double ward = 0.0;

         //first S = 1/2 part
         for(int k = 0;k < L;++k){

            int K_pph = (a + b + k)%L;

            ward += ppharray[K_pph][a + b*L + e*L2 + z*L3 + S*L4 + 2*S_*L4] * ppharray[K_pph][c + d*L + t*L2 + h*L3 + S*L4 + 2*S_*L4]

               + ppharray[K_pph][a + b*L + t*L2 + h*L3 + S*L4 + 2*S_*L4] * ppharray[K_pph][c + d*L + e*L2 + z*L3 + S*L4 + 2*S_*L4];

         }

         (*this)(i,j) = 2.0 / ( (2*S + 1.0)*(2*S_ + 1.0) ) * ward;

         if(S == 1 && S_ == 1){//only then contribution from S = 3/2 part

            ward = 0.0;

            for(int k = 0;k < L;++k){

               int K_pph = (a + b + k)%L;

               ward += ppharray[K_pph + L][a + b*L + e*L2 + z*L3] * ppharray[K_pph + L][c + d*L + t*L2 + h*L3]

                  + ppharray[K_pph + L][a + b*L + t*L2 + h*L3] * ppharray[K_pph + L][c + d*L + e*L2 + z*L3];

            }

            (*this)(i,j) += 4.0 / 9.0  * ward;

         }


      }
   }

}

/**
 * construct a TPTPM by once tracing and once skew-tracing the direct product of two PPHM matrices, already translated to 'array' for for faster access
 */
void TPTPM::dptw(double **ppharray){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int a_,b_;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   int sign;

   double ward;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      sign = 1 - 2*S;

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      a_ = Tools::par(a);
      b_ = Tools::par(b);

      for(int j = 0;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         (*this)(i,j) = 0.0;

         //S'' = 1/2 first
         for(int J = 0;J < 2;++J){

            ward = 0.0;

            for(int k = 0;k < L;++k){

               //(da;bc)
               int K_pph = (k + d + a_)%L;

               ward += ppharray[K_pph][e + z*L + k*L2 + d*L3 + S_*L4 + 2*J*L4] *  ppharray[K_pph][t + h*L + k*L2 + b*L3 + S_*L4 + 2*J*L4]

                  + ppharray[K_pph][e + z*L + k*L2 + b*L3 + S_*L4 + 2*J*L4] *  ppharray[K_pph][t + h*L + k*L2 + d*L3 + S_*L4 + 2*J*L4];

               K_pph = (k + d + b_)%L;

               ward += sign * ( ppharray[K_pph][e + z*L + k*L2 + d*L3 + S_*L4 + 2*J*L4] *  ppharray[K_pph][t + h*L + k*L2 + a*L3 + S_*L4 + 2*J*L4]

                  + ppharray[K_pph][e + z*L + k*L2 + a*L3 + S_*L4 + 2*J*L4] *  ppharray[K_pph][t + h*L + k*L2 + d*L3 + S_*L4 + 2*J*L4] );

               K_pph = (k + c + a_)%L;

               ward += sign * ( ppharray[K_pph][e + z*L + k*L2 + c*L3 + S_*L4 + 2*J*L4] *  ppharray[K_pph][t + h*L + k*L2 + b*L3 + S_*L4 + 2*J*L4]

                  + ppharray[K_pph][e + z*L + k*L2 + b*L3 + S_*L4 + 2*J*L4] *  ppharray[K_pph][t + h*L + k*L2 + c*L3 + S_*L4 + 2*J*L4] );

               K_pph = (k + c + b_)%L;

               ward += ppharray[K_pph][e + z*L + k*L2 + c*L3 + S_*L4 + 2*J*L4] *  ppharray[K_pph][t + h*L + k*L2 + a*L3 + S_*L4 + 2*J*L4]

                  + ppharray[K_pph][e + z*L + k*L2 + a*L3 + S_*L4 + 2*J*L4] *  ppharray[K_pph][t + h*L + k*L2 + c*L3 + S_*L4 + 2*J*L4];

            }

            (*this)(i,j) += ward * (2*J + 1.0) * Tools::g6j(0,0,S,J);

         }

         (*this)(i,j) *= 2.0 / (2*S_ + 1.0);

         //then S'' = 3/2, only
         if(S_ == 1){

            ward = 0.0;

            for(int k = 0;k < L;++k){

               //(da;bc)
               int K_pph = (k + d + a_)%L;

               ward += ppharray[K_pph + L][e + z*L + k*L2 + d*L3] *  ppharray[K_pph + L][t + h*L + k*L2 + b*L3]

                  + ppharray[K_pph + L][e + z*L + k*L2 + b*L3] *  ppharray[K_pph + L][t + h*L + k*L2 + d*L3];

               K_pph = (k + d + b_)%L;

               ward -= ppharray[K_pph + L][e + z*L + k*L2 + d*L3] *  ppharray[K_pph + L][t + h*L + k*L2 + a*L3]

                  + ppharray[K_pph + L][e + z*L + k*L2 + a*L3] *  ppharray[K_pph + L][t + h*L + k*L2 + d*L3];

               K_pph = (k + c + a_)%L;

               ward -= ppharray[K_pph + L][e + z*L + k*L2 + c*L3] *  ppharray[K_pph + L][t + h*L + k*L2 + b*L3]

                  + ppharray[K_pph + L][e + z*L + k*L2 + b*L3] *  ppharray[K_pph + L][t + h*L + k*L2 + c*L3];

               K_pph = (k + c + b_)%L;

               ward += ppharray[K_pph + L][e + z*L + k*L2 + c*L3] *  ppharray[K_pph + L][t + h*L + k*L2 + a*L3]

                  + ppharray[K_pph + L][e + z*L + k*L2 + a*L3] *  ppharray[K_pph + L][t + h*L + k*L2 + c*L3];

            }

            (*this)(i,j) +=  4.0 * ward * Tools::g6j(0,0,S,1);

         }

      }
   }

}

/**
 * construct a TPTPM by twice skew-tracing the direct product of two PPHM matrices, already translated to 'array' for for faster access
 */
void TPTPM::dpw2(double **ppharray){

   int L = Tools::gL();
   int L2 = L*L;
   int L3 = L2*L;
   int L4 = L3*L;

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int a_,b_;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   int sign,sign_;

   double ward;

   int K_pph,m;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      sign = 1 - 2*S;

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      a_ = Tools::par(a);
      b_ = Tools::par(b);

      for(int j = i;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         sign_ = 1 - 2*S_;

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         (*this)(i,j) = 0.0;

         //first S'' = 1/2
         for(int J = 0;J < 2;++J)
            for(int J_ = 0;J_ < 2;++J_){

               ward = 0.0;

               for(int k = 0;k < L;++k){

                  K_pph = (k + d + a_)%L;

                  //1) da;bc he;zt
                  m = (K_pph - h + e + L)%L;

                  ward += ppharray[K_pph][k + d*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + z*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + d*L + m*L2 + z*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + h*L3 + J*L4 + 2*J_*L4];

                  //2) da;bc hz;et
                  m = (K_pph - h + z + L)%L;

                  ward += sign_ * ( ppharray[K_pph][k + d*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + e*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + d*L + m*L2 + e*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] );

                  //3) da;bc te;zh
                  m = (K_pph - t + e + L)%L;

                  ward += sign_ * (ppharray[K_pph][k + d*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + z*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + d*L + m*L2 + z*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] );

                  //4) da;bc tz;eh
                  m = (K_pph - t + z + L)%L;

                  ward += ppharray[K_pph][k + d*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + e*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + d*L + m*L2 + e*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + t*L3 + J*L4 + 2*J_*L4];

                  K_pph = (k + d + b_)%L;

                  //5) db;ac he;zt
                  m = (K_pph - h + e + L)%L;

                  ward += sign * ( ppharray[K_pph][k + d*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + z*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + d*L + m*L2 + z*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] );

                  //6) db;ac hz;et
                  m = (K_pph - h + z + L)%L;

                  ward += sign * sign_ * ( ppharray[K_pph][k + d*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + e*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + d*L + m*L2 + e*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] );

                  //7) db;ac te;zh
                  m = (K_pph - t + e + L)%L;

                  ward += sign * sign_ * ( ppharray[K_pph][k + d*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + z*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + d*L + m*L2 + z*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] );

                  //8) db;ac tz;eh
                  m = (K_pph - t + z + L)%L;

                  ward += sign * ( ppharray[K_pph][k + d*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + e*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + d*L + m*L2 + e*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] );

                  K_pph = (k + c + a_)%L;

                  //9) ca;bd he;zt
                  m = (K_pph - h + e + L)%L;

                  ward += sign * ( ppharray[K_pph][k + c*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + z*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + c*L + m*L2 + z*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] );

                  //10) ca;bd hz;et
                  m = (K_pph - h + z + L)%L;

                  ward += sign * sign_ * ( ppharray[K_pph][k + c*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + e*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + c*L + m*L2 + e*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] );

                  //11) ca;bd te;zh
                  m = (K_pph - t + e + L)%L;

                  ward += sign * sign_ * ( ppharray[K_pph][k + c*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + z*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + c*L + m*L2 + z*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] );

                  //12) ca;bd tz;eh
                  m = (K_pph - t + z + L)%L;

                  ward += sign * ( ppharray[K_pph][k + c*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + e*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + c*L + m*L2 + e*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + b*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] );

                  K_pph = (k + c + b_)%L;

                  //13) cb;ad he;zt
                  m = (K_pph - h + e + L)%L;

                  ward += ppharray[K_pph][k + c*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + z*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + c*L + m*L2 + z*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + h*L3 + J*L4 + 2*J_*L4];

                  //14) cb;ad hz;et
                  m = (K_pph - h + z + L)%L;

                  ward += sign_ * ( ppharray[K_pph][k + c*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + e*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + c*L + m*L2 + e*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + h*L3 + J*L4 + 2*J_*L4] );

                  //15) cb;ad te;zh
                  m = (K_pph - t + e + L)%L;

                  ward += sign_ * ( ppharray[K_pph][k + c*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + z*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + c*L + m*L2 + z*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] );

                  //16) cb;ad tz;eh
                  m = (K_pph - t + z + L)%L;

                  ward += ppharray[K_pph][k + c*L + m*L2 + t*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + e*L3 + J*L4 + 2*J_*L4]

                     + ppharray[K_pph][k + c*L + m*L2 + e*L3 + J*L4 + 2*J_*L4] * ppharray[K_pph][k + a*L + m*L2 + t*L3 + J*L4 + 2*J_*L4];

               }//end of loop over k
               
               (*this)(i,j) += (2*J + 1.0) * (2*J_ + 1.0) * Tools::g6j(0,0,S,J) * Tools::g6j(0,0,S_,J_) * ward;

            }//end of loops over J,J_

         (*this)(i,j) *= 2.0;

         //S = 3/2
         ward = 0.0;

         for(int k = 0;k < L;++k){

            K_pph = (k + d + a_)%L;

            //1) da;bc he;zt
            m = (K_pph - h + e + L)%L;

            ward += ppharray[K_pph + L][k + d*L + m*L2 + h*L3] * ppharray[K_pph + L][k + b*L + m*L2 + z*L3]

               + ppharray[K_pph + L][k + d*L + m*L2 + z*L3] * ppharray[K_pph + L][k + b*L + m*L2 + h*L3];

            //2) da;bc hz;et
            m = (K_pph - h + z + L)%L;

            ward += sign_ * ( ppharray[K_pph + L][k + d*L + m*L2 + h*L3] * ppharray[K_pph + L][k + b*L + m*L2 + e*L3]

               + ppharray[K_pph + L][k + d*L + m*L2 + e*L3] * ppharray[K_pph + L][k + b*L + m*L2 + h*L3] );

            //3) da;bc te;zh
            m = (K_pph - t + e + L)%L;

            ward += sign_ * ( ppharray[K_pph + L][k + d*L + m*L2 + t*L3] * ppharray[K_pph + L][k + b*L + m*L2 + z*L3]

               + ppharray[K_pph + L][k + d*L + m*L2 + z*L3] * ppharray[K_pph + L][k + b*L + m*L2 + t*L3] );

            //4) da;bc tz;eh
            m = (K_pph - t + z + L)%L;

            ward += ppharray[K_pph + L][k + d*L + m*L2 + t*L3] * ppharray[K_pph + L][k + b*L + m*L2 + e*L3]

               + ppharray[K_pph + L][k + d*L + m*L2 + e*L3] * ppharray[K_pph + L][k + b*L + m*L2 + t*L3];

            K_pph = (k + d + b_)%L;

            //5) db;ac he;zt
            m = (K_pph - h + e + L)%L;

            ward += sign * ( ppharray[K_pph + L][k + d*L + m*L2 + h*L3] * ppharray[K_pph + L][k + a*L + m*L2 + z*L3]

               + ppharray[K_pph + L][k + d*L + m*L2 + z*L3] * ppharray[K_pph + L][k + a*L + m*L2 + h*L3] );

            //6) db;ac hz;et
            m = (K_pph - h + z + L)%L;

            ward += sign * sign_ * ( ppharray[K_pph + L][k + d*L + m*L2 + h*L3] * ppharray[K_pph + L][k + a*L + m*L2 + e*L3]

               + ppharray[K_pph + L][k + d*L + m*L2 + e*L3] * ppharray[K_pph + L][k + a*L + m*L2 + h*L3] );

            //7) db;ac te;zh
            m = (K_pph - t + e + L)%L;

            ward += sign * sign_ * ( ppharray[K_pph + L][k + d*L + m*L2 + t*L3] * ppharray[K_pph + L][k + a*L + m*L2 + z*L3]

               + ppharray[K_pph + L][k + d*L + m*L2 + z*L3] * ppharray[K_pph + L][k + a*L + m*L2 + t*L3] );

            //8) db;ac tz;eh
            m = (K_pph - t + z + L)%L;

            ward += sign * ( ppharray[K_pph + L][k + d*L + m*L2 + t*L3] * ppharray[K_pph + L][k + a*L + m*L2 + e*L3]

               + ppharray[K_pph + L][k + d*L + m*L2 + e*L3] * ppharray[K_pph + L][k + a*L + m*L2 + t*L3] );

            K_pph = (k + c + a_)%L;

            //9) ca;bd he;zt
            m = (K_pph - h + e + L)%L;

            ward += sign * ( ppharray[K_pph + L][k + c*L + m*L2 + h*L3] * ppharray[K_pph + L][k + b*L + m*L2 + z*L3]

               + ppharray[K_pph + L][k + c*L + m*L2 + z*L3] * ppharray[K_pph + L][k + b*L + m*L2 + h*L3] );

            //10) ca;bd hz;et
            m = (K_pph - h + z + L)%L;

            ward += sign * sign_ * ( ppharray[K_pph + L][k + c*L + m*L2 + h*L3] * ppharray[K_pph + L][k + b*L + m*L2 + e*L3]

               + ppharray[K_pph + L][k + c*L + m*L2 + e*L3] * ppharray[K_pph + L][k + b*L + m*L2 + h*L3] );

            //11) ca;bd te;zh
            m = (K_pph - t + e + L)%L;

            ward += sign * sign_ * ( ppharray[K_pph + L][k + c*L + m*L2 + t*L3] * ppharray[K_pph + L][k + b*L + m*L2 + z*L3]

               + ppharray[K_pph + L][k + c*L + m*L2 + z*L3] * ppharray[K_pph + L][k + b*L + m*L2 + t*L3] );

            //12) ca;bd tz;eh
            m = (K_pph - t + z + L)%L;

            ward += sign * ( ppharray[K_pph + L][k + c*L + m*L2 + t*L3] * ppharray[K_pph + L][k + b*L + m*L2 + e*L3]

               + ppharray[K_pph + L][k + c*L + m*L2 + e*L3] * ppharray[K_pph + L][k + b*L + m*L2 + t*L3] );

            K_pph = (k + c + b_)%L;

            //13) cb;ad he;zt
            m = (K_pph - h + e + L)%L;

            ward += ppharray[K_pph + L][k + c*L + m*L2 + h*L3] * ppharray[K_pph + L][k + a*L + m*L2 + z*L3]

               + ppharray[K_pph + L][k + c*L + m*L2 + z*L3] * ppharray[K_pph + L][k + a*L + m*L2 + h*L3];

            //14) cb;ad hz;et
            m = (K_pph - h + z + L)%L;

            ward += sign_ * ( ppharray[K_pph + L][k + c*L + m*L2 + h*L3] * ppharray[K_pph + L][k + a*L + m*L2 + e*L3]

               + ppharray[K_pph + L][k + c*L + m*L2 + e*L3] * ppharray[K_pph + L][k + a*L + m*L2 + h*L3] );

            //15) cb;ad te;zh
            m = (K_pph - t + e + L)%L;

            ward += sign_ * ( ppharray[K_pph + L][k + c*L + m*L2 + t*L3] * ppharray[K_pph + L][k + a*L + m*L2 + z*L3]

               + ppharray[K_pph + L][k + c*L + m*L2 + z*L3] * ppharray[K_pph + L][k + a*L + m*L2 + t*L3] );

            //16) cb;ad tz;eh
            m = (K_pph - t + z + L)%L;

            ward += ppharray[K_pph + L][k + c*L + m*L2 + t*L3] * ppharray[K_pph + L][k + a*L + m*L2 + e*L3]

               + ppharray[K_pph + L][k + c*L + m*L2 + e*L3] * ppharray[K_pph + L][k + a*L + m*L2 + t*L3];

         }//end of loop over k

         (*this)(i,j) += 4.0 * 9.0 * Tools::g6j(0,0,S,1) * Tools::g6j(0,0,S_,1) * ward;

      }
   }

}
