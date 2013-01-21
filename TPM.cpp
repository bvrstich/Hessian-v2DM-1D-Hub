#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;
using std::vector;

#include "include.h"

vector< vector<int> > *TPM::t2s;
int ***TPM::s2t;

double **TPM::norm;

int **TPM::block_char;

/**
 * initialize the static lists
 */
void TPM::init(){

   int L = Tools::gL();

   //allocate s2t
   s2t = new int ** [2*L];

   for(int B = 0;B < 2*L;++B){

      s2t[B] = new int * [L];

      for(int i = 0;i < L;++i)
         s2t[B][i] = new int [L];

   }

   norm = new double * [L];

   for(int a = 0;a < L;++a)
      norm[a] = new double [L];

   //allocate t2s
   t2s = new vector< vector<int> > [2*L];

   //initialize the lists
   int tp;
   
   vector<int> v(2);

   //S = 0
   for(int K = 0;K < L;++K){

      tp = 0;

      for(int a = 0;a < L;++a)
         for(int b = a;b < L;++b){

            if( (a + b)%L == K ){

               v[0] = a;
               v[1] = b;

               t2s[K].push_back(v);

               s2t[K][a][b] = tp;
               s2t[K][b][a] = tp;

               if(a == b)
                  norm[a][b] = 1.0/(std::sqrt(2.0));
               else
                  norm[a][b] = 1.0;

               norm[b][a] = norm[a][b];

               ++tp;

            }

         }

   }

   //S = 1
   for(int K = 0;K < L;++K){

      tp = 0;

      for(int a = 0;a < L;++a)
         for(int b = a + 1;b < L;++b){

            if( (a + b)%L == K ){

               v[0] = a;
               v[1] = b;

               t2s[K + L].push_back(v);

               s2t[K + L][a][b] = tp;
               s2t[K + L][b][a] = tp;

               ++tp;

            }

         }

   }

   //allocate the block_char list
   block_char = new int * [2*L];

   for(int B = 0;B < 2*L;++B)
      block_char[B] = new int [2];

   //initialize:
   for(int K = 0;K < L;++K){

      block_char[K][0] = 0;//S
      block_char[K][1] = K;//K

      block_char[K + L][0] = 1;//S
      block_char[K + L][1] = K;//K

   }

}

/**
 * deallocate the static lists and stuff
 */
void TPM::clear(){

   int L = Tools::gL();
   
   for(int B = 0;B < 2*L;++B){

      for(int a = 0;a < L;++a)
         delete [] s2t[B][a];

      delete [] s2t[B];

      delete [] block_char[B];

   }

   delete [] s2t;
   delete [] t2s;

   delete [] block_char;

   for(int a = 0;a < L;++a)
      delete [] norm[a];

   delete [] norm;

}

/**
 * standard constructor for a spinsymmetrical, translationally invariant tp matrix: constructs BlockMatrix object with 2 (M/2) blocks, (M/2) for S = 0 and (M/2) for S = 1,
 */
TPM::TPM() : BlockMatrix(2*Tools::gL()) {

   //set the dimension and degeneracy of the blocks
   for(int K = 0;K < Tools::gL();++K){

      this->setMatrixDim(K,t2s[K].size(),1);//the S = 0 block
      this->setMatrixDim(K + Tools::gL(),t2s[K + Tools::gL()].size(),3);//the S = 1 block

   }

}

/**
 * copy constructor for a spinsymmetrical, translationally invariant tp matrix: constructs BlockMatrix object with 2 (M/2) blocks, (M/2) for S = 0 and (M/2) for S = 1,
 * @param tpm_c The TPM object to be copied into (*this)
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor: if counter == 1 the memory for the static lists t2s en s2t will be deleted.
 * 
 */
TPM::~TPM(){ }

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param S The spin-index
 * @param a first sp momentum index that forms the tp row index i of block B, together with b
 * @param b second sp momentum index that forms the tp row index i of block B, together with a
 * @param c first sp momentum index that forms the tp column index j of block B, together with d
 * @param d second sp momentum index that forms the tp column index j of block B, together with c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int S,int a,int b,int c,int d) const{

   //momentum checks out
   if( (a + b)%Tools::gL() != (c + d)%Tools::gL())
      return 0;

   int K = (a + b)%Tools::gL();

   int B = K + S*Tools::gL();//blockindex is easily deductable from the S and K

   if(S == 0){

      int i = s2t[B][a][b];
      int j = s2t[B][c][d];

      return (*this)(B,i,j);

   }
   else{

      if( (a == b) || (c == d) )
         return 0;
      else{

         int i = s2t[B][a][b];
         int j = s2t[B][c][d];

         int phase = 1;

         if(a > b)
            phase *= -1;
         if(c > d)
            phase *= -1;

         return phase*(*this)(B,i,j);

      }

   }

}

ostream &operator<<(ostream &output,const TPM &tpm_p){

   int S,K;

   for(int B = 0;B < tpm_p.gnr();++B){

      S = tpm_p.block_char[B][0];
      K = tpm_p.block_char[B][1];

      output << "S =\t" << S << "\tK =\t" << K << "\tdimension =\t" << tpm_p.gdim(B) << "\tdegeneracy =\t" << tpm_p.gdeg(B) << std::endl;
      output << std::endl;

      for(int i = 0;i < tpm_p.gdim(B);++i)
         for(int j = 0;j < tpm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << tpm_p.t2s[B][i][0] << "\t" << tpm_p.t2s[B][i][1]

               << "\t" << tpm_p.t2s[B][j][0] << "\t" << tpm_p.t2s[B][j][1] << "\t" << tpm_p(B,i,j) << endl;

         }

      output << std::endl;

   }

   return output;

}

/**
 * construct the spinsymmetrical hubbard hamiltonian in momentum space with on site repulsion U
 * @param U onsite repulsion term
 */
void TPM::hubbard(double U){

   int a,b,c,d;//sp momentum 

   double ward = 1.0/(Tools::gN() - 1.0);

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            //init
            (*this)(B,i,j) = 0;

            //hopping (kinetic energy):
            if(i == j)
               (*this)(B,i,i) = -2.0 * ward * ( cos( 4.0 * a * 3.141592653589793238462 / ( 2.0*Tools::gL() ) )  + cos( 4.0 * b * 3.141592653589793238462 / (2.0*Tools::gL()) ) );

            //on-site repulsion
            if(S == 0)
               (*this)(B,i,j) += norm[a][b] * norm[c][d] * 4.0*U / (2.0*Tools::gL());

         }
      }

   }

   this->symmetrize();

}

/**
 * @return The expectation value of the total spin for the TPM.
 */
double TPM::S_2() const{

   double ward = 0.0;

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      if(S == 0){

         for(int i = 0;i < gdim(B);++i)
            ward += -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) * (*this)(B,i,i);

      }
      else{

         for(int i = 0;i < gdim(B);++i)
            ward += 3.0 * ( -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) + 2.0 ) * (*this)(B,i,i);

      }

   }

   return ward;

}

/**
 * initialize this onto the unitmatrix with trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = Tools::gN()*(Tools::gN() - 1.0)/(Tools::gM()*(Tools::gM() - 1.0));

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         (*this)(B,i,i) = ward;

         for(int j = i + 1;j < gdim(B);++j)
            (*this)(B,i,j) = (*this)(B,j,i) = 0.0;

      }
   }

}

/**
 * access to the TPM norms from outside the class
 */
double TPM::gnorm(int a,int b){

   return norm[a][b];

}

/**
 * access to the lists from outside of the class
 */
int TPM::gt2s(int B,int i,int option){

   return t2s[B][i][option];

}

/**
 * access to the lists from outside of the class
 */
int TPM::gs2t(int block,int a,int b){

   return s2t[block][a][b];

}

/**
 * access to the lists from outside of the class
 */
int TPM::gblock_char(int block,int option){

   return block_char[block][option];

}

/**
 * convert a 2DM object from vector to matrix/TPM form
 * @param grad input Gradient object
 */
void TPM::convert(const Gradient &grad){

   int tpmm;
   int S;

   for(int B = 0;B < 2*Tools::gL();++B){

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i)
         for(int j = i;j < gdim(B);++j){

            tpmm = TPTPM::gt2tpmm(B,i,j);

            (*this)(B,i,j) = grad[tpmm]/ ( 2.0 * Gradient::gnorm(tpmm) * (2.0 * S + 1.0) );

         }

   }

   this->symmetrize();

}

/**
 * @return the dimension of the block corresponding to index B
 */
int TPM::gdim(int B) {

   return t2s[B].size();

}

/**
 * perform a line search what step size in along the Newton direction is ideal.
 * @param t potential scaling factor
 * @param P SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,SUP &P,const TPM &ham){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   P.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta;

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp;

   hulp.L_map(P,S_delta);

   EIG eigen(hulp);

   double a = 0;

   //interval where the bissection method is going to search in: [a=0,b]
   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * perform a line search what step size in along the Newton direction is ideal, this one is used for extrapolation.
 * @param t potential scaling factor
 * @param rdm TPM containing the current approximation of the rdm.
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double TPM::line_search(double t,const TPM &rdm,const TPM &ham){

   SUP P;

   P.fill(rdm);

   P.invert();

   return this->line_search(t,P,ham);

}
