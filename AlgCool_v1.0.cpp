#include <vector>

#include <math.h>
#include <iomanip>
#include <algorithm>

#include "mwArray.h"
#include "MPS.h"
#include "Contractor.h"

#define TMPSTORE 0 // Using temporary storage for intermediate results

#include "IsingHamiltonian.h"
#include "HeisenbergHamiltonian.h"
#include "SpinMPO.h"

#define forn(i,n) for(int i = 0; i < n; ++i)
#define forv(i,v) for(int i = 0; i < v.size(); ++i)

using namespace shrt;

/** \file AlgCool.cpp
 
    testIsing runs the findGroundState routine with the MPO for the
    Ising Hamiltonian \ref <IsingHamiltonian.cpp>, to check convergence
    and then runs imaginary time evolution to check also the
    implementation of the unitary evolution.

    Receives arguments:
    \param <L> (int) length of the chain
    \param <J> (double) parameter \f$J\f$ of the Hamiltonian
    \param <g> (double) parameter \f$g\f$ of the Hamiltonian
    \param <h> (double) parameter \f$h\f$ of the Hamiltonian
    \param <D> (int) maximum bond dimension
    \param <M> (int) maximum number of imaginary time steps 
           (if 0, no imaginary time phase is run)
    \param <delta> (double) step width for imaginary time evolution
    \param <outfname> (char*) name of the output file for the energy
    \param <app> (int) whether the output file is to be kept (app==1)
*/

//Computes exp(i*h*t), where h is given in terms of x*sigma_x + y sigma_y + z sigma_z;
mwArray exp_h(double t, double x, double y, double z){
//complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};   complex_t datay[]={ZERO_c,-I_c,I_c,ZERO_c};   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c}; complex_t data0[]={ONE_c,ZERO_c,ZERO_c,ONE_c};
complex_t data[] ={cos(t)*ONE_c + I_c*z*sin(t), (I_c*x+y*ONE_c)*sin(t), (I_c*x-y*ONE_c)*sin(t),cos(t)*ONE_c - I_c*z*sin(t)};
return mwArray(Indices(2,2), data);
}


//BASIC ARITHMETIC FUNCTIONS
double norm(complex_t z){
 complex_t zz = z*conjugate(z);
 return sqrt(real(zz));
}

//Computes U'*H*U
vector<vector<complex_t> > Conjuga(vector<vector<complex_t> >&H, vector<vector<complex_t> >&U){
	vector<complex_t> zeros(U.size(), ZERO_c);	
	vector<vector<complex_t> > R(U.size(), zeros);
	forv(a,U){
		forv(e,U){
			forv(b,U[0]){
				forv(d,U[0]){
					R[a][e] = R[a][e]+conjugate(U[b][a])*H[b][d]*U[d][e];
				}
			}
		}
	}
	return R;
}

//Computes the trace of A*B
complex_t TrAB(vector<vector<complex_t> >&A, vector<vector<complex_t> >&B){
	complex_t res = ZERO_c;
	forv(i,A){
		forv(j,A[0]){
			res = res + A[i][j]*B[j][i];
		}
	}
	return res;
}


//Computes the Kronecker Product of two matrices
vector<vector<complex_t> > kron(vector<vector<complex_t> >& v1, vector<vector<complex_t> >& v2){
int n1 = v1.size(); int m1 = v1[0].size(); int n2 = v2.size(); int m2 = v2[0].size();
int n = n1*n2, m = m1*m2;
vector<complex_t> zerorow(m,ZERO_c);
vector<vector<complex_t> > res(n,zerorow);
	forn(i1,n1){
		forn(i2,n2){
			forn(j1,m1){
				forn(j2,m2){
					res[i1*n2+i2][j1*m2+j2] = v1[i1][j1]*v2[i2][j2];
				}
			}
		}
	}
return res;
}

//Computes the subtraction of two matrices
vector<vector<complex_t> > matrix_subtract(vector<vector<complex_t> >& v1, vector<vector<complex_t> >& v2){
int n1 = v1.size(); int m1 = v1[0].size(); int n2 = v2.size(); int m2 = v2[0].size();
if ((n1!= n2)||(m2!=n2)){
	cerr << "ERROR. Matrices have different size" << n1 << "x" << m1 << " minus " << n2 << "x" << m2 << endl;
	int xx; cin >> xx;
}
int n = n1, m = m1;
vector<complex_t> zerorow(m,ZERO_c);
vector<vector<complex_t> > res(n,zerorow);
	forn(i1, n){
		forn(j1, m){
				res[i1][j1] = v1[i1][j1]-v2[i1][j1];
		}
	}
return res;
}

//Computes the normal product of two matrices
vector<vector<complex_t> > matrix_multiply(vector<vector<complex_t> >& v1, vector<vector<complex_t> >& v2){
int n1 = v1.size(); int m1 = v1[0].size(); int n2 = v2.size(); int m2 = v2[0].size();
if (m1!= n2){
	cerr << "ERROR. Matrices have different size" << n1 << "x" << m1 << " times " << n2 << "x" << m2 << endl;
	int xx; cin >> xx;
}
int n = n1, m = m2;
vector<complex_t> zerorow(m,ZERO_c);
vector<vector<complex_t> > res(n,zerorow);
	forn(i1, n){
		forn(j1, m){
			forn(k, m1){
				res[i1][j1] = res[i1][j1]+v1[i1][k]*v2[k][j1];
			}
		}
	}
return res;
}

static vector<vector<complex_t> > Id2(){
		vector<complex_t> row(2,ZERO_c);
		vector<vector<complex_t> > res(2,row);
		res[0][0] = ONE_c;
		res[1][1] = ONE_c;
		return res;
}

static vector<vector<complex_t> > SigmaX(){
		vector<complex_t> row(2,ZERO_c);
		vector<vector<complex_t> > res(2,row);
		res[0][1] = ONE_c;
		res[1][0] = ONE_c;
		return res;
}
static vector<vector<complex_t> > SigmaY(){
		vector<complex_t> row(2,ZERO_c);
		vector<vector<complex_t> > res(2,row);
		res[0][1] = -I_c;
		res[1][0] = I_c;
		return res;
}
static vector<vector<complex_t> > SigmaZ(){
		vector<complex_t> row(2,ZERO_c);
		vector<vector<complex_t> > res(2,row);
		res[0][0] = ONE_c;
		res[1][1] = -ONE_c;
		return res;
}

//Prepares a set of small h's that are available
map<int, vector<vector<complex_t> > > prepare_h_bag(int size, int &c){
map<int, vector<vector<complex_t> > > res; c = 0;
	vector<vector<complex_t> > SX = SigmaX();
	vector<vector<complex_t> > SY = SigmaY();
	vector<vector<complex_t> > SZ = SigmaZ();
	vector<vector<complex_t> > I2 = Id2();
	if (size==1){
		res[c++] = SX;
		res[c++] = SY;
		res[c++] = SZ;
	}else if (size==2){
		res[c++] = kron(SX,SX);
		res[c++] = kron(SX,SY);
		res[c++] = kron(SX,SZ);
		res[c++] = kron(SY,SX);
		res[c++] = kron(SY,SY);
		res[c++] = kron(SY,SZ);
		res[c++] = kron(SZ,SX);
		res[c++] = kron(SZ,SY);
		res[c++] = kron(SZ,SZ);
		res[c++] = kron(I2,SX);
		res[c++] = kron(I2,SY);
		res[c++] = kron(I2,SZ);
		res[c++] = kron(SX,I2);
		res[c++] = kron(SY,I2);
		res[c++] = kron(SZ,I2);
	}
return res;
}




//AUX FUNCTIONS
int max(vector<int> & v){
	int res = v[0];
	forv(i,v){ if (v[i]>res) res=v[i];}
	return res;
}
int total(vector<int> & v){
	int res = 0;
	forv(i,v) res+=v[i];
	return res;
}
int left_nonzero(vector<int> & v){
	int res = -1;
	forv(i,v){
		if (v[i]==1){
			res = i; break;
		}
	}
	return res;
}
int right_nonzero(vector<int> & v){
	int res = -1;
	for(int i = v.size()-1; i >= 0; --i){
		if (v[i]==1){
			res = i; break;
		}
	}
	return res;
}

//Computes the AND of two vectors of assumed equal length
vector<int> intersect_masks(vector<int>&v1, vector<int>&v2){
	vector<int> res = v1;
	forv(i,v1){
		if (res[i]!=v2[i]) res[i]=0;
	}
	return res;
}
//Computes the OR of two vectors of assumed equal length
vector<int> union_masks(vector<int>&v1, vector<int>&v2){
	vector<int> res = v1;
	forv(i,v1){
		if (v2[i]==1) res[i]=1;
	}
	return res;
}
//Sorts and removes duplicates, if any
void remove_duplicates(vector<int> & vec){
	sort( vec.begin(), vec.end() );
	vec.erase( unique( vec.begin(), vec.end() ), vec.end() );
}



//DISPLAY FUNCTIONS
void print_vector(vector<int> v){
	forv(i,v) cout << v[i];
	cout << endl;
}


void print_matrix(vector<vector<complex_t> > M){
	forv(i,M){
		forv(j,M[i]) cout << M[i][j] << " ";
		cout << endl;
	} cout << endl;
}




//ALGORITHM-RELATED FUNCTIONS
vector<int> generateLayer(int n, int shifted){
	vector<int> res(n,0); int c = 0;
	forn(i,n){
		res[i] = c;
		if (shifted%2){
			if ((i+2)/2 > (i+1)/2) ++c;
		}else{
			if ((i+1)/2 > i/2) ++c;
		}
	}
	return res;
}
map<int, pair<vector<int>, vector<vector<complex_t > > > > initUnitaries(vector<int>& layer){
	map<int, pair<vector<int>, vector<vector<complex_t > > > > U;
	int L = layer.size();
	forn(i,max(layer)+1){
		pair<vector<int>, vector<vector<complex_t > > > p;
		vector<int> v(L,0);
		forv(j,layer){
			if (layer[j]==i){
				v[j]=1;
			}
		}
		print_vector(v);
		//Start preparing identities
		vector<vector<complex_t > > I2 = Id2(), Ui;
		vector<complex_t> uno(1,ONE_c); Ui.push_back(uno);
		forn(j,total(v)){
			Ui = kron(Ui, I2);
		}
		p.first=v; p.second=Ui;
		U[i]=p;
	}
	return U;
}


map<int, pair<vector<int>, vector<vector<complex_t > > > > representHamiltonian(int L, vector<double>& Jx, vector<double>& Jy, vector<double>& Jz, vector<double>& Jxx, vector<double>& Jyy, vector<double>& Jzz){
	map<int, pair<vector<int>, vector<vector<complex_t > > > > res;
	int c = 0; vector<int> zeros(L,0);
	//Pauli Matrices definition
  //complex_t datax[]={ZERO_c,ONE_c,ONE_c,ZERO_c};   complex_t datay[]={ZERO_c,-I_c,I_c,ZERO_c};   complex_t dataz[]={ONE_c,ZERO_c,ZERO_c,-ONE_c};
	forn(i,L){
		pair<vector<int>, vector<vector<complex_t > > > p;
		vector<int> v = zeros; v[i]=1;
		vector<complex_t> row(2,ZERO_c);
		vector<vector<complex_t> > H(2,row);
		H[0][0] = Jz[i]*ONE_c; H[0][1] = Jx[i]*ONE_c -I_c*Jy[i];
		H[1][0] = Jx[i]*ONE_c +I_c*Jy[i]; H[1][1] = -Jz[i]*ONE_c;
		p.first = v; p.second=H;
		res[c++] = p;
		print_vector(v);
		print_matrix(H);
	}
	//Two-body
	forn(i,L-1){
		pair<vector<int>, vector<vector<complex_t > > > p;
		vector<int> v = zeros; v[i]=1; v[i+1]=1;
		vector<complex_t> row(4,ZERO_c);
		vector<vector<complex_t> > H(4,row);
		H[0][0] = Jzz[i]*ONE_c; H[1][1] = -Jzz[i]*ONE_c; H[2][2] = -Jzz[i]*ONE_c; H[3][3] = Jzz[i]*ONE_c;
		H[0][3] = (Jxx[i]-Jyy[i])*ONE_c; H[1][2] = (Jxx[i]+Jyy[i])*ONE_c; H[2][1] = (Jxx[i]+Jyy[i])*ONE_c; H[3][0] = (Jxx[i]-Jyy[i])*ONE_c;
		p.first = v; p.second=H;
		res[c++] = p;
		print_vector(v);
		print_matrix(H);
	}
	return res;
}



//Takes the support of a unitary operator and finds all terms in the Hamiltonian that overlap with it.
//Then flags all the terms of the Hamiltonian that are needed to update the energy value when we change that unitary
vector<int> relevant_Hamilt_terms_given_U(vector<int>& sup_U, map<int, pair<vector<int>, vector<vector<complex_t > > > >& Hmap){
	vector<int> res;
	for (map<int, pair<vector<int>, vector<vector<complex_t > > > >::iterator it = Hmap.begin(); it != Hmap.end(); ++it){
		vector<int> intersection = intersect_masks(sup_U, it->second.first);
		if (total(intersection) > 0){
			cout << "Overlap with " << it->first << endl;
			//res = union_masks(res, it->second.first);
			res.push_back(it->first);
		}
	}
	return res;
}

map<int, vector<int> > local_terms_for_Unitaries(map<int, pair<vector<int>, vector<vector<complex_t > > > >& Hmap, map<int, pair<vector<int>, vector<vector<complex_t > > > >& U){
	map<int, vector<int> > res;
	for (map<int, pair<vector<int>, vector<vector<complex_t > > > >::iterator it = U.begin(); it != U.end(); ++it){
		res[it->first] = relevant_Hamilt_terms_given_U(it->second.first, Hmap);
	}
	return res;
}











//For each term forming the Hamiltonian returns the unitaries that act on it, directly or indirectly.
map<int, vector<int> > unitaries_for_Hlocal(map<int, pair<vector<int>, vector<vector<complex_t > > > >& Hmap, vector<int>& layer){
	map<int, vector<int> > res;
	for(map<int, pair<vector<int>, vector<vector<complex_t > > > >::iterator it = Hmap.begin(); it != Hmap.end(); ++it){
		cout << it->first << " " << it->second.first << endl;
		vector<int> v;
		forv(i,layer){
			if (it->second.first[i]==1){
				v.push_back(layer[i]);
			}
		}
		remove_duplicates(v);
		res[it->first] = v;
		cout << v << endl;
	}
	return res;
}




map<int, vector<int> > hamilt_to_masks(map<int, pair<vector<int>, vector<vector<complex_t > > > >& Hmap, vector<int>& layer, map<int, vector<int> >& H_to_U){
	map<int, vector<int> > res;
	for(map<int, pair<vector<int>, vector<vector<complex_t > > > >::iterator it = Hmap.begin(); it != Hmap.end(); ++it){
	//	vector<int> mask(it->second.first.size(), 0); //We should not initialize the mask like this. There might be unitaries that do not overlap with any site for instance, then we have a problem. Initialize with the local Hamiltonian support.
		vector<int> mask = it->second.first;
		forv(i, H_to_U[it->first]){
			vector<int> tmp(mask.size(), 0);
			forv(j, mask){
				if (layer[j]==H_to_U[it->first][i]){
					tmp[j]=1;
				}
			}
			mask = union_masks(mask, tmp);
		}
		res[it->first] = mask;
	}
	return res;
}



//Auxiliary function that swaps the order of the last K bits of an integer, to permute the parties in the reduced state... (Tell MC about this)
inline int reverse_bits(int a, int k){
	int res = 0, c = 1;
	forn(i,k){
		res += c*((a/(1<<(k-i-1)))%2);
		c=c*2;
	}
//	return a;//REMOVE ME PLEASE!!!
	return res;
}

//Returns the reduced density matrix of the current MPS at the sites given by the mask, using the contractor.
vector<vector<complex_t> > get_rdms(MPS& current, Contractor& contractor, vector<int>& mask){
	int k = total(mask);
	int inici = -1; int fi = -1;
	forv(i,mask){
		if ((inici==-1)&&mask[i]==1) inici = i;
		if ((fi==-1)&&(inici!=-1)&&mask[i]==0) fi = i-1;
	}
	if (fi==-1) fi = mask.size()-1;
	cout << inici << " to " << fi << "for mask " << mask << endl;
	mwArray rdm = contractor.getRDM(current, inici, fi);//WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!! LEAST SIGNIFICANT DIGIT IS RETURNED ON THE LEFT!!!!!!!!
//	putForMatlab(cout, rdm, "A");
	vector<complex_t>	zeros(rdm.getDimension(1), ZERO_c);
	vector<vector<complex_t> > res(rdm.getDimension(0),zeros);
	forn(a,rdm.getDimension(0)){
		forn(b,rdm.getDimension(1)){
			res[reverse_bits(a,k)][reverse_bits(b,k)] = rdm.getElement(Indices(a,b));
		}
	}
	return res;
}
map<vector<int>, vector<vector<complex_t> > > preset_rdms(MPS& current, Contractor& contractor, map<int,vector<int> >& H_to_maskRho){
	map<vector<int>, vector<vector<complex_t> > > res;
	for(map<int, vector<int> >::iterator it = H_to_maskRho.begin(); it != H_to_maskRho.end(); ++it){
		res[it->second] = get_rdms(current, contractor, it->second);
//		cout << "it->second is" << it->second << endl;
//		print_matrix(res[it->second]);
	}
	return res;
}


//We assume masks come consecutive. Otherwise it would require a more subtle approach with permutations of indices etc.
vector<vector<complex_t > > accomodate_in_mask(vector<int>& mask_master, vector<vector<complex_t > >& op, vector<int>& mask_op){
//	cout << "trying to accomodate " << mask_op << " in " << mask_master << endl;
//	print_matrix(op);
	vector<complex_t> uno(1,ONE_c); vector<vector<complex_t > > res(1,uno);
	vector<vector<complex_t > > I2 = Id2();
	bool flag = false;
	forv(i, mask_master){
		if (mask_master[i]==1){
			if (mask_op[i]==0){
				res = kron(res, I2);
			}else{
				if (flag == false){
					res = kron(res, op);
					flag = true;
				}
			}
		}
	}
//	print_matrix(res);
//	int xx; cin >> xx;
	return res;
}


double total_energy(map<int, vector<int> >& U_to_H, map<int, vector<int> >& H_to_U, map<int, vector<int> >& H_to_maskRho, map<vector<int> ,vector<vector<complex_t> > >& mask_to_Rho, map<int, pair<vector<int>, vector<vector<complex_t > > > >& Hmap, map<int, pair<vector<int>, vector<vector<complex_t > > > >& U){
	double E=0.;
	//For every H_i in the Hamiltonian:
	for(map<int, pair<vector<int>, vector<vector<complex_t > > > >::iterator it = Hmap.begin(); it != Hmap.end(); ++it){
		vector<int> mask = H_to_maskRho[it->first]; //This is the mask on which we do the computation
		vector<vector<complex_t> > H_ext = accomodate_in_mask(mask, it->second.second, it->second.first); //Represent the reduced state there
		//Unitary preparation
		vector<vector<complex_t> > I2 = Id2(), Unitary = I2;
		forn(j, total(mask)-1){Unitary = kron(Unitary, I2);}
		//We now find all the unitaries, and multiply them, if there are more than one.
		forv(j, H_to_U[it->first]){//For every U_j interacting with H_i
			vector<vector<complex_t> > U_ext = accomodate_in_mask(mask, U[H_to_U[it->first][j]].second, U[H_to_U[it->first][j]].first);
			Unitary = matrix_multiply(Unitary, U_ext);//Commutes anyway because of layout
		}
		vector<vector<complex_t> > H_ext2 = Conjuga(H_ext, Unitary);
		double contrib = real(TrAB(mask_to_Rho[mask], H_ext2));
		E=E+contrib;
		cout <<it->first << mask << "  " << contrib << endl;
	}
	return E;
}
double energy_contrib(int Uindex, map<int, vector<int> >& U_to_H, map<int, vector<int> >& H_to_U, map<int, vector<int> >& H_to_maskRho, map<vector<int> ,vector<vector<complex_t> > >& mask_to_Rho, map<int, pair<vector<int>, vector<vector<complex_t > > > >& Hmap, map<int, pair<vector<int>, vector<vector<complex_t > > > >& U){
	double E=0.;
	forv(i, U_to_H[Uindex]){//For every H_i
		vector<int> mask = H_to_maskRho[U_to_H[Uindex][i]]; //This is the mask on which we do the computation
		/*cout << H_to_maskRho[U_to_H[Uindex][i]] << endl;
		print_matrix(mask_to_Rho[H_to_maskRho[U_to_H[Uindex][i]]]);*/
		vector<vector<complex_t> > H_ext = accomodate_in_mask(mask, Hmap[U_to_H[Uindex][i]].second, Hmap[U_to_H[Uindex][i]].first);
		vector<vector<complex_t> > I2 = Id2(), Unitary = I2;
		forn(j, total(mask)-1){Unitary = kron(Unitary, I2);}
		//We now find all the unitaries, and multiply them, if there are more than one.
		forv(j, H_to_U[U_to_H[Uindex][i]]){//For every U_j interacting with H_i
			vector<vector<complex_t> > U_ext = accomodate_in_mask(mask, U[H_to_U[U_to_H[Uindex][i]][j]].second, U[H_to_U[U_to_H[Uindex][i]][j]].first);
			Unitary = matrix_multiply(Unitary, U_ext);//Commutes anyway because of layout
		}
		vector<vector<complex_t> > H_ext2 = Conjuga(H_ext, Unitary);
		E=E+real(TrAB(mask_to_Rho[H_to_maskRho[U_to_H[Uindex][i]]], H_ext2));
	}
	return E;
}

//Optimal difference in energy, given h
double energy_delta_h(int Uindex, map<int, vector<int> >& U_to_H, map<int, vector<int> >& H_to_U, map<int, vector<int> >& H_to_maskRho, map<vector<int> ,vector<vector<complex_t> > >& mask_to_Rho, map<int, pair<vector<int>, vector<vector<complex_t > > > >& Hmap, map<int, pair<vector<int>, vector<vector<complex_t > > > >& U, vector<vector<complex_t> >& h, double & t){
	double E=0.; t=0.;
	complex_t A = ZERO_c, B = ZERO_c, C = ZERO_c;
	forv(i, U_to_H[Uindex]){//For every H_i
		vector<int> mask = H_to_maskRho[U_to_H[Uindex][i]]; //This is the mask on which we do the computation
		/*cout << H_to_maskRho[U_to_H[Uindex][i]] << endl;
		print_matrix(mask_to_Rho[H_to_maskRho[U_to_H[Uindex][i]]]);*/
		vector<vector<complex_t> > H_ext = accomodate_in_mask(mask, Hmap[U_to_H[Uindex][i]].second, Hmap[U_to_H[Uindex][i]].first);
		vector<vector<complex_t> > h_ext = accomodate_in_mask(mask, h, U[Uindex].first); //h acts on the same support as the U we are modifying
		vector<vector<complex_t> > I2 = Id2(), Unitary = I2;
		forn(j, total(mask)-1){Unitary = kron(Unitary, I2);}
		//We now find all the unitaries, and multiply them, if there are more than one.
		forv(j, H_to_U[U_to_H[Uindex][i]]){//For every U_j interacting with H_i
			vector<vector<complex_t> > U_ext = accomodate_in_mask(mask, U[H_to_U[U_to_H[Uindex][i]][j]].second, U[H_to_U[U_to_H[Uindex][i]][j]].first);
			Unitary = matrix_multiply(Unitary, U_ext);//Commutes anyway because of layout
		}
		vector<vector<complex_t> > UHU = Conjuga(H_ext, Unitary);
		vector<vector<complex_t> > tmphH = matrix_multiply(h_ext,H_ext);
		vector<vector<complex_t> > tmpHh = matrix_multiply(H_ext,h_ext);
		vector<vector<complex_t> > comm = matrix_subtract(tmphH, tmpHh);
		vector<vector<complex_t> > tmphHh = Conjuga(H_ext,h_ext);
		vector<vector<complex_t> > UcommU = Conjuga(comm, Unitary);
		vector<vector<complex_t> > UhHhU = Conjuga(tmphHh,Unitary);
		A=A+TrAB(mask_to_Rho[H_to_maskRho[U_to_H[Uindex][i]]], UHU);
		B=B+I_c*TrAB(mask_to_Rho[H_to_maskRho[U_to_H[Uindex][i]]], UcommU);
		C=C+TrAB(mask_to_Rho[H_to_maskRho[U_to_H[Uindex][i]]], UhHhU);
/*		print_matrix(comm);
		print_matrix(UcommU);
		print_matrix(h);
		int xx; cin >> xx;*/
	}
	cout << "A = " << A << "; B = " << B << "; C = " << C << endl;
	double E1 = (real(A) + real(C) + sqrt(real(B)*real(B) + (real(A-C)*real(A-C))))/2;
	double E2 = (real(A) + real(C) - sqrt(real(B)*real(B) + (real(A-C)*real(A-C))))/2;
	cout << E1 << " " << E2 << endl;
	//We also need to compute the optimal t, to later perform the unitary
//	t = .5*atan2(-real(B),real(C-A));
	t = .5*atan2(real(B),real(C-A)); // For some reason, this is the correct t; check calculations again!!
	cout << "t = " << t << endl;
	E = real(A*cos(t)*cos(t)+B*cos(t)*sin(t)+C*sin(t)*sin(t));
	cout << E << endl;
	return E;
}


//returns e^(iht), which is cos(t) + ih sin(t), because h^2 = 1;
vector<vector<complex_t> > exponential_h(vector<vector<complex_t> >& h, double t){
	vector<vector<complex_t> > res = h;
	forv(i, res){
		forv(j, res[0]){
			if (i==j){
				res[i][j] = cos(t)*ONE_c + I_c * h[i][j]*sin(t);
			}else{
				res[i][j] = I_c * h[i][j]*sin(t);
			}
		}
	}
	return res;
}



vector<vector<complex_t> > reconstruct_unitary_from_basis(vector<complex_t>& c){
	static vector<vector<complex_t> > SX = SigmaX();
	static vector<vector<complex_t> > SY = SigmaY();
	static vector<vector<complex_t> > SZ = SigmaZ();
	static vector<vector<complex_t> > I2 = Id2();
	static vector<vector<complex_t> > I2I2 = kron(I2,I2);
	static vector<vector<complex_t> > I2SX = kron(I2,SX);
	static vector<vector<complex_t> > I2SY = kron(I2,SY);
	static vector<vector<complex_t> > I2SZ = kron(I2,SZ);
	static vector<vector<complex_t> > SXI2 = kron(SX,I2);
	static vector<vector<complex_t> > SXSX = kron(SX,SX);
	static vector<vector<complex_t> > SXSY = kron(SX,SY);
	static vector<vector<complex_t> > SXSZ = kron(SX,SZ);
	static vector<vector<complex_t> > SYI2 = kron(SY,I2);
	static vector<vector<complex_t> > SYSX = kron(SY,SX);
	static vector<vector<complex_t> > SYSY = kron(SY,SY);
	static vector<vector<complex_t> > SYSZ = kron(SY,SZ);
	static vector<vector<complex_t> > SZI2 = kron(SZ,I2);
	static vector<vector<complex_t> > SZSX = kron(SZ,SX);
	static vector<vector<complex_t> > SZSY = kron(SZ,SY);
	static vector<vector<complex_t> > SZSZ = kron(SZ,SZ);
	static double sq2 = sqrt(2.);
	vector<vector<complex_t> > U;
	if (c.size() == 4){
		vector<complex_t> zeros(2,ZERO_c); forn(i,2){U.push_back(zeros);}
		forv(i,U){
			forv(j,U[0]){
				U[i][j] = (I2[i][j]*c[0] + SX[i][j]*c[1] + SY[i][j]*c[2] + SZ[i][j]*c[3])/sqrt(2.);
			}
		}
	}else if (c.size() == 16){
		vector<complex_t> zeros(4,ZERO_c); forn(i,4){U.push_back(zeros);}
		forv(i,U){
			forv(j,U[0]){
				U[i][j] = (I2I2[i][j]*c[0] + I2SX[i][j]*c[1] + I2SY[i][j]*c[2] + I2SZ[i][j]*c[3] + SXI2[i][j]*c[4] + SXSX[i][j]*c[5] + SXSY[i][j]*c[6] + SXSZ[i][j]*c[7] + SYI2[i][j]*c[8] + SYSX[i][j]*c[9] + SYSY[i][j]*c[10] + SYSZ[i][j]*c[11] + SZI2[i][j]*c[12] + SZSX[i][j]*c[13] + SZSY[i][j]*c[14] + SZSZ[i][j]*c[15])/2.;
			}
		}
	}
	return U;
}

vector<complex_t> expand_basis(vector<vector<complex_t> >& U){
	static vector<vector<complex_t> > SX = SigmaX();
	static vector<vector<complex_t> > SY = SigmaY();
	static vector<vector<complex_t> > SZ = SigmaZ();
	static vector<vector<complex_t> > I2 = Id2();
	static vector<vector<complex_t> > I2I2 = kron(I2,I2);
	static vector<vector<complex_t> > I2SX = kron(I2,SX);
	static vector<vector<complex_t> > I2SY = kron(I2,SY);
	static vector<vector<complex_t> > I2SZ = kron(I2,SZ);
	static vector<vector<complex_t> > SXI2 = kron(SX,I2);
	static vector<vector<complex_t> > SXSX = kron(SX,SX);
	static vector<vector<complex_t> > SXSY = kron(SX,SY);
	static vector<vector<complex_t> > SXSZ = kron(SX,SZ);
	static vector<vector<complex_t> > SYI2 = kron(SY,I2);
	static vector<vector<complex_t> > SYSX = kron(SY,SX);
	static vector<vector<complex_t> > SYSY = kron(SY,SY);
	static vector<vector<complex_t> > SYSZ = kron(SY,SZ);
	static vector<vector<complex_t> > SZI2 = kron(SZ,I2);
	static vector<vector<complex_t> > SZSX = kron(SZ,SX);
	static vector<vector<complex_t> > SZSY = kron(SZ,SY);
	static vector<vector<complex_t> > SZSZ = kron(SZ,SZ);
	static double sq2 = sqrt(2.);
	vector<complex_t> c;
	if (U.size()==2){
		c.push_back(TrAB(I2,U)/sq2);
		c.push_back(TrAB(SX,U)/sq2);
		c.push_back(TrAB(SY,U)/sq2);
		c.push_back(TrAB(SZ,U)/sq2);
	}else if(U.size() == 4){
		c.push_back(TrAB(I2I2,U)/2);
		c.push_back(TrAB(I2SX,U)/2);
		c.push_back(TrAB(I2SY,U)/2);
		c.push_back(TrAB(I2SZ,U)/2);
		c.push_back(TrAB(SXI2,U)/2);
		c.push_back(TrAB(SXSX,U)/2);
		c.push_back(TrAB(SXSY,U)/2);
		c.push_back(TrAB(SXSZ,U)/2);
		c.push_back(TrAB(SYI2,U)/2);
		c.push_back(TrAB(SYSX,U)/2);
		c.push_back(TrAB(SYSY,U)/2);
		c.push_back(TrAB(SYSZ,U)/2);
		c.push_back(TrAB(SZI2,U)/2);
		c.push_back(TrAB(SZSX,U)/2);
		c.push_back(TrAB(SZSY,U)/2);
		c.push_back(TrAB(SZSZ,U)/2);
	}
/*	print_matrix(U);
	print_matrix(reconstruct_unitary_from_basis(c));
	cout << "Does it check?" << endl;
	int xx; cin >> xx;*/
	return c;
}



MPS update_MPS(MPS current, Contractor& contractor, map<int, pair<vector<int>, vector<vector<complex_t > > > >& U){
	MPS res = current;
	map<int, vector<vector<complex_t> > > Pauli; Pauli[0] = Id2(); Pauli[1] = SigmaX(); Pauli[2] = SigmaY(); Pauli[3] = SigmaZ();
	for(map<int, pair<vector<int>, vector<vector<complex_t > > > >::iterator it = U.begin(); it!= U.end(); ++it){//We apply each unitary to the MPS. Remember that all these unitaries commute as they act on different qubits
		//We need to convert to mwArray for the library to process well.
		if (total(it->second.first) == 1){
			mwArray M = mwArray(Indices(2,2));
			forn(i,2){
				forn(j,2){
					M.setElement(it->second.second[i][j],Indices(i,j));
				}
			}
			res.applyLocalOperator(left_nonzero(it->second.first),M,true);
		}else if(total(it->second.first)==2){
			MPO myMPO(res.getLength());
			mwArray M0 = mwArray(Indices(2,1,2,16));
			mwArray M1 = mwArray(Indices(2,16,2,1));
			vector<complex_t> base = expand_basis(it->second.second);
			forn(alpha,16){
				forn(x,2){forn(y,2){
					M0.setElement(ONE_c*base[alpha]*Pauli[alpha/4][x][y]/2., Indices(x,0,y,alpha));//Continue HERE!!!
					M1.setElement(ONE_c*Pauli[alpha%4][x][y], Indices(x,alpha,y,0));
				}}
			}
			Operator IdOp(reshape(identityMatrix(2),Indices(2,1,2,1)));
			forn(ii,res.getLength()){
				myMPO.setOp(ii, &IdOp, false);
			}
			myMPO.setOp(left_nonzero(it->second.first), new Operator(M0), true);
			myMPO.setOp(right_nonzero(it->second.first), new Operator(M1), true);
			MPS tmp = res;
			MPS tmp_evol(tmp);
			contractor.optimize(myMPO, tmp, tmp_evol);
			res = tmp_evol;
		}
	}
	return res;
}




int main(int argc, const char* argv[]){
	srandom(time(NULL));
	//We begin by opening the input file, which should be input.txt
	ifstream* in;
	in = new ifstream("input.txt");
	//Initialization Contractor
	Contractor& contractor=Contractor::theContractor(); cout<<"Initialized Contractor"<<endl;
	int L, D, d; *in >> L >> D >> d; L = 20;
	vector<double> Jx(L, 0), Jy(L,0), Jz(L,0), Jxx(L,0), Jyy(L,0), Jzz(L,0);
	forn(i, L){*in >> Jx[i] >> Jy[i] >> Jz[i];} 	forn(i, L-1){*in >> Jxx[i] >> Jyy[i] >> Jzz[i];} //Read the parameters of the Hamiltonian

	//Need to do some bit of modifications here to use the prebuilt Hamiltonian;
	forn(i,L){Jx[i] = 0.; Jy[i] = 0; Jz[i] = 0; Jxx[i] = .2; Jyy[i] = 1.; Jzz[i] = -.3;}
	//forn(i,L){Jx[i] = 0*((rand()%100)/50.-1.); Jy[i] = 0*((rand()%100)/50.-1.); Jz[i] = 0*((rand()%100)/50.-1.);}
	Jxx[0] = (rand()%100)/50.-1.; Jyy[0] = (rand()%100)/50.-1.; Jzz[0] = (rand()%100)/50.-1.;
	forn(i,L){Jxx[i]=Jxx[0];Jyy[i]=Jyy[0];Jzz[i]=Jzz[0];}
	MPS ground_state(L,D,d);ground_state.setRandomState(); // the intial state, random
   ground_state.gaugeCond('R',1);
   ground_state.gaugeCond('L',1);
cerr << "burrito" << endl;
	HeisenbergHamiltonian myH(L,Jxx[0]*4, Jyy[0]*4, Jzz[0]*4,0.*2,2);	const MPO& my_hamil=myH.getHMPO(); double EGnd = 0.; contractor.findGroundState(my_hamil,D,&EGnd,ground_state);

	//We save the Hamiltonian in a format we can work on
	map<int, pair<vector<int>, vector<vector<complex_t > > > > Hmap = representHamiltonian(L,Jx,Jy,Jz,Jxx,Jyy,Jzz);

	//We start by making an initial MPS, which shall be random, for now.
	MPS current(L,D,d);
	current.setRandomState();   current.gaugeCond('R',1);   current.gaugeCond('L',1); cout<<"Initialized random state, norm "<<contractor.contract(current,current)<<endl;
	//Let's try to start with a product state instead of all zeroes
	IsingHamiltonian H0(L,d,0.,0.,-1.); const MPO& hamil=H0.getHMPO(0.); double E0 = 0.; contractor.findGroundState(hamil,D,&E0,current);
	/*//Make sure we start with |0>|0>|0>...|0>
	vector<int> tmp_mask(L,0);
	forn(i,L){
		tmp_mask[i]=1;
		print_matrix(get_rdms(current, contractor, tmp_mask));
		tmp_mask[i]=0;
	}*/
	int xpa; cin >> xpa;
	//Load the operators we have available
	int hbag1_count = 0, hbag2_count = 0;
	map<int, vector<vector<complex_t> > > hbag1 = prepare_h_bag(1, hbag1_count);
	map<int, vector<vector<complex_t> > > hbag2 = prepare_h_bag(2, hbag2_count);

	//Main Loop
	//forn(layer_it, 2){
	for(int layer_it = 1; layer_it <100; ++layer_it){
		//Define a Layer of unitaries
		vector<int> layer = generateLayer(L, layer_it);
		cout << "Iteration " << layer_it << ": Current layer is "; print_vector(layer);
		//Step 0 Initialize unitaries to be identity
		map<int, pair<vector<int>, vector<vector<complex_t > > > > U = initUnitaries(layer);

		//We need the following maps:
			//A map that says, for every local unitary, which terms of the Hamiltonian are going to be affected when it is changed. This needs to be changed for each layer
			map<int, vector<int> > U_to_H = local_terms_for_Unitaries(Hmap, U);
			//We need to define, for each Hamiltonian term, which unitaries are required to compute its energy later!!! This needs to be changed for each layer
			map<int, vector<int> > H_to_U = unitaries_for_Hlocal(Hmap, layer);
			//We also need a map that tells us, for every local term of the Hamiltonian, which rdm is required to compute its energy, given the current layer. The output of this map is a mask, to save memory and time.
			map<int, vector<int> > H_to_maskRho = hamilt_to_masks(Hmap, layer, H_to_U);
			//Finally, we need a map that translates a given mask to the actual rdm
			map<vector<int>, vector<vector<complex_t> > > mask_to_Rho = preset_rdms(current, contractor, H_to_maskRho);



/*Some random tests about computing energy*/
/*forn(j,1){
		MPS hahahah(L,D,d);hahahah.setRandomState(); // the intial state, random
   hahahah.gaugeCond('R',1);
   hahahah.gaugeCond('L',1);
	map<vector<int>, vector<vector<complex_t> > > lel = preset_rdms(hahahah, contractor, H_to_maskRho);
	vector<int> mm(L,0);
	mm[1]=1;
	cout << "Hola!!" << endl;
	print_matrix(get_rdms(hahahah, contractor, mm));
	mm[2]=1;
	cout << "Que tal" << endl;
	print_matrix(get_rdms(hahahah, contractor, mm));
	cout << "com va?" << endl;
	mm[3]=1;
	print_matrix(get_rdms(hahahah, contractor, mm));
	
	cout <<j << " energy with rdm:" << total_energy(U_to_H, H_to_U, H_to_maskRho, lel, Hmap, U) << endl;
	cout << "Norm mps: " << contractor.contract(hahahah, hahahah) << endl;
	cout << "Energy mps: " << contractor.contract(hahahah, my_hamil, hahahah) << endl;
	//forn(k,18){
	//	print_matrix(lel[H_to_maskRho[k]]);
	//}
	int haas; cin >> haas;
}*/

/*
cout << "Here are the 1-body RDMS!" << endl;
	vector<int> tmp_mask(L,0);
	forn(i,L){
		tmp_mask[i]=1;
		print_matrix(get_rdms(current, contractor, tmp_mask));
		tmp_mask[i]=0;
	}
int asdg; cin >> asdg;
cout << "Here are the 2-body RDMS!" << endl;
//	tmp_mask(L,0);
	forn(i,L-1){
		tmp_mask[i]=1; tmp_mask[i+1] = 1;
		print_matrix(get_rdms(current, contractor, tmp_mask));
		tmp_mask[i]=0;
	}
int fasdf; cin >> fasdf;
*/

		double E0 = total_energy(U_to_H, H_to_U, H_to_maskRho, mask_to_Rho, Hmap, U);
		cout << "Initial energy: " << E0 << endl;
		cout << "Norm: " << contractor.contract(current, current) << endl;
		cout << "Energy of the current MPS: " << contractor.contract(current, my_hamil, current) << endl;
		cout << "Ground energy: " << EGnd << endl;
		cout << "Norm ground state: " << contractor.contract(ground_state, ground_state) << endl;
		cout << "Energy ground state: " << contractor.contract(ground_state, my_hamil, ground_state) << endl;

		//int xxx; cin >> xxx;

		//Start see-saw iteration
		for(map<int, pair<vector<int>, vector<vector<complex_t > > > >::iterator it = U.begin(); it!= U.end(); ++it){//We pick every unitary in the layer
			cout << "Iteration " << layer_it << ". Optimizing Unitary " << it->first << endl;
			//Let's compute the part of the Hamiltonian affected by the current unitary.
			cout << U_to_H[it->first] << endl;
			double E0 = energy_contrib(it->first, U_to_H, H_to_U, H_to_maskRho, mask_to_Rho, Hmap, U);
			cout << E0 << endl;
			double t = 0.; //Gets overwritten by function
			if (total(it->second.first)==1){
				forn(i,hbag1_count){
					cout << "Now trying with " << endl; print_matrix(hbag1[i]);
					double E1 = energy_delta_h(it->first, U_to_H, H_to_U, H_to_maskRho, mask_to_Rho, Hmap, U, hbag1[i],t);
					//We add the unitary that is given, as it would decrease the energy
					cout << "Adding unitary #" << i << " for time " << t << endl;
					vector<vector<complex_t > > h = hbag1[i], eh = exponential_h(h, t), UU = it->second.second, mm = matrix_multiply(eh, UU);
					it->second.second = mm; //Actual unitary update
					print_matrix(U[it->first].second);
			cout << "Energy is now " << total_energy(U_to_H, H_to_U, H_to_maskRho, mask_to_Rho, Hmap, U) << endl;
			//int aa; cin >> aa;
				}
			}else{
				//double E1 = energy_delta_h(it->first, U_to_H, H_to_U, H_to_maskRho, mask_to_Rho, Hmap, U, hbag2[0],t);
				forn(i,hbag2_count){
					cout << "Now trying with " << endl; print_matrix(hbag2[i]);
					double E1 = energy_delta_h(it->first, U_to_H, H_to_U, H_to_maskRho, mask_to_Rho, Hmap, U, hbag2[i],t);
					//We add the unitary that is given, as it would decrease the energy
					cout << "Adding unitary #" << i << " for time " << t << endl;
					vector<vector<complex_t > > h = hbag2[i], eh = exponential_h(h, t), UU = it->second.second, mm = matrix_multiply(eh, UU);
					it->second.second = mm; //Actual unitary update
					print_matrix(U[it->first].second);
			cout << "Energy is now " << total_energy(U_to_H, H_to_U, H_to_maskRho, mask_to_Rho, Hmap, U) << endl;
			//int aa; cin >> aa;
				}
			}
		}

		cout << "Unitaries: " << endl << endl;
		for(map<int, pair<vector<int>, vector<vector<complex_t > > > >::iterator it = U.begin(); it!= U.end(); ++it){//We pick every unitary in the layer
			cout << "Unitary #" << it->first << endl;
			cout << it->second.first << endl;
			print_matrix(it->second.second);
			cout << expand_basis(it->second.second) << endl;
		}

		double E = total_energy(U_to_H, H_to_U, H_to_maskRho, mask_to_Rho, Hmap, U);
		cout << " Jxx " << Jxx[0] << " Jyy " << Jyy[0] << " Jzz " << Jzz[0] << endl;		
		cout << "Total energy = " << E << endl;
		current = update_MPS(current, contractor, U);
		cout << "Norm of the current MPS: " << contractor.contract(current, current) << endl;
		cout << "Energy of the current MPS: " << contractor.contract(current, my_hamil, current) << endl;
		complex_t overlap = contractor.contract(current, ground_state);
		cout << "OverlapSq:" << overlap*conjugate(overlap) << endl;
		int xx; cin >> xx;
	}
}
