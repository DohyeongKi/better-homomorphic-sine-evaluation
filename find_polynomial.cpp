/*
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

#include <iostream>
#include <fstream>
#include <string>

#include <NTL/ZZ.h>
#include <NTL/RR.h>

using namespace std;
using namespace NTL;

#define M_PIl 	3.141592653589793238462643383279502884L 
#define PI		RR(M_PIl)


//=============================================================
//==Functions==================================================
//=============================================================

//--------------------------------------------------------
//--Find an index corresponding to the maximum value------ 
//--------------------------------------------------------
int max_index(double* array, int length) {
	int max_ind	= 0;
	double max = array[0];
	
	for(int i=1; i<length; i++) {
		if(array[i] > max) {
			max_ind = i;
			max = array[i];
		}
	}
	return max_ind;
}
	
//--------------------------------------------------------
//--Long division between polynomials represented--------- 
//---in Chebyshev basis-----------------------------------
//-------------------------------------------------------- 
void Long_div(RR* p, int n, RR* q, int m, RR* quot, RR* rem) {
//--------------------------------------------------------
//--p : Dividend-------q : Divisor------------------------
//--n : degree of p----m : degree of q--------------------
//--quot : quotient----rem : remainder--------------------
//--------------------------------------------------------
	if(m>0) {
		for(int i=n; i>m; i--) {
			RR ratio = p[i]/q[m];
			for(int j=0; j<=m; j++) {
				p[i-m+j] -= q[j]*ratio;
				p[abs(i-m-j)] -= q[j]*ratio;
			}
			quot[i-m] = RR(2.0)*ratio;
		}
		
		RR ratio = p[m]/q[m];
		for(int j=0; j<=m; j++) 
			p[j] -= q[j]*ratio;
		quot[0] = ratio;
		
	} else {
		for(int i=n; i>=0; i--) {
			quot[i] = p[i]/q[0];
			p[i] = RR(0.0);
		}
	}

	for(int i=0; i<m; i++)
		rem[i] = p[i];
}


//=============================================================
//==Main=======================================================
//=============================================================

int main(int argc, char* argv[]) {

//=============================================================
//==Setting====================================================
//=============================================================
	
	int K = 12;								// I_i = [i-.25-e, i-.25+e] where |i|< K
	
	int deg_bdd = atoi(argv[1]) + 1;		// Bound of the degree +1
	
	int* deg = new int[K];					// deg[i] = The number of nodes in I_i
	for(int i=0; i<K; i++)					// We assume deg[i] = deg[-i]
		deg[i] = 1;							// Initialize all deg[i] to 1	
	int tot_deg = 2*K-1;					// Total number of nodes

	int dev = atoi(argv[2]);				
	double err = 1.0/(1 << atoi(argv[2]));	// Maximum deviation from each i-.25
	
	int sc_num = atoi(argv[3]);				// The number of scaling
	RR sc_fac = conv<RR>(ZZ(1) << sc_num);	// Scaling factor
	
	RR::SetPrecision(1000);

//=============================================================
//==Degree Searching===========================================
//=============================================================
//--------------------------------------------------------
//--Initialize--------------------------------------------
//--------------------------------------------------------

	double* bdd = new double[K];

	double temp = 0;
	for(int i=1; i<=(2*K-1); i++) 
		temp -= log2((double)i);
	temp += (2*K-1)*log2(2*M_PI);
	temp += log2(err);

	for(int i=0; i<K; i++) {
		bdd[i] = temp;
		for(int j=1; j<=K-1-i; j++)
			bdd[i] += log2((double)j + err);
		for(int j=1; j<=K-1+i; j++)
			bdd[i] += log2((double)j + err);
	}

//--------------------------------------------------------
//--Algorithm--------------------------------------------- 
//--1. Find a point that has the largest theoretical error bound.
//--2. Increase degree by one at that point.--------------
//--(If the point is not 0, also increase degree by-------
//---one at negate of that point)------------------------- 
//--3. Check whether total degree is greater than the degree bound
//---If so, end the algorithm. If not, go back to 1.------								
//--------------------------------------------------------

	int max_iter = 200;	// Bound of the number of iteration
	int iter;

	for(iter=0; iter<max_iter; iter++) {
		if(tot_deg >= deg_bdd)
			break;
		int maxi = max_index(bdd, K);	
		
		if(maxi != 0) {
			if((tot_deg+2) > deg_bdd) 
				break; 
	
			for(int i=0; i<K; i++) {
				bdd[i] -= log2(tot_deg+1);
				bdd[i] -= log2(tot_deg+2);
				bdd[i] += 2.0*log2(2.0*M_PI);

				if(i != maxi) {	
					bdd[i] += log2(abs((double)(i-maxi)) + err);
					bdd[i] += log2((double)(i+maxi) + err);
				} else { // i = maxi
					bdd[i] += (log2(err)-1.0);
					bdd[i] += log2(2.0*(double)i + err);
				}
			}
			
			tot_deg += 2;
		} else { // maxi = 0
			bdd[0] -= log2(tot_deg+1);
			bdd[0] += (log2(err)-1.0);
			bdd[0] += log2(2.0*M_PI);
			for(int i=1; i<K; i++) {	
				bdd[i] -= log2(tot_deg+1);
				bdd[i] += log2(2.0*M_PI);
				bdd[i] += log2((double)i + err);	
			}
		
			tot_deg += 1;	
		}
		
		deg[maxi] += 1;
	}
	
	delete[] bdd;

//--------------------------------------------------------
//--Print the Result of Degree Searching------------------
//--------------------------------------------------------

	cout << "==============================================" << endl;
	cout << "==Degree Searching Result=====================" << endl;
	cout << "==============================================" << endl;

	if(iter == max_iter) {
		cout << "More Iteration Needed" << endl;
	} else {
		cout << "Degree of Polynomial : " << tot_deg-1 << endl;
		cout << "Degree : ";
		for(int i=0; i<K; i++) 
			cout << deg[i] << "  ";
		cout << endl;
	}
	cout << "==============================================" << endl;

//=============================================================
//==Find an Interpolation Polynomial===========================
//==Goal : Interpolate cos(2PI x)==============================
//=============================================================	
//--------------------------------------------------------
//--Node Setting------------------------------------------
//--------------------------------------------------------
	
	RR inter_size = RR(1.0)/conv<RR>(((ZZ)(1) << dev)); 
	// Half of the size of each interval
	
	RR* z = new RR[tot_deg];	// Node positions	
	int cnt = 0;
	if((deg[0]%2)!=0)
		z[cnt++] = -RR(0.25);

	for(int i=K-1; i>0; i--) {
		for(int j=1; j<=deg[i]; j++) {
			RR temp = ((RR(2*j-1))*PI)/(RR(2*deg[i]));		
			z[cnt++] = RR(i - 0.25) + inter_size*cos(temp);
			z[cnt++] = RR(-i - 0.25) - inter_size*cos(temp);
		}
	}

	for(int j=1; j<=(deg[0]/2); j++) {
		RR temp = ((RR(2*j-1))*PI)/(RR(2*deg[0]));
		z[cnt++] = RR(-0.25) + inter_size*cos(temp);
		z[cnt++] = RR(-0.25) - inter_size*cos(temp);
	}
	
	for(int i=0; i<tot_deg; i++) 
		z[i] /= sc_fac;
	
	delete[] deg;

//--------------------------------------------------------
//--Algorithm---------------------------------------------
//--------------------------------------------------------

	RR* d = new RR[tot_deg];
	for(int i=0; i<tot_deg; i++) 
		d[i] = cos(RR(2.0)*PI*z[i]);

	for(int j=1; j<tot_deg; j++) {
		for(int l=0; l<tot_deg-j; l++) 
			d[l] = (d[l+1] - d[l]) / (z[l+j] - z[l]);
	}

//=============================================================
//==Compute Chebyshev Coefficients by Solving Matrix Equation==
//==Result Polynomial :     ===================================
//==	c[0]T_0(x) + ... + c[tot_deg-1]T_{tot_deg-1}(x)   =====
//== where T_i(x) is an adjusted Chebyshev polynomial =========
//=============================================================

	tot_deg += 1;
	
	RR* x = new RR[tot_deg];
	for(int i=0; i<tot_deg; i++) 
		x[i] = RR(K)/sc_fac * cos(RR(i)*PI/RR(tot_deg-1));	

	RR* c = new RR[tot_deg];
	RR* p = new RR[tot_deg];
	for(int i=0; i<tot_deg; i++) {
		p[i] = d[0];
		for(int j=1; j<tot_deg-1; j++)
			p[i] = p[i]*(x[i] - z[j]) + d[j];
	}
	
	delete[] z;

	RR** T = new RR*[tot_deg];
	for(int i=0; i<tot_deg; i++)
		T[i] = new RR[tot_deg];
	
	for(int i=0; i<tot_deg; i++) {
		T[i][0] = RR(1.0);
		T[i][1] = x[i]/(RR(K)/sc_fac);
		for(int j=2; j<tot_deg; j++)
			T[i][j] = RR(2.0)*(x[i]/(RR(K)/sc_fac))*T[i][j-1] - T[i][j-2];
	}

	
	for(int i=0; i<tot_deg-1; i++) {
		RR max_abs = abs(T[i][i]);
		int max_index = i;
		for(int j = i+1; j<tot_deg; j++) {
			if(abs(T[j][i]) > max_abs) {
				max_abs = abs(T[j][i]);
				max_index = j;
			}
		}
		
		if(i != max_index) {
			for(int j=i; j<tot_deg; j++) {
				RR temp = T[max_index][j];
				T[max_index][j] = T[i][j];
				T[i][j] = temp;
			}

			RR temp = p[max_index];
			p[max_index] = p[i];
			p[i] = temp;
		}
		
		for(int j=i+1; j<tot_deg; j++)
			T[i][j] /= T[i][i];
		p[i] /= T[i][i];
		T[i][i] = RR(1.0);

		for(int j=i+1; j<tot_deg; j++) {
			p[j] -= T[j][i] * p[i];
			for(int l=i+1; l<tot_deg; l++)
				T[j][l] -= T[j][i] * T[i][l];
			T[j][i] = RR(0.0);
		}	
	}

	c[tot_deg-1] = p[tot_deg-1];
	for(int i=tot_deg-2; i>=0; i--) {
		c[i] = p[i];
		for(int j=i+1; j<tot_deg; j++)
			c[i] -= T[i][j]*c[j];
	}

	tot_deg -= 1;
	
	for(int i=0; i<tot_deg; i++)
		delete[] T[i];
	delete[] T;
	delete[] d;
	delete[] x;
	delete[] p;

//=============================================================
//==Baby Step Giant Step Algorithm=============================
//=============================================================	
//--------------------------------------------------------
//--Parameter Setting-------------------------------------
//--------------------------------------------------------

	int temp_tot_deg = tot_deg; 
	
	int m = 1;					// m = ceil(log2(tot_deg))					
	while (temp_tot_deg > 1) {
		m++;
		temp_tot_deg /= 2;
	}
	
	int* pow2 = new int[m+1];	// pow2[i] = 2^i (i=0, ...,m)
	pow2[0] = 1;
	for(int i=0; i<m; i++)
		pow2[i+1] = 2*pow2[i];
	
	int l;						// l ~ m/2 
	if( m % 2 == 0)	{
		l = m/2;
	} else {						
		int l1 = m/2;
		int l2 = m/2 + 1;
		
		l = (pow2[l1] + pow2[m-l1] - l1 
				<= pow2[l2] + pow2[m-l2] - l2) ? l1 : l2;
		// Choose one that requires less number 
		//		of non-scalar multiplications
	}

//--------------------------------------------------------
//--Algorithm---------------------------------------------
//--Details:    ------------------------------------------
//---1. alg_coef[0][0] represents the interpolation poly--
//---2. alg_coef[i][j] represents polynomial p_{i,j}------
//---3. p_{i,j} = p_{i+1,2j} + p_{i+1,2j+1} T_{2^(m-i-1)}-
//--------------------------------------------------------		
	
	RR*** alg_coef = new RR**[m-l+1];
	for(int i=0; i<m-l+1; i++) {
		alg_coef[i] = new RR*[pow2[i]]; 
		for(int j=0; j<pow2[i]; j++)
			alg_coef[i][j] = new RR[pow2[m-i]];
	}
	
	for(int i=0; i<tot_deg; i++)
		alg_coef[0][0][i] = c[i];
	for(int i=tot_deg; i<pow2[m]; i++)
		alg_coef[0][0][i] = RR(0.0);	

	delete[] c;

	for(int i=0; i<m-l; i++) {
		RR* divisor = new RR[pow2[m-i-1]+1];
		for(int j=0; j<pow2[m-i-1]; j++)
			divisor[j] = RR(0.0);
		divisor[pow2[m-i-1]] = RR(1.0);
		
		for(int j=0; j<pow2[i]; j++) 
			Long_div(alg_coef[i][j], pow2[m-i]-1, divisor, pow2[m-i-1],
				 		alg_coef[i+1][2*j+1], alg_coef[i+1][2*j]);
	}

//--------------------------------------------------------
//--Print Algorithm Coefficients in File------------------
//--------------------------------------------------------

	string path_coef = "./result/coef/Deg";
	path_coef += (to_string(tot_deg-1) 
			+ "Err" + to_string(dev) 
			+ "Scale" + to_string(sc_num) + ".csv");

	ofstream output_coef(path_coef);

	for(int i=0; i<pow2[m-l]; i++) {
		for(int j=0; j<pow2[l]; j++) 	
			output_coef << alg_coef[m-l][i][j] << ", ";
		output_coef << endl;
	}
	
	output_coef.close();

//=============================================================
//==Find Maximum Error and Print Errors in File================
//=============================================================	
//--------------------------------------------------------
//--File Path Setting-------------------------------------
//--------------------------------------------------------

	string path_err = "./result/error/Deg";
	path_err += (to_string(tot_deg-1) 
			+ "Err" + to_string(dev) 
			+ "Scale" + to_string(sc_num) + ".csv");

	ofstream output_err(path_err);
	
	output_err << "Tested Values" << "," << "Real Values" << "," 
			<< "Approximate Values" << "," << "Error" << "," 
			<< "log2(Error)" << endl;
	
//--------------------------------------------------------
//--Test Nodes Setting------------------------------------
//--------------------------------------------------------

	int test_num = 20;		// The number of test points in each interval I_i

	RR** test = new RR*[2*K-1];	
	for(int i=0; i<2*K-1; i++)
		test[i] = new RR[test_num+1];
	
	RR incr = RR(2.0)*inter_size*RR(1.0/test_num);
	test[0][0] = -RR(0.25) - inter_size;
	for(int i=1; i<test_num+1; i++)	
		test[0][i] = test[0][i-1] + incr;
	
	for(int i=1; i<=K-1; i++) {
		test[2*i-1][0] = RR(-i - 0.25) - inter_size;
		for(int j=1; j<test_num+1; j++) 
			test[2*i-1][j] = test[2*i-1][j-1] + incr;
		test[2*i][0] = RR(i - 0.25) - inter_size;
		for(int j=1; j<test_num+1; j++) 
			test[2*i][j] = test[2*i][j-1] + incr;
	}

//--------------------------------------------------------
//--Computation and Print Errors in File------------------
//--------------------------------------------------------

	RR max = RR(-999.0);
	for(int i=0; i<2*K-1; i++) {
		for(int j=0; j<test_num+1; j++) {
			
			RR real = cos(RR(2.0) * PI * test[i][j]); 	// Real value of cos(2PI x)
			RR approx = RR(0.0);						// Approximate value of cos(2PI x)
			
			RR* BS = new RR[pow2[l]]; 		// Baby-step basis  : T_0(x), ... , T_{2^l-1}(x)	
			RR* GS = new RR[m-l];			// Giant-step basis : T_{2^l)(x), ... , T_{2^(m-1)}(x)

			BS[0] = RR(1.0);
			BS[1] = (test[i][j]/RR(K));
		
			for(int k=2; k<pow2[l]; k++)
				BS[k] = RR(2.0)*BS[k/2]*BS[k-k/2] - BS[k-2*(k/2)];

			GS[0] = RR(2.0)*BS[pow2[l-1]]*BS[pow2[l-1]] - RR(1.0);
			for(int k=1; k<m-l; k++)
				GS[k] = RR(2.0)*GS[k-1]*GS[k-1] - RR(1.0);
			
			RR** alg_value = new RR*[m-l+1]; 	// Recall that alg_coef[i][j] represents polynomial p_{i,j}
			for(int k=0; k<m-l+1; k++) 			// alg_value[i][j] : The value of p_{i,j} at x
				alg_value[k] = new RR[pow2[k]]; // alg_value[0][0] : The value of interpolation poly at x  

			for(int k=0; k<pow2[m-l]; k++) {
				RR temp = RR(0.0);
				for(int s=0; s<pow2[l]; s++)
					temp += alg_coef[m-l][k][s]*BS[s];
				alg_value[m-l][k] = temp;
			}
	
			for(int k=m-l-1; k>=0; k--) {
				for(int s=0; s<pow2[k]; s++) 
					alg_value[k][s] = alg_value[k+1][2*s] + GS[m-l-k-1]*alg_value[k+1][2*s+1];
			}
			approx = alg_value[0][0];
			
			for(int k=0; k<sc_num; k++)
				approx = RR(2.0)*approx*approx - RR(1.0);	// double angle formula			

			output_err << test[i][j] << "," << real << "," 
						<< approx << "," << approx-real << ",";

			if(approx-real!=0) {
				if(max < log(abs(approx-real))/log(RR(2.0)))
					max = log(abs(approx-real))/log(RR(2.0));

				output_err << log(abs(approx-real))/log(RR(2.0)) << endl;
			} else {
				output_err << "*" << endl;
			}
		}	
	}

	output_err.close();	

//--------------------------------------------------------
//--Print Maximum Error of the interpolation polynomial---
//--------------------------------------------------------

	cout << "==============================================" << endl;
	cout << "==Baby Step Giant Step Algorithm Result=======" << endl;
	cout << "==============================================" << endl;
	cout << "Max_Error : " << max << endl;
	cout << "==============================================" << endl;
	
	for(int i=0; i<m-l+1; i++) {
		for(int j=0; j<pow2[i]; j++) 
			delete[] alg_coef[i][j];
		delete[] alg_coef[i];
	}
	delete[] alg_coef;
	
	for(int i=0; i<2*K-1; i++)
		delete[] test[i];
	delete[] test;
}
