#include "Simplex.h"
#define NR_END 1
#define EPS 1.0e-6
#define FREEALL             \
	free_ivector(l3, 1, m); \
	free_ivector(l1, 1, n + 1);
#define FREE_ARG char*
using namespace std;

// extern "C" int *ivector(long nl, long nh);
// extern "C" void simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase, int izrov[], int iposv[]);

/*
	* This method is used to create a Simplex object with a
	* given number of variables and constraints.
	* Parameters:
	* - m1: Number of constraints (<=)
	* - m2: Number of constraints (>=)
	* - m3: Number of constraints (=)
*/
Simplex::Simplex(int m1, int m2, int m3){
	// faltan cosas aqui
	this->m1 = m1;
	this->m2 = m2;
	this->m3 = m3;
	this->solved = false;
}
/*
	* This method is used to destroy the Simplex object
	* Parameters:
	* - None.
	* Return:
	* -void.
*/
Simplex::~Simplex(){
}
/*
	* This method is used to get the data from the input file.
	* Parameters:
	* - string filename: the name of the file.
	* Return:
	* - bool: true if the file was loaded correctly, false otherwise.
*/
bool Simplex::loadTxt(string filename){
	string line;
	fstream txt_stream(filename);
	stringstream charStream;
	string number_txt;
	vector<float> *number_float;
	if(!txt_stream.is_open()){
		cout << "Archivo: " << filename << " no encontrado" << endl;
		return false;
	}else{
		cout << "Archivo encontrado " << endl;
	}
	getline(txt_stream, line);
	charStream << line;
	getline(charStream, number_txt, ' ');
	m1 = stoi(number_txt);
	getline(charStream, number_txt, ' ');
	m2 = stoi(number_txt);
	getline(charStream, number_txt, ' ');
	m3 = stoi(number_txt);
	while(getline(txt_stream, line)){
		charStream.clear();
		charStream << line;
		number_float = new vector<float>;
		while(getline(charStream, number_txt, ' ')){
			number_float->push_back(stof(number_txt));
		}
		a.push_back(*number_float);
	}
	N = a[0].size() - 1;
	M = a.size() - 1;
	NP = N + 1;
	MP = M + 2;
	NM1M2 = N + m1 + m2;
	number_float = new vector<float>;
	for (int i = 0; i <= N; i++){
		number_float->push_back(0.0);
	}
	a.push_back(*number_float);
	txt_stream.close();
	return true;
}
/*
	* This method is used to print the structure of the optimization problem
	* to solve and their variables.
	* Parameters:
	* - None.
	* Returns:
	* - void.
*/
void Simplex::print(){
	cout << "m1: " << m1 << " m2:" << m2 << " m3:" << m3 << endl;
	cout << "N:" << N << " M:" << M << endl;
	cout << "NP:" << NP << " MP:" << MP << endl;
	int i = 0, j = 0;
	for (i = 0; i < a.size(); i++){
		for (j = 0; j < a[i].size(); j++){
			cout << a[i][j] << " ";
		}
		cout << endl;
	}
}
// void Simplex::insertConstraint(int b, int var, int type){
// }
/*
	* This method is used to clone a Simplex object
	* Parameters:
	* - None.
	* Return:
	* - Simplex object.
*/
Simplex *Simplex::clone(){
	Simplex *s1 = new Simplex(m1, m2, m3);
	return s1;
}
/*
	* This method is used to solve the problem with the simplex method
	* Parameters:
	* - None.
	* Return:
	* - void.
*/
void Simplex::solve(){
	/*
	#define NP 5        // NP >= N+1
	#define MP 6        // MP >= M+2
	*/
	int i, icase, j, izrov[N + 1], iposv[M + 1];
	// cout << "Iposv tiene tamaño: " << M + 1 << endl;
	float c[MP][NP] =
		{0.0, 1.0, 1.0, 3.0, -0.5,
		740.0, -1.0, 0.0, -2.0, 0.0,
		0.0, 0.0, -2.0, 0.0, 7.0,
		0.5, 0.0, -1.0, 1.0, -2.0,
		9.0, -1.0, -1.0, -1.0, -1.0,
		0.0, 0.0, 0.0, 0.0, 0.0};
	float **a;
	string txt[NM1M2 + 1] =
		{" ", "x1", "x2", "x3", "x4", "y1", "y2", "y3"};
	a = convert_matrix(&c[0][0], 1, MP, 1, NP);
	simplx(a, M, N, m1, m2, m3, &icase, izrov, iposv);
	if (icase == 1){
		cout << endl << "Unbounded objective function" << endl;
	}else if (icase == -1){
		cout << endl << "No solutions satisfy constraints given" << endl;
	}else{
		cout << endl << " ";
		//printf("\n%11s", " ");
		for (i = 1; i <= N; i++){
			if (izrov[i] <= NM1M2){
				cout << txt[izrov[i]];
			}
		}
		cout << endl;
		for (i = 1; i <= M + 1; i++){
			if (i == 1 || iposv[i - 1] <= NM1M2){
				if (i > 1){
					cout << txt[iposv[i - 1]];
				}else{
					cout << "  ";
				}
				cout << fixed << setprecision(2) << setw(10) << a[i][1];
				//printf("%10.10f", a[i][1]);
				for (j = 2; j <= N + 1; j++){
					if (izrov[j - 1] <= NM1M2){
						cout << fixed << setprecision(2) << setw(10) << a[i][j];
						//printf("%10.2f", a[i][j]);
					}
				}
				cout << endl;
			}
		}
	}
	// deben hacer el destructor
	free_convert_matrix(a,1,MP,1,NP);
	free_ivector(iposv,1,M);
	free_ivector(izrov,1,N);
}
/*
	* This method is used to solve an optimization problem using the simplex method.
	* The method is based on the algorithm described in the book
	* Parameters:
	* - a: the matrix of the problem.
	* - m: the number of constraints.
	* - n: the number of variables.
	* - m1_, m2_ m3_: integers obtained from the file with the problem.
	* - icase: an integer initialize to zero.
	* - izrov, iposv: two arrays that will be modified inside the function.
	* Return:
	* -void.
*/
void Simplex::simplx(float **a, int m, int n, int m1_, int m2_, int m3_, int *icase,
			int izrov[], int iposv[]){
	int i, ip, is, k, kh, kp = 0, nl1 = 0, iter = 0;
	int *l1, *l3;
	for(i = 0; i < M + 1; i++){
		for(int jn = 0; jn < N + 1; jn++){
			cout << a[i][jn] << " ";
		}
		cout << endl;
	}
	l3 = (int *)malloc(m2_ * sizeof(int));
	l1 = (int *)malloc(n * sizeof(int));
	float q1, bmax;
	if(m != (m1_ + m2_ + m3_)){
		nrerror("Bad input constraint counts in simplx");
		l1 = ivector(1, n + 1);
		l3 = ivector(1, m);
		nl1 = n;
	}
	for(k = 1; k <= n; k++){
		l1[k-1] = izrov[k-1] = k;
	}
	for(i = 1; i <= m; i++){
		if(i > 2){
			if (a[i-1][1] < 0.0){
			nrerror("Bad input table in simplx");
			}
		}
		iposv[i - 1] = n + i;
	}
	if(m2_ + m3_){
		for(i = 1; i <= m2_; i++){
			l3[i-1] = 1;
			// cout << "l3 es: " << l3[i-1] << endl;
		}
		for(k = 1; k <= (n + 1); k++){
			// cout << "Entramos al for" << endl;
			q1 = 0.0;
			for(i = m1_ + 1; i <= m; i++){
				// cout << "Entramos al for 2" << endl;
				q1 += a[i + 1][k];
			}
			a[m + 2][k] = -q1;
		}
		// Here is the problem, when simplex function tries to use Matrix A.
		for(;;){
			// cout << "Entramos al for 3" << endl;
			simp1(a, m + 1, l1, nl1, 0, &kp, &bmax);
			if(bmax <= EPS && a[m + 2][1] < -EPS){
				// cout << "Entramos al if" << endl;
				*icase = -1;
				FREEALL return;
			}
			else if(bmax <= EPS && a[m + 2][1] <= EPS){
				// cout << "Entramos al else if" << endl;
				for(ip = m1_ + m2_ + 1; ip <= m; ip++){
					// cout << "Entramos al for 4" << endl;
					if(iposv[ip - 1] == (ip + n)){
						// cout << "Entramos al if 2" << endl;
						simp1(a, ip, l1, nl1, 1, &kp, &bmax);
						if (bmax > EPS){
							// cout << "Entramos al if 3" << endl;
							goto One;
						}
					}
				}
				for(i = (m1_ + 1); i <= m1_ + m2_; i++){
					// cout << "Entramos al for 5" << endl;
					// cout << "El valor es: " <<  i - m1_ - 1 << endl;
					if(l3[(i - m1_) - 1] == 1){
						// cout << "Entramos al if 4" << endl;
						for(k = 1; k <= n + 1; k++){
							// cout << "Entramos al for 6" << endl;
							// cout << "A es: " << a[i][k-1] << endl;
							a[i][k-1] = -a[i][k-1];
						}
					}
				}
				//break;
			}
			// cout << "M es: " << m << endl;
			// cout << "N es: " << n << endl;
			// cout << "Kp es: " << kp << endl;
			// cout << "Ip es: " << ip << endl;
			simp2(a, m, n, &ip, kp);
			if (ip == 0){
				*icase = -1;
				FREEALL return;
			}
		One:
			simp3(a, m + 1, n, ip, kp);
			if (iposv[ip] >= (n + m1_ + m2_ + 1)){
				for(k = 1; k <= nl1; k++){
					if(l1[k - 1] == kp){
						break;
					}
				}
				--nl1;
				for (is = k; is <= nl1; is++){
					l1[is] = l1[is + 1];
				}
			}else{
				kh = iposv[ip] - m1_ - n;
				if(kh >= 1 && l3[kh]){
					l3[kh] = 0;
					++a[m + 2][kp + 1];
					for (i = 1; i <= m + 2; i++){
						a[i-1][kp + 1] = -a[i-1][kp + 1];
					}
				}
			}
			is = izrov[kp];
			izrov[kp] = iposv[ip];
			iposv[ip] = is;
		}
	}
	for(;;){
		simp1(a, 0, l1, nl1, 0, &kp, &bmax);
		if (bmax <= EPS){
			*icase = 0;
			FREEALL return;
		}
		simp2(a, m, n, &ip, kp);
		if(ip == 0){
			*icase = 1;
			FREEALL return;
		}
		simp3(a, m, n, ip, kp);
		is = izrov[kp];
		izrov[kp] = iposv[ip];
		iposv[ip] = is;
	}
}
/*
	This method is a part of the Simplex solve method.
	* Parameters:
	* - a: Matrix of the problem.
	* - mm: Number of first row (plus one).
	* - ll: Array of integers.
	* - nll: An integer whose value can be zero or the number of variables.
	* - iabf: An integer number, the value is 0.
	* - kp: An integer number, initial value is 0.
	* - bmax: A float number, initial value is 0.
	* Return:
	* - void	
*/
void Simplex::simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp, float *bmax){
	int k;
	float test = 0.0;
	if (nll <= 0){
		*bmax = 0.0;
	}else{
		*kp = ll[1];
		*bmax = a[mm + 1][*kp + 1];
		for (k = 2; k <= nll; k++){
			if (iabf == 0){
				test = a[mm + 1][ll[k] + 1] - (*bmax);
			}else{
				test = std::fabs(a[mm + 1][ll[k] + 1]) - std::fabs(*bmax);
			}
			if (test > 0.0){
				*bmax = a[mm + 1][ll[k] + 1];
				*kp = ll[k];
			}
		}
	}
}
/*
	This method is a part of the Simplex solve method.
	* Parameters:
	* - a: Matrix of the problem.
	* - m: Number of variables.
	* - n: Number of constraints.
	* - ip: Row of the pivot.
	* - kp: Column of the pivot.
	* Return:
	* - void	
*/
void Simplex::simp2(float **a, int m, int n, int *ip, int kp){
	int k,i;
	float qp,q0,q = 0.0, q1;
	*ip=0;
	for (i = 1; i <= m; i++){
		// cout << "Entramos al for 7" << endl;
		if (a[i-1][kp] < -EPS){
			break;
		}
	}
	if (i > m){
		// cout << "Entramos al if 5" << endl;
		return;
	}
	q1 = -a[i][1]/a[i][kp+1];
	*ip=i;
	for(i = *ip+1; i <= m; i++){
		// cout << "Entramos al for 8" << endl;
		if(a[i+1][kp+1] < -EPS){
			// cout << "Entramos al if 6" << endl;
			q = -a[i+1][1]/a[i+1][kp+1];
			if(q < q1){
				// cout << "Entramos al if 7" << endl;
				*ip=i;
				q1=q;
			}
			else if(q == q1){
				// cout << "Entramos al if 8" << endl;
				for(k=1;k<=n;k++){
					// cout << "Entramos al for 9" << endl;
					qp = -a[*ip+1][k]/a[*ip+1][kp+1];
					q0 = -a[i+1][k]/a[i+1][kp+1];
					if (q0 != qp){
						// cout << "Entramos al if 9" << endl;
						break;
					}
				}
				if (q0 < qp){
					// cout << "Entramos al if 10" << endl;
					*ip=i;
				}
			}
		}
	}
}
/*
	* This method is a part of the Simplex solve method.
	* Parameters:
	* - a: Matrix of the problem.
	* - i1: Number of variables.
	* - i2: Number of constraints.
	* - ip: Row of the pivot.
	* - kp: Column of the pivot.
	* Return:
	* - void
*/
void Simplex::simp3(float **a, int i1, int k1, int ip, int kp){
	int kk,ii;
	float piv;
	piv=1.0/a[ip+1][kp+1];
	for(ii=1;ii<=i1+1;ii++){
		if(ii-1 != ip){
			a[ii][kp+1] *= piv;
			for(kk=1;kk<=k1+1;kk++){
				if (kk-1 != kp){
					a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1];
				}
			}
		}
	}
	for(kk=1;kk<=k1+1;kk++){
		if(kk-1 != kp){
			a[ip+1][kk] *= -piv;
		}
	}
	a[ip+1][kp+1]=piv;
}
/*
	* This method is used to generate a matrix with the given parameters
	* Parameters:
	* - a: First element of the matrix
	* - nrl: number of first row.
	* - nrh: Number of basic variables.
	* - ncl: number of first column.
	* - nch: Number of no-basic variables.
	* Return:
	* - Matrix with the given parameters.
*/
float **Simplex::convert_matrix(float *a, long nrl, long nrh, long ncl, long nch){
	// The parameters are: a = 0.0, nrl = 1, nrh = 6, ncl = 1, nch = 5.
	// nrow = 6, ncol = 5
	long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	float **m;
	m = (float **)malloc((nrow + NR_END) * sizeof(float *));
	if (!m){
		nrerror("Allocation failure in convert_matrix()");
	}
	m += NR_END;
	m -= nrl;
	m[nrl] = a - ncl;
	for (i = 1, j = nrl + 1; i < nrow; i++, j++){
		m[j] = m[j - 1] + ncol;
	}
	return m;
}
void Simplex::nrerror(const char* error_text){
	perror("Numerical Recipes run-time error...\n");
	perror(error_text);
	perror("\n...now exiting to system...\n");
	exit(1);
}
/*
	* This method is used to free the memory allocated for the vector
	* Parameters:
	* - v: vector to free
	* - nrl: first index of the vector
	* - nrh: last index of the vector
	* Return:
	* - void
*/
void Simplex::free_ivector(int *v, long nl, long nh){
	free((FREE_ARG) (v+nl-NR_END));
}
/*
	* This method is used to create a vectors of integers.
	* Parameters:
	* - nl: lower index
	* - nh: higher index
	* Return:
	* - A vector of integers.
*/
int *Simplex::ivector(long nl, long nh){
	int *v;
	v = (int *)malloc((nh-nl+1+NR_END)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}
/*
	* This method is used to free a matrix.
	* Parameters:
	* - b: matrix to free.
	* - nrl: first row.
	* - nrh: last row.
	* - ncl: first column.
	* - nch: last column.
	* Return:
	* - void.
*/
void Simplex::free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch){
	free((FREE_ARG) (b+nrl-NR_END));
}
/*
	* This method is used to get the UpperBound of the problem.
	* Parameters:
	* - None.
	* Return:
	* - UpperBound of the problem.
*/
int Simplex::getUpperBound(){
	return 1;
}