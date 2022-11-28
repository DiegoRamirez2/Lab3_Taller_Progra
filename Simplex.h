#ifndef _Simplex_h_
#define _Simplex_h_
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
using namespace std;
/*
    * This class will be used to model and solve the optimization 
    * problem with simplex method, it will contain the problem 
    * variables and constraints.
*/

class Simplex {
public:
    std::vector<std::vector<float>> a;
    int m1; // numero de restricciones <=
    int m2; // ... >=
    int m3; // igualdad
    int N; // cantidad de variables
    int M; // cantidad de restricciones
    int NP; // cantidad de variable no basicas
    int MP; // cantidad de variables basicas
    bool solved; // si es que esta resuelto
    Simplex(int m1, int m2, int m3);
    ~Simplex();
    void solve(); 
    bool loadTxt(string filename); 
    Simplex *clone(); 
    int getUpperBound(); // This method will be implemented correctly in the future.
    // void insertConstraint(int b, int var, int type);
    void print();
private:
    int NM1M2;
    void simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase,
			int izrov[], int iposv[]);
    void simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp, float *bmax);
    void simp2(float **a, int m, int n, int *ip, int kp);
    void simp3(float **a, int i1, int k1, int ip, int kp);
    void nrerror(const char* error_text);
    float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
    void free_ivector(int *v, long nl, long nh);
    int *ivector(long nl, long nh);
    void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
};
#endif /* _Simplex_h_ */