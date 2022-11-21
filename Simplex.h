#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
using namespace std;
/*
Esta clase resuelve el problema de optimizacion
lineal con variables reales [x1,....,xn]. 
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
    void solve(); // resuelve el simplex
    bool loadTxt(string filename); // carga la matriz de filename
    Simplex *clone(); // copiar el simplex
    void insertConstraint(int b, int var, int type);
    void print(); //
    private:
    int NM1M2;
    void simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase,
			int izrov[], int iposv[]); // este metodo es solo para para de c a c++
    void simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp, float *bmax); // cambiar nombre
    void simp2(float **a, int m, int n, int *ip, int kp); // cambiar nombre
    void simp3(float **a, int i1, int k1, int ip, int kp); // cambiar nombre
    void nrerror(const char* error_text); // Method to print error messages
    float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch); // Method to convert a 1D array to a 2D array
    void free_ivector(int *v, long nl, long nh); // Method to free an integer vector
    int *ivector(long nl, long nh); // Method to allocate an integer vector
    void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch); // Method to free a matrix allocated by convert_matrix
};