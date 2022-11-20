#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
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
    void simplx_tmp(); // este metodo es solo para para de c a c++
    void simp1_tmp(); // cambiar nombre
    void simp2_tmp(); // cambiar nombre
    void simp3_tmp(); // cambiar nombre
};