#include "Simplex.h"


int main(){
    Simplex *s= new Simplex(2,1,1);
    s->loadTxt("sistema.txt");
    s->print();
    // s->solve();
    return 0;
}