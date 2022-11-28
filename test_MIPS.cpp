#include "MIPS.h"


int main(){
    MIPS *m = new MIPS(0);
    Simplex *s= new Simplex(2,1,1);
    s->loadTxt("sistema.txt");
    s->print();
    m->insertSimplex(s);
    return 0;
}