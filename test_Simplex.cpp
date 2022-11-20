#include "Simplex.h"

// ojo esto va a cambiar para mas adelante
// solo funciona con la version de Nov 9
int main(){
    Simplex *s= new Simplex(2,1,1);
    //s->solve();
    s->loadTxt("sistema.txt");
    s->print();
    s->solve();
}