#include "MIPS.h"
using namespace std;

/*
    * This method is used to create a MIPS object with
    * a lower bound.
    * Parameters:
    * - LowerBound: The lower bound of the problem.
    * Return: 
    * - A MIPS object.
*/
MIPS::MIPS(int LowerBound) {
    LB = LowerBound;
    //Simplex_set = new multiset<Simplex, greater<int>>();
}
/*
	* This method is used to destroy the MIPS object
	* Parameters:
	* - None.
	* Return:
	* -void.
*/
MIPS::~MIPS(){
}

/*
    * This method is used to create a key for the multiset
    * Parameters:
    * - s: The simplex object to be used.
    * Return:
    * - An integer key.
*/
int MIPS::makeKey(Simplex *s){
    return s->getUpperBound();
}
/*
    * This method is used to insert a simplex object into the multiset
    * if his key is greater than the lower bound.
    * Parameters:
    * - s: The simplex object to be inserted.
    * Return:
    * - void.
*/
void MIPS::insertSimplex(Simplex *s){
    if(s->getUpperBound() >= LB){
        // Simplex_set->insert(*s);
    }
}
