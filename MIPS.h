#ifndef _MIPS_h_
#define _MIPS_h_
#include "Simplex.h"
#include <set>
#include <iterator>
using namespace std;
/*
    * This class will be used to contain the multiset of simplex objects
    * It will be used to explore the paths of better solutions 
    * through the ordering in the multiset and the upper and lower bounds. 
    * At the moment the multiset will not be included since the function that 
    * allows generating an order in the multiset must be implemented.
*/
class MIPS {
private:
    int LB;
public:
    // multiset<Simplex, greater<int>> *Simplex_set;
    MIPS(int LowerBound);
    ~MIPS();
    int makeKey(Simplex *s);
    void insertSimplex(Simplex *s);
};

#endif /* _MIPS_h_ */