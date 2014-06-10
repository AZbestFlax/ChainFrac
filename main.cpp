#include <iostream>
#include "Chain.h"
#include <vector>

using namespace std;

int main()
{
    //frac x(3,-17,11,-1);
    frac x(11,11,11);
    //frac x(2,3,-2);
    //frac x(7);
    x.print();

    /*vector <int> v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    v.push_back(4);

*/
    chain h(x);
    h.convertToChain();
    h.printChain();
    h.convertToFrac();
    h.printFrac();
    return 0;
}
