#include <iostream>
#include <math.h>
#include <string.h>
#include <time.h>
#include "../CROlib/CRO.hpp"

using namespace std;

/* Before your start using CRO, you shall create your own molecule class
 * inherited from the build-in CROMolecule class. In your own class, the
 * encoding scheme of the problem shall be implemented. The initialization of
 * the initial solution is not done here, but if you plan to have pointer member
 * in your class, do remember to release the memory in the destructor.
 */
class myMol : public CROMolecule
{
public:
    double* x;
    int N;

public:
    virtual ~myMol()
    {
        delete[] x;
    }

    /* This "clone" function is a function that must be override in this new
     * class. The purpose is to tell CROlib how to copy the solution in one
     * molecule to another.
     */
    void clone(const myMol* source)
    {
        memcpy(x, source->x, sizeof(double) * N);
    }
};

/* After you have implement your molecule class, the next thing to do is to
 * implement the fitness function. The function head of this function is defined
 * as the following function. The function name can be changed as you wish.
 */
double molFit(const myMol* m)
{
    return sin(m->x[0]) + cos(m->x[1]);
}

/* With the fitness function we can come to the all kinds of operators of CRO.
 * The function names can be changed as you wish.
 */

/* The molecule initialization function is called whenever a new molecule is to
 * be generated. This function shall return a pointer point to the memory stores
 * a molecule object. So you shall "new" the object by your own. Meanwhile you
 * do not have to worry about the memory management: CROlib will release the
 * memory when appropriate.
 */
myMol* molInit()
{
    myMol* m = new myMol();
    m->N = 2;
    m->x = new double[m->N];
    for (int i=0; i < 2; i++)
        m->x[i] = randDouble(); // Generates random numbers in [0,1)
    return m;
}

// On-wall collision operator.
void molWall(const myMol* m, myMol* t)
{
    t->clone(m);
    if (randDouble()>0.5)
        t->x[1] = t->x[1]+randDouble();
    else
        t->x[0] = t->x[0]+randDouble();
}

// Decomposition operator.
void molDec(const myMol* m1, myMol* t1, myMol* t2)
{
    t1->clone(m1);
    t2->clone(m1);
    t1->x[0] = t1->x[0]+randDouble();
    t2->x[1] = t2->x[1]+randDouble();
}

// Inter-molecular operator.
void molInter(const myMol* m1, const myMol* m2, myMol* t1, myMol* t2)
{
    t1->x[0] = m2->x[0];
    t1->x[1] = m1->x[1];
    t2->x[0] = m1->x[0];
    t2->x[1] = m2->x[1];
}

// Synthesis operator.
void molSyn(const myMol* m1, const myMol* m2, myMol* t1)
{
    t1->clone(m1);
    if (randDouble()>0.5)
        t1->x[1] = m2->x[1];
    else
        t1->x[0] = m2->x[0];
}

// Main function.
int main()
{
    // Initialize the seed.
    srand(time(NULL));
    /* Create a CRO object. In the construction period, the initial environment
     * of the container, including initial population, energy buffer, etc., is
     * created. The input attribute is a CROOperator object which stores the
     * pointers to all operators.
     */
    CROOperator<myMol>* opr = new CROOperator<myMol>(molFit, molInit, molWall,
                                                     molDec, molInter, molSyn);
    CRO<myMol>* c1 = new CRO<myMol>(opr);
    /* The above initialization has no parameter as input, so CROlib will
     * initialize one for you. However, you can change the default parameter
     * after you have create the CRO object as follows.
     */
    c1->param->FELimit = 1E6;
    c1->param->iniKE = 100;
    /* You can also initialize one as input as follows.
     */
    CROParameter* param = new CROParameter();
    param->FELimit = 1E6;
    param->iniKE = 100;
    CRO<myMol>* c2 = new CRO<myMol>(opr, param);
    /* Certainly it is fine if you want to add some parameter to the original
     * CRO parameters by creating a new class inherited from CROParameter class.
     * Again, do not worry about the memory employed by the above CROParameter
     * object. It will be release with the release of CRO object.
     */

    /* After everything is done, run the algorithm! The run() function will
     * return the final global optimum got by CRO.
     */
    cout << c1->run() << endl; // Global optimum.
    cout << c1->optMol->x[0] << "\t" << c1->optMol->x[1] << endl; // Best.

    // Don't forget to release CRO objects!
    delete c1;
    delete c2;
    return 0;
}
