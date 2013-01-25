/* This file is a part of CROlib
 * Copyright (C) 2013 James J. Q. Yu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <math.h>
#include <string.h>
#include "../CROlib/CRO.hpp"

using namespace std;

double molFit(const CROMolecule* m);
void molWall(const CROMolecule* m, CROMolecule* t);
void molDec(const CROMolecule* m1, CROMolecule* t1, CROMolecule* t2);
void molInter(const CROMolecule* m1, const CROMolecule* m2, CROMolecule* t1,
              CROMolecule* t2);
void molSyn(const CROMolecule* m1, const CROMolecule* m2, CROMolecule* t1);

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
    myMol()
    {
        oprWall = molWall;
        oprDec = molDec;
        oprInter = molInter;
        oprSyn = molSyn;
        N = 2;
        x = new double[N];
        for (int i=0; i < 2; i++)
            x[i] = CRORandDouble(); // Generates random numbers in [0,1)
        }

    virtual ~myMol()
    {
        delete[] x;
    }

    /* This "clone" function is a function that must be override in this new
     * class. The purpose is to tell CROlib how to copy the solution in one
     * molecule to another.
     */
    virtual void clone(const myMol* source)
    {
        CROMolecule::clone(source);
        memcpy(x, source->x, sizeof(double) * N);
    }
};

/* After you have implement your molecule class, the next thing to do is to
 * implement the fitness function. The function head of this function is defined
 * as the following function. The function name can be changed as you wish.
 */
double molFit(const CROMolecule* mol)
{
    myMol* m = (myMol*)mol;
    return sin(m->x[0]) + cos(m->x[1]);
}

/* With the fitness function we can come to the all kinds of operators of CRO.
 * The function names can be changed as you wish.
 */

// On-wall collision operator.
void molWall(const CROMolecule* in, CROMolecule* out)
{
    myMol* m = (myMol*)in;
    myMol* t = (myMol*)out;
    t->clone(m);
    if (CRORandDouble()>0.5)
        t->x[1] = t->x[1]+CRORandDouble();
    else
        t->x[0] = t->x[0]+CRORandDouble();
}

// Decomposition operator.
void molDec(const CROMolecule* in1, CROMolecule* out1, CROMolecule* out2)
{
    myMol* m1 = (myMol*)in1;
    myMol* t1 = (myMol*)out1;
    myMol* t2 = (myMol*)out2;
    t1->clone(m1);
    t2->clone(m1);
    t1->x[0] = t1->x[0]+CRORandDouble();
    t2->x[1] = t2->x[1]+CRORandDouble();
}

// Inter-molecular operator.
void molInter(const CROMolecule* in1, const CROMolecule* in2, CROMolecule* out1,
              CROMolecule* out2)
{
    myMol* m1 = (myMol*)in1;
    myMol* m2 = (myMol*)in2;
    myMol* t1 = (myMol*)out1;
    myMol* t2 = (myMol*)out2;
    t1->x[0] = m2->x[0];
    t1->x[1] = m1->x[1];
    t2->x[0] = m1->x[0];
    t2->x[1] = m2->x[1];
}

// Synthesis operator.
void molSyn(const CROMolecule* in1, const CROMolecule* in2, CROMolecule* out1)
{
    myMol* m1 = (myMol*)in1;
    myMol* m2 = (myMol*)in2;
    myMol* t1 = (myMol*)out1;
    t1->clone(m1);
    if (CRORandDouble()>0.5)
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
     * created.
     */
    CRO<myMol>* c1 = new CRO<myMol>(molFit);
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
    CRO<myMol>* c2 = new CRO<myMol>(molFit, param);
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
