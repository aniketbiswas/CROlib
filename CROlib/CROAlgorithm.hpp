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

#pragma once
#include <vector>
#include "CROMolecule.hpp"
#include "CROParameter.hpp"
#include "CRORandom.hpp"

template<class mol>
class CRO
{
public:
    CROParameter* param;
    mol* t1,* t2,* optMol;
    int curFE;
    double optGlobal, eBuffer;
    std::vector<mol*> pop;

public:
    CRO(CROParameter* usrParam = NULL) :
        param(usrParam), optGlobal(1E300)
    {
        if (param == NULL)
            param = new CROParameter();

        t1 = new mol();
        t1->init(param->iniKE);
        t2 = new mol();
        t2->init(param->iniKE);
        optMol = new mol();
        optMol->init(param->iniKE);

        curFE = param->iniPopSize;
        for (int i = 0; i < param->iniPopSize; i++)
        {
            mol* m = new mol();
            m->init(param->iniKE);
            pop.push_back(m);
            update(m);
        }
        eBuffer = param->iniBuffer;
    }

    virtual ~CRO()
    {
        if (param != NULL)
            delete param;
        for (unsigned int i = 0; i < pop.size(); i++)
            delete pop[i];
        pop.clear();
    }

    double run()
    {
        while(curFE < param->FELimit)
        {
            if ((CRORandDouble() > param->collRate)||(pop.size() == 1))
            {
                int pos = CRORandInt(pop.size());
                // Decomposition
                if (pop[pos]->isInactive(param->decThres))
                {
                    mol* m = new mol();
                    m->init(param->iniKE);
                    if (dec(pop[pos], m))
                    {
                        pop.push_back(m);
                        update(pop[pos]);
                        update(m);
                    }
                    else
                    {
                        delete m;
                    }
                    curFE += 2;
                }
                // On-wall
                else
                {
                    wall(pop[pos]);
                    update(pop[pos]);
                    curFE++;
                }
            }
            else
            {
                int pos1 = CRORandInt(pop.size());
                int pos2 = pos1;
                while (pos1 == pos2)
                    pos2 = CRORandInt(pop.size());
                // Synthesis
                if ((pop[pos1]->KE < param->synThres)&&
                        (pop[pos2]->KE < param->synThres))
                {
                    if (syn(pop[pos1], pop[pos2]))
                    {
                        update(pop[pos1]);
                        delete pop[pos2];
                        pop.erase(pop.begin() + pos2);
                    }
                    curFE++;
                }
                // Inter-molecular
                else
                {
                    inter(pop[pos1], pop[pos2]);
                    update(pop[pos1]);
                    update(pop[pos2]);
                    curFE += 2;
                }
            }
        }
        return optGlobal;
    }

private:
    void update(const mol* m)
    {
        if (m->PE < optGlobal)
        {
            optGlobal = m->PE;
            optMol->clone(m);
        }
    }

    void wall(mol* m)
    {
        m->curHitIndex++;
        m->oprWall(m, t1);
        double tempPE = m->oprFit(t1);
        double excessEnergy = m->PE + m->KE - tempPE;
        if (excessEnergy >= 0)
        {
            m->PE = tempPE;
            m->KE = excessEnergy * ((1.0 - param->lossRate) *
                                    CRORandDouble() + param->lossRate);
            m->clone(t1);
            m->update();
            eBuffer += excessEnergy - m->KE;
        }
    }

    bool dec(mol* m1, mol* m2)
    {
        m1->curHitIndex++;
        m1->oprDec(m1, t1, t2);
        double tempPE1 = m1->oprFit(t1);
        double tempPE2 = m2->oprFit(t2);
        double excessEnergy = m1->PE + m1->KE - tempPE1 - tempPE2;
        if ((excessEnergy >= 0) || (excessEnergy + eBuffer >= 0))
        {
            if (excessEnergy >= 0)
            {
                m1->KE = excessEnergy * CRORandDouble();
                m2->KE = excessEnergy - m1->KE;
            }
            else
            {
                eBuffer += excessEnergy;
                m1->KE = eBuffer * CRORandDouble() * CRORandDouble();
                eBuffer -= m1->KE;
                m2->KE = eBuffer * CRORandDouble() * CRORandDouble();
                eBuffer -= m2->KE;
            }
            m1->minHitIndex = 0;
            m1->curHitIndex = 0;
            m1->PE = tempPE1;
            m1->clone(t1);
            m1->update();
            m2->PE = tempPE2;
            m2->clone(t2);
            m2->update();
            return true;
        }
        return false;
    }

    void inter(mol* m1, mol* m2)
    {
        m1->curHitIndex++;
        m2->curHitIndex++;
        m1->oprInter(m1, m2, t1, t2);
        double tempPE1 = m1->oprFit(t1);
        double tempPE2 = m2->oprFit(t2);
        double excessEnergy = m1->PE + m1->KE + m2->PE + m2->KE -
                              tempPE1 - tempPE2;
        if (excessEnergy >= 0)
        {
            m1->PE = tempPE1;
            m1->KE = excessEnergy * CRORandDouble();
            m1->clone(t1);
            m1->update();
            m2->PE = tempPE2;
            m2->KE = excessEnergy - m1->KE;
            m2->clone(t2);
            m2->update();
        }
    }

    bool syn(mol* m1, mol* m2)
    {
        m1->curHitIndex++;
        m1->oprSyn(m1, m2, t1);
        double tempPE = m1->oprFit(t1);
        double excessEnergy = m1->PE + m1->KE + m2->PE + m2->KE - tempPE;
        if (excessEnergy >= 0)
        {
            m1->minHitIndex = 0;
            m1->curHitIndex = 0;
            m1->PE = tempPE;
            m1->KE = excessEnergy;
            m1->clone(t1);
            m1->update();
            return true;
        }
        return false;
    }
};
