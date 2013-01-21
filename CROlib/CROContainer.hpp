#pragma once
#include <vector>
#include "CROMolecule.hpp"
#include "CROParameter.hpp"
#include "CROUtility.hpp"

template<class mol>
class CRO
{
public:
    double (*oprFit)(const mol*);
    mol* (*oprInit)();
    void (*oprWall)(const mol*, mol*);
    void (*oprDec)(const mol*, mol*, mol*);
    void (*oprInter)(const mol*, const mol*, mol*, mol*);
    void (*oprSyn)(const mol*, const mol*, mol*);
    CROParameter* param;
    mol* t1,* t2,* optMol;
    int curFE;
    double optGlobal, eBuffer;
    std::vector<mol*> pop;

public:
    CRO(double (*fit)(const mol*),
        mol* (*init)(),
        void (*wall)(const mol*, mol*),
        void (*dec)(const mol*, mol*, mol*),
        void (*inter)(const mol*, const mol*, mol*, mol*),
        void (*syn)(const mol*, const mol*, mol*),
        CROParameter* usrParam = NULL) :
        oprFit(fit), oprInit(init), oprWall(wall), oprDec(dec),
        oprInter(inter), oprSyn(syn), param(usrParam), optGlobal(1E300)
    {
        if (param == NULL)
            param = new CROParameter();

        t1 = oprInit();
        t2 = oprInit();
        optMol = oprInit();
        initMolValue(t1);
        initMolValue(t2);
        initMolValue(optMol);

        curFE = param->iniPopSize;
        for (int i = 0; i < param->iniPopSize; i++)
        {
            mol* m = oprInit();
            initMolValue(m);
            pop.push_back(m);
            update(m);
        }
        eBuffer = param->iniBuffer;
    };

    virtual ~CRO()
    {
        if (param != NULL)
            delete param;
        for (unsigned int i = 0; i < pop.size(); i++)
            delete pop[i];
        pop.clear();
    };

    double run()
    {
        while(curFE < param->FELimit)
        {
            if ((randDouble() > param->collRate)||(pop.size() == 1))
            {
                int pos = randInt(pop.size());
                // Decomposition
                if (pop[pos]->isInactive(param->decThres))
                {
                    mol* m = oprInit();
                    initMolValue(m);
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
                int pos1 = randInt(pop.size());
                int pos2 = pos1;
                while (pos1 == pos2)
                    pos2 = randInt(pop.size());
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
    };

private:
    void update(const mol* m)
    {
        if (m->PE < optGlobal)
        {
            optGlobal = m->PE;
            optMol->clone(m);
        }
    };

    void initMolValue(mol* m)
    {
        m->KE = param->iniKE;
        m->PE = oprFit(m);
        m->optLocal = m->PE;
    };

    void wall(mol* m)
    {
        m->curHitIndex++;
        oprWall(m, t1);
        double tempPE = oprFit(t1);
        double excessEnergy = m->PE + m->KE - tempPE;
        if (excessEnergy >= 0)
        {
            m->PE = tempPE;
            m->KE = excessEnergy * ((1.0 - param->lossRate) *
                                    randDouble() + param->lossRate);
            m->clone(t1);
            m->update();
            eBuffer += excessEnergy - m->KE;
        }
    };

    bool dec(mol* m1, mol* m2)
    {
        m1->curHitIndex++;
        oprDec(m1, t1, t2);
        double tempPE1 = oprFit(t1);
        double tempPE2 = oprFit(t2);
        double excessEnergy = m1->PE + m1->KE - tempPE1 - tempPE2;
        if ((excessEnergy >= 0) || (excessEnergy + eBuffer >= 0))
        {
            if (excessEnergy >= 0)
            {
                m1->KE = excessEnergy * randDouble();
                m2->KE = excessEnergy - m1->KE;
            }
            else
            {
                eBuffer += excessEnergy;
                m1->KE = eBuffer * randDouble() * randDouble();
                eBuffer -= m1->KE;
                m2->KE = eBuffer * randDouble() * randDouble();
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
    };

    void inter(mol* m1, mol* m2)
    {
        m1->curHitIndex++;
        m2->curHitIndex++;
        oprInter(m1, m2, t1, t2);
        double tempPE1 = oprFit(t1);
        double tempPE2 = oprFit(t2);
        double excessEnergy = m1->PE + m1->KE + m2->PE + m2->KE -
                              tempPE1 - tempPE2;
        if (excessEnergy >= 0)
        {
            m1->PE = tempPE1;
            m1->KE = excessEnergy * randDouble();
            m1->clone(t1);
            m1->update();
            m2->PE = tempPE2;
            m2->KE = excessEnergy - m1->KE;
            m2->clone(t2);
            m2->update();
        }
    };

    bool syn(mol* m1, mol* m2)
    {
        m1->curHitIndex++;
        oprSyn(m1, m2, t1);
        double tempPE = oprFit(t1);
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
    };
};
