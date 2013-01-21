#pragma once

class CROMolecule
{
public:
    double PE, KE, optLocal;
    int minHitIndex, curHitIndex;

public:
    CROMolecule()
    {
        minHitIndex = 0;
        curHitIndex = 0;
    };

    virtual ~CROMolecule() {};

    void update()
    {
        if (PE < optLocal)
        {
            optLocal = PE;
            minHitIndex = curHitIndex;
        }
    };

    bool isInactive(int threshold)
    {
        return (curHitIndex - minHitIndex) > threshold;
    };

    void clone(const CROMolecule* source) {};
};
