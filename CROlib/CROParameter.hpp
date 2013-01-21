#pragma once

class CROParameter
{
public:
    int FELimit, iniPopSize, decThres;
    double iniKE, iniBuffer, collRate, lossRate, synThres;

public:
    CROParameter() : FELimit(100000), iniPopSize(10), decThres(1000),
        iniKE(0), iniBuffer(0), collRate(0.2), lossRate(0.2),
        synThres(10) {};
    virtual ~CROParameter() {};
};
