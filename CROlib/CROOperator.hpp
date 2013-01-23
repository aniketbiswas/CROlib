#pragma once

template<class mol>
class CROOperator
{
public:
    double (*oprFit)(const mol*);
    mol* (*oprInit)();
    void (*oprWall)(const mol*, mol*);
    void (*oprDec)(const mol*, mol*, mol*);
    void (*oprInter)(const mol*, const mol*, mol*, mol*);
    void (*oprSyn)(const mol*, const mol*, mol*);

public:
    CROOperator(double (*fit)(const mol*),
                mol* (*init)(),
                void (*wall)(const mol*, mol*),
                void (*dec)(const mol*, mol*, mol*),
                void (*inter)(const mol*, const mol*, mol*, mol*),
                void (*syn)(const mol*, const mol*, mol*)) :
        oprFit(fit), oprInit(init), oprWall(wall), oprDec(dec),
        oprInter(inter), oprSyn(syn) {}

    virtual ~CROOperator() {}
};

