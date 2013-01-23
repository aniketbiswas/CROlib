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

