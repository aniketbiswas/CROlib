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
#include <fstream>

using namespace std;

class CROParameter
{
public:
    int FELimit, iniPopSize, decThres;
    double iniKE, iniBuffer, collRate, lossRate, synThres;

public:
    CROParameter() : FELimit(100000), iniPopSize(10), decThres(1000),
        iniKE(0), iniBuffer(0), collRate(0.2), lossRate(0.2),
        synThres(10) {}

    CROParameter(const char* iniName) : FELimit(100000), iniPopSize(10),
        decThres(1000), iniKE(0), iniBuffer(0), collRate(0.2), lossRate(0.2),
        synThres(10)
    {
        ifstream iniFile(iniName);
        char name[16];
        double value;
        while (!iniFile.eof())
        {
            iniFile >> name >> value;
            parseConfig(name, value);
        }
        iniFile.close();
    }

    virtual ~CROParameter() {}

    virtual void parseConfig(const char* name, const double value)
    {
        if (strcmp(name, "FEL") == 0)
            FELimit = value;
        else if (strcmp(name, "IPS") == 0)
            iniPopSize = value;
        else if (strcmp(name, "IKE") == 0)
            iniKE = value;
        else if (strcmp(name, "IBF") == 0)
            iniBuffer = value;
        else if (strcmp(name, "CRT") == 0)
            collRate = value;
        else if (strcmp(name, "LRT") == 0)
            lossRate = value;
        else if (strcmp(name, "DTH") == 0)
            decThres = value;
        else if (strcmp(name, "STH") == 0)
            synThres = value;
    }
};
