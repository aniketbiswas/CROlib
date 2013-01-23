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

class CROParameter
{
public:
    int FELimit, iniPopSize, decThres;
    double iniKE, iniBuffer, collRate, lossRate, synThres;

public:
    CROParameter() : FELimit(100000), iniPopSize(10), decThres(1000),
        iniKE(0), iniBuffer(0), collRate(0.2), lossRate(0.2),
        synThres(10) {}
    virtual ~CROParameter() {}
};
