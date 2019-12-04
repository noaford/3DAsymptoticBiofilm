//
//  Simulator.hpp
//  Asymptotic3DBiofilm
//
//  Created by Noah Ford on 1/16/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#ifndef Simulator_hpp
#define Simulator_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "Grid2D.hpp"
#include "Grid.hpp"
#include "inputparams.h"

enum dParamList  {eheight = 0, elength, ewidth, efuildVelocity, einitialSoluteConcentration, einitialBiofilmHeight,
    eqhat, emonodHalfRate, eadet, ebdet, eEndTime, edt
};
enum iParamList  {exGridSize = 0, eyGridSize
};

class Simulator {
private:
    double rhox;//Biomass Concentration
    double rhow;//Inactive material concentration
    double Yxo;//Yield of active biomass due to substrate consumption
    double Ywo;//Yield of EPS due to substrate consumption
    double qhat;//Maximum specific substrate utilization
    double Ko;//Half-maximum-rate concentration for utlization of substrate
    double b;//Endogenous decay rate coefficient
    double Do;//Substrate diffusion coefficient in the biofilm
    double Jo;//from variables in code in mm^2/day 2E-5*10^2*60*60*24;
    double gamma;//%Chemical oxygem demand of VS
    double fD;//Biodegradable fraction of active biomass
    
    double B;
    double W;
    double F;
    double Btilde;
    
    double Lbar;//Scaling Length
    double tbar;//Scaling Factor for time
    
    double Tend; //End time for T in DAYS
    double Tends; //Scaled term
    double dt;
    int Titersmax;
    
    double L;//Width of Biofilm (This variable rescales length variable, so as Ls increases, so will perceived flux)
    double Ls;//Scaled term
    double l;//Width of ENTIRE CHANNEL
    double ls;//Scaled term
    
    double IC;//Initial Concentration
    double ICs; //Scaled term
    
    int N; //x-direction nodes
    double Lx;//x-direction length
    int M; //y-direction nodes
    double Ly;//y-direction length;
    
    double adet;//Erosion Constants
    double bdet;
    double mu;
    
    double Lxs;
    double dx;
    double Lys;
    double dy;
    
    double Doutxs;
    
    double Uinit;
    
    bool GrowBiofilm();
    
    void Solve4Solute();
    
    void CalculateFlow();
    
    Grid* grid;
    
    
public:
    
    Simulator() {
        
        rhox = 1.025;
        rhow = 1.0125;
        Yxo = .583;
        Ywo = .477;
        qhat = 8.;
        Ko = 5.e-7;
        b = .3;
        Do = 146.88;
        Jo = 220;
        gamma = 1.42;
        fD = .8;
        W = rhow*((1.-fD)*b+Ywo*qhat)/(rhox*(Yxo*qhat-b));
        F = 1./(1.+W);
        B = rhox*(qhat+gamma*fD*b)/Ko/(Yxo*qhat-b);
        Btilde = 1.+gamma*fD*b/qhat;
        Lbar = sqrt(Ko*Do/qhat/rhox);
        tbar = 1./(Yxo*qhat-b);
        L = .1;
        Ls = L/Lbar;
        l = .5;
        ls = l/Lbar;
        IC = 6.e-6*3.5;
        ICs = IC/Ko;
        N = 20.;
        Lx = 1.;
        Lxs = Lx/Lbar;
        dx = Lxs/N;
        M = 10.;
        Ly = 1.;
        Lys = Ly/Lbar;
        dy = Lys/M;
        Doutxs = Jo/Lbar/Lbar*tbar;
        Uinit = 1.e1*2./3.*60.*60.*24./Lbar*tbar;
        Titersmax = 1.;
        

    }
    ~Simulator() {
        if(grid)
            delete grid;
    }
    
    void  InitializeVariablesFromVector(const std::vector<int> &iparams,const std::vector<double> &dparams) {
        
        rhox = 1.025;
        rhow = 1.0125;
        Yxo = .583;
        Ywo = .477;
        qhat = dparams[eqhat];
        Ko = dparams[emonodHalfRate];
        b = .3;
        Do = 146.88;
        Jo = 220;
        gamma = 1.42;
        fD = .8;
        W = rhow*((1.-fD)*b+Ywo*qhat)/(rhox*(Yxo*qhat-b));
        F = 1./(1.+W);
        B = rhox*(qhat+gamma*fD*b)/Ko/(Yxo*qhat-b);
        Btilde = 1.+gamma*fD*b/qhat;
        Lbar = sqrt(Ko*Do/qhat/rhox);
        tbar = 1./(Yxo*qhat-b);
        L = dparams[einitialBiofilmHeight];
        Ls = L/Lbar;
        l = dparams[eheight];
        ls = l/Lbar;
        IC = dparams[einitialSoluteConcentration];
        ICs = IC/Ko;
        N = iparams[exGridSize];
        Lx = dparams[elength];
        Lxs = Lx/Lbar;
        dx = Lxs/N;
        M = iparams[eyGridSize];
        Ly = dparams[ewidth];
        Lys = Ly/Lbar;
        dy = Lys/M;
        Doutxs = Jo/Lbar/Lbar*tbar;
        Uinit = dparams[efuildVelocity]*2./3.*60.*60.*24./Lbar*tbar;
        
        adet = dparams[eadet]*tbar/Lbar;//Erosion Constants
        bdet = dparams[ebdet];
        
        Tend = dparams[eEndTime];
        Tends = Tend/tbar;
        dt = dparams[edt]/tbar;
        
    }
    
    void DeleteGrid();
    
    void ReadInputFile();

    void Output();
    
    void Initialize();
        
    void InitializeSoluteandBiomass();
    
    void RunSimulation();
    
};

#endif /* Simulator_hpp */
