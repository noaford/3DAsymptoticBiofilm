//
//  Simulator.cpp
//  Asymptotic3DBiofilm
//
//  Created by Noah Ford on 1/16/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#include "Simulator.hpp"
#include <lapacke.h>

using namespace levelset;

void Simulator::ReadInputFile(){
    
    char infile[] = "params.txt";
    levelset::InputParams* inputFile = new levelset::InputParams(infile);
    double U = inputFile->GetDoubleParam("U"); //in mm/sec. Convert below to mm/day
    N = inputFile->GetIntParam("N");
    M = inputFile->GetIntParam("M");
    
    Tend = inputFile->GetDoubleParam("Tend");
    dt = inputFile->GetDoubleParam("dt");
    
    qhat = inputFile->GetDoubleParam("qhat");
    Ko = inputFile->GetDoubleParam("Ko");
    IC = inputFile->GetDoubleParam("IC");
    
    adet = inputFile->GetDoubleParam("adet");
    bdet = inputFile->GetDoubleParam("bdet");
    
    
    W = rhow*((1-fD)*b+Ywo*qhat)/(rhox*(Yxo*qhat-b));
    F = 1/(1+W);
    Btilde = 1+gamma*fD*b/qhat;
    Lbar = sqrt(Ko*Do/qhat/rhox);
    tbar = 1/(Yxo*qhat-b);
    dt = dt/tbar;
    Ls = L/Lbar;
    ls = l/Lbar;
    ICs = IC/Ko;
    Lxs = Lx/Lbar;
    dx = Lxs/N;
    Lys = Ly/Lbar;
    dy = Lys/M;
    Doutxs = Jo/Lbar/Lbar*tbar;
    Uinit = U*2/3*60*60*24/Lbar*tbar;
    Tends = Tend/tbar;
    Titersmax = (int)round(Tends/dt);
    adet = adet*tbar/Lbar;
    mu = .001/60/60/24/tbar; //Rescale time dimension of mu to match with speed


    
    std::cout << "Size of Grid is " << N << " by " << M << std::endl;
    std::cout << "Initial speed is " << Uinit << std::endl;
    std::cout << "End time is " << Tend << " hours" << std::endl;

}

void Simulator::Initialize(){
    grid = new Grid2D(N,M,dx,dy,Lxs,Lys,ls);
    grid->InitializeChanWidth(Ls);
    grid->InitializeSolute(ICs);
}


void Simulator::DeleteGrid(){
    delete grid;
}


void Simulator::InitializeSoluteandBiomass(){
    grid->InitializeChanWidth(Ls);
    grid->InitializeSolute(ICs);
}

void Simulator::RunSimulation(){
    bool continuesimulation = true;
    for(int i = 0; i<Titersmax;i++){//Change to allow for full iterations
        CalculateFlow();
        Solve4Solute();
        continuesimulation = GrowBiofilm();
        if(!continuesimulation){
            std::cout << "Stopping simulation at time: " << i*dt*tbar << std::endl;
            std::cout << "Because channel is closed" << std::endl;
            break;
        }
    }
}

void Simulator::Output(){
    grid->ScaleUp(Ko,Lbar,tbar);
    grid->OutputFields();
}

void Simulator::CalculateFlow(){
    grid->CalculatePressure(Uinit);
    grid->CalculateVelocity();
}


bool Simulator::GrowBiofilm(){
    return grid->GrowBiofilm(F, Btilde, adet, bdet, mu, dt);
    
}

void Simulator::Solve4Solute(){
    grid->GetSolute(ICs,F,B,Btilde,Jo,Do,Lbar,tbar);
}

