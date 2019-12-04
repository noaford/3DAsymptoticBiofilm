//
//  Grid.hpp
//  Asymptotic3DBiofilm
//
//  Created by Noah Ford on 1/16/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#ifndef Grid2D_hpp
#define Grid2D_hpp

#include <stdio.h>
#include <math.h>
#include <complex>
#include "Grid.hpp"

class Grid2D: public Grid{
private:
    int N;
    int M;
    double Lx;
    double Ly;
    double dx;
    double dy;
    double ls;
    double minwidth;
    double* x;
    double* y;
    double* ChanWidthU;
    double* ChanWidthV;
    double* ChanWidthP;
    double* FieldU;
    double* FieldV;
    double* FieldP;
    double* SoluteMax;
    double* SoluteAvg;
    double* SoluteInterface;
    double* Growth;
    double* Fluxout;
    double* HeightBiofilm;

    
    int inletStart;
    int inletEnd;

public:
        
    Grid2D(int N,int M,double dx,double dy,double Lx,double Ly,double ls):N(N),M(M),dx(dx),dy(dy),Lx(Lx),Ly(Ly),ls(ls){
        x = new double[N+1];
        y = new double[M+1];
        for(int i = 0; i<N+1;i++)
            x[i] = (double)i*dx;
        for(int j = 0; j<M+1;j++)
            y[j] = (double)j*dy;

        ChanWidthU = new double[N*(M+1)];
        FieldU = new double[N*(M+1)];
        ChanWidthV = new double[(N+1)*M];
        FieldV = new double[(N+1)*M];
        ChanWidthP = new double[(N+1)*(M+1)];
        FieldP = new double[(N+1)*(M+1)];
        SoluteMax = new double[(N+1)*(M+1)];
        SoluteAvg = new double[(N+1)*(M+1)];
        SoluteInterface = new double[(N+1)*(M+1)];
        Growth = new double[(N+1)*(M+1)];
        Fluxout = new double[(N+1)*(M+1)];
        HeightBiofilm = new double[(N+1)*(M+1)];

        minwidth = 1.;
        inletStart = floor((M+1)/3);
        inletEnd = floor(2*(M+1)/3);

    }
    
    ~Grid2D(){
        delete[] x;
        delete[] y;
        delete[] ChanWidthU;
        delete[] FieldU;
        delete[] ChanWidthV;
        delete[] FieldV;
        delete[] ChanWidthP;
        delete[] FieldP;
        delete[] SoluteMax;
        delete[] SoluteAvg;
        delete[] SoluteInterface;
        delete[] Growth;
        delete[] Fluxout;
        delete[] HeightBiofilm;
    }
    
    void InitializeChanWidth(double Ls);
    
    void InitializeSolute(double ICs);
    
    bool GrowBiofilm(double F, double Btilde, double adet, double bdet, double mu, double dt);
    
    void GetSolute(double ICs, double F, double B, double Btilde, double Jo, double Do, double Lbar, double tbar);
    
    void CalculatePressure(double U);
    
    void CalculateVelocity();

    void ScaleUp(double Ko, double Lbar, double tbar);
    
    void OutputFields(void);
    
    void HeightExtension(double* Fspeed);
    
    void solve_quartic(const std::complex<double> coefficients[5], std::complex<double> roots[4]);
    
};

#endif /* Grid2D_hpp */
