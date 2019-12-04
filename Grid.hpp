//
//  Grid.hpp
//  Asymptotic3DBiofilm
//
//  Created by Noah Ford on 1/16/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#ifndef Grid_hpp
#define Grid_hpp

#include <stdio.h>
#include <math.h>

class Grid{
private:
 /*   int N;
    double Lx;
    double dx;
    double ls;
    double* x;
    double* ChanWidthU;
    double* ChanWidthP;
    double* FieldU;
    double* FieldP;
    double* Solute;
    double* SoluteInterface;
    double* Growth;
    double* Fluxout;
    double* HeightBiofilm;*/
    
public:
    Grid(){
        
    }
    
    ~Grid(){
        
    }
    virtual void InitializeChanWidth(double Ls){};
    
    virtual void InitializeSolute(double ICs){};
    
    virtual bool GrowBiofilm(double F, double Btilde, double adet, double bdet, double mu, double dt){ return true;};
    
    virtual void GetSolute(double ICs, double F, double B, double Btilde, double Jo, double Do, double Lbar, double tbar){};
    
    virtual void CalculatePressure(double U){};
    
    virtual void CalculateVelocity(){};
    
    virtual void ScaleUp(double Ko, double Lbar, double tbar){};
    
    virtual void OutputFields(void){};
    
    virtual void HeightExtension(void){};
    
};

#endif /* Grid_hpp */
