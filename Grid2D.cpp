//
//  Grid.cpp
//  Asymptotic3DBiofilm
//
//  Created by Noah Ford on 1/16/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#include "Grid2D.hpp"
#include "LinearAlgebra.hpp"
#include "SparseMatrix.hpp"
#include "LinearSolver.hpp"
#include <lapacke.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <queue>
#include <list>

#define ind(I,J) (((I)*maxj)+J)

double pos(double x){
    if (x>=0.)
        return x;
    else
        return 0.;
}

void Grid2D::CalculatePressure(double U){
    std::cout << "Calculating Pressure Field " << std::endl;
    //This Pressure Solver Uses LAPACK
/*
    double* MP = new double[(N+1)*(M+1)*(N+1)*(M+1)](); //Matrix
    double* b = new double[(N+1)*(M+1)]();
    for(int i = 1; i<N; i++){
        for(int j = 1; j<M; j++){
            MP[((N+1)*(j)+i)*(N+1)*(M+1)+((N+1)*(j)+i)] = (-pow(ChanWidthU[j*N+i-1],3)/dx*dy -pow(ChanWidthU[j*N+i],3)/dx*dy
            -pow(ChanWidthV[(j-1)*(N+1)+i],3)/dy*dx -pow(ChanWidthV[j*(N+1)+i],3)/dy*dx);
            MP[((N+1)*(j)+i-1)*(N+1)*(M+1)+((N+1)*(j)+i)] = pow(ChanWidthU[j*N+i-1],3)/dx*dy;
            MP[((N+1)*(j)+i+1)*(N+1)*(M+1)+((N+1)*(j)+i)] = pow(ChanWidthU[j*N+i],3)/dx*dy;
            MP[((N+1)*(j)+i-N-1)*(N+1)*(M+1)+((N+1)*(j)+i)] = pow(ChanWidthV[(j-1)*(N+1)+i],3)/dy*dx;
            MP[((N+1)*(j)+i+N+1)*(N+1)*(M+1)+((N+1)*(j)+i)] = pow(ChanWidthV[j*(N+1)+i],3)/dy*dx;
            b[((N+1)*(j)+i)] = 0.;
        }
    }
    for(int j = 0; j<M+1; j++){
        //Left Boundary
        MP[((N+1)*(j))*(N+1)*(M+1)+((N+1)*(j))] = -1;
        MP[((N+1)*(j)+1)*(N+1)*(M+1)+((N+1)*(j))] = 1;
        b[((N+1)*(j))] = -12/ChanWidthU[j*N]/ChanWidthU[j*N]/ChanWidthU[j*N]*ls*dx;
        //Right Boundary
        MP[((N+1)*(j)+N)*(N+1)*(M+1)+((N+1)*(j)+N)] = 1.;
        b[((N+1)*(j)+N)] = 0.;
    }
    for(int i = 1; i<N; i++){
        //Bottom Boundary
        MP[(i)*(N+1)*(M+1)+(i)] = -1;
        MP[(N+1+i)*(N+1)*(M+1)+(i)] = 1;
        b[i] = 0.;
        //Top Boundary
        MP[((N+1)*(M)+i)*(N+1)*(M+1)+((N+1)*(M)+i)] = -1;
        MP[((N+1)*(M-1)+i)*(N+1)*(M+1)+((N+1)*(M)+i)] = 1;
        b[((N+1)*(M)+i)] = 0;
    }
    
    int numrhs = 1;
    int* pivots = new int[(N+1)*(M+1)]();
    int info = 0;
    int inputN = (N+1)*(M+1);
    int lda = (N+1)*(M+1);
    int ldb = (N+1)*(M+1);
    // dgesv_(&inputN, &numrhs, MP, &lda, pivots, b, &ldb, &info);
    for(int i = 0; i<N+1; i++){
        for(int j = 0; j<M+1; j++){
            //FieldP[j*(N+1)+i] = b[j*(N+1)+i]*U;
        }
    }
    delete[] pivots;
    delete[] b;
    delete[] MP;
    if(info==0)
        std::cout << "Pressure Calculated Successfully" << std::endl;
    else
        std::cout << "Problem in Calculating Pressure with error code: " << info << std::endl;
*/
//    int maxi = N+1;
    
    //This Pressure Solver uses PARDISO
    
    int maxj = M+1;
    bool verboseLogging = false;
    LinearAlgebra::SparseMatrix matJ((N+1)*(M+1));
    LinearAlgebra::SparseVector vecH((N+1)*(M+1));
    LinearAlgebra::SparseVector vecS((N+1)*(M+1));
    
    for(int i = 1; i<N; i++){
        for(int j = 1; j<M; j++){
            matJ(ind(i,j),ind(i,j)) = (-pow(ChanWidthU[j*N+i-1],3)/dx*dy -pow(ChanWidthU[j*N+i],3)/dx*dy
                                                           -pow(ChanWidthV[(j-1)*(N+1)+i],3)/dy*dx -pow(ChanWidthV[j*(N+1)+i],3)/dy*dx);
            matJ(ind(i,j),ind(i-1,j)) = pow(ChanWidthU[j*N+i-1],3)/dx*dy;
            matJ(ind(i,j),ind(i+1,j)) = pow(ChanWidthU[j*N+i],3)/dx*dy;
            matJ(ind(i,j),ind(i,j-1)) = pow(ChanWidthV[(j-1)*(N+1)+i],3)/dy*dx;
            matJ(ind(i,j),ind(i,j+1)) = pow(ChanWidthV[j*(N+1)+i],3)/dy*dx;
            vecH(ind(i,j)) = 0.;
        }
    }
 //Top and Bottom Boundary
    for(int i = 1; i<N; i++){
        //Bottom Boundary
        matJ(ind(i,0),ind(i,0)) = -1.;
        matJ(ind(i,0),ind(i,1))  = 1.;
        vecH(ind(i,0)) = 0.;

        //Top Boundary
        matJ(ind(i,M),ind(i,M)) = -1.;
        matJ(ind(i,M),ind(i,M-1)) = 1.;
        vecH(ind(i,M)) = 0.;
    }
    
    for(int j = 0; j<M+1; j++){
        //Left Boundary
        matJ(ind(0,j),ind(0,j)) = -1.;
        matJ(ind(0,j),ind(1,j))  = 1.;
        
        //Right Boundary
        matJ(ind(N,j),ind(N,j)) = 1.;
        vecH(ind(N,j)) = 0.;
    }
    for(int j = 0; j<inletStart; j++)
        vecH(ind(0,j)) = 0.;
    for(int j = inletStart; j<inletEnd+1; j++)
        vecH(ind(0,j)) = -12./ChanWidthU[j*N]/ChanWidthU[j*N]/ChanWidthU[j*N]*ls*dx;
    for(int j = inletEnd+1; j<M+1; j++)
        vecH(ind(0,j)) = 0.;
    
    vecH(ind(0,inletStart-1)) = -12./ChanWidthU[(inletStart-1)*N]/ChanWidthU[(inletStart-1)*N]/ChanWidthU[(inletStart-1)*N]*ls*dx/2;
    vecH(ind(0,inletEnd+1)) = -12./ChanWidthU[(inletEnd-1)*N]/ChanWidthU[(inletEnd-1)*N]/ChanWidthU[(inletEnd-1)*N]*ls*dx/2;

//Trying to smooth out inlet boundary data
/*    for(int j = 0; j<inletStart-boundary; j++)
        vecH(ind(0,j)) = 0.;
    
    for(int j = inletStart-boundary; j<inletStart; j++)
        vecH(ind(0,j)) = -12./ChanWidthU[j*N]/ChanWidthU[j*N]/ChanWidthU[j*N]*ls*dx*((j-inletStart+boundary)/boundary);
    
    for(int j = inletStart; j<inletEnd+1; j++)
        vecH(ind(0,j)) = -12./ChanWidthU[j*N]/ChanWidthU[j*N]/ChanWidthU[j*N]*ls*dx;
    
    for(int j = inletEnd+1; j<inletEnd+1+boundary; j++)
        vecH(ind(0,j)) = -12./ChanWidthU[j*N]/ChanWidthU[j*N]/ChanWidthU[j*N]*ls*dx*((inletEnd+1+boundary-j)/boundary);
    
    for(int j = inletEnd+1+boundary; j<M+1; j++)
        vecH(ind(0,j)) = -12./ChanWidthU[j*N]/ChanWidthU[j*N]/ChanWidthU[j*N]*ls*dx;
*/


    LinearAlgebra::DirectSolver::Solve( matJ, vecH, vecS, verboseLogging);
    
    
    for(int i = 0; i<N+1; i++){
        for(int j = 0; j<M+1; j++){
            FieldP[j*(N+1)+i] = vecS(ind(i,j))*U;
        }
    }
    
    
}

void Grid2D::CalculateVelocity(){
    for(int i = 0; i<N; i++){
        for(int j = 0; j<M+1; j++){
            FieldU[j*(N)+i] = -1./12./dx*(FieldP[j*(N+1)+i+1]-FieldP[j*(N+1)+i])*ChanWidthU[j*(N)+i]*ChanWidthU[j*(N)+i];
        }
    }
    
    for(int i = 0; i<N+1; i++){
        for(int j = 0; j<M; j++){
            FieldV[j*(N+1)+i] = -1./12./dy*(FieldP[(j+1)*(N+1)+i]-FieldP[j*(N+1)+i])*ChanWidthV[j*(N+1)+i]*ChanWidthV[j*(N+1)+i];
        }
    }
    
}



void Grid2D::InitializeChanWidth(double Ls){
    std::cout << "Initializing Channel Width " << std::endl;

    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++){
            //ChanWidthP[j*(N+1)+i] = ls-2*Ls*(1+sin(x[i]*M_PI/Lx)*sin(y[j]*M_PI/Ly));
            //ChanWidthP[j*(N+1)+i] = std::max(std::min(ls - 2.*sqrt(std::max(Ls*Ls-(x[i]-x[N/2])*(x[i]-x[N/2])-(y[j]-y[M/2])*(y[j]-y[M/2]),0.)),ls),1.);
            ChanWidthP[j*(N+1)+i] = std::max(std::min(ls - 2.*(std::max(Ls -pow(x[i]-x[N/2],2.)/4./Ls - pow(y[j]-y[M/4+4],2.)/4./Ls,0.) +
                                                               std::max(Ls - pow(x[i]-x[N/2],2.)/4./Ls - pow(y[j]-y[3*M/4-4],2.)/4./Ls,0.)),ls),1.);
            //ChanWidthP[j*(N+1)+i] = ls-2*Ls;
            minwidth = std::min(minwidth,ChanWidthP[j*(N+1)+i]);
        }

    for(int i = 0; i<N; i++)
        for(int j = 0; j<M+1; j++)
            ChanWidthU[j*(N)+i] = (ChanWidthP[j*(N+1)+i]+ChanWidthP[j*(N+1)+i+1])/2.;

    for(int i = 0; i<N+1; i++)
        for(int j = 0; j<M; j++)
            ChanWidthV[j*(N+1)+i] = (ChanWidthP[j*(N+1)+i]+ChanWidthP[(j+1)*(N+1)+i])/2.;
    
}

void Grid2D::InitializeSolute(double ICs){
    std::cout << "Initializing Solute " << std::endl;
    for(int i=0;i<N+1;i++)
        for(int j=0;j<M+1;j++){
            SoluteMax[j*(N+1)+i] = ICs*.5;
            SoluteAvg[j*(N+1)+i] = ICs*.5;
            SoluteInterface[j*(N+1)+i] = ICs*.5;}
}

void Grid2D::OutputFields(){
    std::cout << "Outputting Fields to File " << std::endl;
    FILE *file = fopen("data.bin", "wb");
    fwrite(&N,sizeof(int),1,file);
    fwrite(x,sizeof(double),(N+1),file);
    fwrite(&M,sizeof(int),1,file);
    fwrite(y,sizeof(double),(M+1),file);
    fwrite(SoluteMax,sizeof(double),(M+1)*(N+1),file);
    fwrite(SoluteAvg,sizeof(double),(M+1)*(N+1),file);
    fwrite(SoluteInterface,sizeof(double),(M+1)*(N+1),file);
    //fwrite(HeightBiofilm,sizeof(double),(M+1)*(N+1),file);
    fwrite(ChanWidthP,sizeof(double),(M+1)*(N+1),file);
    fwrite(FieldP,sizeof(double),(M+1)*(N+1),file);
    fwrite(FieldU,sizeof(double),(M+1)*(N),file);
    fwrite(FieldV,sizeof(double),(M)*(N+1),file);
    fwrite(Fluxout,sizeof(double),(M+1)*(N+1),file);
    fwrite(Growth,sizeof(double),(M+1)*(N+1),file);

    fclose(file);
    
}

void Grid2D::ScaleUp(double Ko, double Lbar, double tbar){
    for(int i = 0; i<N+1; i++){
        for(int j = 0; j<M+1; j++){
            SoluteMax[j*(N+1)+i] = SoluteMax[j*(N+1)+i]*Ko;
            SoluteAvg[j*(N+1)+i] = SoluteAvg[j*(N+1)+i]*Ko;
            SoluteInterface[j*(N+1)+i] = SoluteInterface[j*(N+1)+i]*Ko;
            Growth[j*(N+1)+i] = Growth[j*(N+1)+i]/tbar*Lbar;
            Fluxout[j*(N+1)+i] = -Fluxout[j*(N+1)+i]/(tbar/Lbar)*Ko;
            ChanWidthP[j*(N+1)+i] = ChanWidthP[j*(N+1)+i]*Lbar;
            HeightBiofilm[j*(N+1)+i] = HeightBiofilm[j*(N+1)+i]*Lbar;
        }
    }
    for(int i = 0; i<N; i++)
        for(int j = 0; j<M+1; j++)
            FieldU[j*N+i] = FieldU[j*N+i] *Lbar/tbar/24/60/60;
    
    for(int i = 0; i<N+1; i++)
        for(int j = 0; j<M; j++)
            FieldU[j*(N+1)+i] = FieldU[j*(N+1)+i] *Lbar/tbar/24/60/60;
    
    for(int i = 1; i<N; i++){
        Fluxout[i] = Fluxout[(N+1)+i];
        Fluxout[M*(N+1)+i] = Fluxout[(M-1)*(N+1)+i];
    }
    for(int i = 0; i<N+1; i++){
        x[i] = x[i] *Lbar;
    }
    for(int j = 0; j<M+1; j++){
        Fluxout[j*(N+1)] = Fluxout[j*(N+1)+1];
        Fluxout[j*(N+1)+N] = Fluxout[j*(N+1)+N-1];
        y[j] = y[j]*Lbar;
    }

    
}


bool Grid2D::GrowBiofilm(double F, double Btilde, double adet, double bdet, double mu, double dt){
    std::cout << "Growing Biofilm " << std::endl;
    bool minwidthlargerthanzero = true;
    double* Fspeed = new double[(N+1)*(M+1)];
    
    for(int i = 0; i<N+1; i++)
        for(int j = 0; j<M+1; j++)
            HeightBiofilm[j*(N+1)+i] = (ls-ChanWidthP[j*(N+1)+i])/2.;

    
    for(int i = 0; i<N+1; i++){
        for(int j = 0; j<M+1; j++){
            double U,V;
            if(i==0)
                U = (3.*FieldU[j*(N)] - FieldU[j*(N)+1])/2.;
            else if (i==N)
                U = (3.*FieldU[j*(N)+N-1] - FieldU[j*(N)+N-2])/2.;
            else
                U = (FieldU[j*(N)+i-1] + FieldU[j*(N)+i])/2.;
            
            if(j==0)
                V = (3.*FieldV[i] - FieldV[1*(N+1)+i])/2.;
            else if(j==M)
                V = (3.*FieldV[(M-1)*(N+1)+i] - FieldV[(M-2)*(N+1)+i])/2.;
            else
                V = (FieldV[(j-1)*(N+1)+i] + FieldV[j*(N+1)+i])/2.;
            
            double rightderiv, leftderiv, upderiv, downderiv, diradjust;
            if(i<N)
                rightderiv = (HeightBiofilm[j*(N+1)+i+1] - HeightBiofilm[j*(N+1)+i])/dx;
            else
                rightderiv = 0.;
            if(i>0)
                leftderiv = (HeightBiofilm[j*(N+1)+i] - HeightBiofilm[j*(N+1)+i-1])/dx;
            else
                leftderiv = 0.;
            if(j<M)
                upderiv = (HeightBiofilm[(j+1)*(N+1)+i] - HeightBiofilm[j*(N+1)+i])/dy;
            else
                upderiv = 0.;
            if(j>0)
                downderiv = (HeightBiofilm[j*(N+1)+i] - HeightBiofilm[(j-1)*(N+1)+i])/dy;
            else
                downderiv = 0.;
            
            double squaredgradient =std::pow(std::max(std::max(rightderiv,-leftderiv), 0.),2.)+std::pow(std::max(std::max(upderiv,-downderiv), 0.),2.);
            diradjust = sqrt(1.+squaredgradient);
            Fspeed[j*(N+1)+i] = 1./sqrt(squaredgradient);
            
            double horizontalheight = 0.;
            if(rightderiv>0. && rightderiv>-leftderiv && i<N)
                horizontalheight = HeightBiofilm[j*(N+1)+i+1];
            else if(-leftderiv>0. && i>0)
                horizontalheight = HeightBiofilm[j*(N+1)+i-1];
            else
                horizontalheight = HeightBiofilm[j*(N+1)+i];
            
            double verticalheight;
            if(upderiv>0. && upderiv>-downderiv && j<M)
                verticalheight = HeightBiofilm[(j+1)*(N+1)+i];
            else if(-downderiv>0. && j>0)
                verticalheight = HeightBiofilm[(j-1)*(N+1)+i];
            else
                verticalheight = HeightBiofilm[j*(N+1)+i];
            
            
            if(HeightBiofilm[j*(N+1)+i]>0.){
                double erosion = adet*pow(6.*mu*sqrt(U*U+V*V)/ChanWidthP[j*(N+1)+i],bdet); //Erosion

                double newheight = (HeightBiofilm[j*(N+1)+i]+verticalheight+horizontalheight)/3.;
                Growth[j*(N+1)+i] = std::max(std::min(1./sqrt(F*Btilde)*(SoluteInterface[j*(N+1)+i]+3)/(SoluteInterface[j*(N+1)+i]+4)
                                                 *log((1+SoluteInterface[j*(N+1)+i]/3.)*(1+SoluteInterface[j*(N+1)+i])),
                                                 /*HeightBiofilm[j*(N+1)+i]*/newheight*SoluteInterface[j*(N+1)+i]/(1.+SoluteInterface[j*(N+1)+i])*diradjust )
                                        - erosion,1.e-16)*diradjust;
            }else
                Growth[j*(N+1)+i] = 0.;
        }
    }
    
    HeightExtension(Fspeed); //Extend Height and growth velocities to rest of biofilm

    
    //Copy new channel width values back into ChanWidthP
    minwidth = 1.;
    for(int i = 0; i<N+1; i++){
        for(int j = 0; j<M+1; j++){
            HeightBiofilm[j*(N+1)+i] += std::max(Growth[j*(N+1)+i]*dt,0.);
            HeightBiofilm[j*(N+1)+i] = std::min(HeightBiofilm[j*(N+1)+i],ls/2.-.5);
            double width = ls - 2.*std::max(HeightBiofilm[j*(N+1)+i],0.);
            if(width<=0.){ minwidthlargerthanzero = false;
                std::cout << "Channel Width is zero at some point" << std::endl;
            }
            ChanWidthP[j*(N+1)+i] =  std::max(width,1.);
            minwidth = std::min(minwidth,width);
        }
    }
    
    //Now use linear interpolation to find values for ChanWidthU & ChanWidthV
    for(int i = 0; i<N; i++){
        for(int j = 0; j<M+1; j++){
            ChanWidthU[j*(N)+i] = (ChanWidthP[j*(N+1)+i]+ChanWidthP[j*(N+1)+i+1])/2.;
        }
    }
    for(int i = 0; i<N+1; i++){
        for(int j = 0; j<M; j++){
            ChanWidthV[j*(N+1)+i] = (ChanWidthP[j*(N+1)+i]+ChanWidthP[(j+1)*(N+1)+i])/2.;
        }
    }
    delete[] Fspeed;
    return minwidthlargerthanzero;

}

void Grid2D::GetSolute(double ICs, double F, double B, double Btilde, double Jo, double Do, double Lbar, double tbar){
    std::cout << "Calculating Solute " << std::endl;
    
    double Fo;
    double omega = .9;//std::min(1.,minwidth);
    int itersmax = 4000;
    double Ds = Jo/Lbar/Lbar*tbar;
    double lastmaxchange = 1.e8;
    double maxchange=0.;
    int iters=0;
    

for(int i = 0; i<N+1; i++)
    for(int j = 0; j<M+1; j++){
        Fo = 4*Jo/Do/ChanWidthP[j*(N+1)+i];
        SoluteInterface[j*(N+1)+i] = (SoluteMax[j*(N+1)+i]-3.)/2. - 3*sqrt(F*Btilde)/(2.*Fo)+
        .5*sqrt((SoluteMax[j*(N+1)+i]+3)*(SoluteMax[j*(N+1)+i]+3)-6*(SoluteMax[j*(N+1)+i]-3.)*sqrt(F*Btilde)/Fo+9*F*Btilde/(Fo*Fo));
    }
    
    for(iters = 0; iters<itersmax; iters++){
        for(int i = 0; i<N+1; i++){
            for(int j = 0; j<M+1; j++){
                Fo = 4*Jo/Do/ChanWidthP[j*(N+1)+i];
                SoluteInterface[j*(N+1)+i] = std::min((SoluteMax[j*(N+1)+i]-3.)/2. - 3*sqrt(F*Btilde)/(2*Fo)+
                .5*sqrt((SoluteMax[j*(N+1)+i]+3)*(SoluteMax[j*(N+1)+i]+3)-6*(SoluteMax[j*(N+1)+i]-3.)*sqrt(F*Btilde)/Fo+9*F*Btilde/(Fo*Fo)),SoluteMax[j*(N+1)+i]);
            }
        }
        
  //      for(int j = 0; j<M+1; j++)
   //         SoluteAvg[j*(N+1)] = ICs;//(2.*SoluteMax[j*(N+1)]+SoluteInterface[j*(N+1)])/3.;
        //Boundary stuff for SOR step
        
        for(int j = 1; j<M; j++){
            for(int i = 1; i<N; i++){
                Fluxout[j*(N+1)+i] = std::min(std::max(-(Jo*4*(SoluteMax[j*(N+1)+i]-SoluteInterface[j*(N+1)+i])/ChanWidthP[j*(N+1)+i]/Lbar)*tbar/Lbar,
                                              -F*B*(ls-ChanWidthP[j*(N+1)+i])/2.*SoluteInterface[j*(N+1)+i]/(1.+SoluteInterface[j*(N+1)+i]) ),0.);
                
            //    if((FieldV[j*(N+1)+i]+FieldV[(j+1)*(N+1)+i])>0)
                SoluteAvg[j*(N+1)+i] = (1.-omega)*SoluteAvg[j*(N+1)+i] +omega/(2*Ds/dx*dy+2*Ds/dy*dx+pos(FieldU[j*(N)+i])*ChanWidthU[j*(N)+i]*dy + pos(-FieldU[j*(N)+i-1])*ChanWidthU[j*(N)+i-1]*dy+ pos(FieldV[j*(N+1)+i])*ChanWidthV[j*(N+1)+i]*dx + pos(-FieldV[(j-1)*(N+1)+i])*ChanWidthV[(j-1)*(N+1)+i]*dx) * (2*Fluxout[j*(N+1)+i]*dx*dy+SoluteAvg[j*(N+1)+i+1]*(Ds/dx*dy + pos(-FieldU[j*(N)+i])*ChanWidthU[j*(N)+i]*dy)+SoluteAvg[j*(N+1)+i-1]*(Ds/dx*dy+pos(FieldU[j*(N)+i-1])*ChanWidthU[j*(N)+i-1]*dy)+SoluteAvg[(j+1)*(N+1)+i]*(Ds/dy*dx+pos(-FieldV[j*(N+1)+i])*ChanWidthV[j*(N+1)+i]*dx)+SoluteAvg[(j-1)*(N+1)+i]*(Ds/dy*dx+pos(FieldV[(j-1)*(N+1)+i])*ChanWidthV[(j-1)*(N+1)+i]*dx));
           //     else
           /*         SoluteAvg[j*(N+1)+i] = (1.-omega)*SoluteAvg[j*(N+1)+i] +omega/(2*Ds/dx*dy+2*Ds/dy*dx+FieldU[j*(N)+i]*ChanWidthU[j*(N)+i]*dy - FieldV[(j-1)*(N+1)+i]*ChanWidthV[(j-1)*(N+1)+i]*dx)*(2*Fluxout[j*(N+1)+i]*dx*dy+SoluteAvg[j*(N+1)+i+1]*Ds/dx*dy+SoluteAvg[j*(N+1)+i-1]*(Ds/dx*dy+FieldU[j*(N)+i-1]*ChanWidthU[j*(N)+i-1]*dy)+SoluteAvg[(j+1)*(N+1)+i]*(Ds/dy*dx-FieldV[j*(N+1)+i]*ChanWidthV[j*(N+1)+i]*dx)+SoluteAvg[(j-1)*(N+1)+i]*(Ds/dy*dx+0*FieldV[(j-1)*(N+1)+i]*ChanWidthV[(j-1)*(N+1)+i]*dx));*/
            }
        }
        
        for(int i = 1; i<N; i++){
            SoluteAvg[i] = SoluteAvg[i+N+1];
            SoluteAvg[M*(N+1)+i] = SoluteAvg[(M-1)*(N+1)+i];
        }
        
        for(int j = 0; j<M+1; j++){
            SoluteAvg[j*(N+1)+N] = SoluteAvg[j*(N+1)+N-1];
        }
        
        for(int j = 0; j<inletStart; j++){
            SoluteAvg[j*(N+1)] = SoluteAvg[j*(N+1)+1];
        }
        for(int j = inletStart; j<inletEnd+1; j++){
            SoluteAvg[j*(N+1)] = ICs;
        }
        for(int j = inletEnd+1; j<M+1; j++){
            SoluteAvg[j*(N+1)] = SoluteAvg[j*(N+1)+1];
        }

        
        maxchange = 0.;
        for(int i = 0; i<N+1; i++)
            for(int j = 0; j<M+1; j++){
                double change;
                double ICadjust = ICs + .01; //Maximum value in integral is ICs. add 1 for stability
                double proposedmax = .5*(3.*SoluteAvg[j*(N+1)+i] - SoluteInterface[j*(N+1)+i]);
                if(0){
                //if((proposedmax>ICadjust)&&(SoluteMax[j*(N+1)+i]>ICadjust) &&(SoluteMax[j*(N+1)+i]>SoluteInterface[j*(N+1)+i]) && iters>100){
                    double sqrtmaxminusIC = sqrt(SoluteMax[j*(N+1)+i]-ICadjust);
                    double sqrtmaxminusinterface = sqrt(SoluteMax[j*(N+1)+i]-SoluteInterface[j*(N+1)+i]);
                    double fnewton = ICadjust*sqrtmaxminusIC/sqrtmaxminusinterface +
                        1./3.*(SoluteInterface[j*(N+1)+i] + 2.*SoluteMax[j*(N+1)+i] - (ICadjust+2.*SoluteMax[j*(N+1)+i])*sqrtmaxminusIC/sqrtmaxminusinterface)
                        -SoluteAvg[j*(N+1)+i];
                    
                    double fprimenewton = -ICadjust*sqrtmaxminusIC/(3.*sqrtmaxminusinterface*(SoluteMax[j*(N+1)+i]-SoluteInterface[j*(N+1)+i]))
                        + ICadjust/(3.*sqrtmaxminusinterface*sqrtmaxminusIC)
                        + 1./3.*(2. + SoluteMax[j*(N+1)+i]*sqrtmaxminusIC/(sqrtmaxminusinterface*(SoluteMax[j*(N+1)+i]-SoluteInterface[j*(N+1)+i]))
                        -SoluteMax[j*(N+1)+i]/(sqrtmaxminusIC*sqrtmaxminusinterface)
                        -2.*sqrtmaxminusIC/sqrtmaxminusinterface);
                    /*double fprimenewton = pow(ICadjust,2.) +3*SoluteInterface[j*(N+1)+i]*SoluteMax[j*(N+1)+i] - 2.*pow(SoluteMax[j*(N+1)+i],2.)
                        +ICadjust*(-3.*SoluteInterface[j*(N+1)+i]+SoluteMax[j*(N+1)+i])
                        -2.*SoluteInterface[j*(N+1)+i]*sqrtmaxminusinterface*sqrtmaxminusIC
                        +2.*SoluteMax[j*(N+1)+i]/(3.*(SoluteMax[j*(N+1)+i]-SoluteInterface[j*(N+1)+i]));*/
                    change = -fnewton/fprimenewton*.005;
                    
                    }else{
                        
                    change = proposedmax - SoluteMax[j*(N+1)+i];
                }
                SoluteMax[j*(N+1)+i] = SoluteMax[j*(N+1)+i]+change;
                maxchange = std::max(maxchange,abs(change));
                
            }
    //    std::cout << "At iterations " << iters << " Max change was " << maxchange << std::endl;
        if((maxchange < 1e-10) && (iters>101))
            break;
        if(maxchange>lastmaxchange)
            omega /= 1.1;// std::max(omega/2.,.01);
        lastmaxchange = maxchange;
    }
    
    
    std::cout << "Solute Solver Ending on iteration " << iters << " with error " << maxchange << std::endl;
    
    
}

void Grid2D::HeightExtension(double* Fspeed){
    
    //Priority Queue to keep max value for negative biofilm height (extrapolation of biofilm surface to adjacent grid points with no biofilm)
    std::priority_queue<std::tuple<double,double,double,int>> q3;

    //Initialize all boundary grid point and put into queue
    
    for(int i=0; i<N+1; i++){
        for(int j=0;j<M+1;j++){
            
            //Check grid below for nonzero Biofilm
            if(HeightBiofilm[j*(N+1)+i] == 0. && i>0){
                
                //Check if grid point to left has biofilm
                if (HeightBiofilm[j*(N+1)+i-1]>0.){
                    double tentativegrowth = Growth[j*(N+1)+i-1];
                    double Ftentative = Fspeed[j*(N+1)+i-1];
                    double tentativeheight = HeightBiofilm[j*(N+1)+i-1]-dx/Ftentative;
                    q3.emplace(std::min(tentativeheight,-1.e-16),Ftentative,tentativegrowth,j*(N+1)+i);
                    
                    //Check if grid below has biofilm too
                    if(j>1)
                        if (HeightBiofilm[(j-1)*(N+1)+i]>0.){
                            double phii = HeightBiofilm[j*(N+1)+i-1];
                            double phij = HeightBiofilm[(j-1)*(N+1)+i];
                            double Fspeedi = Fspeed[(j)*(N+1)+i-1];
                            double Fspeedj = Fspeed[(j-1)*(N+1)+i];
                            
                            const std::complex<double> a = pow(phii-phij,2)*(dx*dx+dy*dy);
                            const std::complex<double> b = -2*pow(phii-phij,2)*(dx*dx*Fspeedj + dy*dy*Fspeedi);
                            const std::complex<double> c = pow(phii-phij,2)*(dx*dx*Fspeedj*Fspeedj + dy*dy*Fspeedi*Fspeedi) - pow(dx*dx+dy*dy,2);
                            const std::complex<double> d = 2*(dx*dx+dy*dy)*(dx*dx*Fspeedj + dy*dy*Fspeedi);
                            const std::complex<double> e = -pow(dx*dx*Fspeedj + dy*dy*Fspeedi,2);
                            
                            const std::complex<double> coefficients[5] = {e, d, c, b, a};
                            std::complex<double> r[4] = {0.,0.,0.,0.};
                            
                            solve_quartic(coefficients, r);
                            
                            double phinew = -1.e10;
                            double Fnew = 0;
                            for(int rootindx = 0; rootindx<4; rootindx++){
                                if(std::imag(r[rootindx]) <1.e-10){
                                    double Ftry = std::real(r[rootindx]);
                                    double phitry = ((dx*dx)*(Fspeedj-Ftry)*phij+dy*dy*(Fspeedi-Ftry)*phii)/(dx*dx*(Fspeedj-Ftry)+dy*dy*(Fspeedi-Ftry));
                                    if(phitry>phinew && phitry < std::min(phii,phij)){
                                        phinew = phitry;
                                        Fnew = Ftry;
                                    }
                                }
                            }
                            tentativeheight = phinew;
                            Ftentative = Fnew;
                            tentativegrowth = (Growth[j*(N+1)+i-1] + Growth[(j-1)*(N+1)+i])/2.;
                            q3.emplace(std::min(tentativeheight,-1.e-16),Ftentative,tentativegrowth,j*(N+1)+i);

                        }
                    //Check if grid above has biofilm too
                    if(j<M-1)
                        if (HeightBiofilm[(j+1)*(N+1)+i]>0.){
                            double phii = HeightBiofilm[j*(N+1)+i-1];
                            double phij = HeightBiofilm[(j+1)*(N+1)+i];
                            double Fspeedi = Fspeed[(j)*(N+1)+i-1];
                            double Fspeedj = Fspeed[(j+1)*(N+1)+i];
                            
                            const std::complex<double> a = pow(phii-phij,2)*(dx*dx+dy*dy);
                            const std::complex<double> b = -2*pow(phii-phij,2)*(dx*dx*Fspeedj + dy*dy*Fspeedi);
                            const std::complex<double> c = pow(phii-phij,2)*(dx*dx*Fspeedj*Fspeedj + dy*dy*Fspeedi*Fspeedi) - pow(dx*dx+dy*dy,2);
                            const std::complex<double> d = 2*(dx*dx+dy*dy)*(dx*dx*Fspeedj + dy*dy*Fspeedi);
                            const std::complex<double> e = -pow(dx*dx*Fspeedj + dy*dy*Fspeedi,2);
                            
                            const std::complex<double> coefficients[5] = {e, d, c, b, a};
                            std::complex<double> r[4] = {0.,0.,0.,0.};
                            
                            solve_quartic(coefficients, r);
                            
                            double phinew = -1.e10;
                            double Fnew = 0;
                            for(int rootindx = 0; rootindx<4; rootindx++){
                                if(std::imag(r[rootindx]) <1.e-10){
                                    double Ftry = std::real(r[rootindx]);
                                    double phitry = ((dx*dx)*(Fspeedj-Ftry)*phij+dy*dy*(Fspeedi-Ftry)*phii)/(dx*dx*(Fspeedj-Ftry)+dy*dy*(Fspeedi-Ftry));
                                    if(phitry>phinew && phitry < std::min(phii,phij)){
                                        phinew = phitry;
                                        Fnew = Ftry;
                                    }
                                }
                            }
                            tentativeheight = phinew;
                            Ftentative = Fnew;
                            tentativegrowth = (Growth[j*(N+1)+i-1] + Growth[(j+1)*(N+1)+i])/2.;
                            q3.emplace(std::min(tentativeheight,-1.e-16),Ftentative,tentativegrowth,j*(N+1)+i);
                        }

                }
            }
            
            
            if(HeightBiofilm[j*(N+1)+i] == 0. && i<N){
                //Check if grid point to right has biofilm
                if (HeightBiofilm[j*(N+1)+i+1]>0.){
                    double tentativegrowth = Growth[j*(N+1)+i+1];
                    double Ftentative = Fspeed[j*(N+1)+i+1];
                    double tentativeheight = HeightBiofilm[j*(N+1)+i+1]-dx/Ftentative;
                    q3.emplace(std::min(tentativeheight,-1.e-16),Ftentative,tentativegrowth,j*(N+1)+i);
                    
                    //Check if grid below has biofilm too
                    if(j>1)
                        if (HeightBiofilm[(j-1)*(N+1)+i]>0.){
                            double phii = HeightBiofilm[j*(N+1)+i+1];
                            double phij = HeightBiofilm[(j-1)*(N+1)+i];
                            double Fspeedi = Fspeed[j*(N+1)+i+1];
                            double Fspeedj = Fspeed[(j-1)*(N+1)+i];
                            
                            const std::complex<double> a = pow(phii-phij,2)*(dx*dx+dy*dy);
                            const std::complex<double> b = -2*pow(phii-phij,2)*(dx*dx*Fspeedj + dy*dy*Fspeedi);
                            const std::complex<double> c = pow(phii-phij,2)*(dx*dx*Fspeedj*Fspeedj + dy*dy*Fspeedi*Fspeedi) - pow(dx*dx+dy*dy,2);
                            const std::complex<double> d = 2*(dx*dx+dy*dy)*(dx*dx*Fspeedj + dy*dy*Fspeedi);
                            const std::complex<double> e = -pow(dx*dx*Fspeedj + dy*dy*Fspeedi,2);
                            
                            const std::complex<double> coefficients[5] = {e, d, c, b, a};
                            std::complex<double> r[4] = {0.,0.,0.,0.};
                            
                            solve_quartic(coefficients, r);
                            
                            double phinew = -1.e10;
                            double Fnew = 0;
                            for(int rootindx = 0; rootindx<4; rootindx++){
                                if(std::imag(r[rootindx]) <1.e-10){
                                    double Ftry = std::real(r[rootindx]);
                                    double phitry = ((dx*dx)*(Fspeedj-Ftry)*phij+dy*dy*(Fspeedi-Ftry)*phii)/(dx*dx*(Fspeedj-Ftry)+dy*dy*(Fspeedi-Ftry));
                                    if(phitry>phinew && phitry < std::min(phii,phij)){
                                        phinew = phitry;
                                        Fnew = Ftry;
                                    }
                                }
                            }
                            
                            tentativeheight = phinew;
                            Ftentative = Fnew;
                            tentativegrowth = (Growth[j*(N+1)+i+1]+ Growth[(j-1)*(N+1)+i])/2.;
                            q3.emplace(std::min(tentativeheight,-1.e-16),Ftentative,tentativegrowth,j*(N+1)+i);
                        }
                    
                    //Check if grid above has biofilm too
                    if(j<M-1)
                        if (HeightBiofilm[(j+1)*(N+1)+i]>0.){
                            double phii = HeightBiofilm[j*(N+1)+i+1];
                            double phij = HeightBiofilm[(j+1)*(N+1)+i];
                            double Fspeedi = Fspeed[j*(N+1)+i+1];
                            double Fspeedj = Fspeed[(j+1)*(N+1)+i];
                            
                            const std::complex<double> a = pow(phii-phij,2)*(dx*dx+dy*dy);
                            const std::complex<double> b = -2*pow(phii-phij,2)*(dx*dx*Fspeedj + dy*dy*Fspeedi);
                            const std::complex<double> c = pow(phii-phij,2)*(dx*dx*Fspeedj*Fspeedj + dy*dy*Fspeedi*Fspeedi) - pow(dx*dx+dy*dy,2);
                            const std::complex<double> d = 2*(dx*dx+dy*dy)*(dx*dx*Fspeedj + dy*dy*Fspeedi);
                            const std::complex<double> e = -pow(dx*dx*Fspeedj + dy*dy*Fspeedi,2);
                            
                            const std::complex<double> coefficients[5] = {e, d, c, b, a};
                            std::complex<double> r[4] = {0.,0.,0.,0.};
                            
                            solve_quartic(coefficients, r);
                            
                            double phinew = -1.e10;
                            double Fnew = 0;
                            for(int rootindx = 0; rootindx<4; rootindx++){
                                if(std::imag(r[rootindx]) <1.e-10){
                                    double Ftry = std::real(r[rootindx]);
                                    double phitry = ((dx*dx)*(Fspeedj-Ftry)*phij+dy*dy*(Fspeedi-Ftry)*phii)/(dx*dx*(Fspeedj-Ftry)+dy*dy*(Fspeedi-Ftry));
                                    if(phitry>phinew && phitry < std::min(phii,phij)){
                                        phinew = phitry;
                                        Fnew = Ftry;
                                    }
                                }
                            }
                            tentativeheight = phinew;
                            Ftentative = Fnew;
                            tentativegrowth = (Growth[j*(N+1)+i+1]+ Growth[(j+1)*(N+1)+i])/2.;
                            q3.emplace(std::min(tentativeheight,-1.e-16),Ftentative,tentativegrowth,j*(N+1)+i);
                        }
                }
            }
            
            //Check grid below for nonzero Biofilm and put in these numbers into queue now
            if(HeightBiofilm[j*(N+1)+i] == 0. && j>0){
                if (HeightBiofilm[(j-1)*(N+1)+i]>0.){
                    double tentativegrowth = Growth[(j-1)*(N+1)+i];
                    double Ftentative = Fspeed[(j-1)*(N+1)+i];
                    double tentativeheight = HeightBiofilm[(j-1)*(N+1)+i]-dy/Ftentative;
                    q3.emplace(tentativeheight,Ftentative,tentativegrowth,j*(N+1)+i);
                }
            }
            
            //Check grid above for nonzero Biofilm and put in these numbers into queue now
            if(HeightBiofilm[j*(N+1)+i] == 0. && j<M){
                if (HeightBiofilm[(j+1)*(N+1)+i]>0.){
                    double tentativegrowth = Growth[(j+1)*(N+1)+i];
                    double Ftentative = Fspeed[(j+1)*(N+1)+i];
                    double tentativeheight = HeightBiofilm[(j+1)*(N+1)+i]-dy/Ftentative;
                    q3.emplace(tentativeheight,Ftentative,tentativegrowth,j*(N+1)+i);
                }
            }
            
        }
    }
    
    //Put them into HeightBiofilm grid
    
    while (!q3.empty()){
        std::tuple<double,double,double,int> point = q3.top(); //Get next largest biofilm heigh
        q3.pop();
        
        int index = std::get<3>(point);
        if(HeightBiofilm[index]>=0.){            //If we haven't already put data into grid here, but height and growth in now
            HeightBiofilm[index] = std::min(std::get<0>(point),-1.e-16);
            Fspeed[index] = std::get<1>(point);
            Growth[index] = std::get<2>(point);
            
 
        }
        
    }
    
    
    
}


// The solve_quartic routine solves the generic quartic equation:
//
//     a * x^4 + b * x^3 + c * x^2 + d * x + e == 0
//
// Usage:
//
//     solve_quartic({e, d, c, b, a}, roots).

static std::complex<double> complex_sqrt(const std::complex<double> & z)
{
    return pow(z, 1. / 2.);
}

static std::complex<double> complex_cbrt(const std::complex<double> & z)
{
    return pow(z, 1. / 3.);
}

void Grid2D::solve_quartic(const std::complex<double> coefficients[5], std::complex<double> roots[4])
{
    // The algorithm below was derived by solving the quartic in Mathematica, and simplifying the resulting expression by hand.
    
    const std::complex<double> a = coefficients[4];
    const std::complex<double> b = coefficients[3] / a;
    const std::complex<double> c = coefficients[2] / a;
    const std::complex<double> d = coefficients[1] / a;
    const std::complex<double> e = coefficients[0] / a;
    
    const std::complex<double> Q1 = c * c - 3. * b * d + 12. * e;
    const std::complex<double> Q2 = 2. * c * c * c - 9. * b * c * d + 27. * d * d + 27. * b * b * e - 72. * c * e;
    const std::complex<double> Q3 = 8. * b * c - 16. * d - 2. * b * b * b;
    const std::complex<double> Q4 = 3. * b * b - 8. * c;
    
    const std::complex<double> Q5 = complex_cbrt(Q2 / 2. + complex_sqrt(Q2 * Q2 / 4. - Q1 * Q1 * Q1));
    const std::complex<double> Q6 = (Q1 / Q5 + Q5) / 3.;
    const std::complex<double> Q7 = 2. * complex_sqrt(Q4 / 12. + Q6);
    
    roots[0] = (-b - Q7 - complex_sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)) / 4.;
    roots[1] = (-b - Q7 + complex_sqrt(4. * Q4 / 6. - 4. * Q6 - Q3 / Q7)) / 4.;
    roots[2] = (-b + Q7 - complex_sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)) / 4.;
    roots[3] = (-b + Q7 + complex_sqrt(4. * Q4 / 6. - 4. * Q6 + Q3 / Q7)) / 4.;
}


