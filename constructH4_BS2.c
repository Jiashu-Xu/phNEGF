//
// Created by xjs on 2022/10/24.
//

#include "constructH4_BS2.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>
struct  Phonon {
    double** pos;
    char* type;
    double latticevector[3][3];
};
int size(double **a){
    return sizeof(a)/sizeof(double)/3;
}
double **zeros(int a,int b){
double **c=(double *)malloc(sizeof(double)*a);
for(int i=0;i<a;i++){
    c[i]=(double *)malloc(sizeof(double)*b);
}
for(int i=0;i<a;i++)
    for(int j=0;j<b;j++)
        c[i][j]=0;
    return c;
}
double **constructH4_BS2(struct Phonon a,struct Phonon b){
    double frc2cm1=521.470833674;
    double massC=sqrt(12.0107);
    double Vasp2cm1=frc2cm1/massC;
    double factor=pow(Vasp2cm1,2);
    double cutoff[4] = {1.4203,1.4203*sqrt(3),1.4203*2,1.4203*sqrt(7)};//1.4203*[1, sqrt(3),2, sqrt(7)];
    int tol=0.1;
    double d1,d2,d3,d4;
    d1=cutoff[0];d2=cutoff[1];d3=cutoff[2];d4=cutoff[3];
    if ((d4-d3 < 2* tol) || (d3-d2 < 2*tol ) || (d2-d1< 2*tol))
        printf("Tolerance is too big\n");
    if (( d1 >= d2)||(d2>= d3)|| (d3>=d4))
        printf("Please put cutoff in ascending order\n");
    double frc[3][4]={{factor*36.5/1.602,factor*8.8/1.602,factor*3/1.602,factor*-1.92/1.602},{factor*9.82/1.602,factor*-0.4/1.602,factor*0.15/1.602,factor*-0.58/1.602},{factor*24.5/1.602,factor*-3.23/1.602,factor*-5.25/1.602,factor*2.29/1.602}};
    double f1,f2,f3,f4,v1,v2,v3,v4,mu1,mu2,mu3,mu4;
    f1=frc[0][0];f2=frc[0][1];f3=frc[0][2];f4=frc[0][3];
    v1=frc[1][0];v2=frc[1][1];v3=frc[1][2];v4=frc[1][3];
    mu1=frc[2][0];mu2=frc[2][1];mu3=frc[2][2];mu4=frc[2][3];
    int natA= size(a.pos),natB= size(b.pos);
    double **Phi= zeros(3*natA,3*natB),**posA=a.pos,**posB=b.pos,*tmp=(double *) malloc(sizeof(double)*3),sum=0,*aa,*bb;
    char *typeA=a.type,*typeB=b.type;
    for(int na=0;na<natA;na++){
        for(int nb=0;nb<natB;nb++){
            for(int k=0;k<3;k++){
                tmp[k]= posB[nb][k]-posA[na][k];
                sum+=tmp[k]*tmp[k];
            }
            for(int k=0;k<3;k++)
                tmp[k]/=sqrt(sum);
            double ex=tmp[0],ey=tmp[1];
            if(abs(tmp[2])>0.01)
                printf("ERROR: Atoms are not in x - y plane\n");
            aa=(double *) malloc(sizeof(double)*(3*na-3*(na-1)-1));
            bb=(double *) malloc(sizeof(double)*(3*nb-3*(nb-1)-1));
            if(abs(sqrt(sum)-d4)<tol){
                double temp[3][3]={{-f4*ex*ex-mu4*ey*ey,-f4*ex*ey+mu4*ex*ey,0},{-f4*ex*ey+mu4*ex*ey,-f4*ey*ey-mu4*ex*ex,0},{0,0,-v4}};
                for(int i=0;i<3;i++)
                    for(int j=0;j<3;j++)
                        Phi[i+3*(na-1)][j+(3*nb-1)]=temp[i][j];
                    for(int i=0;i<3;i++)
                        for(int j=0;j<3;j++){
                            Phi[i+3*(na-1)][j+(3*na-1)]-=temp[i][j]/2;
                            Phi[i+3*(nb-1)][j+(3*nb-1)]-=temp[i][j]/2;
                        }
            }
            else if(abs(sqrt(sum)-d3)<tol){
                double temp[3][3]={{-f3*ex*ex-mu3*ey*ey,-f3*ex*ey+mu3*ex*ey,0},{-f3*ex*ey+mu3*ex*ey,-f3*ey*ey-mu3*ex*ex,0},{0,0,-v3}};
                for(int i=0;i<3;i++)
                    for(int j=0;j<3;j++)
                        Phi[i+3*(na-1)][j+(3*nb-1)]=temp[i][j];
                    for(int i=0;i<3;i++)
                        for(int j=0;j<3;j++){
                            Phi[i+3*(na-1)][j+(3*na-1)]-=temp[i][j]/2;
                            Phi[i+3*(nb-1)][j+(3*nb-1)]-=temp[i][j]/2;
                        }
            }
            else if(abs(sqrt(sum)-d2)<tol){
                double temp[3][3]={{-f2*ex*ex-mu2*ey*ey,-f2*ex*ey+mu2*ex*ey,0},{-f2*ex*ey+mu2*ex*ey,-f2*ey*ey-mu2*ex*ex,0},{0,0,-v2}};
                for(int i=0;i<3;i++)
                    for(int j=0;j<3;j++)
                        Phi[i+3*(na-1)][j+(3*nb-1)]=temp[i][j];
                    for(int i=0;i<3;i++)
                        for(int j=0;j<3;j++){
                            Phi[i+3*(na-1)][j+(3*na-1)]-=temp[i][j]/2;
                            Phi[i+3*(nb-1)][j+(3*nb-1)]-=temp[i][j]/2;
                        }
            }
            else if(abs(sqrt(sum)-d1)<tol){
                double temp[3][3]={{-f1*ex*ex-mu1*ey*ey,-f1*ex*ey+mu1*ex*ey,0},{-f1*ex*ey+mu1*ex*ey,-f1*ey*ey-mu1*ex*ex,0},{0,0,-v1}};
                for(int i=0;i<3;i++)
                    for(int j=0;j<3;j++)
                        Phi[i+3*(na-1)][j+(3*nb-1)]=temp[i][j];
                    for(int i=0;i<3;i++)
                        for(int j=0;j<3;j++){
                            Phi[i+3*(na-1)][j+(3*na-1)]-=temp[i][j]/2;
                            Phi[i+3*(nb-1)][j+(3*nb-1)]-=temp[i][j]/2;
                        }
            }
            else;
        }
    }

return Phi;
}
