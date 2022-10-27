#include <iostream>
#include <graphics.h>
#include "pos2xyz.h"
#include "genRibbon.h"
#include "constructH4_BS2.h"
#include "surfGFph3.h"
#include "plot.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <cmath>
#include <malloc.h>
#define eps pow(2,-52)
//extern "C" int length(double xx[201]);
struct  Phonon{
    double **pos;
    char *type;
    double latticevector[3][3];
};
//void pos2xyz(double pos[3],char type[3],char filename[7]);
void checkLattice(struct Phonon lead, int ll);
double** posmerge(double** a, double** b, int la, int lb);
char* typemerge(char* a, char* b, int la, int lb);
//int size(double** a);
double** zeros(int a, int b);
double** ones(int a, int b);
double** MulMatrix(int m, int n, int k, double** A, double** B);
double** matrix_cal(double** a, double** b, char c, int la, int lb);
double** matrix_trans(double** a, int m, int n);
double** eye(int a);
double trace(double** a, int la);
int main() {
    double cutoff = 1.5;
    int eta = 1e-5, compare = 1, index = 3, withH = 0, nLeadUnitCell = 2, nCenterUnitCell = 2;
    bool xyz = 1;
    int natL, natC, natR, j = 0;
    double xx[201], i, allnat, ** hL0, ** VL01, ** VL, ** hR0, ** VR01, ** VR, * cent, * ind0, * ind1, ** Trans, ** dos, ** ldos, ** allH, ** H0, ** Gr;
    char AZ = 'a';
    struct Phonon leftlead, central, rightlead, all, * Ph;
    char* l = (char*)malloc(sizeof(char) * 5), * c = (char*)malloc(sizeof(char) * 5), * r = (char*)malloc(sizeof(char) * 5), * lcr = (char*)malloc(sizeof(char) * 7);
    for (int i = 0; i < 5; i++) {
        if (i == 0) {
            l[i] = 'L'; c[i] = 'C'; r[i] = 'R';
        }
        else if (i == 1)
            l[i] = c[i] = r[i] = '.';
        else if (i == 2)
            l[i] = c[i] = r[i] = 'x';
        else if (i == 3)
            l[i] = c[i] = r[i] = 'y';
        else if (i == 4)
            l[i] = c[i] = r[i] = 'z';
    }
    lcr[0] = 'L'; lcr[1] = 'C'; lcr[2] = 'R'; lcr[3] = '.'; lcr[4] = 'x'; lcr[5] = 'y'; lcr[6] = 'z';
    i = 0.1;
    while (i <= 2000) {
        xx[j++] = i;
        i += 10;
    }
    Ph = (struct Phonon*)malloc(sizeof(struct Phonon) * 3);
    Ph = genRibbon(AZ, index, withH, nLeadUnitCell, nCenterUnitCell);
    leftlead = Ph[0]; central = Ph[1]; rightlead = Ph[2];
    for (j = 0; j < 3; j++) {
        for (int k = 0; k < 3; k++) {
            leftlead.latticevector[j][k] = -abs(leftlead.latticevector[j][k]);
            rightlead.latticevector[j][k] = abs(rightlead.latticevector[j][k]);
        }
    }
    /*leftlead.latticevector          = - abs( leftlead.latticevector);
    rightlead.latticevector        =  abs(rightlead.latticevector);*/
  /*  if (xyz) {
        pos2xyz(leftlead.pos, leftlead.type, l);
        pos2xyz(central.pos, central.type, c);
        pos2xyz(rightlead.pos, rightlead.type, r);
        pos2xyz(posmerge(posmerge(leftlead.pos, central.pos, 3, 2), rightlead.pos, 5, 2), typemerge(typemerge(leftlead.type, central.type, 3, 2), rightlead.type, 5, 2), lcr);//pos2xyz([leftlead.pos; central.pos;rightlead.pos] ,[leftlead.type; central.type; rightlead.type],'LCR.xyz');
    }*/
    all.pos = posmerge(posmerge(leftlead.pos, central.pos, 3, 2), rightlead.pos, 5, 2);
    all.type = typemerge(typemerge(leftlead.type, central.type, 3, 2), rightlead.type, 5, 2);
    all.latticevector[0][0] = rightlead.pos[0][0] - leftlead.pos[0][0] + rightlead.latticevector[0][0];
    all.latticevector[0][1] = 0; all.latticevector[0][2] = 0;
    all.latticevector[1][0] = 0; all.latticevector[1][1] = 15; all.latticevector[1][2] = 0;
    all.latticevector[2][0] = 0; all.latticevector[2][1] = 0; all.latticevector[2][2] = 15;
    checkLattice(leftlead, 3);
    checkLattice(rightlead, 2);
    natL = size(leftlead.pos);
    natC = size(central.pos);
    natR = size(rightlead.pos);
    allH = constructH4_BS2(all, all);
    allnat = size(all.pos);
    cent = (double*)malloc(sizeof(double) * (natL * 3 + natC * 3 - natL * 3 - 1)); //cent=natL*3+1 : natL*3+natC*3;
    for (int i = 0; i < natL * 3 + natC * 3 - natL * 3 - 1; i++) {
        int j = natL * 3 + 1;
        cent[i] = j;
        j++;
    }
    H0 = (double **)malloc(sizeof(double *) * (natL * 3 + natC * 3 - natL * 3 - 1));
    for (int i = 0; i < natL * 3 + natC * 3 - natL * 3 - 1; i++) {
        H0[i] = (double *)malloc(sizeof(double) * (natL * 3 + natC * 3 - natL * 3 - 1));
    }
    for (int i = 0; i < natL * 3 + natC * 3 - natL * 3 - 1; i++)
        for (int j = 0; j < natL * 3 + natC * 3 - natL * 3 - 1; j++)
            H0[i][j] = allH[i + natL * 3 + 1][j + natL * 3 + 1];
    if (nLeadUnitCell % 2 != 0)
        std::cout << "We suppose there are 2 logical cells in V10. But it is not only 2 ucs now" << std::endl;
    ind0 = (double *)malloc(sizeof(double) * (natL * 3 - natL * 3 / 2 - 1));//ind0=natL*3/2+1:natL*3;
    ind1 = (double *)malloc(sizeof(double) * (natL * 3 / 2 - 1));//ind1=1:natL*3/2;
    for (int i = 0; i < natL * 3 - natL * 3 / 2 - 1; i++) {
        int j = natL * 3 / 2 + 1;
        ind0[i] = j;
        j++;
    }
    for (int i = 0; i < natL * 3 / 2 - 1; i++) {
        ind1[i] = i + 1;
    }
    hL0 = (double **)malloc(sizeof(double *) * (natL * 3 - natL * 3 / 2 - 1));
    for (int i = 0; i < natL * 3 - natL * 3 / 2 - 1; i++) {
        hL0[i] = (double *)malloc(sizeof(double) * (natL * 3 - natL * 3 / 2 - 1));
    }
    VL01 = (double **)malloc(sizeof(double *) * (natL * 3 - natL * 3 / 2 - 1));
    for (int i = 0; i < natL * 3 - natL * 3 / 2 - 1; i++) {
        VL01[i] = (double *)malloc(sizeof(double) * (natL * 3 / 2 - 1));
    }
    VL = (double **)malloc(sizeof(double *) * (natL * 3 - natL * 3 / 2 - 1));
    for (int i = 0; i < natL * 3 - natL * 3 / 2 - 1; i++) {
        VL[i] = (double *)malloc(sizeof(double) * (natL * 3 + natC * 3 - natL * 3 - 1));
    }
    for (int i = 0; i < natL * 3 - natL * 3 / 2 - 1; i++) {
        for (int j = 0; j < natL * 3 - natL * 3 / 2 - 1; j++)
            hL0[i][j] = allH[i + natL * 3 / 2 + 1][j + natL * 3 / 2 + 1];
        for (int j = 0; j < natL * 3 / 2 - 1; j++)
            VL01[i][j] = allH[i + natL * 3 / 2 + 1][j + 1];
        for (int j = 0; j < natL * 3 + natC * 3 - natL * 3 - 1; j++)
            VL[i][j] = allH[i + natL * 3 / 2 + 1][j + natL * 3 + 1];
    }
    //  hL0=allH(ind0, ind0);
    //   VL01=allH(ind0,ind1);
    //  VL=allH(ind0,cent);

    ind0 = (double *)malloc(sizeof(double) * ((natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1));//ind0=(natL+natC)*3+1:(natL+natC)*3+natR*3/2;
    ind1 = (double *)malloc(sizeof(double) * ((natL + natC + natR) * 3 - (natL + natC) * 3 - natR * 3 / 2 - 1));//ind1=(natL+natC)*3+natR*3/2+1:(natL+natC+natR)*3;
    for (int i = 0; i < (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1; i++) {
        int j = (natL + natC) * 3 + 1;
        ind0[i] = j;
        j++;
    }
    for (int i = 0; i < (natL + natC + natR) * 3 - (natL + natC) * 3 + natR * 3 / 2 - 1; i++) {
        int j = (natL + natC) * 3 + natR * 3 / 2 + 1;
        ind0[i] = j;
        j++;
    }
    hR0 = (double **)malloc(sizeof(double *) * ((natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1));
    for (int i = 0; i < (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1; i++) {
        hR0[i] = (double *)malloc(sizeof(double) * ((natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1));
    }
    VR01 = (double **)malloc(sizeof(double *) * ((natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1));
    for (int i = 0; i < (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1; i++) {
        VR01[i] = (double *)malloc(sizeof(double) * ((natL + natC + natR) * 3 - (natL + natC) * 3 - natR * 3 / 2 - 1));
    }
    VR = (double **)malloc(sizeof(double *) * ((natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1));
    for (int i = 0; i < (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1; i++) {
        VR[i] = (double *)malloc(sizeof(double) * (natL * 3 + natC * 3 - natL * 3 - 1));
    }
    for (int i = 0; i < (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1; i++) {
        for (int j = 0; j < natL * 3 - natL * 3 / 2 - 1; j++)
            hR0[i][j] = allH[i + (natL + natC) * 3 + 1][j + (natL + natC) * 3 + 1];
        for (int j = 0; j < (natL + natC + natR) * 3 - (natL + natC) * 3 - natR * 3 / 2 - 1; j++)
            VR01[i][j] = allH[i + (natL + natC) * 3 + 1][j + (natL + natC) * 3 + natR * 3 / 2 + 1];
        for (int j = 0; j < natL * 3 + natC * 3 - natL * 3 - 1; j++)
            VR[i][j] = allH[i + (natL + natC) * 3 + 1][j + natL * 3 + 1];
    }
    //hR0=allH(ind0,ind0);
    //VR01=allH(ind0,ind1);
    //VR=allH(ind0,cent);
    int n = 0;
    double dimC = natC * 3;
    Trans = zeros(length(xx), 1);
    ldos = zeros(length(xx), dimC);
    dos = zeros(length(xx), 1);
    int ene = 0;
    double** gL, ** selfL, ** gR, ** selfR, ** gammaL, ** gammaR;
    gammaL = (double **)malloc(sizeof(double *) * (natL * 3 + natC * 3 - natL * 3 - 1));
    for (int i = 0; i < natL * 3 + natC * 3 - natL * 3 - 1; i++)
        gammaL[i] = (double *)malloc(sizeof(double) * (natL * 3 + natC * 3 - natL * 3 - 1));
    gammaR = (double **)malloc(sizeof(double *) * (natL * 3 + natC * 3 - natL * 3 - 1));
    for (int i = 0; i < natL * 3 + natC * 3 - natL * 3 - 1; i++)
        gammaR[i] = (double *)malloc(sizeof(double) * (natL * 3 + natC * 3 - natL * 3 - 1));
    for (int i = 0; i < length(xx); i++) {
        ene = xx[i] * xx[i];
        gL = surfGFph3(ene, eta, matrix_trans(VL01, natL * 3 - natL * 3 / 2 - 1, natL * 3 / 2 - 1), VL01, hL0); //VL01''ÎªVL01µÄ×ªÖÃ
        selfL = MulMatrix(natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 - natL * 3 / 2 - 1, natL * 3 + natC * 3 - natL * 3 - 1, MulMatrix(natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 - natL * 3 / 2 - 1, natL * 3 - natL * 3 / 2 - 1, matrix_trans(VL, natL * 3 - natL * 3 / 2 - 1, natL * 3 + natC * 3 - natL * 3 - 1), gL), VL);
        gR = surfGFph3(ene,eta, matrix_trans(VR01, (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1, (natL + natC + natR) * 3 - (natL + natC) * 3 - natR * 3 / 2 - 1), VR01, hR0);
        selfR = MulMatrix(natL * 3 + natC * 3 - natL * 3 - 1, (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1, MulMatrix(natL * 3 + natC * 3 - natL * 3 - 1, (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1, (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1, matrix_trans(VR, (natL + natC) * 3 + natR * 3 / 2 - (natL + natC) * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1), gR), VR);
        gammaL = matrix_cal(selfL, matrix_trans(selfL, natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1), '-', natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1);
        gammaR = matrix_cal(selfR, matrix_trans(selfR, natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1), '-', natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1);
        /*for(int j=0;j<natL*3+natC*3-natL*3-1;j++)
            for(int k=0;k<natL*3+natC*3-natL*3-1;k++){
                gammaL[j][k]= static_cast<double>(1i * gammaL[j][k]);
                gammaR[j][k]= static_cast<double>(1i * gammaR[j][k]);
            }*/
            // gammaL = 1i * matrix_cal(selfL,matrix_trans(selfL,natL*3+natC*3-natL*3-1,natL*3+natC*3-natL*3-1),'-',natL*3+natC*3-natL*3-1,natL*3+natC*3-natL*3-1);
            // gammaR = 1i * matrix_cal(selfR,matrix_trans(selfR,natL*3+natC*3-natL*3-1,natL*3+natC*3-natL*3-1),'-',natL*3+natC*3-natL*3-1,natL*3+natC*3-natL*3-1);
        Gr = multiply(eye(size(H0)),inverse(add(add(H0, selfL), selfR), natL * 3 + natC * 3 - natL * 3 - 1));// Gr =eye(size(H0)) / ( (ene + 1i * eta)*eye(size(H0)) - add(add(H0,selfL),selfR);
        Trans[n][0] = trace(MulMatrix(natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1, MulMatrix(natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1, MulMatrix(natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1, gammaL, Gr), gammaR), matrix_trans(Gr, natL * 3 + natC * 3 - natL * 3 - 1, natL * 3 + natC * 3 - natL * 3 - 1)), natL * 3 + natC * 3 - natL * 3 - 1);
        n = n + 1;
    }
    //hold on;
    printf("Ready to plot\n");
    plot(Trans, xx);
    
    //plot(Trans,xx,'b-','DisplayName','w2'); //axis tight;xlabel('Transmission'); ylabel('\omega');
    //save trans.mat xx Trans ;

    return 0;
}
double trace(double** a, int la) {
    double b = 0;
    for (int i = 0; i < la; i++)
        b += a[i][i];
    return b;
}
double** eye(int a) {
    double** b = (double **)malloc(sizeof(double *) * a);
    for (int i = 0; i < a; i++)
        b[i] = (double *)malloc(sizeof(double) * a);
    for (int i = 0; i < a; i++)
        b[i][i] = 1;
    return b;
}
double** matrix_trans(double** a, int m, int n) {
    double** c = (double **)malloc(sizeof(double *) * n);
    for (int i = 0; i < n; i++)
        c[i] = (double *)malloc(sizeof(double) * m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            c[j][i] = a[i][j];
    return c;
}
double** MulMatrix(int m, int n, int k, double** A, double** B) {
    int i, j, p;
    double** C;
    /*C = (double**)malloc(sizeof(double*) * m);
    for (i = 0; i < m; i++) 
        C[i] = (double*)malloc(sizeof(double) * k);*/
    C = zeros(m, k);
    for (i = 0; i < m; i++) {
       // C[i] = (double *)malloc(sizeof(double) * k);
        for (j = 0; j < k; j++) {
            for (p = 0; p < n; p++) {
                C[i][j] += A[i][p] * B[p][j];
            }
        }
    }
    return C;
}
void checkLattice(struct Phonon lead, int ll) {
    double **vvec, re[3], a0 = 0, ** pos0, nat, ** pos1, * dpos, d;
    pos0 = (double **)malloc(sizeof(double *) * ll);
    vvec = (double **)malloc(sizeof(double *) * 1);
    vvec[0] = (double *)malloc(sizeof(double) * 3);
    for (int i = 0; i < ll; i++)
        pos0[i] = (double *)malloc(sizeof(double) * 3);
    for (int i = 0; i < 3; i++) {
        vvec[0][i] = lead.latticevector[0][i];
        /*for (int j = 0; j<3; j++)
            pos0[i][j] = lead.pos[i][j];*/
        vvec[0][i] = vvec[0][i] * vvec[0][i];
        a0 += sqrt(vvec[0][i]);
    }
    for (int i = 0; i < ll; i++)
        for (int j = 0; j < 3; j++)
            pos0[i][j] = lead.pos[i][j];
    if (a0 == 0) {
        printf("Check you lattice vector! It is zero!");
    }
    nat = size(pos0);
    //MulMatrix(nat, 1, 3, ones(nat, 1), reinterpret_cast<double**>(vvec));
    pos1 = add(lead.pos, MulMatrix(nat, 1, 3, ones(nat, 1), vvec));//pos1 = matrix_cal(lead.pos, MulMatrix(nat, 1, 3, ones(nat, 1), vvec), '+', ll, 3);// lead.pos + MulMatrix(nat,1,3,ones(nat,1),vvec);   //ones(nat,1) * vvec;
    dpos = (double*)malloc(sizeof(double) * 3);
    for (int i = 0; i < nat; i++) {
        for (int j = 0; j < nat; j++) {
            for (int k = 0; k < 3; k++) {
                dpos[k] = abs(pos0[i][k] - pos1[j][k]);
                if (dpos[k] < eps) {
                    printf("Lattice constants might be wrong!");
                    break;
                }
            }

        }
    }
}
double** matrix_cal(double** a, double** b, char c, int la, int lb) {
    double** d;
    d = (double **)malloc(sizeof(double *) * la);
    for (int i = 0; i < la; i++)
        d[i] = (double *)malloc(sizeof(double) * lb);
    if (c == '+') {
        for (int i = 0; i < la; i++)
            for (int j = 0; j < lb; j++)
                d[i][j] = a[i][j] + b[i][j];
    }
    else if (c == '-') {
        for (int i = 0; i < la; i++)
            for (int j = 0; j < lb; j++)
                d[i][j] = a[i][j] - b[i][j];
    }
    return d;
}
double** posmerge(double** a, double** b, int la, int lb) {
    double** c = (double **)malloc(sizeof(double *) * (la + lb));
    for (int i = 0; i < la + lb; i++) {
        c[i] = (double *)malloc(sizeof(double) * 3);
    }
    for (int i = 0; i < la; i++)
        for (int j = 0; j < 3; j++)
            c[i][j] = a[i][j];
    for (int i = 0; i < lb; i++)
        for (int j = 0; j < 3; j++)
            c[i + la][j] = b[i][j];
    return c;
}
/*int size(double** a) {
    return sizeof(a) / sizeof(double) / 3;
}*/
char* typemerge(char* a, char* b, int la, int lb) {
    char* c = (char *)malloc(sizeof(char) * (la + lb));
    for (int i = 0; i < la; i++)
        c[i] = a[i];
    for (int i = 0; i < lb; i++)
        c[i + la] = b[i];
    return c;
}
double** zeros(int a, int b) {
    double** C = (double **)malloc(sizeof(double *) * a);
    for (int i = 0; i < a; i++) {
        C[i] = (double *)malloc(sizeof(double) * b);
    }
    for (int i = 0; i < a; i++)
        for (int j = 0; j < b; j++)
            C[i][j] = 0;
    return C;
}
double** ones(int a, int b) {
    double** C = (double **)malloc(sizeof(double *) * a);
    for (int i = 0; i < a; i++) {
        C[i] = (double *)malloc(sizeof(double) * b);
    }
    for (int i = 0; i < a; i++)
        for (int j = 0; j < b; j++)
            C[i][j] = 1;
    return C;
}