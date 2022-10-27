//
// Created by xjs on 2022/10/23.
//

#include "genRibbon.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#define eps pow(2,-52)
struct  Phonon{
    double **pos;
    char *type;
    double latticevector[3][3];
};
struct Phonon *genRibbon(char AZ,int index,int withH,int nCenterSupercell,int nLeadSupercell);
struct Phonon genAGNR(int natC,int natH,double bondCC,double bondCH);
struct Phonon genZGNR(int natC,int natH,double bondCC,double bondCH);
double **posmerge(double **a,double **b,int la,int lb);
char *typemerge(char *a,char *b,int la,int lb);
struct Phonon *genRibbon(char AZ,int index,int withH,int nCenterSupercell,int nLeadSupercell){
    struct Phonon leftlead,central,rightlead;
    struct Phonon Ph[3],uc;
    //Ph[0]=leftlead;Ph[1]=central;Ph[2]=rightlead;
    double bondCC = 1.42020658600; //1.425; %20160622
    double bondCH = 1.08;
    printf("Using bond length aCC=%f\n",&bondCC);
    int natC,natH,nat;
    

    double lattice;
    if ( (AZ=='A') || (AZ =='a')){
        lattice = bondCC*3;
        if(withH){
            natC = 2*index;
            natH = 4;
            nat = natC+natH;
        }
        else{
            natC = 2*index;
            natH = 0;
            nat = natC+natH;
        }
        uc = genAGNR(  natC, natH, bondCC,bondCH);
    }
    else{
        lattice = bondCC*sqrt(3);
        if(withH){
            natC = 2*index;
            natH = 2;
            nat = natC+natH;
        }
        else{
            natC = 2*index;
            natH = 0;
            nat = natC+natH;
        }
        uc = genZGNR(  natC, natH, bondCC,bondCH);
    }
    for (int i = 0; i < 3; i++) {
        Ph[i].pos = (double**)malloc(sizeof(double) * (natC + natH));
        Ph[i].type = (char*)malloc(sizeof(char) * (natC + natH));
        for (int j = 0; j < natC + natH; j++)
            Ph[i].pos[j] = (double*)malloc(sizeof(double) * 3);
    }
    double **pos2=(double **)malloc(sizeof(double*)*(natC+natH));
    for (int i = 0; i < natC + natH; i++)
        pos2[i] = (double*)malloc(sizeof(double) * 3);
    for(int i=0;i<natC+natH;i++)
        uc.pos[i]=(double*)malloc(sizeof(double)*3);
    for(int i=0;i<natC+natH;i++)
        for(int j=0;j<3;j++)
            pos2[i][j]=uc.pos[i][j];

        int* oID=(int*)malloc(sizeof(int)*nat);
        for(int i=0;i<nat;i++)
            *(oID+i)=i+1;
        int dir1=1;
        int dir2=2;
        char* atype2=(char*)malloc(sizeof(char)*(natC+natH));
        for(int i=0;i<natC+natH;i++)
            atype2[i]=uc.type[i];
        for(int n=1;n<=nat;n++){
            for(int i=1;i<=nat-n;i++){
                int j=i+1;
                double a=pos2[i-1][dir1-1];
                double b=pos2[j-1][dir1-1];
                if (a - b > eps){
                    for(int k=1;k<=2;k++){
                        double temp1=pos2[i-1][k-1];
                        pos2[i-1][k-1]=pos2[j-1][k-1];
                        pos2[j-1][k-1]=temp1;
                    }

                    int temp2=oID[i-1];
                    oID[i-1]=oID[j-1];
                    oID[j-1]=temp2;

                    char temp3=atype2[i-1];
                    atype2[i-1]=atype2[j-1];
                    atype2[j-1]=temp3;

                }
                else if( abs(a-b) <= eps){ // almost at the same position along dir1
                    double c=pos2[i-1][dir2-1];
                    double d=pos2[j-1][dir2-1];
                    if ( c - d > eps ){
                        for(int k=1;k<=2;k++){
                            double temp1=pos2[i-1][k-1];
                            pos2[i-1][k-1]=pos2[j-1][k-1];
                            pos2[j-1][k-1]=temp1;
                        }
                        int temp2=oID[i-1];
                        oID[i-1]=oID[j-1];
                        oID[j-1]=temp2;
                        char temp3=atype2[i-1];
                        atype2[i-1]=atype2[j-1];
                        atype2[j-1]=temp3;
                    }
                }
                else ; //return ;//a-b < -eps
            }
        }
        for(int i=0;i<natC+natH;i++)
            for(int j=0;j<3;j++)
                uc.pos[i][j]=pos2[i][j];
            for(int i=0;i<natC+natH;i++)
                uc.type[i]=atype2[i];
            for(int i=0;i<natC+natH;i++)
                for(int j=0;j<3;j++)
                    Ph[0].pos[i][j]=uc.pos[i][j];
                for(int i=0;i<natC+natH;i++)
                    Ph[0].type[i]=uc.type[i];

                double **tmpuc=(double **)malloc(sizeof(double*)*(natC+natH));
                
                for(int i=0;i<natC+natH;i++)
                    tmpuc[i]=(double*)malloc(sizeof(double)*3);
                for(int i=0;i<natC+natH;i++)
                    for(int j=0;j<3;j++)
                        tmpuc[i][j]=uc.pos[i][j];
                    for(int i=1;i<=nLeadSupercell-1;i++){
                        for(int j=0;j<nat;j++)
                            tmpuc[j][0]=uc.pos[j][0];
                        for(int j=0;j<nat;j++)
                            tmpuc[j][0]=tmpuc[j][0]+i*lattice;
                        Ph[0].pos=posmerge(Ph[0].pos,tmpuc,nat*i,nat);    //纵向合并矩阵
                        Ph[0].type=typemerge(Ph[0].type,uc.type,nat*i,nat);
                    }
                    double max=uc.pos[0][1],min=uc.pos[0][1];
                    for(int i=1;i<nat;i++){
                        if(uc.pos[i][1]>max)max=uc.pos[i][1];
                        if(uc.pos[i][1]<min)min=uc.pos[i][1];
                    }
                    double y=max-min;
                    Ph[0].latticevector[0][0] = -lattice*nLeadSupercell;// vacuum space 15 Angstrom
                    Ph[0].latticevector[1][1] = y+15;
                    Ph[0].latticevector[2][2] = 15;

                    for(int i=1;i<=nCenterSupercell;i++){
                        for(int j=0;j<nat;j++)
                            tmpuc[j][0]=uc.pos[j][0];
                        for(int j=0;j<nat;j++)
                            tmpuc[j][0]=tmpuc[j][0]+(i+nLeadSupercell - 1)*lattice;
                        Ph[1].pos=posmerge(Ph[1].pos, tmpuc, nat*(i-1), nat);    //纵向合并矩阵
                        Ph[1].type=typemerge(Ph[1].type,uc.type,nat*(i-1), nat);
                    }
                    Ph[1].latticevector[0][0] = lattice*nCenterSupercell;// vacuum space 15 Angstrom
                    Ph[1].latticevector[1][1] = y+15;
                    Ph[1].latticevector[2][2] = 15;

                    for(int i=1;i<=nLeadSupercell;i++){
                        for(int j=0;j<nat;j++)
                            tmpuc[j][0]=uc.pos[j][0];
                        for(int j=0;j<nat;j++)
                            tmpuc[j][0]=tmpuc[j][0]+(i+nLeadSupercell +nCenterSupercell - 1)*lattice;
                        Ph[2].pos=posmerge(Ph[2].pos,tmpuc,nat*(i-1), nat);    //纵向合并矩阵
                        Ph[2].type=typemerge(Ph[2].type,uc.type,nat*(i-1), nat);
                    }
                    Ph[2].latticevector[0][0] = lattice*nLeadSupercell;// vacuum space 15 Angstrom
                    Ph[2].latticevector[1][1] = y+15;
                    Ph[2].latticevector[2][2] = 15;
                    return Ph;
}



struct Phonon genAGNR(int natC,int natH,double bondCC,double bondCH){
    struct Phonon uc;
    char centerAtomType[2] ={'C','C'};
    char edgeAtomType = 'H';

    int typek0 = 0;
    int typek= typek0;

    double t1 = bondCC/2.0;
    double t2 = bondCC*sqrt(3.0)/2.0;
    double startleftdownpos[3] ={t1, 0.0 , 0.0};

    uc.pos=(double**)malloc(sizeof(double*)*(natC+natH));
    for(int i=0;i<natC+natH;i++)
        uc.pos[i]=(double*)malloc(sizeof(double)*3);
    for(int i=0;i<natC+natH;i++)
        for(int j=0;j<3;j++)
            uc.pos[i][j]=0;
        uc.type=(char*)malloc(sizeof(char)*(natC+natH));
        for(int i=0;i<natC+natH;i++)
            uc.type[i]='0';

        for(int i=1;i<=natC/2;i++){
            if (i == 1){
                for(int j=0;j<3;j++)
                    uc.pos[0][j] = startleftdownpos[j];
                uc.type[0] = centerAtomType[typek];
                typek=1-typek;
            }
            else if( i%2 == 0){
                uc.pos[i-1][0] = uc.pos[i-2][0]-t1;
                uc.pos[i-1][1] = uc.pos[i-2][1]+t2;
                uc.pos[i-1][2] = uc.pos[i-2][2];
                uc.type[i-1] = centerAtomType[typek];
                typek=1-typek;
            }
            else{ //( mod(ii,2) == 1)
                uc.pos[i-1][0] = uc.pos[i-2][0]+t1;
                uc.pos[i-1][1] = uc.pos[i-2][1]+t2;
                uc.pos[i-1][2] = uc.pos[i-2][2];
                uc.type[i-1] = centerAtomType[typek];
                typek=1-typek ;
            }
        }

        if (natH == 4){
            uc.pos[natC][0] =  uc.pos[0][0] - bondCH* 0.5;  // down
            uc.pos[natC][1] =  uc.pos[0][1] - bondCH*sqrt(3.0 )/2.0;
            uc.pos[natC][2] =  uc.pos[0][2];
            if ( natC/2%2 ==0){
                uc.pos[natC+1][0] =  uc.pos[natC/2-1][0] + bondCH* 0.5; // up
                uc.pos[natC+1][1] =  uc.pos[natC/2-1][1] + bondCH* sqrt(3.0 )/2.0;
                uc.pos[natC+1][2] =  uc.pos[natC/2-1][2];
            }
            else{
                uc.pos[natC+1][0] =  uc.pos[natC/2-1][0] - bondCH* 0.5; // up
                uc.pos[natC+1][1] =  uc.pos[natC/2-1][1] + bondCH* sqrt(3.0 )/2.0;
                uc.pos[natC+1][2] =  uc.pos[natC/2-1][2];
            }
            uc.type[natC] = edgeAtomType;
            uc.type[natC+1] = edgeAtomType;
        }

        typek = 1 - typek0;

        for(int i=1;i<=natC/2;i++){
            int j = i + natC/2-1;
            if (i == 1) {
                uc.pos[j][0] = startleftdownpos[0] + bondCC;
                uc.pos[j][1] = startleftdownpos[1];
                uc.pos[j][2] = startleftdownpos[2];
                uc.type[j] = centerAtomType[typek];
                typek=1-typek;
            }
            else if( i%2 == 0){
    uc.pos[j][0] = uc.pos[i-1][0] + 2.0 *bondCC;
    uc.pos[j][1] = uc.pos[i-1][1];
    uc.pos[j][2] = uc.pos[i-1][2];
    uc.type[j] = centerAtomType[typek];
    typek=1-typek;
}
else{ //( mod(ii,2) == 1)
    uc.pos[j][0] = uc.pos[i-1][0] + bondCC;
    uc.pos[j][1] = uc.pos[i-1][1];
    uc.pos[j][2] = uc.pos[i-1][2];
    uc.type[j] = centerAtomType[typek];
    typek=1-typek;
}
        }

        if (natH == 4){
            uc.pos[natC+2][0] =  uc.pos[natC][0] + bondCC+bondCH;
            uc.pos[natC+2][1] =  uc.pos[natC][1];
            uc.pos[natC+2][2] =  uc.pos[natC][2];
            if ( natC/2%2 ==0){
                uc.pos[natC+3][0] =  uc.pos[natC+1][0] + 2.0 *bondCC - bondCH;
                uc.pos[natC+3][1] =  uc.pos[natC+1][1];
                uc.pos[natC+3][2] =  uc.pos[natC+1][2];
            }
            else{
                uc.pos[natC+3][0] =  uc.pos[natC+1][0] + bondCC + bondCH;
                uc.pos[natC+3][1] =  uc.pos[natC+1][1];
                uc.pos[natC+3][2] =  uc.pos[natC+1][2];
                uc.type[natC+2] = edgeAtomType;
                uc.type[natC+3] = edgeAtomType;
            }
        }
        return uc;
}

struct Phonon genZGNR(int natC,int natH,double bondCC,double bondCH){
    struct Phonon uc;
    char centerAtomType[2] ={'C','C'};
    char edgeAtomType = 'H';

    int typek0 = 0;
    int typek= typek0;

    double t1 = bondCC/2.0;
    double t2 = bondCC*sqrt(3.0)/2.0;
    double startleftdownpos[3] ={t1, 0.0 , 0.0};

    uc.pos=(double**)malloc(sizeof(double*)*(natC+natH));
    uc.type = (char*)malloc(sizeof(char) * 3);
    for(int i=0;i<natC+natH;i++){
        uc.pos[i]=(double*)malloc(sizeof(double)*3);
                for(int i=0;i<natC+natH;i++)
                    for(int j=0;j<3;j++)
                        uc.pos[i][j]=0;
                    uc.type[3]="000";}

                    for(int i=1;i<=natC;i++){
                        if (i == 1){
                            uc.pos[i-1][0] = startleftdownpos[0];
                            uc.pos[i-1][1] = startleftdownpos[1];
                            uc.pos[i-1][2] = startleftdownpos[2];
                            uc.type[i-1] = centerAtomType[typek];
                            typek = 1 -  typek;
                        }
                        else if ( i%4 == 2){
                            uc.pos[i-1][0] = uc.pos[i-2][0]-t2;
                            uc.pos[i-1][1] = uc.pos[i-2][1]+t1;
                            uc.pos[i-1][2] = uc.pos[i-2][2];
                            uc.type[i-1] = centerAtomType[typek];
                            typek = 1 -  typek;
                        }
                        else if ( i%4 == 0){
                            uc.pos[i-1][0] = uc.pos[i-2][0]+t2;
                            uc.pos[i-1][1] = uc.pos[i-2][1]+t1;
                            uc.pos[i-1][2] = uc.pos[i-2][2];
                            uc.type[i-1] = centerAtomType[typek];
                            typek = 1 -  typek;
                        }
                        else{ //  mod(ii,4) == 1 or 3
                            uc.pos[i-1][0] = uc.pos[i-2][0];
                            uc.pos[i-1][1] = uc.pos[i-2][1] + bondCC;
                            uc.pos[i-1][2] = uc.pos[i-2][2];
                            uc.type[i-1] = centerAtomType[typek];
                            typek = 1 -  typek;
                        }
                    }

                    if (natH == 2){
                        uc.pos[natC][0] = startleftdownpos[0];
                        uc.pos[natC][1] = startleftdownpos[1] + bondCH;
                        uc.pos[natC][2] = startleftdownpos[2];

                        uc.pos[natC+1][0] =  uc.pos[natC-1][0];
                        uc.pos[natC+1][1] =  uc.pos[natC-1][1] + bondCH;
                        uc.pos[natC+1][2] =  uc.pos[natC-1][2];

                        uc.type[natC] = edgeAtomType;
                        uc.type[natC+1] = edgeAtomType;
                    }
                    return uc;
}


double **posmerge(double **a,double **b,int la,int lb){
    double **c=(double *)malloc(sizeof(double)*(la+lb));
    for(int i=0;i<la+lb;i++){
        c[i]=(double *)malloc(sizeof(double)*3);
    }
    for(int i=0;i<la;i++)
        for(int j=0;j<3;j++)
            c[i][j]=a[i][j];
        for(int i=0;i<lb;i++)
            for(int j=0;j<3;j++)
                c[i+la][j]=b[i][j];
            return c;
}

char *typemerge(char *a,char *b,int la,int lb){
    char *c=(char *)malloc(sizeof(char)*(la+lb));
    for(int i=0;i<la;i++)
        c[i]=a[i];
    for(int i=0;i<lb;i++)
        c[i+la]=b[i];
    return c;
}
