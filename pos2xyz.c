//
// Created by xjs on 2022/10/22.
//

#include "pos2xyz.h"
#include <stdio.h>
#include <math.h>

   /* int size(double** a) {
        return sizeof(a)/sizeof(double)/3;
    }*/
    int length1(char *type){
        if(sizeof(type)/sizeof(char)/3>0)
            return sizeof(type)/sizeof(char)/3<3?3:sizeof(type)/sizeof(char)/3;
        else
            return 0;
    }
    int length(double xx[201]){
        return sizeof(xx)/sizeof(double);

    }
  /*  void pos2xyz(double** pos, char* type, char filename[]) {
        double nat;
        FILE *fid;
        nat=size(pos);
        fid=fopen(filename,'w'); // open file for writing; discard existing contents
        fprintf(fid,'%d\r\n',nat);
        fprintf(fid,' \r\n');

        if (length1(type)>1){
            for (int ii=0;ii<nat;ii++)
                for(int j=0;j<3;j++)
                    fprintf(fid,"%s %12.4f\r\n",type[ii],pos[ii][j]);
        }
        else
            for(int ii=0;ii<nat;ii++)
                for(int j=0;j<3;j++)
                    fprintf(fid," % s % 12.4f % 12.4f % 12.4f\r\n",type[ii],pos[ii][j]);
                fclose(fid);
    }*/

