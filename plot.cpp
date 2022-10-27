//
// Created by xjs on 2022/10/25.
//

#include "plot.h"
#include "pos2xyz.h"
#include <iostream>
#include <easyx.h>
void plot(double **Trans,double xx[]){
    initgraph(1000,1000,NULL);
    for(int i=0;i<length(xx);i++)
        line(Trans[i][0],xx[i],Trans[i+1][0],xx[i+1]);
    getchar();
    closegraph();
}