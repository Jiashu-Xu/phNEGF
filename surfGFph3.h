//
// Created by xjs on 2022/10/23.
//

//#ifndef JIANHUA1_SURFGFPH3_H
//#define JIANHUA1_SURFGFPH3_H
#ifdef __cplusplus
extern "C" {
#endif
double **surfGFph3(double ene,double eta,double **V10,double **V01,double **h0);
double **multiply(double **c,double **d);
double **minus(double **c,double **d);
double **add(double **c,double **d);
double** inverse(double ** mat, int len);    
#ifdef __cplusplus
}
#endif
//#endif //JIANHUA1_SURFGFPH3_H
