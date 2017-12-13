#ifndef GAMMA_TRANSP_H
#define GAMMA_TRANSP_H
/*Author: Ema
This file contains the definition of the gamma parameters
as in:
[1] N. A. Bobrova and P. V. Sasorov, Plasma Phys. Rep. 19, 409 1993
This file contains also some other quantities and parameter useful for computing
quantities like thermal conductivity, electrical resistivity..
*/
/*GAMMA_ilj[i,l,j] is the parameter "capital gamma"
with subscripts i+1,j+1 and superscript l .
You find it in Ref [1]*/
// GAMMA_ilj=np.zeros((10,4,2))
// double Gamma_ijl[10][4][2]=
// GAMMA_ilj[:,:,0] = {[13/4, 2, 0, 0],
//                     [5/2, 0, 0, 0],
//                     [212.7/56, 103.2/56, 0, 0],
//                     [3/2, 0, 0, 0],
//                     [3630.33/(49*16), 123.705/49, 0, 0],
//                     [477/280, 0, 0, 0],
//                     [5866.01/(49*16), 412.69/49, 132.52/49, 0],
//                     [1.2, 1.2, 0, 0],
//                     [1, 0, 0, 0],
//                     [265.32/49, 440.94/49, 793.21/(49*4), 0]]
//
// GAMMA_ilj[:,:,1] = [[939.61/(49*16), 375.74/49, 427.68/49, 129.6/49],
//                     [32079.7/(49*64), 632.025/49, 227.7/49, 0],
//                     [7.161/49, 59.394/49, 127.728/49, 51.84/49],
//                     [42.9675/49, 100.665/49, 70.92/49, 0],
//                     [3.3201/49, 34.1064/49, 95.7888/49, 41.472/49],
//                     [4.608/49, 21.744/49, 36.432/49, 0],
//                     [0.31, 12.08/7, 5.76/7, 0],
//                     [58.752/49, 243.253/49, 294.051/49, 109.47/49],
//                     [149.4/49, 254.46/49, 465.61/(49*4), 0],
//                     [576/700, 1806/700, 1068/700, 0]]

/*GAMMA_ilj[i,l,j] is the parameter "capital gamma"
with subscripts i+1,j+1 and superscript l*/
double Gamma_jil[2][10][4] = {
                    {{13/4, 2, 0, 0},
                    {5/2, 0, 0, 0},
                    {212.7/56, 103.2/56, 0, 0},
                    {3/2, 0, 0, 0},
                    {3630.33/(49*16), 123.705/49, 0, 0},
                    {477/280, 0, 0, 0},
                    {5866.01/(49*16), 412.69/49, 132.52/49, 0},
                    {1.2, 1.2, 0, 0},
                    {1, 0, 0, 0},
                    {265.32/49, 440.94/49, 793.21/(49*4), 0}},
                    {{939.61/(49*16), 375.74/49, 427.68/49, 129.6/49},
                    {32079.7/(49*64), 632.025/49, 227.7/49, 0},
                    {7.161/49, 59.394/49, 127.728/49, 51.84/49},
                    {42.9675/49, 100.665/49, 70.92/49, 0},
                    {3.3201/49, 34.1064/49, 95.7888/49, 41.472/49},
                    {4.608/49, 21.744/49, 36.432/49, 0},
                    {0.31, 12.08/7, 5.76/7, 0},
                    {58.752/49, 243.253/49, 294.051/49, 109.47/49},
                    {149.4/49, 254.46/49, 465.61/(49*4), 0},
                    {576/700, 1806/700, 1068/700, 0}}
                  };

double Gamma_i(int i, double xe, double w);
double Gamma_ij(int i, int j, double w);
double cl_ei(double ne, double T);
double cl_en(double ioniz, double T);
// double cl_ee(double ne, double T);

#endif
