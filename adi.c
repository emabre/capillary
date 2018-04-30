/*Functions and utilities for integration of parabolic (currently resistive
and thermal conduction) terms with the Alternating Directino Implicit algorithm*/

#include "pluto.h"
#include "structs_adi.h"

void ADI(const Data *d, Time_Step *Dts, Grid *grid) {
    static int first_call=1;
    int i,j,k, Ni, Nj;
    static Lines lines[2]; /*I define two of them as they are 1 per direction (r and z)*/

    // print1("\nI am inside the adi function\n");
    if (first_call) {
        // This part is not optimized, I don't care since it runs only once per code run
        // print1("I do the first call loop");
        first_call=0;

        // Count the lines in the domain
        Ni = CountLines(d, grid, IDIR);
        Nj = CountLines(d, grid, JDIR);

        print1("\nNumber of lines in IDIR:%d", Ni);
        print1("\nNumber of lines in JDIR:%d", Nj);
        
        InitializeLines(&lines[IDIR], Ni);
        InitializeLines(&lines[JDIR], Nj);


        // lines[0].dom_line_idx[3] = 4;

        // Find out how the domain is made
        // DOM_LOOP(k,j,i){
        //     if ((int) (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY)) {
        //         interBound[k][j][i] = 33.3; //Just a conventional number
        //     } else {
        //         interBound[k][j][i] = 0.0; //Just a conventional number
        //     }
        // }
    }
    // print1("\nlines[0].dom_line_idx[3]=%d\n", lines[0].dom_line_idx[3]);
    
}


/*Function to initialize lines, I hope this kind of initialization is ok and 
there are no problems with data continuity and similar things*/
void InitializeLines(Lines *lines, int N){

    lines->dom_line_idx = ARRAY_1D(N, int);
    lines->lidx = ARRAY_1D(N, int);
    lines->ridx = ARRAY_1D(N, int);
    lines->N = N;
    lines->lbound = ARRAY_1D(N, Bcs);
    lines->rbound = ARRAY_1D(N, Bcs);
}

int CountLines(Data *d, Grid *grid, int dir) {
    int k,j,i;
    int N=0;
    int now_on_bound, this_on_bound;

    if (dir==IDIR) {
        // I am sweeping i direction
        KDOM_LOOP(k)
            JDOM_LOOP(j){
                now_on_bound = 1;
                IDOM_LOOP(i) {
                    this_on_bound = (int) (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY);
                    if (!now_on_bound && this_on_bound)
                        N++;
                    now_on_bound = this_on_bound;
                }
                if (!now_on_bound) // this means the last line ends with the usual domain ghost cell
                    N++;
            }
    } else if (dir == JDIR) {
        KDOM_LOOP(k)
            IDOM_LOOP(i){
                now_on_bound = 1;
                JDOM_LOOP(j) {
                    this_on_bound = (int) (d->flag[k][j][i] & FLAG_INTERNAL_BOUNDARY);
                    if (!now_on_bound && this_on_bound)
                        N++;
                    now_on_bound = this_on_bound;
                }
            if (!now_on_bound) // this means the last line ends with the usual domain ghost cell
                    N++;
            }
    } else {
        print1("\nWrong choice for dir!");
        QUIT_PLUTO(1);
    }

    return N;    
}