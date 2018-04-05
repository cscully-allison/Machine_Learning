#include "ecos.h"
#include "data.h"

__global__ void ECOSCall(idxint n,idxint m,idxint p,idxint l,idxint ncones,idxint* q, idxint e, pfloat* Gpr ,idxint* Gjc, idxint* Gir, pfloat* Apr, idxint* Ajc, idxint* Air, pfloat* c, pfloat* h, pfloat* b){
    idxint exitflag = ECOS_FATAL;
    pwork* mywork;

    /* set up data */
    mywork = ECOS_setup(n, m, p, l, ncones, q, 0, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);

    if( mywork != NULL ){

        /* solve */
        exitflag = ECOS_solve(mywork);

        /* clean up memory */
        ECOS_cleanup(mywork, 0);

    }

    /* test version number
    ECOS_ver(ver);
    printf("This test has been run on ECOS version %s\n", ver);
     */

    /* explicitly truncate exit code */
    //return (int)exitflag;

}


int main(){
    /*char ver[7];*/




    return 0;
}
