/**********************************************************************/
/*                                                                    */
/* Definitie de tipuri si aritmetica pe numere complexe.              */
/*                                                                    */
/**********************************************************************/

#ifndef _COMPLEX_H
#define _COMPLEX_H

#include "types.h"

typedef struct _complex {
    t_s64       r;
    t_s64       i;
} t_complex;

/* Numarul de biti (inferiori) folositi pentru a reprezenta partea
   subunitara a componentelor unui numar complex */
#define COMPLEX_PRECISION 16

/* Componenta reala a produsului a 2 numere complexe */
#define cmul_r(x,y) (((x).r * (y).r - (x).i * (y).i) / COMPLEX_PRECISION)

/* Componenta imaginara a produsului a 2 numere complexe */
#define cmul_i(x,y) (((x).r * (y).i + (x).i * (y).r) / COMPLEX_PRECISION)

#endif
