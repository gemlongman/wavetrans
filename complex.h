/**********************************************************************/
/*                                                                    */
/* Definitie de tipuri si aritmetica pe numere complexe.              */
/*  定义数据类型和算法复杂。                                                                  */
/**********************************************************************/

#ifndef _COMPLEX_H
#define _COMPLEX_H

#include "types.h"

typedef struct _complex {
    t_s64       r;
    t_s64       i;
} t_complex;

/* Numarul de biti (inferiori) folositi pentru a reprezenta partea
   subunitara a componentelor unui numar complex
   使用的字节数（下游）为代表的部分组成的复杂subunitara号码
*/
#define COMPLEX_PRECISION 16

/* Componenta reala a produsului a 2 numere complexe
真正的产品组件的复杂数字2
 */
#define cmul_r(x,y) (((x).r * (y).r - (x).i * (y).i) / (1 << COMPLEX_PRECISION))

/* Componenta imaginara a produsului a 2 numere `
虚构的数字产品组件2
 */
#define cmul_i(x,y) (((x).r * (y).i + (x).i * (y).r) / (1 << COMPLEX_PRECISION))

/* Modulul unui numar complex
一个复杂的模块
 */
#define cmodulus(x) (((x).r * (x).r + (x).i * (x).i) >> COMPLEX_PRECISION)

#endif
