#ifndef _WAVE_H
#define _WAVE_H

#include "types.h"
#include "complex.h"

/* FIXME urmatoarele sunt valabile doar pe i386 */

/* Obtine un t_s16 din continutul de la adresa indicata de p */
#define WAVE_s16_2_MACHINE(p) (*(t_s16 *)(p))

/* Depune un t_s16 la adresa indicata de p */
#define WAVE_MACHINE_2_s16(p,x) (*(t_s16 *)(p) = (x))

#endif
