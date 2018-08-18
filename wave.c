/*****************************************************************************
 * wavetrans - a simple wave transformation tool                             *
 * Copyright (C) 2005 Radu Rendec <rrendec@yahoo.com>                        *
 *                                                                           *
 * This program is free software; you can redistribute it and/or modify      *
 * it under the terms of the GNU General Public License as published by      *
 * the Free Software Foundation; either version 2 of the License, or         *
 * (at your option) any later version.                                       *
 *                                                                           *
 * Foobar is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 * GNU General Public License for more details.                              *
 *                                                                           *
 * You should have received a copy of the GNU General Public License         *
 * along with Foobar; if not, write to the Free Software                     *
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA *
 *                                                                           *
 *****************************************************************************/

#include <math.h>

#include "wave.h"

#define RANGE(x,min,max) ((x) < (min) ? (min) : ((x) > (max) ? (max) : (x)))

/* Extrage esantioane dintr-un vector de date (obtinut de exemplu dintr-un fisier wave) si le depune intr-un vector de numere complexe.

   n reprezinta numarul de esantioane (si nu numarul de biti al dimensiunii), pentru a putea prelucra si ultima fereastra din fisier, care probabil nu e completa si va trebui umpluta cu 0.
   提取的样品从一个矢量数据（例如，从一个文件得到波），他们将在一个复杂的数字向量。

    n是样本数（而不是字节数的大小），和最后的窗口可以处理的文件，这可能将是不完整的，充满了0。
 */
void complex_promote_s16(int n, int interleave, char *src, t_complex *dst) {
    int i;

    for(i = 0; i < n; i++, src += interleave, dst++) {
        dst->r = (t_s64)WAVE_s16_2_MACHINE(src) * (1 << COMPLEX_PRECISION);
        dst->i = 0;
    }
}

void complex_reduce_s16(int n, int interleave, t_complex *src, char *dst) {
    int i;
    t_s64 reduced;

    for(i = 0; i < n; i++, src++, dst += interleave) {
        reduced = src->r / (1 << COMPLEX_PRECISION);
        /* FIXME: rotunjire atunci cand renunt la bitii de precizie 当我放弃舍入位精度 */
        WAVE_MACHINE_2_s16(dst, (t_s16)RANGE(reduced, -32768, 32767));
    }
}

void complex_promote_u8(int n, int interleave, char *src, t_complex *dst) {
    int i;

    for(i = 0; i < n; i++, src += interleave, dst++) {
        dst->r = ((long)(*(unsigned char *)src) - 128) * (1 << COMPLEX_PRECISION);
        dst->i = 0;
    }
}

void complex_reduce_u8(int n, int interleave, t_complex *src, char *dst) {
    int i;
    t_s64 reduced;

    for(i = 0; i < n; i++, src++, dst += interleave) {
        reduced = src->r / (1 << COMPLEX_PRECISION) + 128;
        /* FIXME: rotunjire atunci cand renunt la bitii de precizie 当我放弃舍入位精度 */
        *(unsigned char *)dst = (unsigned char)RANGE(reduced, 0, 255);
    }
}

void weight_map(int n, t_complex *data) {
    int i, max = 1 << n;

    double precision = 1 << COMPLEX_PRECISION;

    for(i = 0; i < max; i++, data++) {
        data->r = precision * (1 - cos(2.0 * M_PI * i / max)) / 2;
        data->i = 0;
    }
}
