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
#include <stdio.h>
#include <stdlib.h>

#include "complex.h"
#include "fft.h"

#define TEST

#ifdef TEST
t_complex test1[8] = {
    {0x70000, 0}, {0x50000, 0}, {0x60000, 0}, {0x90000, 0},
    {0x70000, 0}, {0x50000, 0}, {0x60000, 0}, {0x90000, 0}};

t_complex test2[4] = {
    {0x70000, 0}, {0x50000, 0}, {0x60000, 0}, {0x90000, 0}};

void dump(int n, t_complex *data) {
    int i;
    double q = (double)(1 << COMPLEX_PRECISION);
    
    for(i = 0; i < n; i++, data++)
        fprintf(stderr, "%d: (%lf, %lf)\n", i, (double)data->r / q,
                (double)data->i / q);
}

void dump_rev_map(int n, int *map) {
    int i;
    
    for(i = 1 << n; i; i--, map++)
        fprintf(stderr, "%d ", *map);
    fprintf(stderr, "\n");
}
#endif

int main(void) {
    int n = 12;
    int size = 1 << n;

    t_complex *w;
    int *rev_map;
    t_complex *data, *aux;

    double fs = 44100;

    double omega = 2 * M_PI * 8000; /* 8000 Hz */

    int i;

    data = malloc(size * sizeof(t_complex));
    aux = malloc(size * sizeof(t_complex));

    w = malloc(size * sizeof(t_complex));
    rev_map = malloc(size * sizeof(int));

    for(i = 0; i < size; i++) {
        data[i].r = cos(omega * ((double)i * 1.0 / fs)) * (1 << COMPLEX_PRECISION);
        data[i].i = 0;
    }

    bit_reverse_map(n, rev_map);
    w_map(n, w);

    bit_reverse(n, rev_map, data, aux);
    dec_time_fft(n, w, aux);

    dump(size, aux);
    return 0;
}
