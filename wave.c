#include <math.h>

#include "wave.h"

#define RANGE(x,min,max) ((x) < (min) ? (min) : ((x) > (max) ? (max) : (x)))

/* Extrage esantioane dintr-un vector de date (obtinut de exemplu dintr-un
   fisier wave) si le depune intr-un vector de numere complexe.

   n reprezinta numarul de esantioane (si nu numarul de biti al
   dimensiunii), pentru a putea prelucra si ultima fereastra din fisier,
   care probabil nu e completa si va trebui umpluta cu 0.
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
        /* FIXME: rotunjire atunci cand renunt la bitii de precizie */
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
        /* FIXME: rotunjire atunci cand renunt la bitii de precizie */
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
