#include "wave.h"

/* Extrage esantioane dintr-un vector de date (obtinut de exemplu dintr-un
   fisier wave) si le depune intr-un vector de numere complexe.

   n reprezinta numarul de esantioane (si nu numarul de biti al
   dimensiunii), pentru a putea prelucra si ultima fereastra din fisier,
   care probabil nu e completa si va trebui umpluta cu 0.
 */
void fft_promote_s16(int n, int interleave, char *src, t_complex *dst) {
    int i;

    for(i = 0; i < n; i++, src += interleave, dst++) {
        dst->r = (t_s64)WAVE_s16_2_MACHINE(src) * (1 << COMPLEX_PRECISION);
        dst->i = 0;
    }
}

void fft_reduce_s16(int n, int interleave, t_complex *src, char *dst) {
    int i;

    for(i = 0; i < n; i++, src++, dst += interleave)
        WAVE_MACHINE_2_s16(dst, (t_s16)(src->r / (1 << COMPLEX_PRECISION)));
    /* FIXME: rotunjire atunci cand renunt la bitii de precizie */
}

void fft_promote_u8(int n, int interleave, char *src, t_complex *dst) {
    int i;

    for(i = 0; i < n; i++, src += interleave, dst++) {
        dst->r = ((long)(*(unsigned char *)src) - 128) * (1 << COMPLEX_PRECISION);
        dst->i = 0;
    }
}

void fft_reduce_u8(int n, int interleave, t_complex *src, char *dst) {
    int i;

    for(i = 0; i < n; i++, src++, dst += interleave)
        *(unsigned char *)src = 0; //FIXME
}
