#include "wave.h"

void fft_promote_s16(int n, int interleave, char *src, complex_t *dst) {
    int i;

    for(i = 0; i < n; i++, src += interleave, dst++) {
        dst->r = (t_s64)WAVE_s16_2_MACHINE(src) * (1 << COMPLEX_PRECISION);
        dst->i = 0;
    }
}
