#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

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
        printf("(%lf, %lf)\n", (double)data->r / q, (double)data->i / q);
}

void dump_rev_map(int n, int *map) {
    int i;
    
    for(i = 1 << n; i; i--, map++)
        printf("%d ", *map);
    printf("\n");
}
#endif

int main(int argc, char **argv) {
    int n, size;
    int *rev_map;
    t_complex *w, *data1, *data2;

    n = 3;

    size = 1 << n;
    rev_map = (int *)malloc(size * sizeof(int));
    assert(rev_map != NULL);
    w = (t_complex *)malloc(size * sizeof(t_complex));
    assert(w != NULL);
    data1 = (t_complex *)malloc(size * sizeof(t_complex));
    assert(data1 != NULL);
    data2 = (t_complex *)malloc(size * sizeof(t_complex));
    assert(data2 != NULL);

    bit_reverse_map(n, rev_map);
    dump_rev_map(n, rev_map);
    w_map(n, w);

    bit_reverse(n, rev_map, test1, data1);
    dec_time_fft(n, w, data1);
    printf("Rezultat:\n");dump(size, data1);

    bit_reverse(n, rev_map, data1, data2);
    dec_time_ifft(n, w, data2);
    printf("Inversa:\n");dump(size,data2);

    return 0;
}
