#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "complex.h"
#include "fft.h"
#include "wave.h"
#include "wavetrans.h"

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

void trans_dummy(void *context, int n, t_complex *data) {
}

int main(int argc, char **argv) {
    /* Variabile pentru transformarea Fourier */
    int n = 4096, size;
    int *rev_map;
    t_complex *w, *data1, *data2;

    /* Variabile pentru formatul wave/pcm */
    t_riff_hdr riff_hdr;
    t_chunk tmp_chunk;
    t_wave_format_ex wave_fmt;

    int channels;
    int interleave;

    unsigned long wave_len = 0;
    t_chunk wave_fmt_chunk;

    t_complex_promote complex_promote = NULL;
    t_complex_reduce complex_reduce = NULL;

    /* Variabile pentru fluxul de date */
    FILE *in = stdin;
    FILE *out = stdout;
    unsigned long data_len = 0;

    /* Variabile de uz general */
    int status;

    /* Variabile pentru transformarile aplicate */
    int out_wav = 1;
    t_trans_f trans_f = trans_dummy;
    void *trans_ctx = NULL;
    int frame_len;

    /* Citire header RIFF */
    status = fread(&riff_hdr, sizeof(t_riff_hdr), 1, in);
    assert(status);

    if(strncmp(riff_hdr.chunk.id, "RIFF", 4)) {
        fprintf(stderr, "Input data is not RIFF\n");
        return 1;
    }

    if(strncmp(riff_hdr.file_type, "WAVE", 4)) {
        fprintf(stderr, "Unknown format '%4s'\n", riff_hdr.file_type);
        return 2;
    }

    data_len = riff_hdr.chunk.size;

    /* Restul datelor sunt chunk-uri specifice wave */
    do {
        assert(data_len >= sizeof(t_chunk));
        status = fread(&tmp_chunk, sizeof(t_chunk), 1, in);
        assert(status);
        data_len -= sizeof(t_chunk);

        if(!strncmp(tmp_chunk.id, "fmt ", 4)) {
            assert(tmp_chunk.size <= sizeof(t_wave_format_ex));
            memset(&wave_fmt, 0, sizeof(t_wave_format_ex));
            assert(data_len >= tmp_chunk.size);
            status = fread(&wave_fmt, tmp_chunk.size, 1, in);
            assert(status);
            data_len -= tmp_chunk.size;

            if(wave_fmt.format_tag != WAVE_FORMAT_PCM) {
                fprintf(stderr, "Unsupported wave format %d\n",
                        (int)wave_fmt.format_tag);
                return 3;
            }

            fprintf(stderr, "WAVE PCM format detected\n");
            fprintf(stderr, "%d bit audio samples\n",
                    (int)wave_fmt.bits_per_sample);
            fprintf(stderr, "%d channel(s), %d Hz\n", (int)wave_fmt.channels,
                    (int)wave_fmt.samples_per_sec);

            switch(wave_fmt.bits_per_sample) {
            case 8:
                complex_promote = complex_promote_u8;
                complex_reduce = complex_reduce_u8;
                // padding value = (t_u8)128
                break;
            case 16:
                complex_promote = complex_promote_s16;
                complex_reduce = complex_reduce_s16;
                // padding value = (t_s16)0
                break;
            default:
                fprintf(stderr, "Unsupported sample type!\n");
                return 4;
            }

            continue;
        }

        if(!strncmp(tmp_chunk.id, "data", 4)) {
            assert(data_len >= tmp_chunk.size);
            wave_len = tmp_chunk.size;
            break;
        }

        /* Sar peste orice alt fel de chunk (pe care nu il inteleg) */
        assert(data_len >= tmp_chunk.size);
        status = fseek(in, tmp_chunk.size, SEEK_CUR);
        assert(!status);
        data_len -= tmp_chunk.size;
    } while(1);

    /* Pregatire memorie pentru tabele si buffere */
    size = 1 << n;
    rev_map = (int *)malloc(size * sizeof(int));
    assert(rev_map != NULL);
    w = (t_complex *)malloc(size * sizeof(t_complex));
    assert(w != NULL);
    data1 = (t_complex *)malloc(size * sizeof(t_complex));
    assert(data1 != NULL);
    data2 = (t_complex *)malloc(size * sizeof(t_complex));
    assert(data2 != NULL);

    /* Calcul interleave */
    channels = wave_fmt.channels;
    interleave = channels * wave_fmt.bits_per_sample / 8;

    /* Initializare tabele */
    bit_reverse_map(n, rev_map);
    w_map(n, w);

    /* Output header riff, format wave si chunk date */
    /* FIXME calcul lungime date */
    status = fwrite(&riff_hdr, sizeof(t_riff_hdr), 1, out);
    assert(status);

    strncpy(wave_fmt_chunk.id, "fmt ", 4);
    wave_fmt_chunk.size = sizeof(t_wave_format_ex);
    status = fwrite(&wave_fmt_chunk, sizeof(t_chunk), 1, out);
    assert(status);
    status = fwrite(&wave_fmt, sizeof(t_wave_format_ex), 1, out);
    assert(status);

    strncpy(tmp_chunk.id, "data", 4);
    tmp_chunk.size = wave_len;
    status = fwrite(&tmp_chunk, sizeof(t_chunk), 1, out);
    assert(status);

    /* Pregatire transformari de framing */
    frame_len = interleave * size;


    /* Aplicarea transformarilor.

       Urmatoarele prelucrari se aplica in paralel pe toate canalele
       (pentru ca sunt intretesute in fisierul wave), dar pentru
       simplitate pseudocodul va fi scris ca pentru un singur canal.
       Toate prelucrarile la nivel de frame sunt in paralel pe canale.

       Algoritmul foloseste 2 buffere raw de dimensiune 2^n / 2
       (jumatate de frame) si 2 buffere de frame (in complex, de
       dimensiune 2^n).

       Algoritmul trebuie sa proceseze cate un frame suplimentar
       (completat cu 0) la inceputul si respectiv sfarsitul fisierului,
       pentru a nu atenua primele, respectiv ultimele 2^n / 2
       esantioane.

       Functia join(frame0, frame1) suprapune a 2-a jumatate din primul
       frame peste prima jumatate din al 2-lea frame, folosind functia
       de pondere.

       buf0 <- 0;
       read(buf1);
       frame1 <- concat(buf0, buf1);
       transform(frame1);
       do {
           buf0 <- buf1;
           read(buf1); // padding cu 0 daca e < 2^n / 2
           frame0 <- frame1;
           frame1 <- concat(buf0, buf1);
           transform(frame1);
           buf0 <- join(frame0, frame1);
           write(buf0);
       } while(length(buf1) == 2^n / 2);
       if(length(buf1) > 0) {
           buf0 <- buf1;
           buf1 <- 0;
           frame0 <- frame1;
           frame1 <- concat(buf0, buf1);
           transform(frame1);
           buf0 <- join(frame0, frame1);
           write(buf0);
       } 
     */
    
    bit_reverse(n, rev_map, test1, data1);
    dec_time_fft(n, w, data1);
    printf("Rezultat:\n");dump(size, data1);

    bit_reverse(n, rev_map, data1, data2);
    dec_time_ifft(n, w, data2);
    printf("Inversa:\n");dump(size,data2);

    free(rev_map);
    free(w);
    free(data1);
    free(data2);

    return 0;
}
