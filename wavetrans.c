#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "complex.h"
#include "fft.h"
#include "wave.h"
#include "wavetrans.h"
#include "types.h"

extern char *optarg;

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
        fprintf(stderr, "(%lf, %lf)\n", (double)data->r / q,
                (double)data->i / q);
}

void dump_rev_map(int n, int *map) {
    int i;
    
    for(i = 1 << n; i; i--, map++)
        fprintf(stderr, "%d ", *map);
    fprintf(stderr, "\n");
}
#endif

struct noise_reduction_context {
    t_complex *noise_profile;
};

void trans_noise_reduction(void *context, int n, t_complex *data) {
}

struct equalizer_context {
};

void trans_equalizer(void *context, int n, t_complex *data) {
}

void trans_dummy(void *context, int n, t_complex *data) {
}

void sample_noise_level(int n, FILE *in, FILE *out) {
}

int main(int argc, char **argv) {
    /* Variabile pentru transformarea Fourier */
    int n = 12, size;
    int *rev_map;
    t_complex *w;
    t_complex *weight;

    /* Variabile pentru formatul wave/pcm */
    t_riff_hdr riff_hdr;
    t_chunk tmp_chunk;
    t_wave_format_ex wave_fmt;

    int channels;
    int interleave;
    int bytes_per_sample;

    unsigned long wave_len = 0;
    t_chunk wave_fmt_chunk;

    t_complex_promote complex_promote = NULL;
    t_complex_reduce complex_reduce = NULL;

    /* Variabile pentru fluxul de date */
    FILE *in = stdin;
    FILE *out = stdout;
    FILE *noise_file = NULL;
    unsigned long data_len = 0;
    int frame_cnt = 1, frames;

    /* Variabile de uz general */
    int status;
    int i, j, opt;

    /* Variabile pentru transformarile aplicate */
    int no_wave_output = 1;
    int do_sample_noise_level = 0;
    t_trans_f trans_f = trans_dummy;
    void *trans_ctx = NULL;
    int frame_len;
    char pad_value[4];

    struct noise_reduction_context noise_reduction_context;
    struct equalizer_context equalizer_context;

    /* Buffer data */
    t_complex **fb0, **fb1, **fb_tmp, **fb_base;
    t_complex *fb_aux;
    char *buf0, *buf1, *buf_pad, *buf_tmp;
    int buf1_len;

    while((opt = getopt(argc, argv, "i:o:sfrne")) != -1) {
        switch(opt) {
        case 'i': /* Input file */
            if((in = fopen(optarg, "r")) == NULL) {
                fprintf(stderr, "Failed opening %s\n", optarg);
                return 15;
            }
            break;
        case 'o': /* Output file */
            if((out = fopen(optarg, "w")) == NULL) {
                fprintf(stderr, "Failed opening %s\n", optarg);
                return 16;
            }
            break;
        case 'f': /* FFT Size */
            if(sscanf(optarg, "%d", &n) != 1) {
                fprintf(stderr, "Invalid number for FFT Size\n");
                return 21;
            }
            if(n < 8 || n > 16) {
                fprintf(stderr, "FFT Size must be between 8 and 16\n");
                return 22;
            }
            break;
        case 's': /* Sample noise level */
            do_sample_noise_level = 1;
            break;
        case 'r': /* Perform noise reduction */
            no_wave_output = 0;
            trans_f = trans_noise_reduction;
            trans_ctx = &noise_reduction_context;
            break;
        case 'n': /* Noise level sample file */
            if((noise_file = fopen(optarg, "r")) == NULL) {
                fprintf(stderr, "Failed opening %s\n", optarg);
                return 17;
            }
            break;
        case 'e': /* Perform Parametric Equalization */
            no_wave_output = 0;
            trans_f = trans_equalizer;
            trans_ctx = &equalizer_context;
            break;
        case '?':
        case ':':
            //FIXME: usage
            return 20;
        }
    }
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

            /* Calcul interleave */
            channels = wave_fmt.channels;
            bytes_per_sample = wave_fmt.bits_per_sample / 8;
            interleave = channels * bytes_per_sample;

            switch(wave_fmt.bits_per_sample) {
            case 8:
                complex_promote = complex_promote_u8;
                complex_reduce = complex_reduce_u8;
                *(unsigned char *)pad_value = 128;
                break;
            case 16:
                complex_promote = complex_promote_s16;
                complex_reduce = complex_reduce_s16;
                *(t_s16 *)pad_value = 0;
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

    size = 1 << n;
    frame_len = interleave * size;

    if(do_sample_noise_level) {
        sample_noise_level(n, in, out);
        fclose(in);
        fclose(out);
        return 0;
    }

    if(trans_f == trans_noise_reduction) {
        // citire & verificare fisier
    }

    if(no_wave_output) {
        fclose(in);
        fclose(out);
        return 0;
    }
    /* Pregatire memorie pentru tabele si buffere */
    rev_map = (int *)malloc(size * sizeof(int));
    assert(rev_map != NULL);
    w = (t_complex *)malloc(size * sizeof(t_complex));
    assert(w != NULL);
    weight = (t_complex *)malloc(size * sizeof(t_complex));
    assert(weight != NULL);

    fb_base = (t_complex **)malloc(2 * channels * sizeof(t_complex *));
    assert(fb_base != NULL);
    for(i = 0; i < 2 * channels; i++) {
        fb_base[i] = (t_complex *)malloc(size * sizeof(t_complex));
        assert(fb0[i] != NULL);
    }
    fb0 = fb_base;
    fb1 = fb_base + channels;
    fb_aux = (t_complex *)malloc(size * sizeof(t_complex));
    assert(fb_aux != NULL);

    buf0 = (char *)malloc(frame_len / 2);
    assert(buf0 != NULL);
    buf1 = (char *)malloc(frame_len / 2);
    assert(buf1 != NULL);
    buf_pad = (char *)malloc(frame_len / 2);
    assert(buf_pad != NULL);

    /* Initializare tabele */
    bit_reverse_map(n, rev_map);
    w_map(n, w);
    weight_map(n, weight);
    for(i = size / 2 * channels, buf_tmp = buf_pad; i;
            i--, buf_tmp += bytes_per_sample)
        memcpy(buf_tmp, pad_value, bytes_per_sample);

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
       while(length(buf1) == 2^n / 2) {
           buf0 <- buf1;
           read(buf1); // padding cu 0 daca e < 2^n / 2
           frame0 <- frame1;
           frame1 <- concat(buf0, buf1);
           transform(frame1);
           buf0 <- join(frame0, frame1);
           write(buf0);
       }
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

    frames = (wave_len + frame_len - 1) / frame_len;

    memcpy(buf0, buf_pad, frame_len / 2);
    buf1_len = fread(buf1, 1, frame_len / 2, in);
    for(i = 0; i < channels; i++) {
        complex_promote(size / 2, interleave,
                buf0 + i * bytes_per_sample, fb1[i]);
        complex_promote(size / 2, interleave,
                buf1 + i * bytes_per_sample, fb1[i] + size / 2);
        
        bit_reverse(n, rev_map, fb1[i], fb_aux);
        dec_time_fft(n, w, fb_aux);
        trans_f(trans_ctx, n, fb_aux);
        bit_reverse(n, rev_map, fb_aux, fb1[i]);
        dec_time_ifft(n, w, fb1[i]);
    }
    while(buf1_len == frame_len / 2) {
        if(!(frame_cnt % 10)) {
            /* De fapt eu numar jumatati de frame-uri */
            fprintf(stderr, "\rFrame %d/%d (%d%%)", frame_cnt / 2, frames,
                    frame_cnt * 50 / frames);
            fflush(stderr);
        }
        frame_cnt++;

        buf_tmp = buf0; buf0 = buf1; buf1 = buf_tmp;
        buf1_len = fread(buf1, 1, frame_len / 2, in);
        if(buf1_len < frame_len / 2)
            memcpy(buf1 + buf1_len, buf_pad, frame_len / 2 - buf1_len);
        fb_tmp = fb0; fb0 = fb1; fb1 = fb_tmp;

        for(i = 0; i < channels; i++) {
            complex_promote(size / 2, interleave,
                    buf0 + i * bytes_per_sample, fb1[i]);
            complex_promote(size / 2, interleave,
                    buf1 + i * bytes_per_sample, fb1[i] + size / 2);

            bit_reverse(n, rev_map, fb1[i], fb_aux);
            dec_time_fft(n, w, fb_aux);
            trans_f(trans_ctx, n, fb_aux);
            bit_reverse(n, rev_map, fb_aux, fb1[i]);
            dec_time_ifft(n, w, fb1[i]);

            for(j = 0; j < size / 2; j++) {
                fb_aux[j].r =
                    (fb0[i][size / 2 + j].r * weight[size / 2 + j].r +
                     fb1[i][j].r * weight[j].r) / (1 << COMPLEX_PRECISION);
            }
            complex_reduce(size / 2, interleave, fb_aux,
                    buf0 + i * bytes_per_sample);
        }
        status = fwrite(buf0, frame_len / 2, 1, out);
        assert(status);
    }
    if(buf1_len) {
        buf_tmp = buf0; buf0 = buf1; buf1 = buf_tmp;
        memcpy(buf1, buf_pad, frame_len / 2);
        fb_tmp = fb0; fb0 = fb1; fb1 = fb_tmp;

        for(i = 0; i < channels; i++) {
            complex_promote(size / 2, interleave,
                    buf0 + i * bytes_per_sample, fb1[i]);
            complex_promote(size / 2, interleave,
                    buf1 + i * bytes_per_sample, fb1[i] + size / 2);

            bit_reverse(n, rev_map, fb1[i], fb_aux);
            dec_time_fft(n, w, fb_aux);
            trans_f(trans_ctx, n, fb_aux);
            bit_reverse(n, rev_map, fb_aux, fb1[i]);
            dec_time_ifft(n, w, fb1[i]);

            for(j = 0; j < size / 2; j++) {
                fb_aux[j].r =
                    (fb0[i][size / 2 + j].r * weight[size / 2 + j].r +
                     fb1[i][j].r * weight[j].r) / (1 << COMPLEX_PRECISION);
            }
            complex_reduce(size / 2, interleave, fb_aux,
                    buf0 + i * bytes_per_sample);
        }
        status = fwrite(buf0, buf1_len, 1, out);
        assert(status);
    }
    fprintf(stderr, "\rFrame %d/%d (100%%)\n", frames, frames);
    fflush(stderr);

    
    /* Eliberare memorie */
    free(rev_map);
    free(w);
    free(weight_map);
    for(i = 0; i < 2 * channels; i++)
        free(fb_base[i]);
    free(fb_base);
    free(fb_aux);
    free(buf0);
    free(buf1);

    fclose(in);
    fclose(out);
    return 0;
}
