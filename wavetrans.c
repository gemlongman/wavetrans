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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "complex.h"
#include "fft.h"
#include "wave.h"
#include "wavetrans.h"
#include "types.h"

extern char *optarg;

struct noise_reduction_context {
    long double **noise;
    long double alfa, beta;
};

void trans_noise_reduction(void *context, int n, int channel, t_complex *data) {
    struct noise_reduction_context *ctx = context;
    int size = 1 << n;
    int i;
    long double mod, est_mod;
    long long coef;
    
    for(i = 0; i < size; i++) {
        mod = sqrtl((long double)data[i].r * data[i].r +
                (long double)data[i].i * data[i].i) /
            (1LL << (2 * COMPLEX_PRECISION));
        if(mod < 0.0000001) /* protectie la underflow 下溢保护*/
            continue;
        if(ctx->alfa == 1) {
            est_mod = mod - ctx->beta * ctx->noise[channel][i];
            est_mod = est_mod > 0 ? est_mod : 0;
        } else {
            est_mod = powl(mod, ctx->alfa) -
                ctx->beta * powl(ctx->noise[channel][i], ctx->alfa);
            est_mod = est_mod > 0 ? powl(est_mod, 1 / ctx->alfa) : 0;
        }
        coef = (est_mod / mod) * (1 << COMPLEX_PRECISION);
        data[i].r = data[i].r * coef / (1 << COMPLEX_PRECISION);
        data[i].i = data[i].i * coef / (1 << COMPLEX_PRECISION);
    }
}

#define MAX_EQ_CELLS 128
struct equalizer_context {
    double f[MAX_EQ_CELLS];
    double a[MAX_EQ_CELLS];
    double q[MAX_EQ_CELLS];
    t_complex *data;

    int ncells;
};

void trans_equalizer(void *context, int n, int channel, t_complex *data) {
    struct equalizer_context *ctx = context;
    int size = 1 << n;
    int i;

    for(i = 0; i < size; i++) {
        data[i].r = data[i].r * ctx->data[i].r / (1 << COMPLEX_PRECISION);
        data[i].i = data[i].i * ctx->data[i].r / (1 << COMPLEX_PRECISION);
    }
}

void trans_dummy(void *context, int n, int channel, t_complex *data) {
}

struct global_context {
    /* Variabile pentru transformarea Fourier 傅里叶变换 */
    int n, size;
    int *rev_map;
    t_complex *w;

    /* Variabile pentru formatul wave/pcm 波/ PCM格式 */
    int channels;
    int interleave;
    int bytes_per_sample;

    t_riff_hdr riff_hdr;
    t_wave_format_ex wave_fmt;

    t_complex_promote complex_promote;
    t_complex_reduce complex_reduce;

    int frame_len;
    char pad_value[4];
    unsigned long wave_len;
};

void global_context_init(struct global_context *c) {
    c->n = 12;
    c->complex_promote = NULL;
    c->complex_reduce = NULL;
}

int riff_read(struct global_context *c, FILE *in) {
    t_chunk tmp_chunk;
    unsigned long data_len = 0;
    int status;

    c->wave_len = 0;

    /* Citire header RIFF 读取头RIFF*/
    status = fread(&c->riff_hdr, sizeof(t_riff_hdr), 1, in);
    assert(status);

    if(strncmp(c->riff_hdr.chunk.id, "RIFF", 4)) {
        fprintf(stderr, "Input data is not RIFF\n");
        return 1;
    }

    if(strncmp(c->riff_hdr.file_type, "WAVE", 4)) {
        fprintf(stderr, "Unknown format '%4s'\n", c->riff_hdr.file_type);
        return 2;
    }

    data_len = c->riff_hdr.chunk.size;

    /* Restul datelor sunt chunk-uri specifice wave 其余的数据块是站点特定波 */
    do {
        assert(data_len >= sizeof(t_chunk));
        status = fread(&tmp_chunk, sizeof(t_chunk), 1, in);
        assert(status);
        data_len -= sizeof(t_chunk);

        if(!strncmp(tmp_chunk.id, "fmt ", 4)) {
            assert(tmp_chunk.size <= sizeof(t_wave_format_ex));
            memset(&c->wave_fmt, 0, sizeof(t_wave_format_ex));
            assert(data_len >= tmp_chunk.size);
            status = fread(&c->wave_fmt, tmp_chunk.size, 1, in);
            assert(status);
            data_len -= tmp_chunk.size;

            if(c->wave_fmt.format_tag != WAVE_FORMAT_PCM) {
                fprintf(stderr, "Unsupported wave format %d\n",
                        (int)c->wave_fmt.format_tag);
                return 3;
            }

            fprintf(stderr, "WAVE PCM format detected\n");
            fprintf(stderr, "%d bit audio samples\n",
                    (int)c->wave_fmt.bits_per_sample);
            fprintf(stderr, "%d channel(s), %d Hz\n", (int)c->wave_fmt.channels,
                    (int)c->wave_fmt.samples_per_sec);

            /* Calcul interleave 计算交织 */
            c->channels = c->wave_fmt.channels;
            c->bytes_per_sample = c->wave_fmt.bits_per_sample / 8;
            c->interleave = c->channels * c->bytes_per_sample;

            switch(c->wave_fmt.bits_per_sample) {
            case 8:
                c->complex_promote = complex_promote_u8;
                c->complex_reduce = complex_reduce_u8;
                *(unsigned char *)(c->pad_value) = 128;
                break;
            case 16:
                c->complex_promote = complex_promote_s16;
                c->complex_reduce = complex_reduce_s16;
                *(t_s16 *)(c->pad_value) = 0;
                break;
            default:
                fprintf(stderr, "Unsupported sample type!\n");
                return 4;
            }

            continue;
        }

        if(!strncmp(tmp_chunk.id, "data", 4)) {
            assert(data_len >= tmp_chunk.size);
            c->wave_len = tmp_chunk.size;
            break;
        }

        /* Sar peste orice alt fel de chunk (pe care nu il inteleg) 跳过任何一块的（我不明白）*/
        assert(data_len >= tmp_chunk.size);
        status = fseek(in, tmp_chunk.size, SEEK_CUR);
        assert(!status);
        data_len -= tmp_chunk.size;
    } while(1);

    c->size = 1 << c->n;
    c->frame_len = c->interleave * c->size;
    return 0;
}

void riff_write(struct global_context *c, FILE *out) {
    t_chunk wave_fmt_chunk;
    t_chunk tmp_chunk;
    int status;

    /* Output header riff, format wave si chunk date 头波，输出数据格式RIFF块， */
    /* FIXME calcul lungime date fixme数据计算长度*/
    status = fwrite(&c->riff_hdr, sizeof(t_riff_hdr), 1, out);
    assert(status);

    strncpy(wave_fmt_chunk.id, "fmt ", 4);
    wave_fmt_chunk.size = sizeof(t_wave_format_ex);
    status = fwrite(&wave_fmt_chunk, sizeof(t_chunk), 1, out);
    assert(status);
    status = fwrite(&c->wave_fmt, sizeof(t_wave_format_ex), 1, out);
    assert(status);

    strncpy(tmp_chunk.id, "data", 4);
    tmp_chunk.size = c->wave_len;
    status = fwrite(&tmp_chunk, sizeof(t_chunk), 1, out);
    assert(status);
}

void init_fft_maps(struct global_context *c) {
    /* Pregatire memorie pentru tabele si buffere 工作表和内存缓冲区 */
    c->rev_map = (int *)malloc(c->size * sizeof(int));
    assert(c->rev_map != NULL);
    c->w = (t_complex *)malloc(c->size * sizeof(t_complex));
    assert(c->w != NULL);

    bit_reverse_map(c->n, c->rev_map);
    w_map(c->n, c->w);
}

void free_fft_maps(struct global_context *c) {
    free(c->rev_map);
    free(c->w);
}

void read_eq_cells(struct equalizer_context *c, FILE *in) {
    char buf[1024], *comment;
    c->ncells = 0;
    while(fgets(buf, sizeof(buf), in) != NULL) {
        assert(strlen(buf) < sizeof(buf) - 1);
        if((comment = strchr(buf, '#')) != NULL)
            *comment = '\0';
        if(sscanf(buf, "%lf%lf%lf", &c->f[c->ncells],
                    &c->a[c->ncells], &c->q[c->ncells]) != 3)
            continue;
        if(++(c->ncells) >= MAX_EQ_CELLS)
            break;
    }
    fprintf(stderr, "Successfully read %d equalizer cells\n", c->ncells);
}

void compute_equalizer_data(struct global_context *c,
        struct equalizer_context *e) {
    /* Calculeaza amplificarea pentru fiecare componenta spectrala in functie de parametrii egalizatorului din vectorii f, a, q din contextul e
    计算每个组件的光谱放大的egalizatorului取决于参数向量F，Q的背景下
     */
    int i, j;
    double *adb;
    double log10 = log(10);
    double k;
    double prod;
    double adbi_k;

    adb = alloca(e->ncells * sizeof(double));
    /* Precalcul amplificare in decibeli pentru fiecare celula 在每一个细胞precalcul分贝的放大器 */
    for(i = 0; i < e->ncells; i++)
        adb[i] = 10 * log(e->a[i]) / log10;
    /*c */
    for(j = 0; j < c->size; j++) {
        /* Calculez frecventa reala k a componentei j in functie de
           frecventa de esantionare 真正的频率分量计算J K取决于采样频率
         */
        k = j <= c->size / 2 ?
            (double)j / c->size * c->wave_fmt.samples_per_sec :
            (double)(c->size - j) / c->size * c->wave_fmt.samples_per_sec;
        /* Calculez amplificarea pe componenta j in functie de influenta
           fiecarei celule 计算放大组件根据每个细胞的影响[J].
         */
        prod = 1;
        for(i = 0; i < e->ncells; i++) {
            adbi_k = adb[i] * exp(-fabs(log(k / e->f[i])) * e->q[i]);
            prod *= pow(10, adbi_k / 10);
        }
        /* fprintf(stderr, "%lf\n", prod); */
        e->data[j].r = prod * (1 << COMPLEX_PRECISION);
        e->data[j].i = 0;
    }
}

void print_help(char *arg) {
    fprintf(stderr, "Usage: %s [options]\n\n"
            "Options:\n"
            "\t-i file\t\tRead input data from file (defaults to stdin)\n"
            "\t-o file\t\tOutput data to file (defaults to stdout)\n"
            "\t-s\t\tPerform noise spectrum sampling\n"
            "\t-f n\t\tFFT size (in power of 2 i.e. 2^n) (defaults to 4096)\n"
            "\t-r file\t\tPerform noise reduction with noise spectrum\n"
            "\t\t\tfrom file\n"
            "\t-e file\t\tPerform parametric equalization with equalizer\n"
            "\t\t\tdata from file\n"
            "\t-d\t\tDummy run (only FFT and iFFT transform)\n"
            "\t-h\t\tThis help\n", arg);
}

int main(int argc, char **argv) {
    struct global_context c;
    /* Variabile pentru fluxul de date 数据流 */
    FILE *in = stdin;
    FILE *out = stdout;
    FILE *noise_file = NULL;
    FILE *eq_file = NULL;
    int frame_cnt = 1, frames;
    t_complex *weight;

    /* Variabile de uz general 通用*/
    int status;
    int i, j, k, opt;

    /* Variabile pentru transformarile aplicate 变化的应用 */
    int no_wave_output = 1;
    int do_sample_noise_level = 0;
    t_trans_f trans_f = NULL;
    void *trans_ctx = NULL;

    struct noise_reduction_context noise_reduction_context;
    struct equalizer_context equalizer_context;

    /* Buffer data */
    t_complex **fb0, **fb1, **fb_tmp, **fb_base;
    t_complex *fb_aux, *noise_fb_aux;
    char *buf0, *buf1, *buf_pad, *buf_tmp;
    int buf1_len;

    noise_reduction_context.alfa = 1;
    noise_reduction_context.beta = 1;

    while((opt = getopt(argc, argv, "i:o:sfr:e:dh")) != -1) {
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
            if(sscanf(optarg, "%d", &c.n) != 1) {
                fprintf(stderr, "Invalid number for FFT Size\n");
                return 21;
            }
            if(c.n < 8 || c.n > 16) {
                fprintf(stderr, "FFT Size must be between 8 and 16\n");
                return 22;
            }
            break;
        case 's': /* Sample noise level */
            do_sample_noise_level = 1;
            break;
        case 'd': /* Dummy run (no transform) 虚拟运行（不改变） */
            trans_f = trans_dummy;
            no_wave_output = 0;
            break;
        case 'r': /* Perform noise reduction 执行降噪*/
            if((noise_file = fopen(optarg, "r")) == NULL) {
                fprintf(stderr, "Failed opening %s\n", optarg);
                return 17;
            }
            no_wave_output = 0;
            trans_f = trans_noise_reduction;
            trans_ctx = &noise_reduction_context;
            break;
        case 'e': /* Perform Parametric Equalization 执行参数均衡 */
            if((eq_file = fopen(optarg, "r")) == NULL) {
                fprintf(stderr, "Failed opening %s\n", optarg);
                return 18;
            }
            read_eq_cells(&equalizer_context, eq_file);
            fclose(eq_file);
            no_wave_output = 0;
            trans_f = trans_equalizer;
            trans_ctx = &equalizer_context;
            break;
        case 'h':
            print_help(argv[0]);
            return 0;
        case '?':
        case ':':
            print_help(argv[0]);
            return 20;
        }
    }

    if(!do_sample_noise_level && trans_f == NULL) {
        print_help(argv[0]);
        return 30;
    }

    global_context_init(&c);
    if((status = riff_read(&c, in)))
        return status;

    if(do_sample_noise_level || trans_f == trans_noise_reduction) {
        /* Alocare buffere pentru mostrele de zgomot 分配的缓冲区的噪声样本*/
        noise_reduction_context.noise =
            (long double **)malloc(c.channels * sizeof(long double *));
        assert(noise_reduction_context.noise != NULL);
        for(i = 0; i < c.channels; i++) {
            noise_reduction_context.noise[i] =
                (long double *)malloc(c.size * sizeof(long double));
            assert(noise_reduction_context.noise[i] != NULL);
        }
        /* noise reprezinta o matrice cu componentele spectrale ale
           zgomotului pentru fiecare canal in parte 噪声是一种噪声谱矩阵为每个通道的组成部分*/
    }

    if(do_sample_noise_level) {
        init_fft_maps(&c);
        for(i = 0; i < c.channels; i++)
            for(k = 0; k < c.size; k++)
                noise_reduction_context.noise[i][k] = 0;
        fb_aux = malloc(c.size * sizeof(t_complex));
        assert(fb_aux != NULL);
        noise_fb_aux = malloc(c.size * sizeof(t_complex));
        assert(noise_fb_aux != NULL);
        buf0 = (char *)malloc(c.frame_len);
        assert(buf0 != NULL);
        for(j = 0; j < 10; j++) {
            status = fread(buf0, c.frame_len, 1, in);
            assert(status);
            for(i = 0; i < c.channels; i++) {
                c.complex_promote(c.size, c.interleave,
                        buf0 + i * c.bytes_per_sample, noise_fb_aux);
                bit_reverse(c.n, c.rev_map, noise_fb_aux, fb_aux);
                dec_time_fft(c.n, c.w, fb_aux);
                for(k = 0; k < c.size; k++) {
                    long double mod;

                    mod = sqrtl((long double)fb_aux[k].r * fb_aux[k].r +
                            (long double)fb_aux[k].i * fb_aux[k].i) /
                        (1LL << (2 * COMPLEX_PRECISION));
                    if(mod > noise_reduction_context.noise[i][k])
                        noise_reduction_context.noise[i][k] = mod;
                }
            }
        }
        for(i = 0; i < c.channels; i++) {
            status = fwrite(noise_reduction_context.noise[i],
                    c.size * sizeof(long double), 1, out);
            assert(status);
            free(noise_reduction_context.noise[i]);
        }
        free(noise_reduction_context.noise);
        free(fb_aux);
        free(noise_fb_aux);
        free_fft_maps(&c);
        fclose(in);
        fclose(out);
        return 0;
    }

    if(trans_f == trans_noise_reduction) {
        for(i = 0; i < c.channels; i++) {
            status = fread(noise_reduction_context.noise[i],
                    c.size * sizeof(long double), 1, noise_file);
            assert(status);
        }
        fclose(noise_file);
        fprintf(stderr, "Noise sample file OK\n");
    }

    if(trans_f == trans_equalizer) {
        equalizer_context.data =
            (t_complex *)malloc(c.size * sizeof(t_complex));
        assert(equalizer_context.data != NULL);
        compute_equalizer_data(&c, &equalizer_context);
    }

    if(no_wave_output) {
        fclose(in);
        fclose(out);
        return 0;
    }

    weight = (t_complex *)malloc(c.size * sizeof(t_complex));
    assert(weight != NULL);

    fb_base = (t_complex **)malloc(2 * c.channels * sizeof(t_complex *));
    assert(fb_base != NULL);
    for(i = 0; i < 2 * c.channels; i++) {
        fb_base[i] = (t_complex *)malloc(c.size * sizeof(t_complex));
        assert(fb_base[i] != NULL);
    }
    fb0 = fb_base;
    fb1 = fb_base + c.channels;
    fb_aux = (t_complex *)malloc(c.size * sizeof(t_complex));
    assert(fb_aux != NULL);

    buf0 = (char *)malloc(c.frame_len / 2);
    assert(buf0 != NULL);
    buf1 = (char *)malloc(c.frame_len / 2);
    assert(buf1 != NULL);
    buf_pad = (char *)malloc(c.frame_len / 2);
    assert(buf_pad != NULL);

    /* Initializare tabele 初始化表*/
    weight_map(c.n, weight);
    init_fft_maps(&c);

    for(i = c.size / 2 * c.channels, buf_tmp = buf_pad; i;
            i--, buf_tmp += c.bytes_per_sample)
        memcpy(buf_tmp, c.pad_value, c.bytes_per_sample);

    riff_write(&c, out);

    /* Aplicarea transformarilor.transformarilor应用。

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
接下来的接触过筛适用于所有通道并行（因为我在文件intretesute波），但pseudocodul会写简单的单通道。所有的prelucrarile框架在水平平行的渠道。

该算法利用原始缓冲区的尺寸2和N / 2 2（半帧和帧缓冲区（2）在复杂的，2 ^ n）的尺寸。

该算法需要处理一个额外的帧（由0）的开始和结束文件分别，为了缓解前不分别过去2 ^ ／2样品。


（frame0 join函数之frame1，2）覆盖了一半的第一在第一个半帧的第二帧，使用函数加权。
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

    frames = (c.wave_len + c.frame_len - 1) / c.frame_len;

    memcpy(buf0, buf_pad, c.frame_len / 2);
    buf1_len = fread(buf1, 1, c.frame_len / 2, in);
    for(i = 0; i < c.channels; i++) {
        c.complex_promote(c.size / 2, c.interleave,
                buf0 + i * c.bytes_per_sample, fb1[i]);
        c.complex_promote(c.size / 2, c.interleave,
                buf1 + i * c.bytes_per_sample, fb1[i] + c.size / 2);
        
        bit_reverse(c.n, c.rev_map, fb1[i], fb_aux);
        dec_time_fft(c.n, c.w, fb_aux);
        trans_f(trans_ctx, c.n, i, fb_aux);
        bit_reverse(c.n, c.rev_map, fb_aux, fb1[i]);
        dec_time_ifft(c.n, c.w, fb1[i]);
    }
    while(buf1_len == c.frame_len / 2) {
        if(!(frame_cnt % 10)) {
            /* De fapt eu numar jumatati de frame-uri 事实上我把一半的帧 */
            fprintf(stderr, "\rFrame %d/%d (%d%%)", frame_cnt / 2, frames,
                    frame_cnt * 50 / frames);
            fflush(stderr);
        }
        frame_cnt++;

        buf_tmp = buf0; buf0 = buf1; buf1 = buf_tmp;
        buf1_len = fread(buf1, 1, c.frame_len / 2, in);
        if(buf1_len < c.frame_len / 2)
            memcpy(buf1 + buf1_len, buf_pad, c.frame_len / 2 - buf1_len);
        fb_tmp = fb0; fb0 = fb1; fb1 = fb_tmp;

        for(i = 0; i < c.channels; i++) {
            c.complex_promote(c.size / 2, c.interleave,
                    buf0 + i * c.bytes_per_sample, fb1[i]);
            c.complex_promote(c.size / 2, c.interleave,
                    buf1 + i * c.bytes_per_sample, fb1[i] + c.size / 2);

            bit_reverse(c.n, c.rev_map, fb1[i], fb_aux);
            dec_time_fft(c.n, c.w, fb_aux);
            trans_f(trans_ctx, c.n, i, fb_aux);
            bit_reverse(c.n, c.rev_map, fb_aux, fb1[i]);
            dec_time_ifft(c.n, c.w, fb1[i]);

            for(j = 0; j < c.size / 2; j++) {
                fb_aux[j].r =
                    (fb0[i][c.size / 2 + j].r * weight[c.size / 2 + j].r +
                     fb1[i][j].r * weight[j].r) / (1 << COMPLEX_PRECISION);
            }
            c.complex_reduce(c.size / 2, c.interleave, fb_aux,
                    buf0 + i * c.bytes_per_sample);
        }
        status = fwrite(buf0, c.frame_len / 2, 1, out);
        assert(status);
    }
    if(buf1_len) {
        buf_tmp = buf0; buf0 = buf1; buf1 = buf_tmp;
        memcpy(buf1, buf_pad, c.frame_len / 2);
        fb_tmp = fb0; fb0 = fb1; fb1 = fb_tmp;

        for(i = 0; i < c.channels; i++) {
            c.complex_promote(c.size / 2, c.interleave,
                    buf0 + i * c.bytes_per_sample, fb1[i]);
            c.complex_promote(c.size / 2, c.interleave,
                    buf1 + i * c.bytes_per_sample, fb1[i] + c.size / 2);

            bit_reverse(c.n, c.rev_map, fb1[i], fb_aux);
            dec_time_fft(c.n, c.w, fb_aux);
            trans_f(trans_ctx, c.n, i, fb_aux);
            bit_reverse(c.n, c.rev_map, fb_aux, fb1[i]);
            dec_time_ifft(c.n, c.w, fb1[i]);

            for(j = 0; j < c.size / 2; j++) {
                fb_aux[j].r =
                    (fb0[i][c.size / 2 + j].r * weight[c.size / 2 + j].r +
                     fb1[i][j].r * weight[j].r) / (1 << COMPLEX_PRECISION);
            }
            c.complex_reduce(c.size / 2, c.interleave, fb_aux,
                    buf0 + i * c.bytes_per_sample);
        }
        status = fwrite(buf0, buf1_len, 1, out);
        assert(status);
    }
    fprintf(stderr, "\rFrame %d/%d (100%%)\n", frames, frames);
    fflush(stderr);

    
    /* Eliberare memorie 释放内存*/
    if(trans_f == trans_noise_reduction) {
        for(i = 0; i < c.channels; i++)
            free(noise_reduction_context.noise[i]);
        free(noise_reduction_context.noise);
    }
    if(trans_f == trans_equalizer) {
        free(equalizer_context.data);
    }
    free_fft_maps(&c);
    free(weight);
    for(i = 0; i < 2 * c.channels; i++)
        free(fb_base[i]);
    free(fb_base);
    free(fb_aux);
    free(buf0);
    free(buf1);

    fclose(in);
    fclose(out);
    return 0;
}
