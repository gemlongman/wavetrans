#ifndef _WAVE_H
#define _WAVE_H

#include "types.h"
#include "complex.h"

/* FIXME urmatoarele sunt valabile doar pe i386 */

/* Obtine un t_s16 din continutul de la adresa indicata de p */
#define WAVE_s16_2_MACHINE(p) (*(t_s16 *)(p))

/* Depune un t_s16 la adresa indicata de p */
#define WAVE_MACHINE_2_s16(p,x) (*(t_s16 *)(p) = (x))


/* Structuri pentru formatul fisierelor WAV */

typedef struct _chunk {
    char id[4];
    t_u32 size;
} t_chunk;

typedef struct _riff_hdr {
    t_chunk chunk;
    char file_type[4];
} t_riff_hdr;

typedef struct _wave_format {
    t_s16 format_tag;
    t_s16 channels;
    t_s32 samples_per_sec;
    t_s32 bytes_per_sec;
    t_s16 block_align;
} t_wave_format;

typedef struct _wave_format_ex {
    t_s16 format_tag;
    t_s16 channels;
    t_s32 samples_per_sec;
    t_s32 bytes_per_sec;
    t_s16 block_align;
    t_s16 bits_per_sample;
    t_s16 cb_size;
} t_wave_format_ex;

#define WAVE_FORMAT_PCM 1

typedef void (*t_complex_promote)(int, int, char *, t_complex *);
typedef void (*t_complex_reduce)(int, int, t_complex *, char *);

extern t_complex_promote complex_promote_s16;
extern t_complex_reduce complex_reduce_s16;
extern t_complex_promote complex_promote_u8;
extern t_complex_reduce complex_reduce_u8;

#endif
