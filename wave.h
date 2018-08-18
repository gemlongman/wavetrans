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

#ifndef _WAVE_H
#define _WAVE_H

#include "types.h"
#include "complex.h"

/* FIXME urmatoarele sunt valabile doar pe i386 接下来的fixme仅在i386 */

/* Obtine un t_s16 din continutul de la adresa indicata de p 得到一个_ T S16的目标地址P含量*/
#define WAVE_s16_2_MACHINE(p) (*(t_s16 *)(p))

/* Depune un t_s16 la adresa indicata de p 提出一个目标地址_ S16 T P */
#define WAVE_MACHINE_2_s16(p,x) (*(t_s16 *)(p) = (x))


/* Structuri pentru formatul fisierelor WAV WAV格式文件的结构*/

#ifdef __GNUC__
#define ATTRIBUTE_PACKED __attribute((packed))
#else
#define ATTRIBUTE_PACKED
#endif

typedef struct _chunk {
    char id[4];
    t_u32 size;
} ATTRIBUTE_PACKED t_chunk;

typedef struct _riff_hdr {
    t_chunk chunk;
    char file_type[4];
} ATTRIBUTE_PACKED t_riff_hdr;

typedef struct _wave_format {
    t_s16 format_tag;
    t_s16 channels;
    t_s32 samples_per_sec;
    t_s32 bytes_per_sec;
    t_s16 block_align;
} ATTRIBUTE_PACKED t_wave_format;

typedef struct _wave_format_ex {
    t_s16 format_tag;
    t_s16 channels;
    t_s32 samples_per_sec;
    t_s32 bytes_per_sec;
    t_s16 block_align;
    t_s16 bits_per_sample;
    /* t_s16 cb_size; */
} ATTRIBUTE_PACKED t_wave_format_ex;

#define WAVE_FORMAT_PCM 1

typedef void (*t_complex_promote)(int, int, char *, t_complex *);
typedef void (*t_complex_reduce)(int, int, t_complex *, char *);

extern void complex_promote_s16(int, int, char *, t_complex *);
extern void complex_reduce_s16(int, int, t_complex *, char *);
extern void complex_promote_u8(int, int, char *, t_complex *);
extern void complex_reduce_u8(int, int, t_complex *, char *);

extern void weight_map(int, t_complex *);
#endif
