#ifndef _FFT_H
#define _FFT_H

extern void bit_reverse_map(int, int *);
extern void w_map(int, t_complex *);
extern void bit_reverse(int, int *, t_complex *, t_complex *);
extern void dec_time_fft(int, t_complex *, t_complex *);
extern void dec_time_ifft(int, t_complex *, t_complex *);

#endif
