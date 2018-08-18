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

/**********************************************************************/
/*                                                                    */
/* Implementarea Transformatei Fourier Rapide si a Transformatei
/* 快速傅里叶变换的实现和测试      */
/* Fourier Rapide Inverse folosind algoritmul de decimare in timp.    */
/* 利用快速傅立叶算法反在提姆的抽取     */
/* Fisierul mai contine si functiile auxiliare necesare implementarii 包含文件和必要的辅助功能，实施 */
/* algoritmului, cum ar fi precalcularea vectorilor W_N si a tabelei  向量的算法，如通过微积分基础课_ N和W tabelei */
/* de mapare pentru ordinea binar inversata.   二进制顺序映射的逆转。                       */
/*                                                                    */
/**********************************************************************/

#include <math.h>
#include "complex.h"

/* Observatii asupra implementarii:
    对实施的观察
   Pentru calculul componentelor spectrale F(k) este suficient un
   vector de numere complexe, care va fi initializat cu functiile
   cele mai granulare. La fiecare iteratie va fi generat nivelul
   urmator, pe care numarul functiilor se injumatateste si domeniul
   lor de definitie se dubleaza.
    构件计算谱f（k）是足够的数字向量的复杂的功能，将初始化最造粒。在每一代如何将产生的水平下一步，它的功能是一半的数量和范围他们定义加倍
    构件计算谱f（k）是一个数字向量足够复杂的功能，将初始化的造粒。在每一代如何将产生一个新的水平，它的功能的数量和范围的定义njumatateste加倍。
   OBS 1
   -----
   
   Este de observat ca pe fiecare nivel avem acelasi numar total
   de componente spectrale, prin urmare valorile functiilor
   auxiliare pot fi reprezentate concatenat. Un exemplu pentru N=8,
   la care pe fiecare nivel am reprezentat functiile auxiliare de
   la un pas al algoritmului:
   是我们注意到，每一层上的相同数量的总频谱分量的值，因此，可以用oncatenat辅助功能。为N = 8，每一层的功能iliare我代表在一步控制算法

   +-------+-------+-------+-------+-------+-------+-------+-------+
   |   0   |   1   |   2   |   3   |   4   |   5   |   6   |   7   |
   +===============================================================+
   | F_000 | F_001 | F_010 | F_011 | F_100 | F_101 | F_110 | F_111 |
   +-------+-------+-------+-------+-------+-------+-------+-------+
   |      F_00     |      F_01     |      F_10     |      F_11     |
   +---------------+---------------+---------------+---------------+
   |              F_0              |              F_1              |
   +-------------------------------+-------------------------------+
   |                               F                               |
   +---------------------------------------------------------------+


   OBS 2
   -----

Functiile F cele mai granulare sunt definite doar in punctul 0,prin urmare suma va avea un singur termen, iar vectorul W fa fi cel unitar. Cu alte cuvinte, valoarea ei va fi egala cu o valoare din domeniul timp.
函数f定义最颗粒只在0点，因此将有一个量的限制，和W是向量做单一。换句话说，它的价值也会与时间域的值。

Mai mult, valoarea corespunzatoare este componenta binar inversata,adica lui F_011 (de pe componenta 011(b) = 3) ii va corespunde esantionul de pe pozitia 110(b) = 6. Prin urmare vectorul de valori complexe poate fi initializat cu datele de intrare in ordine binar inversata.
更多的价值，正确的组件是二进制的，我的意思是他的（F _ 011 011成分（B）= 3）将样本110对应位置（B）= 6。因此，复杂的值初始化向量可以用二进制输入数据的顺序。

   OBS 3
   -----

Vectorul de valori complexe poate fi suprascris la fiecare iteratie.
Notam vectorul cu V. Pentru a calcula F(0) este necesar F_0(0) si
F_1(0), deci conform reprezentarii, in V(0) vom inscrie o valoare
calculata folosind V(0) si V(4). F(4) foloseste tot F_0(0) si F_1(0),
deci pentru V(4) se folosesc tot V(0) si V(4). In plus, V(0) si V(4)
nu mai sunt necesare pentru nici o alta valoare. Rezulta deci ca
vectorul V poate fi suprascris calculand perechi de cate 2
componente.
 */

/* Calculeaza o tabela de mapare care, pentru fiecare indice, contine
   indicele cu valoarea binar inversata. De exemplu daca dimensiunea
   tabelei este N=8 (3 biti), pentru componenta 3 = 011(b) va corespunde
   valoarea 6 = 110(b).

   n este numarul de biti al dimensiunii, iar map este vectorul de
   mapare.

   Dimensiunea zonei de memorie pentru map trebuie sa fie cel putin
   2^n * sizeof(int).
 */
void bit_reverse_map(int n, int *map) {
    int i, j, max = 1 << n;
    register int b, d;

    for(j = 0; j < max; j++, map++) {
        b = j;
        d = 0;
        for(i = 0; i < n; i++) {
            d = (d << 1) | (b & 1);
            b >>= 1;
        }
        *map = d;
    }
}

/* Calculeaza tabela cu vectorii W, unde W_N ^ k = exp(2 * PI * k / N).

   n este numarul de biti al dimensiunii, iar map este vectorul de
   valori complexe care va fi generat.

   Dimensiunea zonei de memorie pentru map trebuie sa fie cel putin
   2^n * sizeof(t_complex).
 */
void w_map(int n, t_complex *map) {
    int i, max = 1 << n;
    double N = max;
    double p = (double)((int)1 << COMPLEX_PRECISION);
    double arg;

    for(i = 0; i < max; i++) {
        arg = 2.0 * M_PI * i / N;
        map[i].r = p * cos(arg);
        map[i].i = p * sin(arg);
    }
}

/* Copiaza datele dintr-un vector in altul, in ordinea binar inversa.

   n este numarul de biti al dimensiunii, iar map este o harta de indici
   pentru ordinea binar inversata (eventual generata cu bit_reverse_map).
 */
void bit_reverse(int n, int *map, t_complex *src, t_complex *dst) {
    int i, max = 1 << n;
    for (i = 0; i < max; i++, src++, map++)
        dst[*map] = *src;
}

/* Datele trebuie sa fie in ordine binar inversata.
 */
void dec_time_fft(int n, t_complex *w, t_complex *data) {
    int max = 1 << n;
    int size = 1;
    int mul_k = 1 << (n - 1);
    int w0_exp, w1_exp, w1_exp_init = mul_k;
    int i, j, k;

    /* Masca va avea primii n biti setati si va fi folosita pentru
       calculul rapid al restului impartirii la 2^n, necesar normalizarii
       indexului in tabela de puteri ale lui W_N (puterile lui W_N sunt
       periodice, de perioada 2^n)
     */
    int w_exp_mask = max - 1;
    t_complex q0, q1, w_pow;

    for(i = 0; i < n; i++) {
        /* Bucla pentru nivelul de "recursivitate" */
        for(j = 0; j < max; j += size << 1) {
            /* Bucla pentru functiile F_xx de pe nivelul curent; fiecare
               pas reprezinta calculul unei singure functii F_xx de pe
               nivelul curent
             */
            w0_exp = 0;
            /* w1_exp = k_mul * size; */ /* produsul e constant */
            w1_exp = w1_exp_init;

            for(k = 0; k < size; k++) {
                /* Indicele pentru componenta din functia F_xx curenta;
                   coincide cu notatia din documentatie 
                 */

                /* Calcul exponent pentru W_N si extragerea valorii
                   pentru puterea lui W_N din tabela precalculata.
                   
                   Valoarea e de fapt negativa, asa ca fac o corectie dupa
                   ce calculez restul impartirii.
                 */
                w_pow = w[max - (w0_exp & w_exp_mask ? : max)];

                q0.r = data[j + k].r + cmul_r(data[j + size + k], w_pow);
                q0.i = data[j + k].i + cmul_i(data[j + size + k], w_pow);

                w_pow = w[max - (w1_exp & w_exp_mask ? : max)];

                q1.r = data[j + k].r + cmul_r(data[j + size + k], w_pow);
                q1.i = data[j + k].i + cmul_i(data[j + size + k], w_pow);

                data[j + k] = q0;
                data[j + size + k] = q1;

                /* Actualizare exponenti w */
                w0_exp += mul_k;
                w1_exp += mul_k;
            }
        }
        size <<= 1;
        mul_k >>= 1;
    }
}

/* Calculeaza Transformata Fourier Discreta Inversa folosind algoritmul
   de decimare in timp.
   
   Datele trebuie sa fie in ordine binar inversata.

   Avand in vedere simetria relatiilor, algoritmul este aproape identic
   cu cel pentru transformarea directa. Singurele diferente sunt
   urmatoarele:
   - Nu mai este necesara corectia de semn pentru exponentii vectorilor
     W_N.
   - Apare un factor de 1/2^n, care se distribuie in cate un factor 1/2
     pentru fiecare iteratie.
 */
void dec_time_ifft(int n, t_complex *w, t_complex *data) {
    int max = 1 << n;
    int size = 1;
    int mul_k = 1 << (n - 1);
    int w0_exp, w1_exp, w1_exp_init = mul_k;
    int i, j, k;

    int w_exp_mask = max - 1;
    t_complex q0, q1, w_pow;

    for(i = 0; i < n; i++) {
        /* Bucla pentru nivelul de "recursivitate" */
        for(j = 0; j < max; j += size << 1) {
            /* Bucla pentru functiile F_xx de pe nivelul curent; fiecare
               pas reprezinta calculul unei singure functii F_xx de pe
               nivelul curent
             */
            w0_exp = 0;
            /* w1_exp = k_mul * size; */ /* produsul e constant */
            w1_exp = w1_exp_init;

            for(k = 0; k < size; k++) {
                /* Indicele pentru componenta din functia F_xx curenta;
                   coincide cu notatia din documentatie 
                 */
                w_pow = w[w0_exp & w_exp_mask];

                q0.r = (data[j + k].r + cmul_r(data[j + size + k], w_pow)) / 2;
                q0.i = (data[j + k].i + cmul_i(data[j + size + k], w_pow)) / 2;

                w_pow = w[w1_exp & w_exp_mask];

                q1.r = (data[j + k].r + cmul_r(data[j + size + k], w_pow)) / 2;
                q1.i = (data[j + k].i + cmul_i(data[j + size + k], w_pow)) / 2;

                data[j + k] = q0;
                data[j + size + k] = q1;

                /* Actualizare exponenti w */
                w0_exp += mul_k;
                w1_exp += mul_k;
            }
        }
        size <<= 1;
        mul_k >>= 1;
    }
}
