/**********************************************************************/
/*                                                                    */
/* Implementarea Transformatei Fourier Rapide si a Transformatei      */
/* Fourier Rapide Inverse folosind algoritmul de decimare in timp.    */
/*                                                                    */
/* Fisierul mai contine si functiile auxiliare necesare implementarii */
/* algoritmului, cum ar fi precalcularea vectorilor W_N si a tabelei  */
/* de mapare pentru ordinea binar inversata.                          */
/*                                                                    */
/**********************************************************************/

#include <math.h>
#include "complex.h"

/* Observatii asupra implementarii:

   Pentru calculul componentelor spectrale F(k) este suficient un
   vector de numere complexe, care va fi initializat cu functiile
   cele mai granulare. La fiecare iteratie va fi generat nivelul
   urmator, pe care numarul functiilor se injumatateste si domeniul
   lor de definitie se dubleaza.

   OBS 1
   -----
   
   Este de observat ca pe fiecare nivel avem acelasi numar total
   de componente spectrale, prin urmare valorile functiilor
   auxiliare pot fi reprezentate concatenat. Un exemplu pentru N=8,
   la care pe fiecare nivel am reprezentat functiile auxiliare de
   la un pas al algoritmului:

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

   Functiile F cele mai granulare sunt definite doar in punctul 0,
   prin urmare suma va avea un singur termen, iar vectorul W fa fi
   cel unitar. Cu alte cuvinte, valoarea ei va fi egala cu o valoare
   din domeniul timp.

   Mai mult, valoarea corespunzatoare este componenta binar inversata,
   adica lui F_011 (de pe componenta 011(b) = 3) ii va corespunde
   esantionul de pe pozitia 110(b) = 6. Prin urmare vectorul de valori
   complexe poate fi initializat cu datele de intrare in ordine binar
   inversata.


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
    int i, j, max = (1 << n) - 1;
    register int b, d;

    for(j = 0; j < max; j--, map++) {
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
    int w_exp;
    int i, j, k;

    /* Masca va avea primii n biti setati si va fi folosita pentru
       calculul rapid al restului impartirii la 2^n, necesar normalizarii
       indexului in tabela de puteri ale lui W_N (puterile lui W_N sunt
       periodice, de perioada 2^n)
     */
    int w_exp_mask = max - 1;
    t_complex w_f1, w_pow;

    for(i = 0; i < n; i++) {
        /* Bucla pentru nivelul de "recursivitate" */
        for(j = 0; j < max; j += size << 1) {
            /* Bucla pentru functiile F_xx de pe nivelul curent; fiecare
               pas reprezinta calculul unei singure functii F_xx de pe
               nivelul curent
             */
            w_exp = 0;
            for(k = 0; k < size; k++) {
                /* Indicele pentru componenta din functia F_xx curenta;
                   coincide cu notatia din documentatie 
                 */

                /* Calcul exponent pentru W_N si extragerea valorii
                   pentru puterea lui W_N din tabela precalculata.
                   
                   Valoarea e de fapt negativa, asa ca fac o corectie dupa
                   ce calculez restul impartirii.
                 */
                w_pow = w[max - (w_exp & w_exp_mask ? : max)];

                /* Calculez F(k). Avand in vedere ca se calculeaza in paralel
                   cele 2 ramuri, produsul W_N * F1 (identic pentru cele 2
                   ramuri) se poate precalcula. Cum acest produs e singurul
                   dependent de datele din a 2-a jumatate, ordinea este
                   urmatoarea:
                   - calcul produs W_N * F1
                   - actualizare valoare din a 2-a jumatate
                   - actualizare valoare din prima jumatate
                 */
                w_f1.r = cmul_r(data[j + k + size], w_pow);
                w_f1.i = cmul_i(data[j + k + size], w_pow);

                data[j + k + size].r = data[j + k].r - w_f1.r;
                data[j + k + size].i = data[j + k].i - w_f1.i;
                
                data[j + k].r += w_f1.r;
                data[j + k].i += w_f1.i;
                
                w_exp += mul_k;
            }
        }
        size <<= 1;
        mul_k >>= 1;
    }
}
