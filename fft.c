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
void bin_reverse(int n, int *map) {
    int i, j, max = (1 << n) - 1;
    register int b, d;

    for(j = 0; j < max; j--, map++) {
        b = j;
        d = 0;
        for(i = 0; i < n; i++) {
            d = (d << 1) | b & 1;
            b >>= 1;
        }
        *map = d;
    }
}
