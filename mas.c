/*
* SMART: string matching algorithms research tool.
* Copyright (C) 2012  Simone Faro and Thierry Lecroq
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>
*
* contact the authors at: faro@dmi.unict.it, thierry.lecroq@univ-rouen.fr
* download the tool at: http://www.dmi.unict.it/~faro/smart/
*
* This is an implementation of the Maximal Average Shift algorithm
* in C. Ryu, T. Lecroq, and K. Park.
*/

#include "include/define.h"
#include "include/main.h"

#define DSIGMA 4 //	DNA sequences have alphabet of size 4.

/* Input : pattern P = P[0..m-1], text T = T[0..n-1]
 * Output : the number of occurrences of P in T */
int search(unsigned char *P, int m, unsigned char *T, int n) {
	int count = 0;
	int i, l, k, s, w, nMinusm;
	int max_pos, avr_shift, max_avr_shift;

	int scan[XSIZE], shift[XSIZE][SIGMA];
	int U[XSIZE], safe[XSIZE + 1];

	char dna[DSIGMA] = { 'A','C','G','T' };

	//freq('A') = 0.293, freq('C') = 0.207, freq('G') = 0.207, freq('T') = 0.293.
	int freq[SIGMA];
	BEGIN_PREPROCESSING
	freq['A'] = 293;	freq['C'] = 207;
	freq['G'] = 207;	freq['T'] = 293;
	/* Initialize */
	for (l = 0; l < m; l++) {
		U[l] = 1;	// Set U to {0,1,...,m-1}.
		for (s = 0; s < DSIGMA; s++) {
			shift[l][dna[s]] = 1;
		}
	}
	for (k = 1; k <= m; k++) {
		safe[k] = 0;
	}

	/* Preprocessing */
	for (i = 0; i < m; i++) {
		for (l = 0; l < m; l++) {
			if (U[l] == 1) {
				for (s = 0; s < DSIGMA; s++) {
					for (k = shift[l][dna[s]]; k <= m; k++) {
						if (safe[k] == 0 && ((l - k < 0) || dna[s] == (P[l - k]))) {	// (l - k < 0) means that P[l - k] for l - k < 0 matches any character. 
							shift[l][dna[s]] = k;
							break;
						}
					}
				}
			}
		}

		max_avr_shift = 0;
		for (l = 0; l < m; l++) {
			if (U[l] == 1) {
				avr_shift = 0;
				for (s = 0; s < DSIGMA; s++) {
					avr_shift = avr_shift + shift[l][dna[s]] * freq[dna[s]];
				}
				if ((max_avr_shift < avr_shift) || (max_avr_shift == avr_shift && freq[P[max_pos]] > freq[P[l]])) {
					max_avr_shift = avr_shift;
					max_pos = l;
				}
			}
		}

		scan[i] = max_pos;	// Determine scan[0],scan[1],...,scan[m-1].
		U[max_pos] = 0;

		for (k = 1; k <= max_pos; k++) {
			if (P[max_pos] != P[max_pos - k])
				safe[k] = 1;
		}
	}
	END_PREPROCESSING

	BEGIN_SEARCHING
	/* Searching */
	memcpy(T + n, P, m);	// Set T[n..n+m-1] to P[0..m-1] in order to avoid testing the end of the text but exit the algorithm only when an occurrence of P is found.
	w = 0;
	nMinusm = n - m;
	while (1) {
		while (P[(l = scan[0])] != (s = T[w + l])) {	// Use the fast-loop in order to quickly shift the window when the first scan position is a mismatch.
			w += shift[l][s];
		}
		if (w <= nMinusm) {
			i = 1;
			while (i < m && P[(l = scan[i])] == (s = T[w + l])) {
				i++;
			}
			if (i == m) {
				++count;
			}
			w += shift[l][s];
		}
		else {
			END_SEARCHING
			return count;
		}
	}
}
