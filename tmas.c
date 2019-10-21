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
* This is an implementation of the Tuned Maximal Average Shift algorithm
* in C. Ryu, T. Lecroq, and K. Park.
*/

#include "include/define.h"
#include "include/main.h"

#define DSIGMA 4 //	DNA sequences have alphabet of size 4.

/* Input : pattern P = P[0..m-1], text T = T[0..n-1]
 * Output : the number of occurrences of P in T */
int search(unsigned char *P, int m, unsigned char *T, int n) {
	int count = 0, compare = 0;
	int i, l, k, s, w;
	int max_pos, avr_shift, max_avr_shift;
	int f, sh;

	int scan[XSIZE][XSIZE], shift[XSIZE][XSIZE][SIGMA];
	int U[XSIZE], safe[XSIZE + 1];

	char dna[DSIGMA] = { 'A','C','G','T' };

	//freq('A') = 0.293, freq('C') = 0.207, freq('G') = 0.207, freq('T') = 0.293.
	int freq[SIGMA];
	BEGIN_PREPROCESSING
	freq['A'] = 293;      freq['C'] = 207;
	freq['G'] = 207;      freq['T'] = 293;

	/* Initialize */
	for (l = 0; l < m; l++) {
		U[l] = 1;	// Set U to {0,1,...,m-1}.
	}
	for (f = 0; f < m; f++) {
		for (l = 0; l < m; l++) {
			for (s = 0; s < DSIGMA; s++) {
				shift[f][l][dna[s]] = 1;
			}
		}
	}
	for (k = 1; k <= m; k++) {
		safe[k] = 0;
	}

	/* Preprocessing */
	// The first text character scanned in the previous window is to the left of the current window : set f to position m-1.
	for (i = 0; i < m; i++) {
		for (l = 0; l < m; l++) {
			if (U[l] == 1) {
				for (s = 0; s < DSIGMA; s++) {
					for (k = shift[m - 1][l][dna[s]]; k <= m; k++) {
						if (safe[k] == 0 && ((l - k < 0) || dna[s] == (P[l - k]))) {	// (l - k < 0) means that P[l - k] for l - k < 0 matches any character. 
							shift[m - 1][l][dna[s]] = k;
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
					avr_shift = avr_shift + shift[m - 1][l][dna[s]] * freq[dna[s]];
				}
				if ((max_avr_shift < avr_shift) || (max_avr_shift == avr_shift && freq[P[max_pos]] > freq[P[l]])) {
					max_avr_shift = avr_shift;
					max_pos = l;
				}
			}
		}

		scan[m - 1][i] = max_pos;	// Determine scan[m-1][0],scan[m-1][1],...,scan[m-1][m-1].
		U[max_pos] = 0;

		for (k = 1; k <= max_pos; k++) {
			if (P[max_pos] != P[max_pos - k])
				safe[k] = 1;
		}
	}

	// The first text character scanned in the previous window is within the current window : f can be position 0 to m-2.
	for (f = 0; f < m - 1; f++) {
		for (l = 0; l < m; l++) {
			U[l] = 1;
		}
		U[f] = 0;
		for (k = 1; k <= m; k++) {
			safe[k] = 0;
		}
		for (k = 1; k <= f; k++) {
			if (P[f] != P[f - k])
				safe[k] = 1;
		}

		for (i = 0; i < m - 1; i++) {
			for (l = 0; l < m; l++) {
				if (U[l] == 1) {
					for (s = 0; s < DSIGMA; s++) {
						for (k = shift[f][l][dna[s]]; k <= m; k++) {
							if (safe[k] == 0 && ((l - k < 0) || dna[s] == (P[l - k]))) {	// (l - k < 0) means that P[l - k] for l - k < 0 matches any text character. 
								shift[f][l][dna[s]] = k;
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
						avr_shift = avr_shift + shift[f][l][dna[s]] * freq[dna[s]];
					}
					if ((max_avr_shift < avr_shift) || (max_avr_shift == avr_shift && freq[P[max_pos]] > freq[P[l]])) {
						max_avr_shift = avr_shift;
						max_pos = l;
					}
				}
			}

			scan[f][i] = max_pos;	// Determine scan[f][0],scan[f][1],...,scan[f][m-1].
			U[max_pos] = 0;

			for (k = 1; k <= max_pos; k++) {
				if (P[max_pos] != P[max_pos - k])
					safe[k] = 1;
			}
		}

		scan[f][m - 1] = f;
		shift[f][f][P[f]] = shift[m - 1][scan[m - 1][m - 1]][P[scan[m - 1][m - 1]]];
	}
	END_PREPROCESSING

	BEGIN_SEARCHING
	w = 0;
	f = m - 1;
	while (w <= n - m) {
		while (w <= n - m && P[(l = scan[f][0])] != (s = T[w + scan[f][0]])) {
			w += (sh = shift[f][l][s]);
			f = scan[f][0] - sh;
			if (f < 0)
				f = m - 1;
		}
		i = 1;
		while (i < m && P[(l = scan[f][i])] == (s = T[w + l])) {
			i++;
		}
		if (i == m && w <= n - m) {
			++count;
		}
		w = w + (sh = shift[f][l][s]);
		f = scan[f][0] - sh;
		if (f < 0)
			f = m - 1;
	}

	END_SEARCHING
	return count;
}
