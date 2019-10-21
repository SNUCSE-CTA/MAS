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
* This is an implementation of the q-gram Maximal Average Shift algorithm
* in C. Ryu, T. Lecroq, and K. Park.
*/

#include "include/define.h"
#include "include/main.h"

#define Q 4	// 4-gram
#define DSIGMA 256	// Fingerprint of 4-gram DNA sequences have alphabet of size 256.

/* DNA sequences encoding
 * A = 01000001 ¡æ 00
 * C = 01000011 ¡æ 01
 * G = 01000111 ¡æ 11
 * T = 01010100 ¡æ 10 */

#define CountCG(s) (((s&64)>>6)+ ((s&16)>>4)+((s&4)>>2)+((s&1))) // Count the number of characters C and G in fingerprint.
#define HP(i) ((((P[i]&6))<<5)|(((P[i+1]&6))<<3)|(((P[i+2]&6))<<1)|((P[i+3]&6)>>1))	// h(P[i..i+3])
#define HT(i) ((((T[i]&6))<<5)|(((T[i+1]&6))<<3)|(((T[i+2]&6))<<1)|((T[i+3]&6)>>1)) // h(T[i..i+3])

 /* Input : pattern P = P[0..m-1], text T = T[0..n-1]
  * Output : the number of occurrences of P in T */
int search(unsigned char *P, int m, unsigned char *T, int n) {
	int count = 0;
	int i, l, k, s, w, wPlusr;
	int r = m % Q, mOverQ = m / Q, nMinusmPlusr;
	int matchshift, fscan, sh;
	int max_pos, avr_shift, max_avr_shift;

	int scan[XSIZE / Q], shift[XSIZE / Q][DSIGMA];
	int U[XSIZE / Q], safe[XSIZE + 1];
	int fshift[DSIGMA];

	//freq('A') = 0.293, freq('C') = 0.207, freq('G') = 0.207, freq('T') = 0.293.
	int rfreq[Q + 1], sfreq[DSIGMA];
	//rfreq : reference frequency according to the number of characters C and G in q-gram
	BEGIN_PREPROCESSING
	rfreq[0] = 74;	// (0.293)^4*(0.207)^0 = 0.0074
	rfreq[1] = 52;	// (0.293)^3*(0.207)^1 = 0.0052
	rfreq[2] = 37;	// (0.293)^2*(0.207)^2 = 0.0037
	rfreq[3] = 26;	// (0.293)^1*(0.207)^3 = 0.0026
	rfreq[4] = 18;	// (0.293)^0*(0.207)^4 = 0.0018
	for (s = 0; s < DSIGMA; s++) {
		sfreq[s] = rfreq[CountCG(s)];	// sfreq : frequency of fingerprint s
	}

	/* Initialize */
	for (l = 0; l < mOverQ; l++) {
		U[l] = 1;	// Set U to {0,1,...,m/Q-1}.
		for (s = 0; s < DSIGMA; s++) {
			shift[l][s] = 1;
		}
	}
	for (k = 1; k <= m; k++) {
		safe[k] = 0;
	}

	/* Preprocessing */
	for (i = 0; i < mOverQ; i++) {
		for (l = 0; l < mOverQ; l++) {
			if (U[l] == 1) {
				for (s = 0; s < DSIGMA; s++) {
					for (k = shift[l][s]; k <= m; k++) {
						if (safe[k] == 0 && ((Q*l - k + r < 0) || (s == HP(Q*l - k + r)))) {	// (Q*l - k + r < 0) means that h(P[Q*l-k+r..Q*l-k+r+q-1]) for Q*l - k + r < 0 matches any fingerprint. 
							shift[l][s] = k;
							break;
						}
					}
				}
			}
		}

		max_avr_shift = 0;
		for (l = 0; l < mOverQ; l++) {
			if (U[l] == 1) {
				avr_shift = 0;
				for (s = 0; s < DSIGMA; s++) {
					avr_shift = avr_shift + shift[l][s] * sfreq[s];
				}
				if ((max_avr_shift < avr_shift) || (max_avr_shift == avr_shift && sfreq[HP(Q*max_pos - k + r)] > sfreq[HP(Q*l - k + r)])) {
					max_avr_shift = avr_shift;
					max_pos = l;
				}
			}
		}

		scan[i] = max_pos;	// Determine scan[0],scan[1],...,scan[m/Q-1].
		U[max_pos] = 0;

		for (k = 1; k <= Q * max_pos + r; k++) {
			if (HP(Q*max_pos + r) != HP(Q*max_pos - k + r))
				safe[k] = 1;
		}

	}
	END_PREPROCESSING
	BEGIN_SEARCHING

	matchshift = shift[scan[mOverQ - 1]][HP(Q*scan[mOverQ - 1] + r)]; // Shift length when the whole pattern matches. 
	for (i = 0; i < mOverQ; i++)
		shift[scan[i]][HP(Q*scan[i] + r)] = 0; // Make shift length to 0 when the comparision has a match in order to simplify the text search.
	for (s = 0; s < DSIGMA; s++)
		fshift[s] = shift[scan[0]][s];	// Copy the shift length when the first comparsion. 
	fscan = Q * scan[0];	// Copy the scan position when the first comparsion. 
	/* Searching */
	memcpy(T + n, P, m);	// Set T[n..n+m-1] to P[0..m-1] in order to avoid testing the end of the text but exit the algorithm only when an occurrence of P is found.
	w = 0;	wPlusr = w + r; 	nMinusmPlusr = n - m + r;
	while (1) {
		while (sh = fshift[HT(wPlusr + fscan)]) {	// Use the fast-loop in order to quickly shift the window when the first scan position is a mismatch.
			wPlusr += sh;
		}
		if (wPlusr <= nMinusmPlusr) {
			i = 1;
			while (i < mOverQ && !(shift[l = scan[i]][s = HT(wPlusr + Q * l)])) {
				i++;
			}
			if (i == mOverQ) {
				k = 0;
				w = wPlusr - r;
				while (k < r && P[k] == T[w + k]) {
					k++;
				}
				if (k == r) {
					count++;
				}
				wPlusr += matchshift;
			}
			else
				wPlusr += shift[l][s];
		}
		else {
			END_SEARCHING
			return count;
		}
	}
}
