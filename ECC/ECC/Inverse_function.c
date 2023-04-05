#include "Type.h"

int check_value_Binary(word* A) {

	int i = 0;
	for (i = 0; i < 7; i++) {
		if (A[i] != 0)
			return -1;
	}
	if (A[7] != 1)
		return -1;
	return 1;
}

void Right_Shift(word* A) {

	word i = 0;
	int cnt = 0;
	for (cnt = 7; cnt > 0; cnt--) {
		i = (A[cnt - 1] & 0x1) << 31;
		A[cnt] = (A[cnt] >> 1) + i;
	}

	A[0] = A[0] >> 1;

}

int Binary_Inverse(word* P, word* A, word* output) {

	word u[WORD_LEN];
	word v[WORD_LEN];
	word x1[WORD_LEN];
	word x2[WORD_LEN];
	word mid_value[WORD_LEN];

	word check_u = 0;
	word check_v = 0;

	word even_odd_value = 0;
	word value = 0;
	word carry = 0;
	word borrow = 0;

	memset(x1, 0, sizeof(word) * WORD_LEN);  x1[7] = 1;
	memset(x2, 0, sizeof(word) * WORD_LEN);
	memcpy(u, A, sizeof(word) * WORD_LEN);
	memcpy(v, P, sizeof(word) * WORD_LEN);

	int i = 0;

	check_u = check_value_Binary(u);
	check_v = check_value_Binary(v);

	while ((check_u == -1) && (check_v == -1)) {

		even_odd_value = u[7] % 2;
		if (even_odd_value == 0) {

			Right_Shift(u);

			if ((x1[7] % 2) == 0)
				Right_Shift(x1);
			else {
				carry = addition(x1, P, mid_value);
				memcpy(x1, mid_value, sizeof(word) * WORD_LEN);
				Right_Shift(x1);
				x1[0] = x1[0] + (carry << 31);
			}

		}

		//3.2
		even_odd_value = v[7] % 2;
		if (even_odd_value == 0) {

			Right_Shift(v);

			if ((x2[7] % 2) == 0)
				Right_Shift(x2);

			else {
				carry = addition(x2, P, mid_value);
				memcpy(x2, mid_value, sizeof(word) * WORD_LEN);
				Right_Shift(x2);
				x2[0] = x2[0] + (carry << 31);

			}
		}

		//3.3
		value = compare(u, v);

		if (value == -1) {
			subtract_mod_P(v, u, P, mid_value);
			memcpy(v, mid_value, sizeof(word) * WORD_LEN);


			subtract_mod_P(x2, x1, P, mid_value);
			memcpy(x2, mid_value, sizeof(word) * WORD_LEN);



		}
		else {
			subtract_mod_P(u, v, P, mid_value);
			memcpy(u, mid_value, sizeof(word) * WORD_LEN);

			subtract_mod_P(x1, x2, P, mid_value);
			memcpy(x1, mid_value, sizeof(word) * WORD_LEN);
		}

		check_u = check_value_Binary(u);
		check_v = check_value_Binary(v);

	}

	if (check_u == 1) {
		memcpy(output, x1, sizeof(word) * WORD_LEN);
		return 0;
	}

	else {
		memcpy(output, x2, sizeof(word) * WORD_LEN);
		return 0;
	}

}

int Fermat_based_inversion(word* z, word* P, word* output) {

	word z3[WORD_LEN];
	word z15[WORD_LEN];

	word mid_value[WORD_LEN];
	word mid_output[WORD_LEN * 2];
	word t0[WORD_LEN];
	word t1[WORD_LEN];
	word t2[WORD_LEN];
	word t3[WORD_LEN];
	word t4[WORD_LEN];
	word t5[WORD_LEN];

	int i = 0;

	//z3
	OS_64bit_version(z, z, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	OS_64bit_version(mid_value, z, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, z3);

	//z^15
	OS_64bit_version(z3, z3, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	OS_64bit_version(mid_value, z3, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, z15);

	//t0
	OS_64bit_version(z15, z15, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	OS_64bit_version(mid_value, z3, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, t0);

	//t1
	OS_64bit_version(t0, t0, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value); //t0 ^ 2

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value); //t0 ^ 4

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value); //t0 ^ 8

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value); //t0 ^ 16

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value); //t0 ^ 32

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value); //t0 ^ 64

	OS_64bit_version(mid_value, t0, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, t1);

	//t2
	OS_64bit_version(t1, t1, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	for (i = 1; i < 12; i++) {
		OS_64bit_version(mid_value, mid_value, mid_output);
		Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	}
	OS_64bit_version(mid_value, t1, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	for (i = 0; i < 6; i++) {
		OS_64bit_version(mid_value, mid_value, mid_output);
		Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	}
	OS_64bit_version(mid_value, t0, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, t2);

	//t3
	OS_64bit_version(t2, t2, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	OS_64bit_version(mid_value, z3, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, t3);

	//t4
	OS_64bit_version(t3, t3, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	for (i = 1; i < 32; i++) {
		OS_64bit_version(mid_value, mid_value, mid_output);
		Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	}
	OS_64bit_version(mid_value, z, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	for (i = 1; i < 96; i++) {
		OS_64bit_version(mid_value, mid_value, mid_output);
		Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	}
	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, t4);

	//t5
	OS_64bit_version(t4, t4, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	for (i = 1; i < 32; i++) {
		OS_64bit_version(mid_value, mid_value, mid_output);
		Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	}
	OS_64bit_version(mid_value, t3, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	for (i = 0; i < 32; i++) {
		OS_64bit_version(mid_value, mid_value, mid_output);
		Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	}
	OS_64bit_version(mid_value, t3, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, t5);

	//t
	OS_64bit_version(t5, t5, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	for (i = 1; i < 30; i++) {
		OS_64bit_version(mid_value, mid_value, mid_output);
		Fast_Reduction_NIST_P256(mid_output, P, mid_value);
	}
	OS_64bit_version(mid_value, t2, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	OS_64bit_version(mid_value, mid_value, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, mid_value);

	OS_64bit_version(mid_value, z, mid_output);
	Fast_Reduction_NIST_P256(mid_output, P, output);

	return 0;

}
