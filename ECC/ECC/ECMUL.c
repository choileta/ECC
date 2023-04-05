#include "Type.h"

int OS_64bit_version(word* A, word* B, word* output) {

	int i = 0, j = 0;		//loop data
	Long UV = 0;			//64bit data
	word U = 0, V = 0;		//32bit data

	//Initialize output
	memset(output, 0, 2 * WORD_LEN * sizeof(word));

	//Main_Part
	for (i = WORD_LEN - 1; i >= 0; i--) {
		U = 0;
		for (j = WORD_LEN - 1; j >= 0; j--) {
			UV = ((Long)A[i] * (Long)B[j]);
			UV += output[i + j + 1];
			UV += U;
			U = UV >> 32;
			output[i + j + 1] = (UV & 0xffffffff);
		}
		output[i] = U;
	}
	return 0;
}


int Fast_Reduction_NIST_P256(word* A, word* P, word* output) {

	int i = 0;		//loop data
	int time = 0;	//Reduction-time

	//Reduction-value
	word s1[8] = { 0, };
	word s2[8] = { 0, };
	word s3[8] = { 0, };
	word s4[8] = { 0, };
	word s5[8] = { 0, };
	word s6[8] = { 0, };
	word s7[8] = { 0, };
	word s8[8] = { 0, };
	word s9[8] = { 0, };

	word mid_value1[8];
	word mid_value2[8];
	word mid_value3[8];
	word mid_value4[8];
	word mid_value5[8];

	//s1
	s1[0] = A[8];
	s1[1] = A[9];
	s1[2] = A[10];
	s1[3] = A[11];
	s1[4] = A[12];
	s1[5] = A[13];
	s1[6] = A[14];
	s1[7] = A[15];

	//s2
	s2[0] = A[0];
	s2[1] = A[1];
	s2[2] = A[2];
	s2[3] = A[3];
	s2[4] = A[4];

	//s3
	s3[1] = A[0];
	s3[2] = A[1];
	s3[3] = A[2];
	s3[4] = A[3];

	//s4
	s4[0] = A[0];
	s4[1] = A[1];
	s4[5] = A[5];
	s4[6] = A[6];
	s4[7] = A[7];

	//s5
	s5[0] = A[7];
	s5[1] = A[2];
	s5[2] = A[0];
	s5[3] = A[1];
	s5[4] = A[2];
	s5[5] = A[4];
	s5[6] = A[5];
	s5[7] = A[6];

	//s6
	s6[0] = A[5];
	s6[1] = A[7];
	s6[5] = A[2];
	s6[6] = A[3];
	s6[7] = A[4];

	//s7
	s7[0] = A[4];
	s7[1] = A[6];
	s7[4] = A[0];
	s7[5] = A[1];
	s7[6] = A[2];
	s7[7] = A[3];

	//s8
	s8[0] = A[3];
	s8[2] = A[5];
	s8[3] = A[6];
	s8[4] = A[7];
	s8[5] = A[0];
	s8[6] = A[1];
	s8[7] = A[2];

	//s9
	s9[0] = A[2];
	s9[2] = A[4];
	s9[3] = A[5];
	s9[4] = A[6];
	s9[6] = A[0];
	s9[7] = A[1];

	time -= subtract(s4, s6, mid_value1);
	time -= subtract(s5, s9, mid_value2);
	time -= subtract(s2, s7, mid_value3);
	time -= subtract(s2, s8, mid_value4);

	time += addition(mid_value1, mid_value2, mid_value5);
	time += addition(mid_value3, mid_value5, mid_value1);
	time += addition(mid_value1, mid_value4, mid_value2);

	time += addition(mid_value2, s1, mid_value1);
	time += addition(mid_value1, s3, mid_value2);
	time += addition(mid_value2, s3, output);

	if (time > 0) {
		for (i = 0; i < time; i++) subtract(output, P, output);
		return 0;
	}
	else if (time == 0)
		return 0;
	else
		for (i = 0; i > time; i--) addition(P, output, output);
	return 0;

}
