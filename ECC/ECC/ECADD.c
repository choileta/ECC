#include "Type.h"

/**
 * @brief		덧셈 함수
 * @detail		Output = A + B
 * @param		word* A			값을 더할 첫 번째 데이터
 * @param		word* B			값을 더할 두 번째 데이터
 * @param		word* output	값을 저장할 배열
 * @return		carry			최상위비트 carry Data
*/
int addition(word* A, word* B, word* output) {
	int carry = 0;	//carry
	int subcarry = 0;
	int i = 0;		//loop data

	for (i = 7; i >= 0; i--) {
		output[i] = A[i] + B[i];
		if (output[i] < A[i])
			subcarry = 1;
		else
			subcarry = 0;

		output[i] += carry;

		if (output[i] < carry)
			subcarry = 1;

		carry = subcarry;
	}

	return carry;
}



/**
 * @brief		뺄셈 함수
 * @detail		Output = A - B
 * @param		word* A			값을 뺄 첫 번째 데이터
 * @param		word* B			값을 뺄 두 번째 데이터
 * @param		word* output	값을 저장할 배열
 * @return		borrow			최상위비트 borrow Data
*/
int subtract(word* A, word* B, word* output) {

	int borrow = 0;	//borrow
	int i = 0;		//loop data
	word value = 0;

	for (i = 7; i >= 0; i--) {
		output[i] = A[i] - B[i] - borrow;
		if (A[i] < B[i] || ((A[i] == B[i]) && borrow == 1))
			borrow = 1;
		else
			borrow = 0;
	}
	return borrow;
}

/**
 * @brief		크기 비교
 * @detail		input의 2개의 배열 대소관계 비교
 * @param		word* A		비교할 첫 번째 데이터
 * @param		word* B		비교할 첫 번째 데이터
 * @return		1			if(A>B)
 * @return		0			if(A=B)
 * @return		-1			if(A<B)
*/
int compare(word* A, word* B) {

	int i = 0;//loop counter
	for (i = 0; i < WORD_LEN; i++) {
		if (A[i] > B[i])
			return 1;
		else if (A[i] < B[i])
			return -1;
	}
	return 0;

}


/**
 * @brief		모듈러 덧셈
 * @detail		Output = (A + B) mod P
 * @param		word* A			값을 더할 첫 번째 데이터
 * @param		word* B			값을 더할 두 번째 데이터
 * @param		word* P			모듈러 P
 * @param		word* output	값을 저장할 배열
 * @return		int				함수 성공 여부
*/
int addition_mod_P(word* A, word* B, word* P, word* output) {

	int carry = 0;
	int comp = 0;
	int i = 0;
	word value[WORD_LEN];
	carry = addition(A, B, value);
	comp = compare(value, P);

	if (carry == 1 || comp != -1)
		subtract(value, P, output);
	else {
		memcpy(output, value, WORD_LEN * sizeof(word));
	}
	return 0;
}

/**
 * @brief		모듈러 뺄셈
 * @detail		Output = (A - B) mod P
 * @param		word* A			값을 뺄 첫 번째 데이터
 * @param		word* B			값을 뺄 두 번째 데이터
 * @param		word* P			모듈러 P
 * @param		word* output	값을 저장할 배열
 * @return		int				함수 성공 여부
*/
int subtract_mod_P(word* A, word* B, word* P, word* output) {

	int borrow = 0;		//borrow data
	int i = 0;			//loop data
	word value[WORD_LEN];

	borrow = subtract(A, B, value);

	if (borrow) {
		addition(value, P, output);
	}
	else
		memcpy(output, value, 8 * sizeof(word));

	return 0;
}
