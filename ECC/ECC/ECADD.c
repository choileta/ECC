#include "Type.h"

/**
 * @brief		���� �Լ�
 * @detail		Output = A + B
 * @param		word* A			���� ���� ù ��° ������
 * @param		word* B			���� ���� �� ��° ������
 * @param		word* output	���� ������ �迭
 * @return		carry			�ֻ�����Ʈ carry Data
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
 * @brief		���� �Լ�
 * @detail		Output = A - B
 * @param		word* A			���� �� ù ��° ������
 * @param		word* B			���� �� �� ��° ������
 * @param		word* output	���� ������ �迭
 * @return		borrow			�ֻ�����Ʈ borrow Data
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
 * @brief		ũ�� ��
 * @detail		input�� 2���� �迭 ��Ұ��� ��
 * @param		word* A		���� ù ��° ������
 * @param		word* B		���� ù ��° ������
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
 * @brief		��ⷯ ����
 * @detail		Output = (A + B) mod P
 * @param		word* A			���� ���� ù ��° ������
 * @param		word* B			���� ���� �� ��° ������
 * @param		word* P			��ⷯ P
 * @param		word* output	���� ������ �迭
 * @return		int				�Լ� ���� ����
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
 * @brief		��ⷯ ����
 * @detail		Output = (A - B) mod P
 * @param		word* A			���� �� ù ��° ������
 * @param		word* B			���� �� �� ��° ������
 * @param		word* P			��ⷯ P
 * @param		word* output	���� ������ �迭
 * @return		int				�Լ� ���� ����
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
