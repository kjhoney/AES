#include <stdio.h>
#include <utility>
#define IP 0x0169					// Irreducible Polynomials  : x^8 + x^6 + x^5 + x^3 + x^1
#define CONSTANT 0x15
using namespace std;

int key_matrix[4][4][11];			// 11������ key
int state[4][4];					// �� ������ �Լ��� ������ �� text ���¸� �����ϴ� ���
int temp[4];
int RC[11] = {0, };					// RC

void RC_init() {					// RC �ʱ�ȭ �Լ�
	printf("RC: ");
	RC[1] = 1;
	for (int i = 2; i < 11; i++) {	// 1���� left shift �ϸ鼭 �迭�� ����
		RC[i] = (RC[i - 1] << 1);
		if (RC[i] > 0xff)			// 8��Ʈ�� �Ѿ�� mod M
			RC[i] ^= IP;
	}
	for (int i = 1; i < 11; i++)
		printf("%02X ", RC[i]);
	puts("");
}

void print_state() {				// ���� 4x4 text ��� ���� ��� �Լ�
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			printf("%02X ", state[j][i]);
	}
	puts("");
}

void print_key(int round) {			// �� ���庰 4x4 key ��� ���� ��� �Լ�
	printf("ROUND %d: ", round);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++)
			printf("%02X ", key_matrix[j][i][round]);
	}
	puts("");
}

int gf_add(int a, int b) {			// Galois Field �� ���ϱ� ����
	return a ^ b;
}

int gf_mul(int a, int b) {			// Galois Field �� ���ϱ� ����
	int tmp = 0;
	int hi_bit_set;					// Carry �Ǻ� ��Ʈ
	for (int i = 0; i < 8; i++) {
		if (b & 1)					// b�� �� �Ʒ� �� ��Ʈ �˻�
			tmp ^= a;				// 1�� �� a�� ������
		hi_bit_set = (a & 0x80);	// a�� �� �� ��Ʈ�� 1�̸� Carry �߻�
		a <<= 1;					// a�� left shift
		if (hi_bit_set)				// Carry �߻��� mod M ����
			a ^= IP;
		b >>= 1;					// b�� right shift
	}
	return tmp;
}

int gf_inv(int a, int M) {			// Galois Field �� ������ ���� ������ ã�� �Լ�
	if (!a) return 0;
	for (int i = 0; i < 256; i++) {	// 1���� 255 ���� a * x = 1 (mod M) �� x�� ã���� return
		if (gf_mul(a, i) == 1)
			return i;
	}
}

/*int Extended_Euclid(int m, int b) {
	int A[3] = { 1, 0, m };
	int B[3] = { 0, 1, b };
	while (1) {
		if (B[2] == 0)
			return 0;
		else if (B[2] == 1)
			return B[1];
		int Q = A[2] / B[2];
		printf("Q : %d\n", Q);
		int T[3] = { A[0] ^ gf_mul(Q,B[0]), A[1] ^ gf_mul(Q,B[1]), A[2] ^ gf_mul(Q ,B[2]) };
		for (int i = 0; i < 3; i++) {
			A[i] = B[i];
			B[i] = T[i];
			printf("%02X ", B[i]);
		}
		puts("");
	}
}*/

int SBOX(int input) {						// S BOX ���� �Լ�
	int tmp, inv, ret = 0;	
	int matrix[8][8] = { 1,0,0,0,1,1,1,1,	// Affine ��ȯ ���
						1,1,0,0,0,1,1,1,
						1,1,1,0,0,0,1,1,
						1,1,1,1,0,0,0,1,
						1,1,1,1,1,0,0,0,
						0,1,1,1,1,1,0,0,
						0,0,1,1,1,1,1,0,
						0,0,0,1,1,1,1,1 };
	inv = gf_inv(input, IP);				// inv = input�� ������ ���� ����
	for (int k = 0; k < 8; k++) {			// SBOX ��İ�
		tmp = 0;
		for (int l = 0; l < 8; l++) {
			if (inv & (1 << l))
				tmp ^= matrix[k][l];
		}
		ret |= (tmp << k);					// 8��Ʈ �� �ڸ����� 8x1����� �� �ڸ���� ����
	}
	ret ^= CONSTANT;						// ����� XOR����
	return ret;								// ����� ��ȯ
}

int Inverse_SBOX(int input) {
	int tmp, tmp2, ret = 0;
	int newc = 0;
	int inv_matrix[8][8] = { 0,1,0,1,0,0,1,0,	// Affine ��ȯ ���
							0,0,1,0,1,0,0,1,
							1,0,0,1,0,1,0,0,
							0,1,0,0,1,0,1,0,
							0,0,1,0,0,1,0,1,
							1,0,0,1,0,0,1,0,
							0,1,0,0,1,0,0,1,
							1,0,1,0,0,1,0,0 };
	tmp = input % 16;
	input = (tmp << 4) + (input / 16);
	for (int k = 0; k < 8; k++) {			// Inverse SBOX ��İ� ( ������� �Է°��� ����� ��� ������ )
		tmp = 0; tmp2 = 0;
		for (int l = 0; l < 8; l++) {
			if (input & (1 << l))
				tmp ^= inv_matrix[k][l];
			if (CONSTANT & (1 << (7 - l)))
				tmp2 ^= inv_matrix[k][l];
		}
		ret |= (tmp << k);					// 8��Ʈ �� �ڸ����� 8x1����� �� �ڸ���� ����
		newc |= (tmp2 << (7 - k));			// ���� ������ ���
	}
	ret ^= newc;							// ����� XOR ����
	ret = gf_inv(ret, IP);					// ��Ŀ����� ���� �ڿ� ������ ã�� ����
	return ret;								// ����� ��ȯ
}

int substitute_byte(int byte,int flag) {	// 8��Ʈ ������ �޾� ġȯ�ϴ� �Լ�   Encryption flag : 0    Decryption flag : 1
	int row, column;
	row = byte >> 4;						// SBOX ���� ������ �� 4��Ʈ
	column = byte % 16;						// SBOX ���� ������ �� 4��Ʈ
	if (!flag)
		byte = SBOX(byte);					// Encryption�� SBOX ����
	else 
		byte = Inverse_SBOX(byte);			// Decryption�� INVERSE SBOX ����
	return byte;
}

void substitute_bytes(int flag) {			// Substitute byte �ܰ�				 Encryption flag : 0    Decryption flag : 1
	for (int i = 0; i < 4; i++) {			// 4x4 ����� �� ����Ʈ�� �� �Լ��� �̿��� ġȯ
		for (int j = 0; j < 4; j++)
			state[i][j] = substitute_byte(state[i][j],flag);
	}
	printf("SB: ");							// ġȯ�� state �迭 ���
	print_state();
}

void R_function(int round) {				// Key Expansion �� �� ���� �� R Function
	temp[0] = key_matrix[1][3][round - 1];	// �� Round�� �������� ���� shift
	temp[1] = key_matrix[2][3][round - 1];
	temp[2] = key_matrix[3][3][round - 1];
	temp[3] = key_matrix[0][3][round - 1];
	for (int i = 0; i < 4; i++)
		temp[i] = substitute_byte(temp[i], 0); // �� Round�� ������ �� ġȯ
	temp[0] ^= RC[round];					// �� �� ����Ʈ�� �ش� Round�� RC ���� XOR ����
}

void key_expansion() {						// Key Expansion
	print_key(0);							// ���� Key ��� ( ROUND 0 )
	for (int i = 1; i < 11; i++) {			
		R_function(i);						// ���� �� R Function ����
		for (int j = 0; j < 4; j++) {		// �� ���� ù��° ���� R Function ���� ���� �� Round�� ù��°�� XOR ����
			key_matrix[j][0][i] = temp[j] ^ key_matrix[j][0][i - 1];
			for (int k = 1; k < 4; k++)		// ������ ���� �� Round�� ���� �ش� Round�� �ٷ� �� ���� XOR ����
				key_matrix[j][k][i] = key_matrix[j][k - 1][i] ^ key_matrix[j][k][i - 1];

		}
		print_key(i);						// �ش� Round Key ���
	}
}

void add_round_key(int round) {				// Add Round Key �ܰ� ( Encryption�� Decryption�� �������� ���� )
	for (int i = 0; i < 4; i++) {			// �ش� Round�� Key �� Text�� XOR ����
		for (int j = 0; j < 4; j++) {
			state[i][j] = state[i][j] ^ key_matrix[i][j][round];
		}
	}
	printf("AR: ");
	print_state();							// Text ���
}

void shift_row(int flag) {					// Shift Row �ܰ�    Encryption flag : 0    Decryption flag : 1
	int tmp;
	if (!flag) {							// Encryption �� �� �ึ�� left shift
		for (int i = 1; i < 4; i++) {
			for (int j = 1; j <= i; j++) {	// 2��° ���� 1��, 3��° ���� 2��, 4��° ���� 3�� shift
				tmp = state[i][0];
				state[i][0] = state[i][1];
				state[i][1] = state[i][2];
				state[i][2] = state[i][3];
				state[i][3] = tmp;
			}
		}
	}
	else {									// Decryption �� �� �ึ�� right shift
		for (int i = 1; i < 4; i++) {
			for (int j = 1; j <= i; j++) {	// 2��° ���� 1��, 3��° ���� 2��, 4��° ���� 3�� shift
				tmp = state[i][0];
				state[i][0] = state[i][3];
				state[i][3] = state[i][2];
				state[i][2] = state[i][1];
				state[i][1] = tmp;
			}
		}
	}
	printf("SR: ");
	print_state();							// Text ���
}

void mix_column(int flag) {					// Mix Column �ܰ�    Encryption flag : 0    Decryption flag : 1
	int matrix[4][4] = {					// Encryption�� ������ ���
		0x02, 0x03, 0x01, 0x01,
		0x01, 0x02, 0x03, 0x01,
		0x01, 0x01, 0x02, 0x03,
		0x03, 0x01, 0x01, 0x02
	};
	int inv_matrix[4][4] = {				// Decryption�� ������ ���
		0x0E, 0x0B, 0x0D, 0x09,
		0x09, 0x0E, 0x0B, 0x0D,
		0x0D, 0x09, 0x0E, 0x0B,
		0x0B, 0x0D, 0x09, 0x0E,
	};
	for (int i = 0; i < 4; i++) {			// ��İ�
		int temp[4] = { 0, };
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				if (!flag)
					temp[j] ^= gf_mul(matrix[j][k], state[k][i]);
				else
					temp[j] ^= gf_mul(inv_matrix[j][k], state[k][i]);
			}
		}
		state[0][i] = temp[0];
		state[1][i] = temp[1];
		state[2][i] = temp[2];
		state[3][i] = temp[3];
	}
	printf("MC: ");
	print_state();							// Text ���
}

int main() {
	int read;
	FILE * plain;
	FILE * key;
	FILE * cipher;
	FILE * decrypt;
	plain = fopen("plain.bin", "rb");
	key = fopen("key.bin", "rb");			// plain.bin �� key.bin ���� �б���� ����
	cipher = fopen("cipher.bin", "wb");
	decrypt = fopen("decrypt.bin", "wb");

	RC_init();								// RC �ʱ�ȭ

	printf("PLAIN: ");
	for (int i = 0; i < 4; i++) {			// plain.bin ������ 1����Ʈ�� �о� state�� ���� �� ���
		for (int j = 0; j < 4; j++) {
			read = fgetc(plain);
			state[j][i] = read;
			printf("%02X ", read);
		}
	}

	printf("\nKEY: ");
	for (int i = 0; i < 4; i++) {			// key.bin ������ 1����Ʈ�� �о� key ��Ŀ� ���� �� ���
		for (int j = 0; j < 4; j++) {
			read = fgetc(key);
			key_matrix[j][i][0] = read;
			printf("%02X ", read);
		}
	}

	printf("\n\n<------ENCRYPTION------>\n\n"); // ENCRYPTION
	printf("KEY EXPANSION\n");

	key_expansion();						// Key expansion ����

	for (int i = 0; i < 11; i++) {			// �� Round ���� �ؾߵǴ� �ܰ�� ����
		printf("\nROUND %d\n", i);
		if (i) {							// ROUND 0 ������ Add round key�� ����
			substitute_bytes(0);
			shift_row(0);
			if (i != 10)					// ROUND 10 ������ Mix column �ܰ踦 �������� ����
				mix_column(0);
		}
		add_round_key(i);
	}

	printf("\nCIPHER: ");
	print_state();							// ���� Ciphertext ���

	for (int i = 0; i < 4; i++)				// cipher.bin ���Ͽ� state ����� 1Byte�� ��
		for (int j = 0; j < 4; j++)
			fputc(state[j][i], cipher);

	printf("\n\n<------DECRYPTION------>\n\n"); // DECRYPTION

	for (int i = 0; i < 11; i++) {			// �� Round ���� �ؾߵǴ� �ܰ踦 ����
		printf("\nROUND %d\n", i);
		if (i) {							// ROUND 0 ������ Add round key�� ����
			shift_row(1);
			substitute_bytes(1);
		}
		add_round_key(10 - i);				// Add round key�� Key�� ROUND 10���� �Ųٷ� ���
		if (i && i != 10)					// ROUND 10 ������ Mix column �ܰ踦 �������� ����
			mix_column(1);
	}

	printf("\nDECRPYTED: ");
	print_state();							// ���� Plaintext ���

	for (int i = 0; i < 4; i++)				// decrypt.bin ���Ͽ� state ����� 1Byte�� ��
		for (int j = 0; j < 4; j++)
			fputc(state[j][i], decrypt);

	fclose(plain);
	fclose(key);
	fclose(cipher);
	fclose(decrypt);
}