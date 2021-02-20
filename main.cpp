#include "LowMC.h"
#include <iostream>
#include <ctime>
using namespace std;

void testLayerWithOneSBox() {
	LowMC lowmc(20, 20, 1, 23);
	bool pdiff[20] = { 0,0,0,0,1,1,1,1,0,0,0,1,0,1,0,0,1,1,0,1 };//choose it
	bool cdiff[20];//difference in the ciphertext
	bool p0[20], k[20], c0[20], p1[20], c1[20];
	bool eachRoundOut0[23][20], eachRoundOut1[23][20], eachRoundDiff[23][20];
	bool eachRoundSBoxOut0[23][20], eachRoundSBoxOut1[23][20], eachRoundSBoxDiff[23][20];
	int psize = 20, ksize = 20;
	srand(time(NULL));
	bool check[20] = { 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0 };

	int testTimes = 1000;
	cout << "There are in total "<<testTimes << " experiments!" << endl;
	UINT64 result[3];
	UINT64 sum[3] = {0,0,0};
	for (int test = 0; test < testTimes; test++) {
		cout << endl << "Experiments: Times " << test + 1 <<" ("<<testTimes<<" tests in total)" << endl;
		for (int i = 0; i < psize; i++) {
			p0[i] = rand() % 2;
			k[i] = rand()%2;
			p1[i] = p0[i] ^ pdiff[i];
		}
		lowmc.encrypt(p0, k, psize, ksize, c0, 23, eachRoundOut0, eachRoundSBoxOut0);//r0=6
		lowmc.encrypt(p1, k, psize, ksize, c1, 23, eachRoundOut1, eachRoundSBoxOut1);
		for (int i = 0; i < 20; i++) {
			cdiff[i] = c0[i] ^ c1[i];
		}

		//output the difference after 6 rounds
		for (int i = 0; i < 23; i++) {
			//cout << "The difference after "<<i+1<<" rounds: ";
			for (int j = 0; j < 20; j++) {
				eachRoundDiff[i][j] = eachRoundOut0[i][j] ^ eachRoundOut1[i][j];
				eachRoundSBoxDiff[i][j] = eachRoundSBoxOut0[i][j] ^ eachRoundSBoxOut1[i][j];
				//cout << eachRoundSBoxDiff[i][j];
			}
			//cout << endl;
		}
		lowmc.setAccurateRoundDiff(eachRoundDiff, eachRoundSBoxDiff);
		lowmc.constructForwardDiffR1(check, 7, 6);
		lowmc.enumerateBackwardDiffR2(cdiff, 6, 7, 10,result);
		sum[0] += result[0];
		sum[1] += result[1];
		sum[2] += result[2];
		cout << "Average time to enumerate trails:" << sum[1] / (test+1) << endl;
		cout << "Average number of valid trails:" << sum[0] / (test+1) << endl;
		cout << "Average time to recover the key:" << sum[2] / (test+1) << endl;
	}
	//cout << "Average time to enumerate trails:" << sum[1]/testTimes << endl;
	//cout << "Average number of valid trails:" << sum[0] / testTimes << endl;
	//cout << "Average time to recover the key:" << sum[2] / testTimes << endl;
	//lowmc.computeR0();
}

void exhaustiveSearch() {
	LowMC lowmc(20, 20, 1, 23);
	bool p[20],c[20],k[20];
	int psize = 20, ksize = 20;
	for (int i = 0; i < psize; i++) {
		p[i] = rand() % 2;
	}
	for (int i = 0; i < 0x100000; i++) {
		for (int j = 0; j < 20; j++) {
			k[j] = (i >> j) & 0x1;
		}
		lowmc.encrypt(p, k, psize, ksize, c, 23);
	}
}

void testLayerWithFullSBox() {
	LowMC lowmc(21, 21, 7, 4);
	bool bestDiff[21];
	lowmc.chooseBestInputDiff(bestDiff);
	cout << "best diff:" << endl;
	for (int i = 0; i < 7; i++) {
		cout << bestDiff[3 * i];
		cout << bestDiff[3 * i+1];
		cout << bestDiff[3 * i+2];
		cout << " ";
	}
	cout << endl;

	bool pdiff[21] = { 0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };//choose it
	bool cdiff[21];//difference in the ciphertext
	bool p0[21], k[21], c0[21], p1[21], c1[21];
	bool eachRoundOut0[4][21], eachRoundOut1[4][21], eachRoundDiff[4][21];
	bool eachRoundSBoxOut0[4][21], eachRoundSBoxOut1[4][21], eachRoundSBoxDiff[4][21];
	int psize = 21, ksize = 21;
	int smallerTimes = 0, largerTimes = 0,testTimes=10000;
	for (int test = 0; test < testTimes; test++) {
		cout << endl << "Experiments: Times " <<dec<< test + 1 << endl;
		
		bool isdesired = false;
		while (!isdesired) {
			for (int i = 0; i < psize; i++) {
				p0[i] = rand() % 2;
				k[i] = rand() % 2;
				p1[i] = p0[i] ^ pdiff[i];
			}

			lowmc.encryptFull(p0, k, psize, ksize, c0, 4, eachRoundOut0, eachRoundSBoxOut0);
			lowmc.encryptFull(p1, k, psize, ksize, c1, 4, eachRoundOut1, eachRoundSBoxOut1);
			for (int i = 0; i < psize; i++) {
				cdiff[i] = c0[i] ^ c1[i];
			}

			//the internal state diff
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 21; j++) {
					eachRoundDiff[i][j] = eachRoundOut0[i][j] ^ eachRoundOut1[i][j];
					eachRoundSBoxDiff[i][j] = eachRoundSBoxOut0[i][j] ^ eachRoundSBoxOut1[i][j];
				}
			}

			//choose a desired diff
			bool find = true;
			for (int i = 0; i < psize; i++) {
				if (eachRoundSBoxDiff[0][i] != bestDiff[i]) {
					find = false;
					break;
				}
			}
			isdesired = find;
		}
		cout << "Found a desired diff" << endl;
		cout << "The input difference of the 2nd round:";
		for (int j = 0; j < 7; j++) {
			cout << eachRoundDiff[0][j * 3];
			cout << eachRoundDiff[0][j * 3+1];
			cout << eachRoundDiff[0][j * 3+2];
			cout << " ";
		}
		cout << endl;

		cout << "The output difference of the Sbox in the 4th round:";
		for (int j = 0; j < 7; j++) {
			cout << eachRoundSBoxDiff[3][j * 3];
			cout << eachRoundSBoxDiff[3][j * 3 + 1];
			cout << eachRoundSBoxDiff[3][j * 3 + 2];
			cout << " ";
		}
		cout << endl;

		int requiredNum = 0;
		//cout << "required number of active sboxes:";
		int inac = 0;
		for (int j = 0; j < 7; j++) {
			if (eachRoundSBoxDiff[3][3 * j] == 0
				&& eachRoundSBoxDiff[3][3 * j + 1] == 0
				&& eachRoundSBoxDiff[3][3 * j + 2]==0) {
				inac++;
			}
		}
		requiredNum = (21 + 3 * inac - 2) - (7 - inac) * 2;
		//cout << "required extra equations:" << dec<<requiredNum << endl;

		int secondNum = 7-lowmc.getInactiveNum(eachRoundDiff[0],7);
		//cout << "secondNum:" << secondNum << endl;

		int q = 7 - secondNum;//#(inactive Sboxes in the 2nd round)
		int t = inac;//#(inactive sboxes in the last round)
		cout << "q:" << dec <<q<< "; t:" << t << endl;
		
		int f=lowmc.startTestingFullSBoxLayer(cdiff, eachRoundDiff[0], eachRoundDiff[2],requiredNum,secondNum,q,t);
		if (f == 1) {
			smallerTimes++;
		}
		else {
			largerTimes++;
		}
	}
	cout << "testTimes: " <<dec<< testTimes << endl;
	cout << "the times that #diffs is smaller than expected: " << dec << smallerTimes << endl;
	cout << "the times that #diffs is larger than expected: " << dec << largerTimes << endl;
}

int main() {
	cout << "1 -> test the construction with a partial S-box layer (1000 tests in a few minutes)" << endl;
	cout << "2 -> test the construction with a full S-box layer (10000 tests in about 2 minute)" << endl;
	int cmd;
	cout << endl << "input command (1/2):";
	cin >> cmd;
	if (cmd == 1) {
		testLayerWithOneSBox();
	}
	else {
		testLayerWithFullSBox();
	}
	return 0;
}