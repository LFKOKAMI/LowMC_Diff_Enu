#include "LowMC.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
using namespace std;

LowMC::LowMC(int BS, int KS, int M, int R) {
	bs = BS;
	ks = KS;
	m = M;
	r = R;

	linear = new matrix [r];
	invLinear = new matrix[r];
	keyMa = new matrix[r + 1];
	
	if (bs != 3 * m) {
		loadFileNonFull();
	}
	else {
		loadFileFull();
	}

	srand(time(NULL));
}

void LowMC::loadFileNonFull() {
	ifstream linearFile("linear.txt");
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			for (int k = 0; k < bs; k++) {
				linearFile >> linear[i].ma[j][k];
			}
		}
		linear[i].r = bs;
		linear[i].c = bs;
	}
	linearFile.close();
	//inverse
	ifstream invLinearFile("inverse.txt");
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			for (int k = 0; k < bs; k++) {
				invLinearFile >> invLinear[i].ma[j][k];
			}
		}
		invLinear[i].r = bs;
		invLinear[i].c = bs;
	}
	invLinearFile.close();
	//constant
	ifstream constantFile("constant.txt");
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			constantFile >> constant.ma[i][j];
		}
	}
	constant.r = r;
	constant.c = bs;
	constantFile.close();
	//key
	ifstream keyFile("key.txt");
	for (int i = 0; i < r + 1; i++) {
		for (int j = 0; j < bs; j++) {
			for (int k = 0; k < bs; k++) {
				keyFile >> keyMa[i].ma[j][k];
			}
		}
		keyMa[i].r = bs;
		keyMa[i].c = ks;
	}
	keyFile.close();
}

void LowMC::loadFileFull() {
	ifstream linearFile("linearFull.txt");
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			for (int k = 0; k < bs; k++) {
				linearFile >> linear[i].ma[j][k];
			}
		}
		linear[i].r = bs;
		linear[i].c = bs;
	}
	linearFile.close();
	//inverse
	ifstream invLinearFile("inverseFull.txt");
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			for (int k = 0; k < bs; k++) {
				invLinearFile >> invLinear[i].ma[j][k];
			}
		}
		invLinear[i].r = bs;
		invLinear[i].c = bs;
	}
	invLinearFile.close();
	//constant
	ifstream constantFile("constantFull.txt");
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			constantFile >> constant.ma[i][j];
		}
	}
	constant.r = r;
	constant.c = bs;
	constantFile.close();
	//key
	ifstream keyFile("keyFull.txt");
	for (int i = 0; i < r + 1; i++) {
		for (int j = 0; j < bs; j++) {
			for (int k = 0; k < bs; k++) {
				keyFile >> keyMa[i].ma[j][k];
			}
		}
		keyMa[i].r = bs;
		keyMa[i].c = ks;
	}
	keyFile.close();
}

void LowMC::setAccurateRoundDiff(bool eachRoundDiff[][20], bool eachRoundSBoxDiff[][20]) {
	for (int i = 0; i < r; i++) {
		for (int j = 0; j < bs; j++) {
			accurateRoundDiff[i][j] = eachRoundDiff[i][j];
			accurateRoundSBoxDiff[i][j] = eachRoundSBoxDiff[i][j];
		}
	}

}

void LowMC::encrypt(bool p0[], bool k[], int psize, int ksize,bool c[],int rounds,bool eachRoundOut[][20], bool eachRoundSBoxOut[][20]) {
	matrix rk;
	rk.r = r + 1;
	rk.c = bs;

	bool* p;
	p = new bool[psize];
	for (int i = 0; i < psize; i++) {
		p[i] = p0[i];
	}

	for (int i = 0; i < r + 1; i++) {
		matrixMul(keyMa[i], k, rk.ma[i]);
	}
	//whitening key
	for (int i = 0; i < psize; i++) {
		p[i] = p[i] ^ rk.ma[0][i];
	}

	int t[3];
	bool* pt = new bool[psize];
	for (int i = 0; i < rounds; i++) {
		//s-box (only works for the first 3m bits)
		for (int j = 0; j < m; j++) {
			t[0] = p[3*j] ^ (p[1+ 3 * j] & p[2+ 3 * j]);
			t[1] = p[3 * j] ^ p[1+ 3 * j] ^ (p[3 * j] & p[2+ 3 * j]);
			t[2] = p[3 * j] ^ p[1+ 3 * j] ^ p[2+ 3 * j] ^ (p[3 * j] & p[1+ 3 * j]);
			p[3 * j] = t[0];
			p[3 * j+1] = t[1];
			p[3 * j+2] = t[2];
		}
		//the value after SBox
		for (int j = 0; j < psize; j++) {
			eachRoundSBoxOut[i][j] = p[j];
		}

		//linear matrix
		matrixMul(linear[i], p, pt);

		//constant addition
		for (int j = 0; j < psize; j++) {
			p[j] = pt[j] ^ constant.ma[i][j];
		}

		//key addition
		for (int j = 0; j < psize; j++) {
			p[j] = p[j] ^ rk.ma[i+1][j];
		}

		for(int j=0;j<psize;j++){
			eachRoundOut[i][j] = p[j];
		}
	}
	for (int i = 0; i < psize; i++) {
		c[i] = p[i];
	}
	delete[]pt;
	delete[]p;
}

void LowMC::encryptFull(bool p0[], bool k[], int psize, int ksize, bool c[], int rounds, bool eachRoundOut[][21], bool eachRoundSBoxOut[][21]) {
	matrix rk;
	rk.r = r + 1;
	rk.c = bs;

	bool* p;
	p = new bool[psize];
	for (int i = 0; i < psize; i++) {
		p[i] = p0[i];
	}

	for (int i = 0; i < r + 1; i++) {
		matrixMul(keyMa[i], k, rk.ma[i]);
	}
	//whitening key
	for (int i = 0; i < psize; i++) {
		p[i] = p[i] ^ rk.ma[0][i];
	}

	int t[3];
	bool* pt = new bool[psize];
	for (int i = 0; i < rounds; i++) {
		//s-box (only works for the first 3m bits)
		for (int j = 0; j < m; j++) {
			t[0] = p[3 * j] ^ (p[1 + 3 * j] & p[2 + 3 * j]);
			t[1] = p[3 * j] ^ p[1 + 3 * j] ^ (p[3 * j] & p[2 + 3 * j]);
			t[2] = p[3 * j] ^ p[1 + 3 * j] ^ p[2 + 3 * j] ^ (p[3 * j] & p[1 + 3 * j]);
			p[3 * j] = t[0];
			p[3 * j + 1] = t[1];
			p[3 * j + 2] = t[2];
		}
		//the value after SBox
		for (int j = 0; j < psize; j++) {
			eachRoundSBoxOut[i][j] = p[j];
		}

		//linear matrix
		matrixMul(linear[i], p, pt);

		//constant addition
		for (int j = 0; j < psize; j++) {
			p[j] = pt[j] ^ constant.ma[i][j];
		}

		//key addition
		for (int j = 0; j < psize; j++) {
			p[j] = p[j] ^ rk.ma[i + 1][j];
		}

		for (int j = 0; j < psize; j++) {
			eachRoundOut[i][j] = p[j];
		}
	}
	for (int i = 0; i < psize; i++) {
		c[i] = p[i];
	}
	delete[]pt;
	delete[]p;
}

void LowMC::encrypt(bool p0[], bool k[], int psize, int ksize, bool c[], int rounds) {
	matrix rk;
	rk.r = r + 1;
	rk.c = bs;

	bool* p;
	p = new bool[psize];
	for (int i = 0; i < psize; i++) {
		p[i] = p0[i];
	}

	for (int i = 0; i < r + 1; i++) {
		matrixMul(keyMa[i], k, rk.ma[i]);
	}
	//whitening key
	for (int i = 0; i < psize; i++) {
		p[i] = p[i] ^ rk.ma[0][i];
	}

	int t[3];
	bool* pt = new bool[psize];
	for (int i = 0; i < rounds; i++) {
		//s-box (only works for the first 3m bits)
		for (int j = 0; j < m; j++) {
			t[0] = p[3 * j] ^ (p[1 + 3 * j] & p[2 + 3 * j]);
			t[1] = p[3 * j] ^ p[1 + 3 * j] ^ (p[3 * j] & p[2 + 3 * j]);
			t[2] = p[3 * j] ^ p[1 + 3 * j] ^ p[2 + 3 * j] ^ (p[3 * j] & p[1 + 3 * j]);
			p[3 * j] = t[0];
			p[3 * j + 1] = t[1];
			p[3 * j + 2] = t[2];
		}

		//linear matrix
		matrixMul(linear[i], p, pt);

		//constant addition
		for (int j = 0; j < psize; j++) {
			p[j] = pt[j] ^ constant.ma[i][j];
		}

		//key addition
		for (int j = 0; j < psize; j++) {
			p[j] = p[j] ^ rk.ma[i + 1][j];
		}
	}
	for (int i = 0; i < psize; i++) {
		c[i] = p[i];
	}
	delete[]pt;
	delete[]p;
}

void LowMC::decrypt(bool p[], bool k[], int psize, int ksize, bool c[],int rounds) {
	matrix rk;
	rk.r = r + 1;
	rk.c = bs;
	for (int i = 0; i < r + 1; i++) {
		matrixMul(keyMa[i], k, rk.ma[i]);
	}

	int t[3];
	bool* pt = new bool[psize];
	for (int i = rounds-1; i >= 0; i--) {
		//key addition
		for (int j = 0; j < psize; j++) {
			p[j] = p[j] ^ rk.ma[i + 1][j];
		}

		//constant addition
		for (int j = 0; j < psize; j++) {
			pt[j] = p[j] ^ constant.ma[i][j];
		}

		//linear matrix
		matrixMul(invLinear[i], pt, p);

		//inverse of s-box (only works for the first 3m bits)
		for (int j = 0; j < m; j++) {
			t[0] = p[3 * j] ^ p[1 + 3 * j] ^ (p[1 + 3 * j] & p[2 + 3 * j]);
			t[1] = p[1+3 * j] ^ (p[3 * j] & p[2 + 3 * j]);
			t[2] = p[3 * j] ^ p[1 + 3 * j] ^ p[2 + 3 * j] ^ (p[3 * j] & p[1 + 3 * j]);
			p[3 * j] = t[0];
			p[3 * j + 1] = t[1];
			p[3 * j + 2] = t[2];
		}
	}
	//whitening key
	for (int i = 0; i < psize; i++) {
		p[i] = p[i] ^ rk.ma[0][i];
	}

	for (int i = 0; i < psize; i++) {
		cout << p[i];
	}
	cout << endl;

	delete[]pt;
}

void LowMC::computeR0() {
	// r0=bs/3m
	int r0 = bs / (3 * m);
	cout << "r0:" << r0 << endl;

	matrix eq;
	eq.r = r0 * 3 * m;
	eq.c = bs+1;
	clearMatrix(eq);
	int cnt = 0;

	matrix initState;
	initState.r = bs;
	initState.c = bs;
	clearMatrix(initState);
	for (int i = 0; i < bs; i++) {
		initState.ma[i][i] = 1;
	}

	for (int i = 0; i < r0; i++) {
		//derive equations
		for (int j = 0; j < 3 * m; j++) {
			for (int k = 0; k < initState.c; k++) {
				eq.ma[cnt][k] = initState.ma[j][k];
			}
			eq.ma[cnt][bs] = 0;
			cnt++;
		}

		for (int j = 0; j < 3 * m; j++) {
			for (int k = 0; k < initState.c; k++) {
				initState.ma[j][k] = 0;
			}
		}

		matrixMul(linear[i], initState);
	}
	//gauss elimination
	gauss(eq);

	vector<vector<bool>> sol;
	sol.clear();
	int solNum = 0;
	storeSolutions(sol, eq, solNum);
	cout << "solutions the input difference:"<< endl;
	for (int i = 0; i < solNum; i++) {
		for (int j = 0; j < eq.c - 1; j++) {
			cout << sol[i][j] << ",";
		}
		cout << endl;
	}
}

void LowMC::constructForwardDiffR1(bool diff[],int r1,int r0) {
	//check the first 3m bits of inputDiff: inputDiff[0,1,...,3m-1]
	//introudce at most 2m variables in the first round
	//introduce 3*(r1-2)*m variables in the next r1-2 rounds
	//connect the output difference in the last round
	matrix eq;
	eq.r = bs;
	eq.c = 2 * m + 3 * (r1 - 2) * m + 1;//eq.c=2+15+1=18 for r1=7,m=1
	clearMatrix(eq);
	//matrix start;
	startPart.r = bs;
	startPart.c= 2 * m + 3 * (r1 - 2) * m + 1;
	clearMatrix(startPart);

	//the first round
	//the variables are chosen accorinding to the input difference
	int cnt = 0;//the indices of variables
	for (int i = 0; i < m; i++) {
		returnLinearRelationForwards(diff[3 * i], diff[3 * i + 1], diff[3 * i + 2], eq, i, cnt);
	}
	for (int i = m * 3; i < bs; i++) {
		eq.ma[i][eq.c - 1] = diff[i];
	}
	//copy eq to start for check
	for (int i = 0; i < 3*m; i++) {
		for (int j = 0; j < startPart.c; j++) {
			startPart.ma[i][j] = eq.ma[i][j];
		}
	}
	for (int i = 3 * m; i < bs; i++) {
		startPart.ma[i][startPart.c - 1] = diff[i];
	}
	//output eq
	//cout << "eq:" << endl;
	//outputMatrix(eq);

	//apply the linear transform
	matrixMul(linear[r0], eq);
	//cout << "variables for the first sbox layer:" << cnt << endl;
	/*bool var[17],result[20];
	var[0] = accurateRoundSBoxDiff[r0][1];
	var[1] = accurateRoundSBoxDiff[r0][2];
	for (int i = 2; i < 17; i++) {
		var[i] = 0;
	}
	for (int i = 0; i < bs; i++) {
		result[i] = eq.ma[i][eq.c - 1];
		for (int j = 0; j < eq.c-1; j++) {
			if (eq.ma[i][j] == 1) {
				result[i] ^= var[j];
			}
		}
	}
	//compare result with accurateRoundDiff[r0]
	cout << "diff:" << endl;
	for (int i = 0; i < bs; i++) {
		cout << result[i];
	}
	cout << endl;
	for (int i = 0; i < bs; i++) {
		cout << accurateRoundDiff[r0][i];
	}
	cout << endl;
	system("pause");
	*/

	//the next r1-2 (5) rounds
	for (int i = 1; i < r1-1; i++) {
		//replace equations with new variables
		for (int j = 0; j < 3 * m; j++) {
			for (int k = 0; k < cnt; k++) {//first clear
				eq.ma[j][k] = 0;
			}
			eq.ma[j][cnt] = 1;//replace with a new variable v[cnt];
			eq.ma[j][eq.c - 1] = 0;
			cnt++;
		}
		//linear transform
		matrixMul(linear[i+r0], eq);
	}
	//cout << "cnt:" << cnt << endl;
	//cout << "original cnt:" << eq.c << endl;
	//cout << "Expressions for the difference:" << endl;
	//outputMatrix(eq);
	/*//test the correctness of eq by assigning a correct value for vars
	var[0] = accurateRoundSBoxDiff[r0][1];
	var[1] = accurateRoundSBoxDiff[r0][2];
	for (int i = 0; i < 5; i++) {
		var[3 * i + 2] = accurateRoundSBoxDiff[r0 + 1 + i][0];
		var[3 * i + 2 + 1] = accurateRoundSBoxDiff[r0 + 1 + i][1];
		var[3 * i + 2 + 2] = accurateRoundSBoxDiff[r0 + 1 + i][2];
	}
	for (int i = 0; i < bs; i++) {
		result[i] = eq.ma[i][eq.c - 1];
		for (int j = 0; j < eq.c - 1; j++) {
			if (eq.ma[i][j] == 1) {
				result[i] ^= var[j];
			}
		}
	}
	//compare result with accurateRoundDiff[r0+r1-2]
	cout << "diff:" << endl;
	for (int i = 0; i < bs; i++) {
		cout << result[i];
	}
	cout << endl;
	for (int i = 0; i < bs; i++) {
		cout << accurateRoundDiff[r0+r1-2][i];
	}
	cout << endl;
	system("pause");*/

	//matrix connectPart;
	connectPart.r = bs - 3 * m;
	connectPart.c = eq.c+connectPart.r;//store an identity matrix
	int varNum = eq.c - 1;
	clearMatrix(connectPart);
	for (int i = 3 * m; i < bs; i++) {
		for (int j = 0; j < eq.c; j++) {
			connectPart.ma[i - 3 * m][j] = eq.ma[i][j];
		}
		connectPart.ma[i - 3 * m][eq.c + i - 3 * m] = 1;
	}
	//cout << "connect part:" << endl;
	//outputMatrix(connectPart);

	//store original expression
	//matrix original;
	original.r = eq.r;
	original.c = eq.c;
	for (int i = 0; i < eq.r; i++) {
		for (int j = 0; j < eq.c; j++) {
			original.ma[i][j] = eq.ma[i][j];
		}
	}

	gauss(connectPart,varNum);
	//cout << "gauss on connectPart:" << endl;
	//outputMatrix(connectPart);

	//further simplify
	int index = 0;
	bool find = false;
	for (int i = 0; i < connectPart.r; i++) {
		find = false;
		for (int j = index; j < varNum; j++) {
			if (connectPart.ma[i][j]) {
				index = j;
				find = true;
				break;
			}
		}
		if (find) {//it is 1 in connectPart[i][index], eliminate 1 above
			for (int k = 0; k < i; k++) {
				if (connectPart.ma[k][index]) {//add row[i] to row [k]
					for (int t = i; t < connectPart.c; t++) {
						connectPart.ma[k][t] = connectPart.ma[k][t] ^ connectPart.ma[i][t];
					}
				}
			}
			index++;
		}
	}
	//cout << "further simplify connectPart:" << endl;
	//outputMatrix(connectPart);
	//cout << "row: " << connectPart.r << "; col:" << connectPart.c << endl;
	//system("pause");
	//change the size of connectPart
	connectPart.c = varNum + 1;

	//record transform A making original to connectPart.
	//matrix A;
	A.r = connectPart.r;
	A.c = connectPart.r;
	for (int i = 0; i < A.r; i++) {
		for (int j = 0; j < A.r; j++) {
			A.ma[i][j] = connectPart.ma[i][varNum + 1 + j];
		}
	}

	/*//test the method to compute the solution
	bool iscorrect = true;
	//output the correct solution
	cout << "corret sol:" << endl;
	for (int i = 0; i < 17; i++) {
		cout << var[i];
	}
	cout << endl;

	bool* tested = new bool[17];
	for (int i = 0; i < 17; i++) {
		tested[i] = accurateRoundDiff[r0 + r1 - 2][i + 3];
	}
	matrixMul(A, tested);
	for (int i = 0; i < 17; i++) {
		connectPart.ma[i][connectPart.c - 1] ^= tested[i];
	}
	//solve the equations
	vector<vector<bool> >sol;
	sol.clear();
	int solNum = 0;
	storeSolutions(sol, connectPart, solNum);
	if (solNum == 0) {
		cout << "no solution" << endl;
	}
	else {
		bool initialCheck = true;
		for (int i = 0; i < solNum; i++) {
			iscorrect = true;
			for (int j = 0; j < sol[i].size(); j++) {
				cout << sol[i][j];
				if (sol[i][j] != var[j]) {
					iscorrect = false;
				}
			}
			if (iscorrect) {
				cout << " true" << endl;
				//check sol[i]
				if (checkSolutionsForward(sol[i], r0, r1)) {
					cout << "final correct" << endl;
				}
			}
			else {
				cout << "wrong" << endl;
			}
		}
	}
	for (int i = 0; i < solNum; i++) {
		sol[i].clear();
	}
	sol.clear();
	delete[]tested;
	system("pause");*/
}

bool LowMC::checkSolutionsForward(vector<bool>& sol, int r0, int r1) {
	bool stateDiff[20];
	bool isCorrect = true;
	for (int k = 0; k < 3 * m; k++) {
		stateDiff[k] = startPart.ma[k][startPart.c - 1];
		for (int t = 0; t < startPart.c - 1; t++) {
			if (startPart.ma[k][t]) {
				stateDiff[k] ^= sol[t];
			}
		}
	}
	for (int k = 3 * m; k < bs; k++) {
		stateDiff[k] = startPart.ma[k][startPart.c - 1];
	}

	//apply linear transform
	int varIndex = 0;
	for (int r = r0; r < r0 + r1 - 2; r++) {
		matrixMul(linear[r], stateDiff);
		int a = stateDiff[0];
		int b = stateDiff[1];
		int c = stateDiff[2];
		a = a * 4 + b * 2 + c;
		varIndex = 2 + 3 * (r - r0);
		int d = sol[varIndex];
		int e = sol[varIndex + 1];
		int f = sol[varIndex + 2];
		d = d * 4 + e * 2 + f;
		if (table[a][d] == 0) {
			isCorrect = false;
			//cout << "the diff after Sbox after " << r+1 << " rounds is wrong" << endl;
			break;
		}
		else {//next round (update stateDiff)
			stateDiff[0] = sol[varIndex];
			stateDiff[1] = sol[varIndex+1];
			stateDiff[2] = sol[varIndex + 2];
		}
	}
	return isCorrect;
}

int LowMC::connectDiffs(bool diff[],int r0,int r1) {
	//use connectPart to connect diff[3m,...,bs]
	bool *truncated,*result;
	truncated = new bool[bs - 3 * m];
	result = new bool[bs - 3 * m];
	//apply A to truncated;
	for (int i = 3 * m; i < bs; i++) {
		truncated[i - 3 * m] = diff[i];
	}

	matrixMul(A, truncated, result);
	//check result 
	int validNum=0;
	bool *current;
	current = new bool[bs-3*m];
	for (int i = 0; i < bs - 3 * m; i++) {
		current[i] = connectPart.ma[i][connectPart.c - 1];
		connectPart.ma[i][connectPart.c - 1] ^= result[i];//connect diff
	}
	bool* stateDiff;
	stateDiff = new bool[bs];

	vector<vector<bool> >sol;
	sol.clear();
	int solNum = 0;
	storeSolutions(sol, connectPart, solNum);
	bool check = true;
	if (solNum == 0) {
		validNum=0;
	}

	else {
		for (int i = 0; i < solNum; i++) {
			check = true;
			for (int k = 0; k < 3; k++) {
				stateDiff[k] = original.ma[k][original.c-1];
				for (int t = 0; t < original.c-1; t++) {
					if (original.ma[k][t]) {
						stateDiff[k] ^= sol[i][t];
					}
				}
			}

			//check stateDiff[0,1,2] with diff[0,1,2]
			//we may also use additional equations
			//and this phase will be skipped
			int a = stateDiff[0];
			int b = stateDiff[1];
			int c = stateDiff[2];
			a = a * 4 + b * 2 + c;
			int d = diff[0];
			int e = diff[1];
			int f = diff[2];
			d = d * 4 + e * 2 + f;
			if (table[a][d]==0) {
				check = false;
			}
			
			if (check) {//further check the forward diff
				check = checkSolutionsForward(sol[i], r0, r1);
				if (check) {
					validNum++;
				}
			}
		}
	}
	for (int i = 0; i < solNum; i++) {
		sol[i].clear();
	}
	sol.clear();

	//return to the initial state
	for (int i = 0; i < bs - 3 * m; i++) {
		connectPart.ma[i][connectPart.c - 1]=current[i];
	}

	delete[]truncated;
	delete[]result;
	delete[]current;
	delete[]stateDiff;
	return validNum;
}

void LowMC::enumerateBackwardDiffR2(bool diff[], int r0, int r1,int r2,UINT64 result[]) {
	bool diffOut[256][256];
	bool diffIn[256][256];
	bool record[256][256];
	//matrixMul(invLinear[r - 1], diff, diffRound[0]);//reverse linear
	bool output[256];
	for (int i = 0; i < bs; i++) {
		output[i] = diff[i];
	}
	int cnt[3] = {0,0,0};
	int inacNum[10] = { 0 };

	computeNext(cnt,inacNum,r0,r1,r2, output, 0, r-1, diffOut,diffIn);

	cout << "Total enumerated trails:" << cnt[1] << endl;
	cout << "Total valid trails:" << cnt[0] << endl;
	cout << "The distribution of #(inactive Sboxes) in all the valid trails:" << endl;
	int sum = 0;
	for (int i = 0; i < 10; i++) {
		cout << i << " : " << inacNum[i] << endl;
	}
	cout << "The expected time to recover the key:" << cnt[2] << endl;

	result[0] = cnt[0];//valid diffs
	result[1] = cnt[1];//total diffs
	result[2] = cnt[2];//key
}

void LowMC::computeNext(int cnt[],int inacNum[], int r0, int r1,int r2,bool output[], int index, int roundNum, bool diffOut[][256],bool diffIn[][256]) {
	if (index == r2) {
		//output and check
		int inactive = 0;
		for (int i = 0; i < 10; i++) {
			if (diffIn[i][0] == 0 && diffIn[i][1] == 0 && diffIn[i][2] == 0) {
				inactive++;
			}
		}
		cnt[1]++;//count the number of trails enumerated backwards

		bool isCorrect = true;
		//compare diffIn[index-1] with accurateEachRound 
		//to check whether the correct one is skipped
		for (int i = 0; i < 10; i++) {
			//compare diffIn[i] with accurateRoundDiff[r-2-i]
			for (int j = 0; j < 20; j++) {
				if (diffIn[i][j] != accurateRoundDiff[r - 2 - i][j]) {
					isCorrect = false;
					break;
				}
			}
			if (!isCorrect) {
				break;
			}
		}
		bool tested[256];
		if (isCorrect) {//it is the accurate one
			matrixMul(invLinear[roundNum], diffIn[index - 1], tested);
			int find = connectDiffs(tested, r0, r1);
			if (find > 0) {
				cout << "The solution for the accurate trail is found!"<<endl;
				cnt[0] += find;
				inacNum[inactive] += find;
				int cost = (1 << (2*inactive));
				cnt[2] = cnt[2] + cost*find;
			}
		}
		else {//even if it is not the accurate one, do the same to filter
			matrixMul(invLinear[roundNum], diffIn[index - 1], tested);
			int find = connectDiffs(tested, r0, r1);
			if (find>0) {
				cnt[0] += find;
				inacNum[inactive] += find;
				int cost = (1 << (2 * inactive));
				cnt[2] = cnt[2] + cost*find;
			}
		}
		return;
	}

	matrixMul(invLinear[roundNum], output, diffOut[index]);//reverse linear
	if (diffOut[index][0] == 0 && diffOut[index][1] == 0 && diffOut[index][2] == 0) {
		diffIn[index][0] = 0;
		diffIn[index][1] = 0;
		diffIn[index][2] = 0;
		for (int i = 3; i < 20; i++) {
			diffIn[index][i] = diffOut[index][i];
		}
		computeNext(cnt,inacNum,r0,r1,r2,diffIn[index], index + 1, roundNum - 1, diffOut,diffIn);
	}
	else {
		int a = diffOut[index][0];
		int b = diffOut[index][1];
		int c = diffOut[index][2];
		a = a * 4 + b * 2 + c;
		//exhaust all possible input difference
		for (int i = 0; i < 4; i++) {
			diffIn[index][0] = (ddt[a][i] >> 2) & 0x1;
			diffIn[index][1] = (ddt[a][i] >> 1) & 0x1;
			diffIn[index][2] = ddt[a][i] & 0x1;
			for (int i = 3; i < 20; i++) {
				diffIn[index][i] = diffOut[index][i];
			}
			computeNext(cnt,inacNum,r0,r1,r2, diffIn[index], index + 1, roundNum - 1, diffOut, diffIn);
		}
	}
}

void LowMC::returnLinearRelationForwards(bool x0, bool x1, bool x2,matrix &matrix,int sNum,int &cnt) {
	//(z0,z1,z2) = s(x0,x1,x2)
	if (x0 == 0 && x1 == 0 && x2 == 1) {
		matrix.ma[sNum * 3][cnt] = 1;
		matrix.ma[sNum * 3+1][cnt+1] = 1;
		matrix.ma[sNum * 3 + 2][matrix.c - 1] = 1;
		cnt = cnt + 2;//two variables are added
	}
	else if (x0 == 0 && x1 == 1 && x2 == 0) {
		matrix.ma[sNum * 3][cnt] = 1;
		matrix.ma[sNum * 3 + 1][matrix.c-1] = 1;
		matrix.ma[sNum * 3 + 2][cnt+1] = 1;
		cnt = cnt + 2;//two variables are added
	}
	else if (x0 == 1 && x1 == 0 && x2 == 0) {
		matrix.ma[sNum * 3][matrix.c-1] = 1;
		matrix.ma[sNum * 3 + 1][cnt] = 1;
		matrix.ma[sNum * 3 + 2][cnt+1] = 1;
		cnt = cnt + 2;//two variables are added
	}
	else if (x0 == 0 && x1 == 1 && x2 == 1) {
		matrix.ma[sNum * 3][cnt] = 1;

		matrix.ma[sNum * 3+1][cnt+1] = 1;
		matrix.ma[sNum * 3 + 2][cnt + 1] = 1;
		matrix.ma[sNum * 3 + 2][matrix.c - 1] = 1;
		cnt = cnt + 2;
	}
	else if (x0 == 1 && x1 == 1 && x2 == 0) {
		matrix.ma[sNum * 3][cnt] = 1;
		matrix.ma[sNum * 3 + 1][cnt] = 1;
		matrix.ma[sNum * 3 + 1][matrix.c - 1] = 1;

		matrix.ma[sNum * 3 + 2][cnt + 1] = 1;
		cnt = cnt+2;
	}
	else if (x0 == 1 && x1 == 0 && x2 == 1) {
		matrix.ma[sNum * 3][cnt] = 1;
		matrix.ma[sNum * 3 + 1][cnt + 1] = 1;
		matrix.ma[sNum * 3 + 2][cnt] = 1;
		matrix.ma[sNum * 3 + 2][matrix.c - 1] = 1;
		cnt = cnt + 2;
	}
	else {//1,1,1
		matrix.ma[sNum * 3][cnt] = 1;
		matrix.ma[sNum * 3 + 1][cnt + 1] = 1;

		matrix.ma[sNum * 3 + 2][cnt] = 1;
		matrix.ma[sNum * 3 + 2][cnt+1] = 1;
		matrix.ma[sNum * 3 + 2][matrix.c-1] = 1;
		cnt = cnt + 2;
	}
}

//the full s-box layer
void LowMC::chooseBestInputDiff(bool bestDiff[]) {
	//check matrix linear[0]
	bool* output;
	output = new bool[bs];
	int max = 0, maxPos = 0, maxVal = 0;
	for (int i = 0; i < m; i++) {
		for (int t = 0; t < bs; t++) {
			bestDiff[t] = 0;
		}
		for (int val = 1; val < 8; val++) {
			//only best[3*i,3*i+1,3*i+2] is active
			bestDiff[3 * i] = (val >> 2) & 0x1;
			bestDiff[3 * i+1] = (val >> 1) & 0x1;
			bestDiff[3 * i+2] = val & 0x1;

			matrixMul(linear[0], bestDiff, output);
			//count #(Inactive Sbox) in output
			int cnt = 0;
			for (int t = 0; t < m; t++) {
				if (output[3 * t] == 0
					&& output[3 * t + 1] == 0
					&& output[3 * t + 2] == 0) {
					cnt++;
				}
			}
			if (cnt > max) {
				max = cnt;
				maxPos = i;
				maxVal = val;
			}
		}
	}
	cout << "max:" << max << endl;
	cout << "maxPos:" << maxPos << endl;
	cout << "maxVal:" << maxVal << endl;
	//assign value to bestDiff
	for (int t = 0; t < bs; t++) {
		bestDiff[t] = 0;
	}
	//only best[3*i,3*i+1,3*i+2] is active
	bestDiff[3 * maxPos] = (maxVal >> 2) & 0x1;
	bestDiff[3 * maxPos + 1] = (maxVal >> 1) & 0x1;
	bestDiff[3 * maxPos + 2] = maxVal & 0x1;

	delete[]output;
}

int LowMC::getInactiveNum(bool state[], int sboxCnt) {
	int cnt = 0;
	for (int i = 0; i < sboxCnt; i++) {
		if (state[i * 3] == 0
			&& state[i * 3 + 1] == 0
			&& state[i * 3 + 2] == 0) {
			cnt++;
		}
	}
	return cnt;
}

void LowMC::findNext(int requiredNum,int secondAcNum,int &total, int &totalTime,int boxIndex, bool input[], bool result[], matrix& expression, bool correctOne[]) {
	if (boxIndex == m) {
		//check result
		//total++;
		//enumerate diffs via solving equations
		bool test[21];
		matrixMul(linear[1], result, test);
		int inactive = getInactiveNum(test, m);
		int cnt =constructDiffEquations(test, expression, correctOne);
		//compute the complexity
		//count the active Sboxes
		int time = 0;
		//requiredNum (eqs)
		int eqNum = 2 * m + 2 * secondAcNum;
		time = 1 << (2 * inactive);//in actual attacks, this can be smaller as we do not need to fully linearize it
		if (eqNum < requiredNum) {
			int b = 1 << (requiredNum - eqNum);
			time = time * b;
		}
		//cout << time << endl;
		time = time * cnt;
		totalTime += time;
		total+=cnt;
		return;
	}
	else {
		int t = boxIndex * 3;
		if (input[t] == 0 && input[t + 1] == 0 && input[t + 2] == 0) {
			result[t] = 0;
			result[t + 1] = 0;
			result[t + 2] = 0;
			findNext(requiredNum, secondAcNum, total, totalTime, boxIndex + 1, input, result, expression, correctOne);
		}
		else {
			int a = input[t];
			int b = input[t + 1];
			int c = input[t + 2];
			a = a * 4 + b * 2 + c;
			//exhaust all possible input difference
			for (int i = 0; i < 4; i++) {
				result[t] = (ddt[a][i] >> 2) & 0x1;
				result[t+1] = (ddt[a][i] >> 1) & 0x1;
				result[t+2] = ddt[a][i] & 0x1;
				findNext(requiredNum, secondAcNum, total, totalTime, boxIndex + 1, input, result, expression, correctOne);
			}
		}
	}
}

int LowMC::startTestingFullSBoxLayer(bool finalDiff[],bool inputDiff[], bool correctOne[],int requiredNum,int secondNum,int q,int t) {
	matrix expression;
	constructExpressions(finalDiff, expression);
	//enumerating the input diffs of the third round
	bool result[21];
	int total = 0;
	int boxIndex = 0;
	int totalTime = 0;

	int exp0 = pow(2, 0.858 * m - 2 * t);
	exp0 = exp0 * pow(2, 2 * (m-q))+ pow(2, 2 * (m - q));
	int exp1 = 0;
	if (5 * t <= 6 * m + ks + 2 - 2 * q) {
		exp1 = pow(2, 3 * m - 2 * q - 2 * t);
	}
	else {
		exp1 = pow(2, 3 * t - 2);
	}

	int isSmaller = 1;
	findNext(requiredNum,secondNum,total,totalTime, boxIndex, inputDiff, result, expression, correctOne);
	cout << "Time to enumerate diffs: 0x" << hex<<total<<", expected: 0x"<<exp0;
	if (total <= exp0) {
		cout << ", smaller than expected" << endl;
	}
	else {
		isSmaller = 0;
		cout << ", larger than expected" << endl;
		//system("pause");
	}
	cout << "Time to recover key (upper bound): 0x" <<hex<< totalTime<< ", expected:0x" << exp1 << endl;
	return isSmaller;
}

void LowMC::constructExpressions(bool diff[], matrix& eq) {
	int varNum = 0;
	eq.r = bs;
	eq.c = 2 * m;
	bool* con = new bool [bs];//the constant part in the expression
	for (int i = 0; i < bs; i++) {
		con[i] = 0;
	}
	clearMatrix(eq);

	//apply inverse to diff -> diff becomes the output diff of s-box
	matrixMul(invLinear[r - 1], diff);

	int inactive = 0;
	//assign variables based on diff
	for (int i = 0; i < m; i++) {//deal with each s-box
		int a = diff[3 * i];
		int b = diff[3 * i + 1];
		int c = diff[3 * i + 2];
		a = a * 4 + b * 2 + c;
		if (a == 0) {//(0,0,0)
			//do nothing
			inactive++;
		}
		else if (a == 1) {//(0,0,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i+1][varNum+1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;//two variables are added
		}
		else if (a == 2) {//(0,1,0)
			eq.ma[3 * i][varNum] = 1;
			con[3 * i + 1] = 1;
			eq.ma[3 * i+2][varNum+1] = 1;
			varNum += 2;
		}
		else if (a == 4) {//(1,0,0)
			con[3 * i] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			eq.ma[3 * i + 2][varNum+1] = 1;
			varNum += 2;
		}
		else if (a == 3) {//(0,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i+1][varNum + 1] = 1;
			eq.ma[3 * i+2][varNum + 1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
		else if (a == 6) {//(1,1,0)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum] = 1;
			con[3 * i + 1] = 1;
			eq.ma[3 * i + 2][varNum + 1] = 1;
			varNum += 2;
		}
		else if (a == 5) {//(1,0,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
		else {//(1,1,1)
			eq.ma[3 * i][varNum] = 1;
			eq.ma[3 * i + 1][varNum + 1] = 1;
			eq.ma[3 * i + 2][varNum] = 1;
			eq.ma[3 * i + 2][varNum+1] = 1;
			con[3 * i + 2] = 1;
			varNum += 2;
		}
	}
	//cout << "varNum: " << varNum << endl;
	//cout << "inactive:" << inactive << endl;
	//change eq
	eq.c = varNum + 1;
	for (int i = 0; i < bs; i++) {
		eq.ma[i][eq.c - 1] = con[i];
	}
	original.r = eq.r;
	original.c = eq.c;
	for (int i = 0; i < original.r; i++) {
		for (int j = 0; j < original.c; j++) {
			original.ma[i][j] = eq.ma[i][j];
		}
	}
	//outputMatrix(original);

	//apply inverse
	matrixMul(invLinear[r - 2], eq);

	delete[]con;
}

int LowMC::constructDiffEquations(bool diff[], matrix& expression,bool correctOne[]) {
	matrix eq;
	int currentRow = 0;
	eq.c = expression.c;
	eq.r = bs;
	clearMatrix(eq);
	bool iszero = true;
	for (int i = 0; i < m; i++) {//deal with each s-box
		iszero = true;
		for (int j = 3 * i; j < 3 * i + 3; j++) {
			//if diff[j]=1, add expression at row j to the current row of eq
			if (diff[j]) {
				iszero = false;
				for (int t = 0; t < expression.c; t++) {
					eq.ma[currentRow][t] ^= expression.ma[j][t];
				}
			}
		}
		if (iszero) {//the input diff is zero
			//add three rows of expression to eq
			for (int j = 3 * i; j < 3 * i + 3; j++) {
				for (int t = 0; t < expression.c; t++) {
					eq.ma[currentRow][t] ^= expression.ma[j][t];
				}
				currentRow++;
			}
		}
		else {//the input diff is not zero
			eq.ma[currentRow][eq.c - 1] ^= 1;//add the constant
			currentRow++;//one equation
		}
	}
	eq.r = currentRow;
	//cout << "The eqs to connect diffs:" << endl;
	//outputMatrix(eq);
	//cout << "number of eqs:" << currentRow << endl;
	//cout << "number of unknowns:" << eq.c-1 << endl;

	gauss(eq);
	//outputMatrix(eq);
	vector<vector<bool> >sol;
	sol.clear();
	int solNum = 0;
	storeSolutions(sol, eq, solNum);
	if (solNum == 0) {
		//cout << "infesible" << endl;
	}
	else {
		//cout << "number of solutions:" << solNum << endl;
	}

	//recover the actual diff, use expression to compute the diff
	bool diffInBox[21];
	bool isCorrect = true;
	for (int i = 0; i < solNum; i++) {
		for (int k = 0; k < bs; k++) {
			diffInBox[k] = original.ma[k][original.c - 1];
			for (int t = 0; t < original.c - 1; t++) {
				if (original.ma[k][t]) {
					diffInBox[k] = diffInBox[k]^sol[i][t];
				}
			}
		}
		//check whether the correct is in the solution set
		isCorrect = true;
		for (int k = 0; k < bs; k++) {
			if (diffInBox[k] != correctOne[k]) {
				isCorrect = false;
			}
		}
		if (isCorrect) {
			cout << "The correct one is being enumerated!" << endl;
		}
	}

	for (int i = 0; i < solNum; i++) {
		sol[i].clear();
	}
	sol.clear();
	return solNum;
}

void LowMC::outputMatrix(matrix ma) {
	for (int i = 0; i < ma.r; i++) {
		for (int j = 0; j < ma.c; j++) {
			cout << ma.ma[i][j] << " ";
		}
		cout << endl;
	}
}

void LowMC::gauss(matrix &eqSys) {
	int variableNum = eqSys.c - 1;
	bool isFirst = false;
	int targetRow = 0;

	for (int i = 0; i < variableNum; i++) {
		isFirst = true;
		for (int j = targetRow; j < eqSys.r; j++) {
			if (isFirst && eqSys.ma[j][i]) {
				isFirst = false;
				swap(eqSys.ma[j], eqSys.ma[targetRow]);
				targetRow++;
			}
			else {
				if (eqSys.ma[j][i]) {//apply Gauss
					for (int k = i; k < eqSys.c; k++) {
						eqSys.ma[j][k] ^= eqSys.ma[targetRow - 1][k];
					}
				}
			}
		}
	}
}

void LowMC::gauss(matrix& eqSys,int col) {
	int variableNum = col;
	bool isFirst = false;
	int targetRow = 0;

	for (int i = 0; i < variableNum; i++) {
		isFirst = true;
		for (int j = targetRow; j < eqSys.r; j++) {
			if (isFirst && eqSys.ma[j][i]) {
				isFirst = false;
				swap(eqSys.ma[j], eqSys.ma[targetRow]);
				targetRow++;
			}
			else {
				if (eqSys.ma[j][i]) {//apply Gauss
					for (int k = i; k < eqSys.c; k++) {
						eqSys.ma[j][k] ^= eqSys.ma[targetRow - 1][k];
					}
				}
			}
		}
	}
}

void LowMC::storeSolutions(vector<vector<bool> >& sol,matrix &eqSys, int &solNum) {
	vector<int> lead;
	vector<int> freebits;
	freebits.clear();
	lead.clear();
	bool* isFree;
	isFree = new bool[eqSys.c - 1];
	memset(isFree, 1, eqSys.c - 1);

	int start = 0;
	for (int r = 0; r < eqSys.r; r++) {
		while (start < eqSys.c - 1 && eqSys.ma[r][start] == 0) {
			start++;
		}
		if (start == eqSys.c - 1) {
			break;
		}
		lead.push_back(start);
		isFree[start] = false;
		start++;
	}

	if (lead.size() < eqSys.r) {
		for (int j = lead.size(); j < eqSys.r; j++) {
			if (eqSys.ma[j][eqSys.c - 1] != 0) {
				solNum = 0;
				return;
			}
		}
	}

	for (int i = 0; i < eqSys.c - 1; i++) {
		if (isFree[i]) {
			freebits.push_back(i);
		}
	}
	//cout << "free size:" << freebits.size() << endl;
	//cout << "lead size:" << lead.size() << endl;*/

	vector<bool> eachsol;
	eachsol.clear();
	eachsol.resize(eqSys.c - 1);
	int solSize = 1 << freebits.size();
	for (int i = 0; i < solSize; i++) {
		for (int j = 0; j < freebits.size(); j++) {
			eachsol[freebits[j]] = (i >> j) & 0x1;
		}
		for (int k = lead.size()-1; k>=0; k--) {
			//compute eachsol[lead[k]] use row= k
			eachsol[lead[k]] = eqSys.ma[k][eqSys.c-1];
			for (int j = lead[k] + 1; j < eqSys.c-1; j++) {
				if (eqSys.ma[k][j] == 1) {
					eachsol[lead[k]] = eachsol[lead[k]]^eachsol[j];
				}
			}
		}
		solNum++;
		sol.push_back(eachsol);
	}

	delete[]isFree;
	freebits.clear();
	lead.clear();
	eachsol.clear();
}

void LowMC::clearMatrix(matrix &ma) {
	for (int i = 0; i<ma.r; i++) {
		for (int j = 0; j<ma.c; j++) {
			ma.ma[i][j] = 0;
		}
	}
}

void LowMC::matrixMul(matrix &m, bool x[], bool y[]) {
	for (int i = 0; i < m.r; i++) {
		y[i] = 0;
		for (int j = 0; j < m.c; j++) {
			y[i] = y[i] ^ (m.ma[i][j] & x[j]);
		}
	}
}

void LowMC::matrixMul(matrix& m, bool x[]) {
	bool y[256];
	for (int i = 0; i < m.r; i++) {
		y[i] = 0;
		for (int j = 0; j < m.c; j++) {
			y[i] = y[i] ^ (m.ma[i][j] & x[j]);
		}
	}
	for (int i = 0; i < m.r; i++) {
		x[i] = y[i];
	}
}

void LowMC::matrixMul(matrix& m1, matrix& m2) {
	matrix m3;
	//cout << m1.r << " " << m1.c <<" "<<m2.c<< endl;
	for (int i = 0; i < m1.r; i++) {
		for (int j = 0; j < m2.c; j++) {
			m3.ma[i][j] = 0;
			for (int k = 0; k < m1.c; k++) {
				m3.ma[i][j] = m3.ma[i][j] ^ (m1.ma[i][k] & m2.ma[k][j]);
			}
		}
	}
	//outputMatrix(m3);
	for (int i = 0; i < m2.r; i++) {
		for (int j = 0; j < m2.c; j++) {
			m2.ma[i][j] = m3.ma[i][j];
		}
	}
}