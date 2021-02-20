#ifndef _LOWMC_H_
#define _LOWMC_H_

#include <vector>
using namespace std;
typedef unsigned long long UINT64;

const int ddt[8][4] = {
{0,0,0,0},
{1,3,5,7},
{2,3,6,7},
{1,2,5,6},
{4,6,5,7},
{1,4,3,6},
{4,2,5,3},
{1,2,4,7}
};

const bool table[8][8] = {
	1,0,0,0,0,0,0,0,
	0,1,0,1,0,1,0,1,
	0,0,1,1,0,0,1,1,
	0,1,1,0,0,1,1,0,
	0,0,0,0,1,1,1,1,
	0,1,0,1,1,0,1,0,
	0,0,1,1,1,1,0,0,
	0,1,1,0,1,0,0,1
};

struct matrix {
	int r;
	int c;
	bool ma[256][256];
};

class LowMC {
private:
	int bs;//blocksize;
	int ks;//keysize
	int m;//# sbox
	int r;//rounds
	matrix *linear;
	matrix *invLinear;
	matrix constant;
	matrix *keyMa;
	bool **diffEq;
	bool **keyEq;
	bool *outputDiff;
	matrix A;//used for connection
	matrix connectPart;
	int zeroRowStart = 0;
	matrix startPart;
	matrix original;
	bool accurateRoundDiff[256][20];
	bool accurateRoundSBoxDiff[256][20];
public:
	LowMC(int bs,int ks,int m,int r);
	void loadFileNonFull();
	void loadFileFull();
	void matrixMul(matrix &m, bool x[], bool y[]);
	void matrixMul(matrix& m, bool x[]);
	void matrixMul(matrix& m1, matrix& m2);
	void encrypt(bool p[], bool k[],int psize, int ksize,bool c[],int rounds,bool eachRoundOut[][20], bool eachRoundSBoxOut[][20]);
	void encryptFull(bool p[], bool k[], int psize, int ksize, bool c[], int rounds, bool eachRoundOut[][21], bool eachRoundSBoxOut[][21]);
	void encrypt(bool p[], bool k[], int psize, int ksize, bool c[], int rounds);
	void decrypt(bool p[], bool k[], int psize, int ksize, bool c[],int rounds);
	void setAccurateRoundDiff(bool eachRoundDiff[][20], bool eachRoundSBoxDiff[][20]);
	//verify the case when there is only 1 sbox
	void computeR0();//r0=6
	void constructForwardDiffR1(bool diff[],int r1,int r0);//r1=7,r0=6
	void enumerateBackwardDiffR2(bool diff[], int r0,int r1, int r2,UINT64 result[]);//r2=10;
	void computeNext(int cnt[],int inacNum[], int r0,int r1,int r2, bool output[], int index, int roundNum, bool diffOut[][256], bool diffIn[][256]);
	void returnLinearRelationForwards(bool x0, bool x1, bool x2,matrix &matrix,int sNum,int &cnt);
	int connectDiffs(bool diff[],int r0,int r1);
	void gauss(matrix &eqSys);
	void gauss(matrix& eqSys,int col);
	void storeSolutions(vector<vector<bool> >&sol,matrix &ma,int &solNum);
	bool checkSolutionsForward(vector<bool>& sol,int r0,int r1);

	//verify the case when a full sbox layer is used
	void chooseBestInputDiff(bool bestDiff[]);
	void findNext(int requiredNum,int secondAcNum,int &total,int &totalTime, int boxIndex, bool input[], bool result[], matrix& expression, bool correctOne[]);
	int startTestingFullSBoxLayer(bool finalDiff[],bool inputDiff[],bool correctOne[], int requiredNum, int secondNum,int q,int t);
	void constructExpressions(bool diff[],matrix &eq);
	int constructDiffEquations(bool diff[],matrix &expression,bool correctOne[]);

	//output
	void clearMatrix(matrix &ma);
	void outputMatrix(matrix ma);
	int getInactiveNum(bool state[], int sboxCnt);
};

#endif
