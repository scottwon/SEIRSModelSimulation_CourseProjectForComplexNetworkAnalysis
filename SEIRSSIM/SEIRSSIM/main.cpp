#include<iostream>
#include<fstream>
#include<cstdlib>
#include<ctime>
#include<cmath>
using namespace std;
const int nPerson = 10000;
const double pAdj = 0.08;
const int totalTime = 1000;
const double pS = 0.8, pE = 0.1, pI = 0.05, pR = 0.05;
const int timeLim = 2;
const int NInfo1 = 100;
const int NInfo2 = 300;
const double resistence = 0.03;
const double forgetRate = 0.04;
bool AdjMat[nPerson][nPerson];
double RelMat[nPerson][nPerson];
double SimMat[nPerson][nPerson];
double Info1[nPerson], Info2[nPerson];
int LastTime[nPerson];
int NTimes[nPerson];
int Status[nPerson][totalTime];
int degree[nPerson];
double rnd[nPerson][totalTime][3];
ofstream o("curve.txt");
void init() {
	memset(AdjMat, 0, sizeof(bool)*nPerson*nPerson);
	memset(degree, 0, sizeof(int)*nPerson);
	srand((int)time(0));
	for (int i = 0;i < nPerson;i++) {
		double r1 = ((double)rand()) / RAND_MAX;
		double r2 = ((double)rand()) / RAND_MAX;
		Info1[i] = r1*2-1;
		Info2[i] = r2*2-1;
		LastTime[i] = -1;
		NTimes[i] = 0;
	}
	srand((int)time(0));
	for (int i = 0;i < nPerson-1;i++) {
		for (int j = i + 1;j < nPerson;j++) {
			double r3 = ((double)rand()) / RAND_MAX;
			double r4 = ((double)rand()) / RAND_MAX;
			RelMat[i][j] = r4;
			RelMat[j][i] = r4;
			double tmp = (Info1[i] * Info1[j] + Info2[i] * Info2[j]) / (sqrt(Info1[i]*Info1[i]+Info2[i]*Info2[i])*sqrt(Info1[j]*Info1[j]+Info2[j]*Info2[j]));
			SimMat[i][j] = tmp;
			SimMat[j][i] = tmp;
			if (r3 <= pAdj) {
				AdjMat[i][j] = 1;
				AdjMat[j][i] = 1;
				degree[i]++;
				degree[j]++;
			}
		}
	}
	srand((int)time(0));
	for (int i = 0;i < nPerson;i++) {
		for (int j = 0;j < totalTime;j++) {
			double rn1= ((double)rand()) / RAND_MAX;
			double rn2= ((double)rand()) / RAND_MAX;
			double rn3 = ((double)rand()) / RAND_MAX;
			rnd[i][j][0] = rn1;
			rnd[i][j][1] = rn2;
			rnd[i][j][2] = rn3;
		}
	}
	for (int i = 0;i < nPerson*pS;i++) {
		Status[i][0] = 1;
	}
	for (int i = nPerson*pS;i < nPerson*(pS+pE);i++) {
		Status[i][0] = 2;
	}
	for (int i = nPerson*(pS+pE);i < nPerson*(pS+pE+pI);i++) {
		Status[i][0] = 3;
	}
	for (int i = nPerson*(pS+pE+pI);i < nPerson;i++) {
		Status[i][0] = 4;
	}
}
void output(int t) {
	for (int i = 0;i < nPerson;i++) {
		cout << Status[i][t] << " ";
	}
	cout << endl;
}
double simSR(int t) {
	int cnt = 0;
	double avg = 0.0;
	for (int i = 0;i < nPerson;i++) {
		for (int j = 0;j < nPerson;j++) {
			if (i == j)continue;
			if (Status[i][t] == 1 && Status[j][t] == 4) {
				avg += SimMat[i][j];
				cnt++;
			}
		}
	}
	avg /= cnt;
	return avg;
}
double avgDegree(int t) {
	int cnt = 0;
	double avg = 0.0;
	for (int i = 0;i < nPerson;i++) {
		if (Status[i][t] != 1)continue;
		cnt++;
		for (int j = 0;j < nPerson;j++) {
			if (i == j)continue;
			if (Status[j][t] == 1 && AdjMat[i][j]==1) {
				avg++;
			}
		}
	}
	avg /= cnt;
	return avg;
}
void stat(int t) {
	int nS = 0, nE = 0, nI = 0, nR = 0;
	for (int i = 0;i < nPerson;i++) {
		if (Status[i][t] == 1) {
			nS++;
		}
		else if (Status[i][t] == 2) {
			nE++;
		}
		else if (Status[i][t] == 3) {
			nI++;
		}
		else if (Status[i][t] == 4) {
			nR++;
		}
	}
	o<<t<<" "<< nS << " " << nE << " " << nI << " " << nR << endl;
}
void calc() {
	cout << "CALC"<<endl;
	stat(0);
	for (int t = 1;t < totalTime;t++) {
		cout << "t=" << t << endl;
		for (int i = 0;i < nPerson;i++) {
			if (Status[i][t - 1] == 1) {
				int potentialE = 0, potentialI = 0, potentialR = 0,totalI=0;
				for (int j = 0;j < nPerson;j++) {
					if (AdjMat[i][j]==0)continue;
					if (Status[j][t - 1] != 3)continue;
					totalI++;
					if (RelMat[i][j] > 0.5 && SimMat[i][j] > 0.5) {
						potentialI++;
					}
					else if (RelMat[i][j] > 0.5) {
						potentialE++;
					}
					else if (SimMat[i][j] <= 0) {
						potentialR++;
					}
				}
				double pot = rnd[i][t][0];
				double pe, pi, pr;
				pe = (double)potentialE / (double)degree[i];
				pi = (double)potentialI / (double)degree[i];
				pr = (double)potentialR / (double)degree[i];
				if (pot < pe) {
					Status[i][t] = 2;
					LastTime[i] = t;
					NTimes[i] = 1;
				}
				else if (pot < (pe + pi)) {
					Status[i][t] = 3;
				}
				else if (pot < (pe + pi + pr)) {
					Status[i][t] = 4;
				}
				else {
					Status[i][t] = 1;
				}
			}
			else if (Status[i][t - 1] == 2) {
				for (int j = 0;j < nPerson;j++) {
					if (AdjMat[i][j] == 0)continue;
					if (Status[j][t - 1] == 3) {
						NTimes[i]++;
					}
				}
				if (t - LastTime[i] >= timeLim && NTimes[i]<NInfo1) {
					Status[i][t] = 4;
					NTimes[i] = 0;
					LastTime[i] = -1;
				}
				else if (t-LastTime[i]<=timeLim && NTimes[i] >= NInfo2) {
					Status[i][t] = 3;
					NTimes[i] = 0;
					LastTime[i] = -1;
				}
				else {
					Status[i][t] = 2;
					if (t - LastTime[i] == timeLim) {
						NTimes[i] = 0;
						LastTime[i] = t;
					}
				}
			}
			else if (Status[i][t - 1] == 3) {
				int numI = 0;
				int numR = 0;
				for (int j = 0;j < nPerson;j++) {
					if (AdjMat[i][j] == 0)continue;
					if (Status[j][t - 1] == 3)numI++;
					if (Status[j][t - 1] == 4)numR++;
				}
				if (numI>0.95*degree[i] ||numR>0.95*degree[i]) {
					Status[i][t] = 4;
				}
				else if(rnd[i][t][2]<resistence) {
					Status[i][t] = 4;
				}
				else {
					Status[i][t] = 3;
				}
			}
			else {
				double temp = rnd[i][t][1];
				if (temp <= forgetRate) {
					Status[i][t] = 1;
				}
				else {
					Status[i][t] = 4;
				}
			}
		}
		stat(t);
		//cout << simSR(t) << endl;
		//cout << avgDegree(t) << endl;
/*		if ((t + 1) % 100 == 0) {
			output(t);
		}*/
	}
}

int main() {
	init();
	calc();
	system("pause");
	return 0;
}