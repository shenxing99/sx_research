#include<iostream>
#include<stdlib.h>
#include<fstream>
#include<time.h>
#include<math.h>
#include<vector>
#include<array> 
#include<set>
#include<map>
#include<string>

using namespace std;

#define N 300  //����ڵ����� 
#define maximum 25 //���Ž� 
#define Itermax 100000   //���������� 
#define RunTime 30  //�������еĴ��� 

int total_iteration_times;

string filename = "p-hat300-2 (25).txt";

array< array<int, N>, N > matrix;//�ڽӾ��� 
array< array<int, N>, N > matrix_temp;
vector< vector<int> > adjList;
vector<int> adjList_line;
int deg[N];
int maxdegree;
int mindegree;
double avedegree;
int edgenumber;

array< int, N > globalSolution;//ȫ�����Ž� 
int globalSize;//ȫ�����Ž��нڵ�ĸ��� 
array< int, N > solution;
array<int, N> km;
array<int, N> ke;
array<int, N> kd;

vector<int> NS0;
vector<int> NS1;
vector<int> NS2;
vector<int> NSmore;
vector<int> NSk;
array< int, N > tabuList;//���tabuList[i]=0����˵���ڵ�iû�б����� 

vector<int> newNS0;
vector<int> newNS1;
vector<int> newNS2;
vector<int> newNSmore;

vector<int> kVertexOutFromS;
vector<int> S;
vector<int> VcutS;

array< array<int, N>, RunTime > solu; // RunTime�����У�ÿһ�εõ��Ľ⼯
int totalCandidate[RunTime]; // RunTime�����У�ÿһ�εõ��Ľ⼯��Ԫ�ظ��������õ��Ķ�������Ԫ�ظ���
double totalTime[RunTime]; // RunTime�����У�ÿһ�ε�����ʱ��
double averageCandidate;// RunTime�����У���������С��ƽ��ֵ 
double stdCandidate;// RunTime�����У���������С�ı�׼��
int maxCandidate;// RunTime��������candidate�����ֵ 
int minCandidate;// RunTime��������candidate����Сֵ 
int maxCandidateNumber;// RunTime�������еõ����candidate�Ĵ���
double averageTime;// RunTime�����е�ƽ������ʱ��
double runningTime; // ÿһ�����е�����ʱ�� 
int run_times; // ��¼���д��� run_times = run+1 


//********************************************  ��������  **********************************************
void globalInitialVariable();
void initialSettings();
void generateInitialSolution();
void initialTabuList();
void obtainVcutS();
void computeKM();
void computeKE();
void computeKD();
void computeNSset();
void getKvalue_NSk();
void intensificationStep();
void zero_one_swap(int vertex);
int computeVerNeighborKEvalue(int vertex);
void checkTabuTime_ObtainNewNSset();
void one_one_swap(int vertex);
void updateTabuList_oneOneSwap(int k, int vertex);
void diversificationStep();
void strongPerturbation_1();
void strongPerturbation_2();
void k_one_swap(int vertex);
void updateTabuList_Kmore();
void weakPerturbation();

void read();
void read_without_e();
void readMatrix();
void complement_graph();
void obtainAdjList();
void obtainVerDegree();
void saveResults();
void printResults();
double computeAverage();
double computeStandardDeviation();
int computeMaxCandidate();
int computeMinCandidate();
int computeMaxTimes();
double computeAverageTime();
int computeSolutionSize();
void printSolution();
void printS();

//********************************************  main  **********************************************
int main(int argc, char* argv[]) {
	int i;
	int j;
	clock_t start, end;
	clock_t start1, end1;
	globalInitialVariable();

	read();
	complement_graph();
	//read_without_e();
	//readMatrix();

	obtainAdjList();
	obtainVerDegree();

	int run;
	for (run = 0; run<RunTime; run++) {
		srand((unsigned int)time(NULL));
		run_times = run + 1;
		start = clock();
		initialSettings();
		generateInitialSolution();//������ʼ��
		int initialSolutionSize = computeSolutionSize();
		globalSolution = solution;//����ȫ�����Ž�
		globalSize = computeSolutionSize();//ȫ�����Ž�Ľڵ����
		initialTabuList();
		getKvalue_NSk();

		for (int iter = 0; iter < Itermax; iter++) {
			checkTabuTime_ObtainNewNSset();//�õ�ֻ����û�б����ɵĽڵ��NS����
			int newNS0size = newNS0.size();
			int newNS1size = newNS1.size();
			int newNS01size = newNS0size + newNS1size;
			if (newNS01size != 0) {
				//�����Խ�����ǿ����
				intensificationStep();//Kֵ��NSk�ͽ��ɱ�ĸ��¶����Ӻ�������
				int newSolutionSize = computeSolutionSize();
				if (newSolutionSize >= globalSize) {
					globalSolution = solution;//����ȫ�����Ž�
					globalSize = computeSolutionSize();//ȫ�����Ž�Ľڵ����
				}
			}
			else {
				diversificationStep();//Kֵ��NSk�ͽ��ɱ�ĸ��¶����Ӻ�������
			}
			int globalSolutionSize = computeSolutionSize();
			if(globalSolutionSize==maximum){
				cout<<"iter = "<<iter<<endl;
				total_iteration_times = iter;
				iter = Itermax;
			}
		}
		//������ϣ��õ�Itermax��ѭ�������������Ž�globalSolution�����СglobalSize
		end = clock();
		runningTime = (double)(end - start) / CLOCKS_PER_SEC;
		solution = globalSolution; 
		solu[run] = solution; // ��run�����еõ��Ľ⼯ 
		int finalSolutionSize = computeSolutionSize();
		totalCandidate[run] = finalSolutionSize; // ��run�����е����Ž� 
		totalTime[run] = runningTime; // ��run�����е�����ʱ�� 
		printf("run = %d, initialSolutionSize = %d, finalSolutionSize= %d \n", run_times, initialSolutionSize, finalSolutionSize);
	}
	printf("!!!code end!!!");
	saveResults();
	printResults();

	system("pause");
	return 0;
}
/********************  �Ӻ���  ****************************/

void globalInitialVariable() {
	int i;
	int j;
	for (i = 0; i<N; i++) {
		for (j = 0; j<N; j++) {
			matrix[i][j] = 0;
			matrix_temp[i][j] = 0;
		}
	}
	for (i = 0; i<N; i++) {
		adjList.push_back(adjList_line);
		deg[i] = 0;
	}
	maxdegree = 0;
	mindegree = 0;
	avedegree = 0;
	edgenumber = 0;
	for (i = 0; i<RunTime; i++) {
		totalCandidate[i] = 0;
		for (j = 0; j<N; j++) {
			solu[i][j] = 0;
		}
	}
	averageCandidate = 0;
	stdCandidate = 0;
	maxCandidate = 0;
	minCandidate = 0;
	maxCandidateNumber = 0;
	averageTime = 0;
	runningTime = 0;
}

void initialSettings() {
	int i;
	for (i = 0; i<N; i++) {
		km[i] = N;
		ke[i] = N;
		kd[i] = N;
		tabuList[i] = 0;//��ʼʱ��ÿ���ڵ㶼�������� 
		solution[i] = 0;
		globalSolution[i] = 0;
	}
	NS0.clear();
	NS1.clear();
	NS2.clear();
	NSmore.clear();
	NSk.clear();
	S.clear();//��ʼʱ���⼯Ϊ�� 
	for (auto j = 0; j<N; j++) {
		VcutS.push_back(j);//��ʼʱ��V\SΪ����ڵ��ȫ�� 
	}
}

void generateInitialSolution() {
	S.clear();
	vector<int> candidateV;//���Խ���S�Ľڵ�ļ��� 
	array<int, N> flag;
	for (auto k = 0; k<N; k++) {
		candidateV.push_back(k);
	}
	for (int i = 0; i<N; i++) {
		flag[i] = 1;
	}
	while (candidateV.size() != 0) {
		int r = (int)(candidateV.size()*rand()) / (RAND_MAX + 1);
		int vertex = candidateV[r];
		S.push_back(vertex);
		flag[vertex] = 0;
		for (auto j = 0; j<adjList[vertex].size(); j++) {
			int vertexNeighbor = adjList[vertex][j];
			flag[vertexNeighbor] = 0;
		}
		candidateV.clear();
		for (int i = 0; i<N; i++) {
			if (flag[i] == 1) {
				candidateV.push_back(i);
			}
		}
	}
	for (auto i = 0; i < S.size(); i++) {
		int vertex = S[i];
		solution[vertex] = 1;
	}
}

void initialTabuList() {
	for (int i = 0; i < N; i++) {
		tabuList[i] = 0;//��ʼʱ��ÿ���ڵ㶼û�н���ʱ��
	}
}
void obtainVcutS() {
	VcutS.clear();
	for (int i = 0; i < N; i++) {
		if (solution[i] == 0) {
			VcutS.push_back(i);
		}
	}
}

void computeKM() {
	obtainVcutS();
	for (int i = 0; i < N; i++) {
		km[i] = N;
	}
	for (auto i = 0; i < VcutS.size(); i++) {
		int vertex = VcutS[i];
		int kmValue = 0;
		for (auto j = 0; j < adjList[vertex].size(); j++) {
			int vertexNeighbor = adjList[vertex][j];
			if (solution[vertexNeighbor] == 1) {
				kmValue = kmValue + 1;
			}
		}
		km[vertex] = kmValue;
	}
}

void computeKE() {
	obtainVcutS();
	for (int i = 0; i < N; i++) {
		ke[i] = N;
	}
	for (auto i = 0; i < S.size(); i++) {
		int vertex = S[i];
		int keValue = 0;
		for (auto j = 0; j < adjList[vertex].size(); j++) {
			int vertexNeighbor = adjList[vertex][j];
			if ((solution[vertexNeighbor] == 0) && (km[vertexNeighbor] == 1)) {
				keValue = keValue + 1;
			}
		}
		ke[vertex] = keValue;
	}
}

void computeKD() {
	obtainVcutS();
	for (int i = 0; i < N; i++) {
		kd[i] = N;
	}
	for (auto i = 0; i < VcutS.size(); i++) {
		int vertex = VcutS[i];
		int kdValue = 0;
		for (auto j = 0; j < adjList[vertex].size(); j++) {
			int vertexNeighbor = adjList[vertex][j];
			if (solution[vertexNeighbor] == 0) {
				kdValue = kdValue + 1;
			}
		}
		kd[vertex] = kdValue;
	}
}

void computeNSset() {
	NS0.clear();
	NS1.clear();
	NS2.clear();
	NSmore.clear();
	for (int i = 0; i < N; i++) {
		if (km[i] == 0) {
			NS0.push_back(i);
		}
		else if (km[i] == 1) {
			NS1.push_back(i);
		}
		else if (km[i] == 2) {
			NS2.push_back(i);
		}
		else if ((km[i] > 2) && (km[i] != N)) {
			NSmore.push_back(i);
		}
	}
}

void getKvalue_NSk() {
	obtainVcutS();
	computeKM();
	computeKE();
	computeKD();
	computeNSset();
}

void intensificationStep() {
	int choosedVertex;
	auto newNS0size = newNS0.size();
	if (newNS0size != 0) {
		//���Խ���0-1swap
		int r = (int)(newNS0size*rand()) / (RAND_MAX + 1);
		int vertex = newNS0[r];
		zero_one_swap(vertex);
	}
	else {
		//����1-1swap
		auto NS1size = NS1.size();
		auto NS2size = NS2.size();
		auto NSmoresize = NSmore.size();
		if (NS1size > NS2size + NSmoresize) {
			for (vector<int>::iterator iter=NS1.begin(); iter < NS1.end(); iter++) {
				int vertex = *iter;
				for (auto j = 0; j < adjList[vertex].size(); j++) {
					int vertexNeighbor = adjList[vertex][j];
					if ((solution[vertexNeighbor] == 1)&&(ke[vertexNeighbor]==1)) {
						NS1.erase(iter);
						break;
					}
				}
				//�����ɾ���п������﷨����
			}
		}
		//else {
			//����ѡ����򣬴�NS1(Ӧ�ô�newNS1��ѡ����Ϊ������Ĳ���û�б����ɵĽڵ�)��ѡ��һ���ڵ����1-1swap
			vector<int> NS_;
			NS_.clear();
			//�ҵ�NS1�еĽڵ���S�е��ھӽڵ��keֵ���Ľڵ㣬������NS_��
			int maxVer = newNS1[0];
			int max = computeVerNeighborKEvalue(maxVer);
			for (auto j = 0; j < newNS1.size(); j++) {
				int vertex = newNS1[j];
				int vertexNeighborKEvalue = computeVerNeighborKEvalue(vertex);
				if (vertexNeighborKEvalue > max) {
					max = vertexNeighborKEvalue;
					maxVer = vertex;
				}
			}
			//�õ�NS_
			for (auto j = 0; j < newNS1.size(); j++) {
				int vertex = newNS1[j];
				int vertexNeighborKEvalue = computeVerNeighborKEvalue(vertex);
				if (vertexNeighborKEvalue == max) {
					NS_.push_back(vertex);
				}
			}
			//��NS_��ѡ��һ���ڵ����1-1swap
			if (NS_.size() == 1) {
				choosedVertex = NS_[0];
				one_one_swap(choosedVertex);
			}
			else {
				//NS_�в�ֹһ���ڵ㣬���NS_�еĽڵ���ѡ��kd���Ľڵ�
				int maxKDver = NS_[0];
				int maxKD = kd[maxKDver];
				for (auto i = 0; i < NS_.size(); i++) {
					int vertex = NS_[i];
					if (kd[vertex] > maxKD) {
						maxKD = kd[vertex];
						maxKDver = vertex;
					}
				}
				choosedVertex = maxKDver;
				one_one_swap(choosedVertex);
			}
		
	}
}

void zero_one_swap(int vertex) {
	int solutionSize1 = computeSolutionSize();
	//���ڵ�vertex����⼯������Ҫ�ӽ⼯���Ƴ��ڵ�
	solution[vertex] = 1;
	S.push_back(vertex);
	obtainVcutS();
	//���¸���Kֵ��NSk����
	getKvalue_NSk();
	int solutionSize2 = computeSolutionSize();
}

int computeVerNeighborKEvalue(int vertex) {
	//����NS1�еĽڵ�vertex��solution�е��ھӽڵ��KEֵ
	int verNeighborKEvalue;
	for (auto i = 0; i < adjList[vertex].size(); i++) {
		int vertexNeighbor = adjList[vertex][i];
		if (solution[vertexNeighbor] == 1) {
			verNeighborKEvalue = ke[vertexNeighbor];
		}
	}
	return verNeighborKEvalue;
}

void checkTabuTime_ObtainNewNSset() {
	newNS0.clear();
	newNS1.clear();
	newNS2.clear();
	newNSmore.clear();
	for (auto i = 0; i < NS0.size(); i++) {
		int vertex = NS0[i];
		if (tabuList[vertex] == 0) {
			//�ýڵ�Ľ��ɳ���Ϊ0��˵��û�б�����
			newNS0.push_back(vertex);
		}
	}
	for (auto i = 0; i < NS1.size(); i++) {
		int vertex = NS1[i];
		if (tabuList[vertex] == 0) {
			//�ýڵ�Ľ��ɳ���Ϊ0��˵��û�б�����
			newNS1.push_back(vertex);
		}
	}
	for (auto i = 0; i < NS2.size(); i++) {
		int vertex = NS2[i];
		if (tabuList[vertex] == 0) {
			//�ýڵ�Ľ��ɳ���Ϊ0��˵��û�б�����
			newNS2.push_back(vertex);
		}
	}
	for (auto i = 0; i < NSmore.size(); i++) {
		int vertex = NSmore[i];
		if (tabuList[vertex] == 0) {
			//�ýڵ�Ľ��ɳ���Ϊ0��˵��û�б�����
			newNSmore.push_back(vertex);
		}
	}
}

void one_one_swap(int vertex) {
	//vertex��NS1�еĽڵ㣬���ڵ�vertex���������S��ͬʱ��S�н�vertex��һ���ھ��Ƴ�
	//���ȣ��ҵ�vertex��S�е�һ���ھӽڵ㣬������ڵ��S���Ƴ�
	int vertexNeighbor_inS;
	for (auto i = 0; i < adjList[vertex].size(); i++) {
		int vertexNeighbor = adjList[vertex][i];
		if (solution[vertexNeighbor] == 1) {
			vertexNeighbor_inS = vertexNeighbor;
			break;
		}
	}
	//���½⼯solution S VcutS
	solution[vertex] = 1;
	solution[vertexNeighbor_inS] = 0;
	S.clear();
	for (int i = 0; i < N; i++) {
		if (solution[i] == 1) {
			S.push_back(i);
		}
	}
	obtainVcutS();
	//���¸���Kֵ��NSk����
	getKvalue_NSk();
	//���½��ɱ�
	updateTabuList_oneOneSwap(1,vertexNeighbor_inS);
}

void updateTabuList_oneOneSwap(int k,int vertex) {
	//���ȣ���ԭ�����ڽ���ʱ��Ľڵ㣬�����ʱ���1
	for (int i = 0; i < N; i++) {
		if (tabuList[i] != 0) {
			tabuList[i] = tabuList[i] - 1;
		}
	}
	//Ȼ�󣬸��ݽ��н��������ࣨkֵ��������S���Ƴ��Ľڵ�vertex�Ľ���ʱ������Ϊtt
	if (k == 1) {
		int NS1size = NS1.size();
		int NS2size = NS2.size();
		int NSmoresize = NSmore.size();
		if (NS1size < NS2size + NSmoresize) {
			int r = (int)(NS1size*rand()) / (RAND_MAX + 1);
			int tt = 10 + r;
			tabuList[vertex] = tt;
		}
	}
	else if (k > 1) {
		int tt = 7;
		tabuList[vertex] = tt;
	}
}

void diversificationStep() {
	int NS1size = NS1.size();
	int NS2size = NS2.size();
	int NSmoresize = NSmore.size();
	if (NS1size > NS2size + NSmoresize) {
		strongPerturbation_1();//��newNSmore
	}
	else {
		double r = (double)(rand()) / (RAND_MAX + 1);
		if (r >= 0.5) {
			strongPerturbation_2();//��NSmore,�����ǽڵ��Ƿ񱻽���
		}
		else {
			weakPerturbation();//��newNS2
		}
	}
}

void strongPerturbation_1() {
	//��newNSmore
	int maxKDver = newNSmore[0];
	int maxKD = kd[maxKDver];
	for (auto i = 0; i < NSmore.size(); i++) {
		int vertex = NSmore[i];
		if (kd[vertex] > maxKD) {
			maxKD = kd[vertex];
			maxKDver = vertex;
		}
	}
	int choosedVertex = maxKDver;
	int k = km[maxKDver];
	k_one_swap(choosedVertex);
}

void strongPerturbation_2() {
	//�����ǽڵ�Ľ���ʱ�䣬��NSmore
	int r = (int)(NSmore.size()*rand()) / (RAND_MAX + 1);
	int choosedVertex = NSmore[r];
	k_one_swap(choosedVertex);
}

void k_one_swap(int vertex) {
	//��vertex����⼯����Ҫ�ӽ⼯���Ƴ�k���ڵ�
	//�ӽ⼯�н��ڵ�vertex���ھӽڵ��Ƴ�
	kVertexOutFromS.clear();//��¼�ӽ⼯���Ƴ�ȥ��k���ڵ㣬���ڸ��½��ɱ�
	for (auto i = 0; i < adjList[vertex].size(); i++) {
		int vertexNeighbor = adjList[vertex][i];
		if (solution[vertexNeighbor] == 1) {
			solution[vertexNeighbor] = 0;
			kVertexOutFromS.push_back(vertexNeighbor);
		}
	}
	//���ڵ�vertex����⼯
	solution[vertex] = 1;
	//����S VcutS
	S.clear();
	for (int i = 0; i < N; i++) {
		if (solution[i] == 1) {
			S.push_back(i);
		}
	}
	obtainVcutS();
	//���¸���Kֵ��NSk����
	getKvalue_NSk();
	//���½��ɱ�
	updateTabuList_Kmore();
}

void updateTabuList_Kmore() {
	//���ȣ���ԭ�����ڽ���ʱ��Ľڵ㣬�����ʱ���1
	for (int i = 0; i < N; i++) {
		if (tabuList[i] != 0) {
			tabuList[i] = tabuList[i] - 1;
		}
	}
	//Ȼ�󣬽�����k-1 swap������S���Ƴ���k���ڵ�Ľ���ʱ������Ϊtt
	int tt = 7;
	for (auto i = 0; i < kVertexOutFromS.size(); i++) {
		int vertex = kVertexOutFromS[i];
		tabuList[vertex] = tt;
	}
}

void weakPerturbation() {
	//��newNS2
	int maxKDver = newNS2[0];
	int maxKD = kd[maxKDver];
	for (auto i = 0; i < NS2.size(); i++) {
		int vertex = NS2[i];
		if (kd[vertex] > maxKD) {
			maxKD = kd[vertex];
			maxKDver = vertex;
		}
	}
	int choosedVertex = maxKDver;
	int k = km[maxKDver];
	k_one_swap(choosedVertex);
}



//********************************************  �Ӻ�������  **********************************************
//***************  1  ��ȡ�ļ�������ıߵı�ʾ��  **************************
void read() // �õ��������ڽӾ��� �� ����¼�˱ߵľ���edgeMatrix[M][2] 
{
	int i;
	FILE* fp = fopen(filename.c_str(), "r");
	if (!fp)
	{
		printf("�ļ��򿪴���");
	}

	char e;
	int startnode;
	int endnode;

	i = 0;
	while (!feof(fp))
	{
		fscanf(fp, "%c %d %d \n", &e, &startnode, &endnode);
		matrix[startnode - 1][endnode - 1] = 1;
		matrix[endnode - 1][startnode - 1] = 1;
		i = i + 1;
	}

	fclose(fp);
}

//***************  1  ��ȡ�ļ�������ıߵı�ʾ��  **************************
void read_without_e() // �õ��������ڽӾ��� �� ����¼�˱ߵľ���edgeMatrix[M][2] 
{
	int i;
	FILE* fp = fopen(filename.c_str(), "r");
	if (!fp)
	{
		printf("�ļ��򿪴���");
	}

	int startnode;
	int endnode;

	i = 0;
	while (!feof(fp))
	{
		fscanf(fp, "%d %d \n", &startnode, &endnode);
		matrix[startnode - 1][endnode - 1] = 1;
		matrix[endnode - 1][startnode - 1] = 1;
		i = i + 1;
	}

	fclose(fp);
}

//***************  2  ��ȡ�ڽӾ����ʾ������  **************************
void readMatrix()
{
	FILE* fp = fopen(filename.c_str(), "r");
	if (fp == NULL)
	{
		printf("�ļ���ʧ��");
		exit(1);
	}
	char C;
	int row = 0;
	int col = 0;
	int i1 = 0;//����matrix���±���Ҫÿ�ε���1����������������������������i1��j1
	int j1 = 0;
	while ((C = fgetc(fp)) != EOF)
	{
		if (C == '0')
		{
			matrix[i1][j1] = 0;
			col++;
			j1++;
		}
		else if (C == '1')
		{
			matrix[i1][j1] = 1;
			col++;
			j1++;
		}
		else if (C == ' ')
		{
			col++;
		}
		else if (C == '\n')
		{
			col++;
			row++;
			i1++;
			j1 = 0;
		}
	}
	fclose(fp);
}

//***************************  5  ������Ĳ�ͼ  *****************************
void complement_graph()
{
	int i;
	int j;
	for (i = 0; i<N; i++)
	{
		for (j = 0; j<N; j++)
		{
			matrix_temp[i][j] = matrix[i][j];
		}
	}
	for (i = 0; i<N; i++)
	{
		for (j = 0; j<N; j++)
		{
			if (matrix_temp[i][j] == 1)
			{
				matrix[i][j] = 0;
			}
			else if (matrix_temp[i][j] == 0)
			{
				matrix[i][j] = 1;
			}
		}
	}
	for (i = 0; i<N; i++)
	{
		matrix[i][i] = 0;
	}
}

//************************************  6  �õ������ڽӱ�  ************************************
void obtainAdjList()
{
	int i;
	int j;

	for (i = 0; i<N; i++)
	{
		adjList[i].clear();
	}
	for (i = 0; i<N; i++)
	{
		for (j = 0; j<N; j++)
		{
			if (matrix[i][j] == 1)
			{
				adjList[i].push_back(j); //  �õ����ڽӱ��У�ÿһ���ڵ���ھӽڵ㶼�Ա�ŵ��������� 
			}
		}
	}
}

//********************************  7  �õ�����Ķ���Ϣ  *****************************
void obtainVerDegree()
{
	int i;

	for (i = 0; i<N; i++)
	{
		deg[i] = adjList[i].size();
	}

	int max = deg[0];
	for (i = 0; i<N; i++)
	{
		if (max<deg[i])
		{
			max = deg[i];
		}
	}
	maxdegree = max; // ���� 

	int min = deg[0];
	for (i = 0; i<N; i++)
	{
		if (min>deg[i])
		{
			min = deg[i];
		}
	}
	mindegree = min; // ��С��  

	double sumdegree = 0;
	for (i = 0; i<N; i++)
	{
		sumdegree = sumdegree + deg[i];
	}
	avedegree = (double)((double)sumdegree / (double)N); // ƽ���� 

	edgenumber = (int)(sumdegree / 2); // �ܱ��� 
}


//**********************************  26  ��ʵ�������ı���ʽ���д洢  ************************************
void saveResults()
{
	int i;
	int j;
	averageCandidate = computeAverage();
	stdCandidate = computeStandardDeviation();
	maxCandidate = computeMaxCandidate();
	minCandidate = computeMinCandidate();
	maxCandidateNumber = computeMaxTimes();
	averageTime = computeAverageTime();

	string name = "result_" + filename;
	ofstream outfile_tm(name, std::ios::trunc);

	//ofstream outfile_tm("newIGLS_result.txt");
	if (!outfile_tm)
	{
		cerr << "open error!" << endl;
		exit(1);
	}

	outfile_tm << "The results: " << endl << endl;
	outfile_tm << "1" << endl;
	outfile_tm << "Network information:" << endl;
	outfile_tm << "network name    : " << filename.c_str() << endl;
	outfile_tm << "network vertices:" << N << endl;
	outfile_tm << "total run times : " << RunTime << endl;
	outfile_tm << "edge number     : " << edgenumber << endl;
	outfile_tm << "maximum degree  : " << maxdegree << endl;
	outfile_tm << "minimum degree  : " << mindegree << endl;
	outfile_tm << "average degree  : " << avedegree << endl << endl;

	outfile_tm << "2" << endl;
	outfile_tm << "The Parameters:" << endl;
	outfile_tm << "Itermax = " << Itermax << endl;

	outfile_tm << "3" << endl;
	outfile_tm << "The MIS:" << endl;
	outfile_tm << "averageCandidate  : " << averageCandidate << endl;
	outfile_tm << "stdCandidate      : " << stdCandidate << endl;
	outfile_tm << "maxCandidate      : " << maxCandidate << endl;
	outfile_tm << "minCandidate      : " << minCandidate << endl;
	outfile_tm << "maxCandidateNumber: " << maxCandidateNumber << endl;
	outfile_tm << "averageTime       : " << averageTime << endl;
	outfile_tm << "total_iteration_times��" << total_iteration_times <<endl<<endl; 

	outfile_tm << "4" << endl;
	outfile_tm << "The detailed results:" << endl;
	outfile_tm << "total candidate of RunTime: " << endl;
	for (i = 0; i<RunTime; i++)
	{
		outfile_tm << totalCandidate[i] << " ";
		if ((i + 1) % 10 == 0) {
			outfile_tm << endl;
		}
	}
	outfile_tm << endl;
	outfile_tm << "the solution of RunTime: " << endl;
	for (i = 0; i<RunTime; i++)
	{
		for (j = 0; j<N; j++)
		{
			if (solu[i][j] == 1)
			{
				outfile_tm << j << " ";
			}
		}
		outfile_tm << endl;
	}
	outfile_tm << endl;


}

//*******************************  27  ������  ***************************************** 
void printResults()
{
	int i;

	averageCandidate = computeAverage();
	stdCandidate = computeStandardDeviation();
	maxCandidate = computeMaxCandidate();
	minCandidate = computeMinCandidate();
	maxCandidateNumber = computeMaxTimes();
	averageTime = computeAverageTime();

	printf("\n*********************  ������ %d ��  *************************��\n\n", RunTime);
	cout << "Network name is " << filename << endl;
	printf("Network vertices = %d, run times = %d\n\n", N, RunTime);
	printf("1��Parameters:\n");
	printf("Itermax = %d\n", Itermax);


	printf("2���õ��Ľ⼯��С�ֱ�Ϊ��\n");
	for (i = 0; i<RunTime; i++)
	{
		printf("%d ", totalCandidate[i]);
		if ((i + 1) % 10 == 0) {
			cout << endl;
		}
	}
	printf("\n\n");

	printf("3���⼯��С��ƽ��ֵΪ��%f����׼��Ϊ��%f \n\n", averageCandidate, stdCandidate);

	printf("4���⼯��С���ֵΪ��%d����СֵΪ��%d, ȡ�����ֵ�Ĵ���Ϊ��%d\n\n", maxCandidate, minCandidate, maxCandidateNumber);

	printf("5��ƽ������ʱ��Ϊ��%f \n\n", averageTime);
	
	printf("6���ܵ���������%d\n",total_iteration_times);

	/*printf("1���õ��Ľ⼯�е�Ԫ�طֱ�Ϊ��\n");
	for (i = 0; i<RunTime; i++)
	{
	for (j = 0; j<N; j++)
	{
	if (solu[i][j] == 1)
	{
	printf("%d ", j);
	}
	}
	printf("\n");
	}
	printf("\n");*/

	/*printf("6��ÿһ�����е���Ӧ�Ⱥ���ֵ��");
	int f[RunTime];
	for( i=0; i<RunTime; i++ )
	{
	f[i] = computeCandidate( i );
	}
	for( i=0; i<RunTime; i++ )
	{
	printf("%d ", f[i]);
	}
	printf("\n\n");*/

}

//*******************************  20  ����RunTime�����еĽ⼯��Сƽ��ֵ  *****************************
double computeAverage() {
	int i;
	double average = 0;
	int sum = 0;
	for (i = 0; i<RunTime; i++) {
		sum = sum + totalCandidate[i];
	}
	average = (double)((double)sum / (double)RunTime);

	return average;
}

//******************************  21  ����RunTime�����еĽ⼯��С�ı�׼��  ******************************* 
double computeStandardDeviation() {
	int i;
	double average = computeAverage();
	double standardDeviation = 0;
	double sum = 0;
	for (i = 0; i<RunTime; i++) {
		sum = sum + (totalCandidate[i] - average)*(totalCandidate[i] - average);
	}
	standardDeviation = (double)(sqrt((double)sum / (double)RunTime));

	return standardDeviation;
}

//*****************************  22 ����RunTime�����У��õ�������  *************************************
int computeMaxCandidate()
{
	int i;
	int max = totalCandidate[0];
	for (i = 0; i<RunTime; i++) {
		if (max<totalCandidate[i]) {
			max = totalCandidate[i];
		}
	}

	return max;
}

//*****************************  23  ����RunTime�����У��õ�����С��  ***************************************
int computeMinCandidate()
{
	int i;
	int min = totalCandidate[0];
	for (i = 0; i<RunTime; i++) {
		if (min>totalCandidate[i]) {
			min = totalCandidate[i];
		}
	}
	return min;
}

//***************************  24  ����RunTime�����У����Ž���ֵĴ���  **************************************
int computeMaxTimes()
{
	int i;
	int max = computeMaxCandidate();
	int times = 0;
	for (i = 0; i<RunTime; i++) {
		if (totalCandidate[i] == max) {
			times = times + 1;
		}
	}
	return times;
}

//****************************  25  ����RunTime�����е�ƽ������ʱ��  *******************************************
double computeAverageTime()
{
	int i;
	double averageTime = 0;
	double sum = 0;
	for (i = 0; i<RunTime; i++) {
		sum = sum + totalTime[i];
	}
	averageTime = (double)((double)sum / (double)RunTime);

	return averageTime;
}

// *****************  15  ����solution�Ĵ�С   ************************
int computeSolutionSize()
{
	int i;
	int sizeSum = 0;
	for (i = 0; i<N; i++) {
		if (solution[i] == 1) {
			sizeSum = sizeSum + 1;
		}
	}
	return sizeSum;
}


void printSolution() {
	int i;
	for (i = 0; i < N; i++) {
		cout << solution[i] << ' ';
	}
	cout << endl;
}

void printS() {
	cout << "the vertices in S: ";
	for (auto i = 0; i < S.size(); i++) {
		cout << S[i] << " ";
	}
	cout << endl;
}


