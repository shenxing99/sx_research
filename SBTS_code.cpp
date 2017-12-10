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

#define N 300  //网络节点总数 
#define maximum 25 //最优解 
#define Itermax 100000   //最大迭代次数 
#define RunTime 30  //独立运行的次数 

int total_iteration_times;

string filename = "p-hat300-2 (25).txt";

array< array<int, N>, N > matrix;//邻接矩阵 
array< array<int, N>, N > matrix_temp;
vector< vector<int> > adjList;
vector<int> adjList_line;
int deg[N];
int maxdegree;
int mindegree;
double avedegree;
int edgenumber;

array< int, N > globalSolution;//全局最优解 
int globalSize;//全局最优解中节点的个数 
array< int, N > solution;
array<int, N> km;
array<int, N> ke;
array<int, N> kd;

vector<int> NS0;
vector<int> NS1;
vector<int> NS2;
vector<int> NSmore;
vector<int> NSk;
array< int, N > tabuList;//如果tabuList[i]=0，则说明节点i没有被禁忌 

vector<int> newNS0;
vector<int> newNS1;
vector<int> newNS2;
vector<int> newNSmore;

vector<int> kVertexOutFromS;
vector<int> S;
vector<int> VcutS;

array< array<int, N>, RunTime > solu; // RunTime次运行，每一次得到的解集
int totalCandidate[RunTime]; // RunTime次运行，每一次得到的解集的元素个数，即得到的独立集的元素个数
double totalTime[RunTime]; // RunTime次运行，每一次的运行时间
double averageCandidate;// RunTime次运行，独立集大小的平均值 
double stdCandidate;// RunTime次运行，独立集大小的标准差
int maxCandidate;// RunTime次运行中candidate的最大值 
int minCandidate;// RunTime次运行中candidate的最小值 
int maxCandidateNumber;// RunTime次运行中得到最大candidate的次数
double averageTime;// RunTime次运行的平均运行时间
double runningTime; // 每一次运行的运行时间 
int run_times; // 记录运行次数 run_times = run+1 


//********************************************  函数声明  **********************************************
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
		generateInitialSolution();//产生初始解
		int initialSolutionSize = computeSolutionSize();
		globalSolution = solution;//更新全局最优解
		globalSize = computeSolutionSize();//全局最优解的节点个数
		initialTabuList();
		getKvalue_NSk();

		for (int iter = 0; iter < Itermax; iter++) {
			checkTabuTime_ObtainNewNSset();//得到只考察没有被禁忌的节点的NS集合
			int newNS0size = newNS0.size();
			int newNS1size = newNS1.size();
			int newNS01size = newNS0size + newNS1size;
			if (newNS01size != 0) {
				//即可以进行增强操作
				intensificationStep();//K值、NSk和禁忌表的更新都在子函数里面
				int newSolutionSize = computeSolutionSize();
				if (newSolutionSize >= globalSize) {
					globalSolution = solution;//更新全局最优解
					globalSize = computeSolutionSize();//全局最优解的节点个数
				}
			}
			else {
				diversificationStep();//K值、NSk和禁忌表的更新都在子函数里面
			}
			int globalSolutionSize = computeSolutionSize();
			if(globalSolutionSize==maximum){
				cout<<"iter = "<<iter<<endl;
				total_iteration_times = iter;
				iter = Itermax;
			}
		}
		//搜索完毕，得到Itermax次循环搜索到的最优解globalSolution及其大小globalSize
		end = clock();
		runningTime = (double)(end - start) / CLOCKS_PER_SEC;
		solution = globalSolution; 
		solu[run] = solution; // 第run次运行得到的解集 
		int finalSolutionSize = computeSolutionSize();
		totalCandidate[run] = finalSolutionSize; // 第run次运行的最优解 
		totalTime[run] = runningTime; // 第run次运行的运行时间 
		printf("run = %d, initialSolutionSize = %d, finalSolutionSize= %d \n", run_times, initialSolutionSize, finalSolutionSize);
	}
	printf("!!!code end!!!");
	saveResults();
	printResults();

	system("pause");
	return 0;
}
/********************  子函数  ****************************/

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
		tabuList[i] = 0;//初始时，每个节点都不被禁忌 
		solution[i] = 0;
		globalSolution[i] = 0;
	}
	NS0.clear();
	NS1.clear();
	NS2.clear();
	NSmore.clear();
	NSk.clear();
	S.clear();//初始时，解集为空 
	for (auto j = 0; j<N; j++) {
		VcutS.push_back(j);//初始时，V\S为网络节点的全集 
	}
}

void generateInitialSolution() {
	S.clear();
	vector<int> candidateV;//可以进入S的节点的集合 
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
		tabuList[i] = 0;//初始时，每个节点都没有禁忌时间
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
		//可以进行0-1swap
		int r = (int)(newNS0size*rand()) / (RAND_MAX + 1);
		int vertex = newNS0[r];
		zero_one_swap(vertex);
	}
	else {
		//进行1-1swap
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
				//这里的删除有可能有语法问题
			}
		}
		//else {
			//按照选择规则，从NS1(应该从newNS1中选择，因为这里面的才是没有被禁忌的节点)中选择一个节点进行1-1swap
			vector<int> NS_;
			NS_.clear();
			//找到NS1中的节点在S中的邻居节点的ke值最大的节点，保存在NS_中
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
			//得到NS_
			for (auto j = 0; j < newNS1.size(); j++) {
				int vertex = newNS1[j];
				int vertexNeighborKEvalue = computeVerNeighborKEvalue(vertex);
				if (vertexNeighborKEvalue == max) {
					NS_.push_back(vertex);
				}
			}
			//从NS_中选择一个节点进行1-1swap
			if (NS_.size() == 1) {
				choosedVertex = NS_[0];
				one_one_swap(choosedVertex);
			}
			else {
				//NS_中不止一个节点，则从NS_中的节点中选择kd最大的节点
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
	//将节点vertex让入解集，不需要从解集中移除节点
	solution[vertex] = 1;
	S.push_back(vertex);
	obtainVcutS();
	//更新各种K值和NSk集合
	getKvalue_NSk();
	int solutionSize2 = computeSolutionSize();
}

int computeVerNeighborKEvalue(int vertex) {
	//计算NS1中的节点vertex在solution中的邻居节点的KE值
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
			//该节点的禁忌长度为0，说明没有被禁忌
			newNS0.push_back(vertex);
		}
	}
	for (auto i = 0; i < NS1.size(); i++) {
		int vertex = NS1[i];
		if (tabuList[vertex] == 0) {
			//该节点的禁忌长度为0，说明没有被禁忌
			newNS1.push_back(vertex);
		}
	}
	for (auto i = 0; i < NS2.size(); i++) {
		int vertex = NS2[i];
		if (tabuList[vertex] == 0) {
			//该节点的禁忌长度为0，说明没有被禁忌
			newNS2.push_back(vertex);
		}
	}
	for (auto i = 0; i < NSmore.size(); i++) {
		int vertex = NSmore[i];
		if (tabuList[vertex] == 0) {
			//该节点的禁忌长度为0，说明没有被禁忌
			newNSmore.push_back(vertex);
		}
	}
}

void one_one_swap(int vertex) {
	//vertex是NS1中的节点，将节点vertex放入独立集S，同时从S中将vertex的一个邻居移除
	//首先，找到vertex在S中的一个邻居节点，将这个节点从S中移除
	int vertexNeighbor_inS;
	for (auto i = 0; i < adjList[vertex].size(); i++) {
		int vertexNeighbor = adjList[vertex][i];
		if (solution[vertexNeighbor] == 1) {
			vertexNeighbor_inS = vertexNeighbor;
			break;
		}
	}
	//更新解集solution S VcutS
	solution[vertex] = 1;
	solution[vertexNeighbor_inS] = 0;
	S.clear();
	for (int i = 0; i < N; i++) {
		if (solution[i] == 1) {
			S.push_back(i);
		}
	}
	obtainVcutS();
	//更新各种K值和NSk集合
	getKvalue_NSk();
	//更新禁忌表
	updateTabuList_oneOneSwap(1,vertexNeighbor_inS);
}

void updateTabuList_oneOneSwap(int k,int vertex) {
	//首先，对原来存在禁忌时间的节点，其禁忌时间减1
	for (int i = 0; i < N; i++) {
		if (tabuList[i] != 0) {
			tabuList[i] = tabuList[i] - 1;
		}
	}
	//然后，根据进行交换的种类（k值），将从S中移除的节点vertex的禁忌时间设置为tt
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
		strongPerturbation_1();//用newNSmore
	}
	else {
		double r = (double)(rand()) / (RAND_MAX + 1);
		if (r >= 0.5) {
			strongPerturbation_2();//用NSmore,不考虑节点是否被禁忌
		}
		else {
			weakPerturbation();//用newNS2
		}
	}
}

void strongPerturbation_1() {
	//用newNSmore
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
	//不考虑节点的禁忌时间，用NSmore
	int r = (int)(NSmore.size()*rand()) / (RAND_MAX + 1);
	int choosedVertex = NSmore[r];
	k_one_swap(choosedVertex);
}

void k_one_swap(int vertex) {
	//将vertex让入解集，需要从解集中移除k个节点
	//从解集中将节点vertex的邻居节点移除
	kVertexOutFromS.clear();//记录从解集中移除去的k个节点，用于更新禁忌表
	for (auto i = 0; i < adjList[vertex].size(); i++) {
		int vertexNeighbor = adjList[vertex][i];
		if (solution[vertexNeighbor] == 1) {
			solution[vertexNeighbor] = 0;
			kVertexOutFromS.push_back(vertexNeighbor);
		}
	}
	//将节点vertex放入解集
	solution[vertex] = 1;
	//更新S VcutS
	S.clear();
	for (int i = 0; i < N; i++) {
		if (solution[i] == 1) {
			S.push_back(i);
		}
	}
	obtainVcutS();
	//更新各种K值和NSk集合
	getKvalue_NSk();
	//更新禁忌表
	updateTabuList_Kmore();
}

void updateTabuList_Kmore() {
	//首先，对原来存在禁忌时间的节点，其禁忌时间减1
	for (int i = 0; i < N; i++) {
		if (tabuList[i] != 0) {
			tabuList[i] = tabuList[i] - 1;
		}
	}
	//然后，进行了k-1 swap，将从S中移除的k个节点的禁忌时间设置为tt
	int tt = 7;
	for (auto i = 0; i < kVertexOutFromS.size(); i++) {
		int vertex = kVertexOutFromS[i];
		tabuList[vertex] = tt;
	}
}

void weakPerturbation() {
	//用newNS2
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



//********************************************  子函数定义  **********************************************
//***************  1  读取文件（网络的边的表示）  **************************
void read() // 得到了网络邻接矩阵 ， 并记录了边的矩阵edgeMatrix[M][2] 
{
	int i;
	FILE* fp = fopen(filename.c_str(), "r");
	if (!fp)
	{
		printf("文件打开错误");
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

//***************  1  读取文件（网络的边的表示）  **************************
void read_without_e() // 得到了网络邻接矩阵 ， 并记录了边的矩阵edgeMatrix[M][2] 
{
	int i;
	FILE* fp = fopen(filename.c_str(), "r");
	if (!fp)
	{
		printf("文件打开错误");
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

//***************  2  读取邻接矩阵表示的网络  **************************
void readMatrix()
{
	FILE* fp = fopen(filename.c_str(), "r");
	if (fp == NULL)
	{
		printf("文件打开失败");
		exit(1);
	}
	char C;
	int row = 0;
	int col = 0;
	int i1 = 0;//由于matrix的下标需要每次递增1，所以这里另外设置两个迭代器i1和j1
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

//***************************  5  求网络的补图  *****************************
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

//************************************  6  得到网络邻接表  ************************************
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
				adjList[i].push_back(j); //  得到的邻接表中，每一个节点的邻居节点都以标号的升序排列 
			}
		}
	}
}

//********************************  7  得到网络的度信息  *****************************
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
	maxdegree = max; // 最大度 

	int min = deg[0];
	for (i = 0; i<N; i++)
	{
		if (min>deg[i])
		{
			min = deg[i];
		}
	}
	mindegree = min; // 最小度  

	double sumdegree = 0;
	for (i = 0; i<N; i++)
	{
		sumdegree = sumdegree + deg[i];
	}
	avedegree = (double)((double)sumdegree / (double)N); // 平均度 

	edgenumber = (int)(sumdegree / 2); // 总边数 
}


//**********************************  26  将实验结果以文本格式进行存储  ************************************
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
	outfile_tm << "total_iteration_times：" << total_iteration_times <<endl<<endl; 

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

//*******************************  27  结果输出  ***************************************** 
void printResults()
{
	int i;

	averageCandidate = computeAverage();
	stdCandidate = computeStandardDeviation();
	maxCandidate = computeMaxCandidate();
	minCandidate = computeMinCandidate();
	maxCandidateNumber = computeMaxTimes();
	averageTime = computeAverageTime();

	printf("\n*********************  共运行 %d 次  *************************：\n\n", RunTime);
	cout << "Network name is " << filename << endl;
	printf("Network vertices = %d, run times = %d\n\n", N, RunTime);
	printf("1、Parameters:\n");
	printf("Itermax = %d\n", Itermax);


	printf("2、得到的解集大小分别为：\n");
	for (i = 0; i<RunTime; i++)
	{
		printf("%d ", totalCandidate[i]);
		if ((i + 1) % 10 == 0) {
			cout << endl;
		}
	}
	printf("\n\n");

	printf("3、解集大小的平均值为：%f，标准差为：%f \n\n", averageCandidate, stdCandidate);

	printf("4、解集大小最大值为：%d，最小值为：%d, 取得最大值的次数为：%d\n\n", maxCandidate, minCandidate, maxCandidateNumber);

	printf("5、平均运行时间为：%f \n\n", averageTime);
	
	printf("6、总迭代次数：%d\n",total_iteration_times);

	/*printf("1、得到的解集中的元素分别为：\n");
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

	/*printf("6、每一次运行的适应度函数值：");
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

//*******************************  20  计算RunTime次运行的解集大小平均值  *****************************
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

//******************************  21  计算RunTime次运行的解集大小的标准差  ******************************* 
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

//*****************************  22 计算RunTime次运行，得到的最大解  *************************************
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

//*****************************  23  计算RunTime次运行，得到的最小解  ***************************************
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

//***************************  24  计算RunTime次运行，最优解出现的次数  **************************************
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

//****************************  25  计算RunTime次运行的平均运行时间  *******************************************
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

// *****************  15  计算solution的大小   ************************
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


