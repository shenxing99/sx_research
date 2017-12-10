/*2017 05 06 update IGLS*/

#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include<iostream>
#include<math.h>
#include<time.h>
#include<vector>
#include<array> 
#include<string>
using namespace std;

#define N 4000 // 网络节点数 
string filename = "frb100-40_4000_572774.txt";

#define b 1.5 // PDG参数  
#define RunTime 5 // 统计RunTime次的运行结果，求平均

#define T 1000// 最大迭代次数
#define tmax 1000

//#define ml 20//记忆长度
int ml = 0; // ml = initialSolutionSize*ratio; 
double ratio = 0.90;

#define Tep_ratio 0.75
#define temperature 100
double Tep = temperature;

//double r = 0.75;
//double Tep = 100; //Tep更新 Tep = Tep*0.75
//double temperature = Tep;



//string filename;
//********************************************  变量定义  *******************************************
array< array<int, N>, N > matrix; // 网络的邻接矩阵 
array< array<int, N>, N > matrix_temp; // 用于求网络的补图
vector< vector<int> > adjList; //网络邻接表，存储节点的邻居节点的标号 
vector< vector<int> > adjList_temp;
vector<int> line;
int deg[N]; // 网络中每个节点的度
int maxdegree; // 最大度 
int mindegree; // 最小度 
double avedegree; // 平均度 
int edgenumber; // 网络总边数

array< int, N > solution; // 解集
array< int, N > solution_temp;

int tightness[N]; // 每一个节点的紧密度 
array< int, N > perSolution; // 三段结构表示的解集
int brock1Num;
int brock2Num;
int brock3Num; // perSolution[N]中，每一块中的节点个数，随着perSolution更新而更新 

double payoff1[N];
double payoff0[N];

array< array<int, N>, RunTime > solu; // RunTime次运行，每一次得到的解集
int totalCandidate[RunTime]; // RunTime次运行，每一次得到的解集的元素个数，即得到的独立集的元素个数
int totalIterationTimes[RunTime];
double totalTime[RunTime]; // RunTime次运行，每一次的运行时间
double averageCandidate;// RunTime次运行，独立集大小的平均值 
double stdCandidate;// RunTime次运行，独立集大小的标准差
int maxCandidate;// RunTime次运行中candidate的最大值 
int minCandidate;// RunTime次运行中candidate的最小值 
int maxCandidateNumber;// RunTime次运行中得到最大candidate的次数
double averageTime;// RunTime次运行的平均运行时间
double runningTime; // 每一次运行的运行时间 
int run_times; // 记录运行次数 run_times = run+1 

vector<int> candidateSet;//对于每一个极大独立集，极大集S中可以执行博弈改进的节点集合 

vector<vector<int>> memory;//memory.push_back(line) 最初的时候，每个节点的memory都为空
vector<int> vertexPDG;
int flag;//判断是否是移除去的ver又重新进入独立集，如果flag=0；则说明ver没有重新进入独立集

array<int, N> formerSolution;
array<int, N> latterSolution;
int formerSize;
int latterSize;

array<int, N> bianhao;

//****************** 子函数 **********************
void read();
void read_without_e();
void readMatrix();
void initialVariable();
void initialVariable_RunTime();
void complement_graph();
void obtainAdjList();
void obtainVerDegree();

void perturbSolution();
void initialPop();
void PDG();
void sortDegree();
int computeSolutionSize();
void initialCandidateSet();
int rselectFromCandidateSet();
void updateCandidateSet();
void updateSolution();
void updateMemory();
void obtainVertexPDGorder(int ver);
void PDGimprove();

double computeAverage();
double computeStandardDeviation();
int computeMaxCandidate();
int computeMinCandidate();
int computeMaxTimes();
double computeAverageTime();

void saveResults();
void printResults();

void removeVer(int vertex);

array<int, N> bestSolution;
int bestSize;

void updateSolution_Tep();


//********************************************  main  *******************************************
int main(int argc, char* argv[]) {
	clock_t start, end;
	clock_t start_perturb, end_perturb;
	srand((unsigned)time(NULL));
	initialVariable();
	//filename = "gen200_p0.9_44 (44).txt";

	read();
	//complement_graph();
	//read_without_e(); 
	//readMatrix(); 

	obtainAdjList();
	obtainVerDegree();
	sortDegree();//按照节点的度从大到小排序 
	initialVariable_RunTime();

	for (int run = 0; run<RunTime; run++) {
		run_times = run + 1;
		srand((unsigned)time(NULL));
		candidateSet.clear();
		initialPop();
		PDG();
		int initialSolutionSize = computeSolutionSize();
		ml = (int)(initialSolutionSize*ratio);
		start = clock();
		double time_perturb = 0;
		for (int iter = 0; iter<T; iter++) {
			formerSolution = solution;//扰动和局部搜索前的解 
			formerSize = computeSolutionSize();
			
			bestSolution = solution;
			bestSize = formerSize;
			
			start_perturb = clock();
			perturbSolution();
			end_perturb = clock();
			time_perturb = time_perturb + (double)(end_perturb - start_perturb) / CLOCKS_PER_SEC;
			
			
			for( int i=0; i<N; i++ ){
				memory[i].clear();//清空所有节点的记忆 
			}
			initialCandidateSet();
			int t = 0;//记录while循环进行的次数 
			int notImprove = 0;//while循环中，解集没得到改进的次数 
			int improve = 0;//while循环中，解集得到改进的次数 
			while (candidateSet.size() >0 ) {
				t = t + 1;
				//cout << "while循环: " << "t=" << t << "," << "candidateSet.size()=" << candidateSet.size()<<endl;
							
				/*从candidateSet中随机选择一个节点，将这个节点从S中移除，同时为这个节点设置记忆,
				，然后，把这个节点放在其邻居节点的末尾，让这些节点进行博弈，得到纳什均衡，
				将纳什均衡中策略为C的节点放入独立集，形成新的独立集S。
				如果S中的节点v执行完这一操作后，只能自身进入独立集，那么，这个节点就被重新放入独立集，
				同时，记忆里又增加一个C。
				对于S中的节点，那些没有记忆、同时与在非独立集中有记忆的节点不相邻的节点， 才可以被选中，
				被放入candidateSet中 。
				记忆更新：记忆长度为ml，如果所有记忆都为C或者都为D，那么就清空记忆，由有记忆变为无记忆
				*/
				
				int size1 = computeSolutionSize();
				int ver = rselectFromCandidateSet();//从candidateSet中随机选择一个节点
				removeVer(ver);
				updateMemory();//对于已经有记忆的节点，使其记忆增加1；
				               //如果记忆已满，则清空记忆，且直接加入candidateSet 
				solution[ver] = 0;
				obtainVertexPDGorder(ver);
				PDGimprove();
				int size2 = computeSolutionSize();

				if (size2 > bestSize) {
					bestSize = size2;
					bestSolution = solution;
				}

				if (solution[ver] == 0) {//被移除的节点没有重新进入独立集 
					memory[ver].push_back(0);
					if(size2==size1){
						//notImprove = notImprove+1;
						//1换1，则新进入独立集的节点要携带记忆，即不能进入candidateSet 
						for (auto k = 0; k < adjList[ver].size(); k++) {
							int n = adjList[ver][k];
							if (solution[n] == 1) {
								memory[n].push_back(1);
							}
						} 
					}
					else if(size2>size1){
						for (auto k = 0; k < adjList[ver].size(); k++) {
							int n = adjList[ver][k];
							if (solution[n] == 1) {
								candidateSet.push_back(n);
								//break;
							}
						} 
					}
				}
				else if (solution[ver] == 1) {
					notImprove = notImprove+1;
					memory[ver].push_back(1);
				}
				/*if(notImprove==500){
					candidateSet.clear();
				}*/
				if(t == tmax) {
					candidateSet.clear();
				}
			}
			
			latterSolution = solution;
			latterSize = computeSolutionSize();
			
			updateSolution_Tep();
			
			//solution = bestSolution;
			
			//cout << "while循环: " << "t=" << t << "," << "candidateSet.size()=" << candidateSet.size()<<endl;
				
			//solution = bestSolution;
		}
		//**********************************************
		solution = bestSolution;
		int finalSolutionSize = computeSolutionSize();
		end = clock();
		runningTime = (double)(end - start) / CLOCKS_PER_SEC;
		runningTime = runningTime - time_perturb;
		solu[run] = solution; // 第run次运行得到的解集 
		totalCandidate[run] = finalSolutionSize; // 第run次运行的最优解 
		totalTime[run] = runningTime; // 第run次运行的运行时间 
		printf("run = %d , initialSolutionSize = %d , finalSolutionSize= %d \n", run_times, initialSolutionSize, finalSolutionSize);
	}
	saveResults();
	printResults();
	system("pause");
	return 0;
}

//***************  1  读取文件（网络的边的表示）  **************************
void read() // 得到了网络邻接矩阵 ， 并记录了边的矩阵edgeMatrix[M][2] 
{
	int i;
	FILE* fp = fopen(filename.c_str(),"r"); 
	
	//FILE* fp = fopen("brock200_2 (12).txt","r");
	//FILE* fp = fopen("vertex6.txt", "r");
	//FILE* fp = fopen("C125.9(34).txt","r"); 
	//FILE* fp = fopen("p_hat1500-3_847244 (94).txt", "r");

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
	FILE* fp = fopen(filename.c_str(),"r"); 
	//FILE* fp = fopen("1dc.1024_24063.txt","r");
	//FILE* fp = fopen("vertex6.txt","r");
	//FILE* fp = fopen("LFR100 15 50（1）.txt", "r");

	if (!fp)
	{
		printf("文件打开错误");
	}

	//char e;
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
	FILE* fp;
	fp = fopen(filename.c_str(),"r"); 
	//fp = fopen( "ws100.txt" , "r" );
	//fp = fopen("er2000(10).txt", "r");
	//fp = fopen( "s1000b.txt" , "r" );

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

//***********************  3  变量初始化  **********************
void initialVariable()
{
	int i;
	int j;
	run_times = 0;
	for (i = 0; i<N; i++)
	{
		for (j = 0; j<N; j++)
		{
			matrix[i][j] = 0;
			matrix_temp[i][j] = 0;
		}
	}
	for (i = 0; i<N; i++)
	{
		adjList.push_back(line);
		adjList_temp.push_back(line);
		memory.push_back(line);
	}
	for (i = 0; i<N; i++)
	{
		deg[i] = 0;
		solution[i] = 0;
	}
	maxdegree = 0;
	mindegree = 0;
	avedegree = 0;
	edgenumber = 0;
}

//*******************  4  初始化与每一次运行有关的变量  ************************* 
void initialVariable_RunTime()
{
	int i;
	int j;
	for (i = 0; i<RunTime; i++)
	{
		totalCandidate[i] = 0;
		totalTime[i] = 0;
	}
	for (i = 0; i<RunTime; i++)
	{
		for (j = 0; j<N; j++)
		{
			solu[i][j] = 0;
		}
	}
	averageCandidate = 0;
	stdCandidate = 0;
	maxCandidate = 0;
	minCandidate = 0;
	maxCandidateNumber = 0;
	averageTime = 0;
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

//*******************************  8  随机产生初始种群 ****************************************
void initialPop()
{
	int i;

	int random[N];
	for (i = 0; i<N; i++){
		random[i] = 0;
	}
	for (i = 0; i<N; i++){
		random[i] = (int)(100 * rand()) / (RAND_MAX + 1);
	}

	for (i = 0; i<N; i++){
		if (random[i] >= 50){
			solution[i] = 1;
		}
		else{
			solution[i] = 0;
		}
	}
}

//****************************  9  PDG（不打乱节点的博弈顺序，从节点0开始，依次博弈）  **************************
void PDG() // 普通的PDG，用于产生最初的解solution 
{
	int i;
	//int j;
	int t = 1;
	int strategy_temp[N];
	for (i = 0; i<N; i++){
		strategy_temp[i] = 0;
	}
	for (i = 0; i<N; i++){
		payoff1[i] = 0;
		payoff0[i] = 0;
	}

	while (t){
		t = 0;
		for (i = 0; i<N; i++){
			payoff1[i] = 0;
			payoff0[i] = 0;
		}
		for (i = 0; i<N; i++){
			strategy_temp[i] = solution[i];
			for (auto j = 0; j<adjList[i].size(); j++){
				int neibor = adjList[i][j];
				if (solution[neibor] == 1){
					payoff1[i] = payoff1[i] + 1;
					payoff0[i] = payoff0[i] + b;
				}
				else if (solution[neibor] == 0){
					payoff1[i] = payoff1[i] + 0;
					payoff0[i] = payoff0[i] + 0;
				}
			}
			if (payoff1[i] >= payoff0[i]){
				solution[i] = 1;
			}
			else if (payoff1[i]<payoff0[i]){
				solution[i] = 0;
			}

			if (t == 0){
				if (strategy_temp[i] != solution[i]){
					t = 1;
				}
				else if (strategy_temp[i] == solution[i]){
					t = 0;
				}
			}
		}

	}
}

//****************************  10  根据节点的度的大小进行排序，更新邻接矩阵   *****************
void sortDegree() {
	int i;
	int j;
	for (i = 0; i<N; i++) {
		bianhao[i] = i;
	}

	//根据节点的度值，从大到小排序 
	for (i = 0; i <= N - 2; i++) {
		for (j = i + 1; j <= N - 1; j++) {
			if (deg[i]<deg[j]) {
				int temp = deg[j];
				deg[j] = deg[i];
				deg[i] = temp;

				int temp1 = bianhao[j];
				bianhao[j] = bianhao[i];
				bianhao[i] = temp1;
			}
		}
	}

	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			adjList_temp[i].clear();
		}
	}
	for (i = 0; i<N; i++){
		int oldBianhao = bianhao[i];
		for (auto j = 0; j<adjList[oldBianhao].size(); j++){
			for (int x = 0; x<N; x++) {
				if (bianhao[x] == adjList[oldBianhao][j]){
					adjList_temp[i].push_back(x);
				}
			}
		}
	}

	//根据邻接表，重写邻接矩阵
	for (i = 0; i<N; i++){
		adjList[i].clear();
	}
	for (i = 0; i<N; i++){
		for (j = 0; j<adjList_temp[i].size(); j++){
			adjList[i].push_back(adjList_temp[i][j]);
		}
	}

	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			matrix[i][j] = 0;
		}
	}
	for (i = 0; i<N; i++){
		for (j = 0; j<adjList[i].size(); j++){
			matrix[i][adjList[i][j]] = 1;
		}
	}
}


// *****************  15  计算solution的大小   ************************
int computeSolutionSize()
{
	int i;
	int sizeSum = 0;
	for (i = 0; i<N; i++){
		if (solution[i] == 1){
			sizeSum = sizeSum + 1;
		}
	}
	return sizeSum;
}

//***************  16 初始化candidateSet ***************** 
void initialCandidateSet(){
	int i;
	//所有独立集中的节点都可以初始化微candidateSet中的节点 
	for (i = 0; i<N; i++){
		if (solution[i] == 1){
			candidateSet.push_back(i);
		}
	}
}

//*************** 17 从candidate中随机选择一个节点，将其移出solution ************** 
int rselectFromCandidateSet(){
	int size = candidateSet.size();
	int r = (int)(size*rand()) / (RAND_MAX + 1);
	int selectedVertex = candidateSet[r];
	return selectedVertex;
}

//************* 18 更新candidateSet ********************** 
void updateCandidateSet(){
	int i;
	int j;
	int f[N];

	candidateSet.clear();
	for (i = 0; i<N; i++){
		f[i] = 1;
		if (solution[i] == 1){
			if((memory[i].size()==ml)||(memory[i].size()==0)){
				candidateSet.push_back(ml);
			}
			/*for (j = 0; j<adjList[i].size(); j++){
				int neighbor = adjList[i][j];
				if ((memory[neighbor].size()>0) && (memory[neighbor].size()<ml)){
					f[i] = 0;
				}
			}
			if (f[i] == 1){
				if ((memory[i].size() == 0)){
					candidateSet.push_back(i);
				}
			}*/
		}
	}
}

//*********************  19 一轮局部搜索结束后，更新解集  ****************** 
void updateSolution() {
	if (latterSize>formerSize) {
		solution = latterSolution;
	}
	else {
		int r = (int)(100 * rand()) / (RAND_MAX + 1);
		if (r<5) {
			solution = latterSolution;
		}
		else {
			solution = formerSolution;
		}
	}
}

//****************  20 更新节点的memory ****************** 
void updateMemory() {
	int i;
	vector<int> vertex;
	vertex.clear();
	for (i = 0; i<N; i++) {
		if ((memory[i].size()>0) && (memory[i].size()<ml)) {
			if (solution[i] == 0) {
				memory[i].push_back(0);
			}
			else if (solution[i] == 1) {
				memory[i].push_back(1);
			}
		}
		if (memory[i].size() == ml) {
			memory[i].clear();
			if(solution[i]==1){
				candidateSet.push_back(i); //直接更新candidateSet 
				//vertex.push_back(i);
			}
		}
	}
	
}

//***************  21 ver放在其邻居节点的后面 进行博弈 ******************* 
void obtainVertexPDGorder(int ver) {
	vertexPDG.clear();
	for (auto i = 0; i<adjList[ver].size(); i++) {
		int verNeighbor = adjList[ver][i];
		vertexPDG.push_back(verNeighbor);
	}
	vertexPDG.push_back(ver);

}

//***********************  22 ver邻居+ver 博弈达到纳什均衡 ************* 
void PDGimprove() {
	int i;
	//vertexPDG中的节点进行博弈
	int t = 1;
	vector<int> strategy_temp;
	for (auto k = 0; k < vertexPDG.size(); k++){
		strategy_temp.push_back(0);
	}
	for (i = 0; i<vertexPDG.size(); i++){
		payoff1[i] = 0;
		payoff0[i] = 0;
	}
	while (t){
		t = 0;
		for (i = 0; i<vertexPDG.size(); i++) {
			payoff1[i] = 0;
			payoff0[i] = 0;
		}
		for (i = 0; i<vertexPDG.size(); i++) {
			int v = vertexPDG[i];
			strategy_temp[i] = solution[v];
			for (int j = 0; j<adjList[v].size(); j++) {
				int neighbor = adjList[v][j];
				if (solution[neighbor] == 1) {
					payoff1[i] = payoff1[i] + 1;
					payoff0[i] = payoff0[i] + b;
				}
				else if (solution[neighbor] == 0) {
					payoff1[i] = payoff1[i] + 0;
					payoff0[i] = payoff0[i] + 0;
				}
			}
			if (payoff1[i] >= payoff0[i]) {
				solution[v] = 1;
			}
			else if (payoff1[i]<payoff0[i]) {
				solution[v] = 0;
			}
			if (t == 0) {
				if (strategy_temp[i] != solution[v]) {
					t = 1;
				}
				else if (strategy_temp[i] == solution[v]) {
					t = 0;
				}
			}
		}
	}
}

//*******************************  20  计算RunTime次运行的解集大小平均值  *****************************
double computeAverage(){
	int i;
	int j;
	double average = 0;
	int sum = 0;
	for (i = 0; i<RunTime; i++)	{
		sum = sum + totalCandidate[i];
	}
	average = (double)((double)sum / (double)RunTime);

	return average;
}

//******************************  21  计算RunTime次运行的解集大小的标准差  ******************************* 
double computeStandardDeviation(){
	int i;
	double average = computeAverage();
	double standardDeviation = 0;
	double sum = 0;
	for (i = 0; i<RunTime; i++)	{
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
	for (i = 0; i<RunTime; i++)	{
		if (max<totalCandidate[i]){
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
	for (i = 0; i<RunTime; i++)	{
		if (min>totalCandidate[i]){
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
	for (i = 0; i<RunTime; i++){
		if (totalCandidate[i] == max){
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
	for (i = 0; i<RunTime; i++)	{
		sum = sum + totalTime[i];
	}
	averageTime = (double)((double)sum / (double)RunTime);

	return averageTime;
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

	string name = "result_"+filename;
	ofstream outfile_tm(name,std::ios::trunc);
	
	//ofstream outfile_tm("newIGLS_result.txt");
	if (!outfile_tm)
	{
		cerr << "open error!" << endl;
		exit(1);
	}
	
	outfile_tm << "The results: " << endl << endl;
	outfile_tm <<"1"<<endl; 
	outfile_tm <<"Network information:"<< endl;	
	outfile_tm << "network name    : "<< filename.c_str()<<endl;
	outfile_tm << "network vertices: " << N <<endl;
	outfile_tm << "total run times : " << RunTime << endl;
	outfile_tm << "edge number     : " << edgenumber <<endl;
	outfile_tm << "maximum degree  : " << maxdegree <<endl;
	outfile_tm << "minimum degree  : " << mindegree <<endl;
	outfile_tm << "average degree  : " << avedegree <<endl<<endl;
	
	outfile_tm <<"2"<<endl; 
	outfile_tm <<"The Parameters:"<<endl;
	outfile_tm <<"PDG b = "<<b<<endl;
	outfile_tm <<"Tmax  = "<<T<<endl;
	outfile_tm <<"tmax  = "<<tmax<<endl;
	outfile_tm <<"ml ratio    = "<<ratio<<endl;
	outfile_tm <<"initial Tep = "<<temperature<<endl;
	outfile_tm <<"Tep ratio   = "<<Tep_ratio<<endl<<endl;
	
	outfile_tm <<"3"<<endl; 
	outfile_tm << "The MIS:"<<endl;
	outfile_tm << "averageCandidate  : " << averageCandidate <<endl;
	outfile_tm << "stdCandidate      : " << stdCandidate <<endl;
	outfile_tm << "maxCandidate      : " << maxCandidate <<endl;
	outfile_tm << "minCandidate      : " << minCandidate <<endl;
	outfile_tm << "maxCandidateNumber: " << maxCandidateNumber <<endl;
	outfile_tm << "averageTime       : " << averageTime << endl << endl;
	
	outfile_tm <<"4"<<endl; 
	outfile_tm<<"The detailed results:"<<endl;
	outfile_tm << "total candidate of RunTime: " <<endl;
	for (i = 0; i<RunTime; i++)
	{
		outfile_tm << totalCandidate[i] << " ";
		if((i+1)%10==0){
			outfile_tm<<endl;
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
	int j;

	averageCandidate = computeAverage();
	stdCandidate = computeStandardDeviation();
	maxCandidate = computeMaxCandidate();
	minCandidate = computeMinCandidate();
	maxCandidateNumber = computeMaxTimes();
	averageTime = computeAverageTime();

	printf("\n*********************  共运行 %d 次  *************************：\n\n", RunTime);

	printf("Network vertices = %d, run times = %d\n\n", N, RunTime);
	printf("1、Parameters:\n");
	printf("Tmax = %d, tmax = %d,\nml ratio = %f,\ninitial Tep = %d, r = %f,\nPDG b = %f\n\n", T, tmax, ratio, temperature, Tep_ratio, b);
	
	
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

	printf("2、得到的解集大小分别为：\n ");
	for (i = 0; i<RunTime; i++)
	{
		printf("%d ", totalCandidate[i]);
		if((i+1)%10==0){
			cout<<endl;
		}
	}
	printf("\n\n");

	printf("3、解集大小的平均值为：%f，标准差为：%f \n\n", averageCandidate, stdCandidate);

	printf("4、解集大小最大值为：%d，最小值为：%d, 取得最大值的次数为：%d\n\n", maxCandidate, minCandidate, maxCandidateNumber);

	printf("5、平均运行时间为：%f \n\n", averageTime);


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

void perturbSolution() {
	int i;
	int j;

	array<int, N> solutionTemp;
	solutionTemp = solution;

	int s = computeSolutionSize();
	s = 2 * s;
	int r = (int)(s*rand()) / (RAND_MAX + 1);
	int k;
	if (r != 1) {
		k = 1;
		vector<int> noneSolution;
		for (i = 0; i<N; i++) {
			if ((solution[i] == 0) && memory[i].size() == 0) {
				noneSolution.push_back(i);
			}
		}
		int r1 = (int)(noneSolution.size()*rand()) / (RAND_MAX + 1);
		int vertex = noneSolution[r1];
		solution[vertex] = 1;
		for (j = 0; j<adjList[vertex].size(); j++) {
			solution[adjList[vertex][j]] = 0;
		}
	}
	else {
		k = 2;
	}
}

void removeVer(int vertex){
	for( vector<int>::iterator iter=candidateSet.begin(); iter<candidateSet.end(); iter++ ){
		if(vertex==*iter){
			candidateSet.erase(iter);
			break;
		}
	}
}

void updateSolution_Tep(){
	int deltaSize = latterSize-formerSize;
	if(deltaSize>0){
		solution = latterSolution;
	}
	else{
		double p = (double)(rand())/(RAND_MAX + 1);
		if( exp((double)(deltaSize/Tep))>p ){
			solution = latterSolution;
		}
		else{
			solution = formerSolution;
		}
	}
	//Tep = (double) ( 0.75*(double)Tep ) ;
	Tep = 0.75*Tep;
}






