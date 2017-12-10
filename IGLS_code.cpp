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

#define N 4000 // ����ڵ��� 
string filename = "frb100-40_4000_572774.txt";

#define b 1.5 // PDG����  
#define RunTime 5 // ͳ��RunTime�ε����н������ƽ��

#define T 1000// ����������
#define tmax 1000

//#define ml 20//���䳤��
int ml = 0; // ml = initialSolutionSize*ratio; 
double ratio = 0.90;

#define Tep_ratio 0.75
#define temperature 100
double Tep = temperature;

//double r = 0.75;
//double Tep = 100; //Tep���� Tep = Tep*0.75
//double temperature = Tep;



//string filename;
//********************************************  ��������  *******************************************
array< array<int, N>, N > matrix; // ������ڽӾ��� 
array< array<int, N>, N > matrix_temp; // ����������Ĳ�ͼ
vector< vector<int> > adjList; //�����ڽӱ��洢�ڵ���ھӽڵ�ı�� 
vector< vector<int> > adjList_temp;
vector<int> line;
int deg[N]; // ������ÿ���ڵ�Ķ�
int maxdegree; // ���� 
int mindegree; // ��С�� 
double avedegree; // ƽ���� 
int edgenumber; // �����ܱ���

array< int, N > solution; // �⼯
array< int, N > solution_temp;

int tightness[N]; // ÿһ���ڵ�Ľ��ܶ� 
array< int, N > perSolution; // ���νṹ��ʾ�Ľ⼯
int brock1Num;
int brock2Num;
int brock3Num; // perSolution[N]�У�ÿһ���еĽڵ����������perSolution���¶����� 

double payoff1[N];
double payoff0[N];

array< array<int, N>, RunTime > solu; // RunTime�����У�ÿһ�εõ��Ľ⼯
int totalCandidate[RunTime]; // RunTime�����У�ÿһ�εõ��Ľ⼯��Ԫ�ظ��������õ��Ķ�������Ԫ�ظ���
int totalIterationTimes[RunTime];
double totalTime[RunTime]; // RunTime�����У�ÿһ�ε�����ʱ��
double averageCandidate;// RunTime�����У���������С��ƽ��ֵ 
double stdCandidate;// RunTime�����У���������С�ı�׼��
int maxCandidate;// RunTime��������candidate�����ֵ 
int minCandidate;// RunTime��������candidate����Сֵ 
int maxCandidateNumber;// RunTime�������еõ����candidate�Ĵ���
double averageTime;// RunTime�����е�ƽ������ʱ��
double runningTime; // ÿһ�����е�����ʱ�� 
int run_times; // ��¼���д��� run_times = run+1 

vector<int> candidateSet;//����ÿһ�����������������S�п���ִ�в��ĸĽ��Ľڵ㼯�� 

vector<vector<int>> memory;//memory.push_back(line) �����ʱ��ÿ���ڵ��memory��Ϊ��
vector<int> vertexPDG;
int flag;//�ж��Ƿ����Ƴ�ȥ��ver�����½�������������flag=0����˵��verû�����½��������

array<int, N> formerSolution;
array<int, N> latterSolution;
int formerSize;
int latterSize;

array<int, N> bianhao;

//****************** �Ӻ��� **********************
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
	sortDegree();//���սڵ�ĶȴӴ�С���� 
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
			formerSolution = solution;//�Ŷ��;ֲ�����ǰ�Ľ� 
			formerSize = computeSolutionSize();
			
			bestSolution = solution;
			bestSize = formerSize;
			
			start_perturb = clock();
			perturbSolution();
			end_perturb = clock();
			time_perturb = time_perturb + (double)(end_perturb - start_perturb) / CLOCKS_PER_SEC;
			
			
			for( int i=0; i<N; i++ ){
				memory[i].clear();//������нڵ�ļ��� 
			}
			initialCandidateSet();
			int t = 0;//��¼whileѭ�����еĴ��� 
			int notImprove = 0;//whileѭ���У��⼯û�õ��Ľ��Ĵ��� 
			int improve = 0;//whileѭ���У��⼯�õ��Ľ��Ĵ��� 
			while (candidateSet.size() >0 ) {
				t = t + 1;
				//cout << "whileѭ��: " << "t=" << t << "," << "candidateSet.size()=" << candidateSet.size()<<endl;
							
				/*��candidateSet�����ѡ��һ���ڵ㣬������ڵ��S���Ƴ���ͬʱΪ����ڵ����ü���,
				��Ȼ�󣬰�����ڵ�������ھӽڵ��ĩβ������Щ�ڵ���в��ģ��õ���ʲ���⣬
				����ʲ�����в���ΪC�Ľڵ������������γ��µĶ�����S��
				���S�еĽڵ�vִ������һ������ֻ������������������ô������ڵ�ͱ����·����������
				ͬʱ��������������һ��C��
				����S�еĽڵ㣬��Щû�м��䡢ͬʱ���ڷǶ��������м���Ľڵ㲻���ڵĽڵ㣬 �ſ��Ա�ѡ�У�
				������candidateSet�� ��
				������£����䳤��Ϊml��������м��䶼ΪC���߶�ΪD����ô����ռ��䣬���м����Ϊ�޼���
				*/
				
				int size1 = computeSolutionSize();
				int ver = rselectFromCandidateSet();//��candidateSet�����ѡ��һ���ڵ�
				removeVer(ver);
				updateMemory();//�����Ѿ��м���Ľڵ㣬ʹ���������1��
				               //�����������������ռ��䣬��ֱ�Ӽ���candidateSet 
				solution[ver] = 0;
				obtainVertexPDGorder(ver);
				PDGimprove();
				int size2 = computeSolutionSize();

				if (size2 > bestSize) {
					bestSize = size2;
					bestSolution = solution;
				}

				if (solution[ver] == 0) {//���Ƴ��Ľڵ�û�����½�������� 
					memory[ver].push_back(0);
					if(size2==size1){
						//notImprove = notImprove+1;
						//1��1�����½���������Ľڵ�ҪЯ�����䣬�����ܽ���candidateSet 
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
			
			//cout << "whileѭ��: " << "t=" << t << "," << "candidateSet.size()=" << candidateSet.size()<<endl;
				
			//solution = bestSolution;
		}
		//**********************************************
		solution = bestSolution;
		int finalSolutionSize = computeSolutionSize();
		end = clock();
		runningTime = (double)(end - start) / CLOCKS_PER_SEC;
		runningTime = runningTime - time_perturb;
		solu[run] = solution; // ��run�����еõ��Ľ⼯ 
		totalCandidate[run] = finalSolutionSize; // ��run�����е����Ž� 
		totalTime[run] = runningTime; // ��run�����е�����ʱ�� 
		printf("run = %d , initialSolutionSize = %d , finalSolutionSize= %d \n", run_times, initialSolutionSize, finalSolutionSize);
	}
	saveResults();
	printResults();
	system("pause");
	return 0;
}

//***************  1  ��ȡ�ļ�������ıߵı�ʾ��  **************************
void read() // �õ��������ڽӾ��� �� ����¼�˱ߵľ���edgeMatrix[M][2] 
{
	int i;
	FILE* fp = fopen(filename.c_str(),"r"); 
	
	//FILE* fp = fopen("brock200_2 (12).txt","r");
	//FILE* fp = fopen("vertex6.txt", "r");
	//FILE* fp = fopen("C125.9(34).txt","r"); 
	//FILE* fp = fopen("p_hat1500-3_847244 (94).txt", "r");

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
	FILE* fp = fopen(filename.c_str(),"r"); 
	//FILE* fp = fopen("1dc.1024_24063.txt","r");
	//FILE* fp = fopen("vertex6.txt","r");
	//FILE* fp = fopen("LFR100 15 50��1��.txt", "r");

	if (!fp)
	{
		printf("�ļ��򿪴���");
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

//***************  2  ��ȡ�ڽӾ����ʾ������  **************************
void readMatrix()
{
	FILE* fp;
	fp = fopen(filename.c_str(),"r"); 
	//fp = fopen( "ws100.txt" , "r" );
	//fp = fopen("er2000(10).txt", "r");
	//fp = fopen( "s1000b.txt" , "r" );

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

//***********************  3  ������ʼ��  **********************
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

//*******************  4  ��ʼ����ÿһ�������йصı���  ************************* 
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

//*******************************  8  ���������ʼ��Ⱥ ****************************************
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

//****************************  9  PDG�������ҽڵ�Ĳ���˳�򣬴ӽڵ�0��ʼ�����β��ģ�  **************************
void PDG() // ��ͨ��PDG�����ڲ�������Ľ�solution 
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

//****************************  10  ���ݽڵ�ĶȵĴ�С�������򣬸����ڽӾ���   *****************
void sortDegree() {
	int i;
	int j;
	for (i = 0; i<N; i++) {
		bianhao[i] = i;
	}

	//���ݽڵ�Ķ�ֵ���Ӵ�С���� 
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

	//�����ڽӱ���д�ڽӾ���
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


// *****************  15  ����solution�Ĵ�С   ************************
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

//***************  16 ��ʼ��candidateSet ***************** 
void initialCandidateSet(){
	int i;
	//���ж������еĽڵ㶼���Գ�ʼ��΢candidateSet�еĽڵ� 
	for (i = 0; i<N; i++){
		if (solution[i] == 1){
			candidateSet.push_back(i);
		}
	}
}

//*************** 17 ��candidate�����ѡ��һ���ڵ㣬�����Ƴ�solution ************** 
int rselectFromCandidateSet(){
	int size = candidateSet.size();
	int r = (int)(size*rand()) / (RAND_MAX + 1);
	int selectedVertex = candidateSet[r];
	return selectedVertex;
}

//************* 18 ����candidateSet ********************** 
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

//*********************  19 һ�־ֲ����������󣬸��½⼯  ****************** 
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

//****************  20 ���½ڵ��memory ****************** 
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
				candidateSet.push_back(i); //ֱ�Ӹ���candidateSet 
				//vertex.push_back(i);
			}
		}
	}
	
}

//***************  21 ver�������ھӽڵ�ĺ��� ���в��� ******************* 
void obtainVertexPDGorder(int ver) {
	vertexPDG.clear();
	for (auto i = 0; i<adjList[ver].size(); i++) {
		int verNeighbor = adjList[ver][i];
		vertexPDG.push_back(verNeighbor);
	}
	vertexPDG.push_back(ver);

}

//***********************  22 ver�ھ�+ver ���Ĵﵽ��ʲ���� ************* 
void PDGimprove() {
	int i;
	//vertexPDG�еĽڵ���в���
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

//*******************************  20  ����RunTime�����еĽ⼯��Сƽ��ֵ  *****************************
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

//******************************  21  ����RunTime�����еĽ⼯��С�ı�׼��  ******************************* 
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

//*****************************  22 ����RunTime�����У��õ�������  *************************************
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

//*****************************  23  ����RunTime�����У��õ�����С��  ***************************************
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

//***************************  24  ����RunTime�����У����Ž���ֵĴ���  **************************************
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

//****************************  25  ����RunTime�����е�ƽ������ʱ��  *******************************************
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

//*******************************  27  ������  ***************************************** 
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

	printf("\n*********************  ������ %d ��  *************************��\n\n", RunTime);

	printf("Network vertices = %d, run times = %d\n\n", N, RunTime);
	printf("1��Parameters:\n");
	printf("Tmax = %d, tmax = %d,\nml ratio = %f,\ninitial Tep = %d, r = %f,\nPDG b = %f\n\n", T, tmax, ratio, temperature, Tep_ratio, b);
	
	
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

	printf("2���õ��Ľ⼯��С�ֱ�Ϊ��\n ");
	for (i = 0; i<RunTime; i++)
	{
		printf("%d ", totalCandidate[i]);
		if((i+1)%10==0){
			cout<<endl;
		}
	}
	printf("\n\n");

	printf("3���⼯��С��ƽ��ֵΪ��%f����׼��Ϊ��%f \n\n", averageCandidate, stdCandidate);

	printf("4���⼯��С���ֵΪ��%d����СֵΪ��%d, ȡ�����ֵ�Ĵ���Ϊ��%d\n\n", maxCandidate, minCandidate, maxCandidateNumber);

	printf("5��ƽ������ʱ��Ϊ��%f \n\n", averageTime);


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






