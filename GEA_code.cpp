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
#include<algorithm>

using namespace std;

#define N 450 //网络节点总数 
#define M 200 //种群规模
#define r 0.002 //博弈收益比率 
#define RunTime 10  //独立运行的次数 
#define Gen 100 //种群迭代的次数 
#define length 3//交叉的长度 
#define mlsTime 500
#define Proba 1//变异率 
//string filename = "vertex6.txt";
//string filename = "C125.9(34).txt";
//string filename = "brock200_4 (17).txt";
string filename = "frb450.txt";

array< array<int, N>, N > matrix;//邻接矩阵 
array< array<int, N>, N > matrix_temp;
vector< vector<int> > adjList;
vector<int> line;
int deg[N];
int maxdegree;
int mindegree;
double degSum;
double avedegree;
int edgenumber;

array< array<int, N>, M> popu;//M行N列 
array< array<int, N>, M> newPopu;
array< array<int, N>, 4> p1;//4行N列 

array<int, M> popuSize;//每一个个体表示的覆盖集的规模 
array<int, N> finalIndividual;//最终的个体 

vector<int> s1;//解集中，与非解集集合中的节点只有一个连边的节点的集合 
vector<int> s2;//解集中，与非解集集合中的节点有多条连边的节点的集合 
vector<int> notResSet;//非解集集合 
vector<int> nodeNei;//选出进行局部搜索的节点，重新进行博弈的邻居节点

//double r; 

double start, end, runT;

array< array<int, N>, RunTime > solu;
int totalCandidate[RunTime];
double totalTime[RunTime];

double averageCandidate;
double stdCandidate;
int maxCandidate;
int minCandidate;
int maxCandidateNumber;
double averageTime;

//********************* 子函数  *************************** 

void globalInitial();
void read(); 
void read_without_e();
void readMatrix();
void complement_graph();
void getAdjList_degree();
void initialPopulation();
void snowGame(int num);
void popuSnowGame();
void p1_SnowGame(int num);
void partitionSet(int n);
int findNode(int n);
void neiGame(int n);
void gameAgain(int n, int node);
void localSearch(int n);
void cross();
void mutation();
int computeSize(array<int, N> vec);
int computeFitness(array<int, N> vec);
void choose(int t);
void newPopu_to_popu();

void individualSwap();
void mutation(int n);

double computeAverage();
double computeStandardDeviation();
int computeMaxCandidate();
int computeMinCandidate();
int computeMaxTimes();
double computeAverageTime();

void saveResults();
void printResults();




int main(){
	cout<<filename<<":"<<endl;
	cout<<"M = "<<M<<endl;
	cout<<"N = "<<N<<endl;
	
	//r = 1/maxdegree-1/N;
	
	int i;
	int j;
	clock_t start, end;
	clock_t start1, end1;
	globalInitial();
	
	//读取网络信息 
	//read();
	//complement_graph();
	//read_without_e();
	readMatrix();
	/*cout<<"邻接矩阵为："<<endl; 
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			cout<<matrix[i][j]<<" ";
		}
		cout<<endl;
	} 
	cout<<endl;*/

	getAdjList_degree();//得到网络中每个节点的邻居、度值信息 
	/*cout<<"maxdegree = "<<maxdegree<<endl;
	cout<<"mindegree = "<<mindegree<<endl;
	cout<<"degSum = "<<degSum<<endl;
	cout<<"avedegree = "<<avedegree<<endl;
	cout<<"edgenumber = "<<edgenumber<<endl;*/
	
	
	int run;
	srand( (unsigned)time(NULL) );
	//initialPopulation();
	for(run=0; run<RunTime; run++){//共独立运行RunTime次 
		start = clock();
		initialPopulation();//初始化种群，是M*N的矩阵 
		
		/*cout<<"初始种群为："<<endl; 
		for(int i=0; i<M; i++){
			for(int j=0; j<N; j++){
				cout<<popu[i][j]<<" ";
			}
			cout<<endl;
		} 
		cout<<endl;*/
		
		popuSnowGame();//种群中每一个个体进行博弈达到纳什均衡 
		
		/*cout<<"第一次博弈后，种群为："<<endl; 
		for(int i=0; i<M; i++){
			for(int j=0; j<N; j++){
				cout<<popu[i][j]<<" ";
			}
			cout<<endl;
		} 
		cout<<endl;*/
		
		int g;
		for(g=0; g<Gen; g++){//整个种群更新Gen代 
			if(g>=1){
				newPopu_to_popu();
			}
			for(int i=0; i<M/4; i++){
				individualSwap();//随机打乱种群 
			} 
			/*每次选择两个个体，这两个个体进行局部搜索、交叉、变异，
			产生两个新的个体，然后在四个个体中选择最优异的两个个体，进入新种群*/ 
			for(int t=0; t<M; t = t+2){ 
				/*p1[0] = popu[t];
				p1[1] = popu[t+M/2];
				p1[2] = popu[t];//交叉、变异都是针对p1[2]、p1[3]进行的 
				p1[3] = popu[t+M/2];*/
				p1[0] = popu[t];
				p1[1] = popu[t+1];
				p1[2] = popu[t];//交叉、变异都是针对p1[2]、p1[3]进行的 
				p1[3] = popu[t+1];
				cross();
				mutation(2);
				mutation(3);
				p1_SnowGame(2);
				p1_SnowGame(3);
				for(int i=0; i<=3; i++){
					localSearch(i);
				}
				/*localSearch(0);
				localSearch(1);
				localSearch(2);
				localSearch(3);*/
				choose(t);
			}
		}
		for(int i=0; i<M; i++){
			popuSize[i] = computeSize(popu[i]);
		}
		int minSize = popuSize[0];
		finalIndividual = popu[0];
		for(int i=0; i<M; i++){
			if(popuSize[i]<minSize){
				minSize = popuSize[i];
				finalIndividual = popu[i];
			}
		}
		
		solu[run] = finalIndividual;
		totalCandidate[run] = minSize;
		
		cout<<"run = "<<run<<","<<"solution size = "<<minSize<<endl;
		
		end = clock();
		runT = (long double)(end-start)/CLOCKS_PER_SEC;
		totalTime[run] = runT;
	}
	saveResults();
	printResults();
	system("pause");
	return 0;
}


//******************************  子函数  ********************************** 
//***************  1  变量初始化  **************************
void globalInitial(){
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			matrix[i][j] = 0;
			matrix_temp[i][j] = 0; 
		}
		deg[i] = 0;
	} 
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			popu[i][j] = 0;
		}
	}
	for(int i=0; i<N; i++){
		adjList.push_back(line);
	}
	maxdegree = 0;
	mindegree = 0;
	degSum = 0;
	avedegree = 0;
	edgenumber = 0;
} 

//***************  2 读取文件（网络的边的表示）  **************************
void read() // 得到了网络邻接矩阵 ， 并记录了边的矩阵edgeMatrix[M][2] 
{
	int i;
	FILE* fp = fopen(filename.c_str(), "r");
	if (!fp){
		printf("文件打开错误");
	}

	char e;
	int startnode;
	int endnode;

	i = 0;
	while (!feof(fp)){
		fscanf(fp, "%c %d %d \n", &e, &startnode, &endnode);
		matrix[startnode - 1][endnode - 1] = 1;
		matrix[endnode - 1][startnode - 1] = 1;
		i = i + 1;
	}

	fclose(fp);
}

//***************  3  读取文件（网络的边的表示）  **************************
void read_without_e() // 得到了网络邻接矩阵 ， 并记录了边的矩阵edgeMatrix[M][2] 
{
	int i;
	FILE* fp = fopen(filename.c_str(), "r");
	if (!fp){
		printf("文件打开错误");
	}

	int startnode;
	int endnode;

	i = 0;
	while (!feof(fp)){
		fscanf(fp, "%d %d \n", &startnode, &endnode);
		matrix[startnode - 1][endnode - 1] = 1;
		matrix[endnode - 1][startnode - 1] = 1;
		i = i + 1;
	}

	fclose(fp);
}

//***************  4 读取邻接矩阵表示的网络  **************************
void readMatrix()
{
	FILE* fp = fopen(filename.c_str(), "r");
	if (fp == NULL){
		printf("文件打开失败");
		exit(1);
	}
	char C;
	int row = 0;
	int col = 0;
	int i1 = 0;//由于matrix的下标需要每次递增1，所以这里另外设置两个迭代器i1和j1
	int j1 = 0;
	while ((C = fgetc(fp)) != EOF){
		if (C == '0'){
			matrix[i1][j1] = 0;
			col++;
			j1++;
		}
		else if (C == '1'){
			matrix[i1][j1] = 1;
			col++;
			j1++;
		}
		else if (C == ' '){
			col++;
		}
		else if (C == '\n'){
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
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			matrix_temp[i][j] = matrix[i][j];
		}
	}
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			if (matrix_temp[i][j] == 1){
				matrix[i][j] = 0;
			}
			else if (matrix_temp[i][j] == 0){
				matrix[i][j] = 1;
			}
		}
	}
	for (i = 0; i<N; i++){
		matrix[i][i] = 0;
	}
}

//***************************   6 得到邻居节点和节点度  *****************************
void getAdjList_degree(){
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			if(matrix[i][j]==1){
				adjList[i].push_back(j);
			}
		}
	}
	for(int i=0; i<N; i++){
		deg[i] = adjList[i].size();
	}
	maxdegree = deg[0];
	mindegree = deg[0];
	degSum = 0;
	for(int i=0; i<N; i++){
		if(maxdegree<deg[i]){
			maxdegree = deg[i];
		}
		if(mindegree>deg[i]){
			mindegree = deg[i];
		}
		degSum = degSum+deg[i];
	}
	avedegree = degSum / (double)N;
	edgenumber = (int)degSum / 2; 
}

//*******************  7 按照节点的度，初始化种群  *********************************
void initialPopulation(){
	/*for(int i=0; i<N; i++){
		cout<<deg[i]<<" ";
	}
	cout<<endl;*/
	for(int m=0; m<M; m++){
		for(int i=0; i<N; i++){
			int littleDegNum = 0;
			for(int j=0; j<N; j++){
				if(deg[j]<deg[i]){
					littleDegNum+=deg[j];
				}
			}
			//littleDegNum-=deg[i];
			int randNum = (int)(rand()%(int)degSum);
			//cout<<"randNum = "<<randNum<<" "<<"littleDegNum = "<<littleDegNum<<" ";
			if(randNum<littleDegNum){
				popu[m][i] = 1;
			}
			else{
				popu[m][i] = 0;
			}
		}
		//cout<<endl;
	}
} 

//***************  8  种群中每一个个体进行博弈  **************************
void snowGame(int num){
	//种群中的第num个个体进行博弈 
	int i;
	int t = 1;
	int strategy_temp[N];
	double payoff1[N];
	double payoff0[N];
	
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
			strategy_temp[i] = popu[num][i];
			for (auto j = 0; j<adjList[i].size(); j++){
				int neibor = adjList[i][j];
				if (popu[num][neibor] == 1){
					payoff1[i] = payoff1[i] + 1;
					payoff0[i] = payoff0[i] + 1+r;
				}
				else if (popu[num][neibor] == 0){
					payoff1[i] = payoff1[i] + 1-r;
					payoff0[i] = payoff0[i] + 0;
				}
			}
			if (payoff1[i] > payoff0[i]){
				popu[num][i] = 1;
			}
			else if (payoff1[i]<=payoff0[i]){
				popu[num][i] = 0;
			}

			if (t == 0){
				if (strategy_temp[i] != popu[num][i]){
					t = 1;
				}
				else if (strategy_temp[i] == popu[num][i]){
					t = 0;
				}
			}
		}

	}
}


//***************  9  调用子函数8，种群进行博弈  **************************
void popuSnowGame(){
	for(int i=0; i<M; i++){
		snowGame(i); 
	}
}

//***************  10  p1中第num个个体进行博弈  **************************
void p1_SnowGame(int num){
	int i;
	int t = 1;
	int strategy_temp[N];
	double payoff1[N];
	double payoff0[N];
	
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
			strategy_temp[i] = p1[num][i];
			for (auto j = 0; j<adjList[i].size(); j++){
				int neibor = adjList[i][j];
				if (p1[num][neibor] == 1){
					payoff1[i] = payoff1[i] + 1;
					payoff0[i] = payoff0[i] + 1+r;
				}
				else if (p1[num][neibor] == 0){
					payoff1[i] = payoff1[i] + 1-r;
					payoff0[i] = payoff0[i] + 0;
				}
			}
			if (payoff1[i] >= payoff0[i]){
				p1[num][i] = 1;
			}
			else if (payoff1[i]<payoff0[i]){
				p1[num][i] = 0;
			}

			if (t == 0){
				if (strategy_temp[i] != p1[num][i]){
					t = 1;
				}
				else if (strategy_temp[i] == p1[num][i]){
					t = 0;
				}
			}
		}

	}
}

//***************  11  将p1[n]表示的个体，划分出s1 s2 notResSet  **************************
void partitionSet(int n){
	s1.clear();
	s2.clear();
	notResSet.clear(); 
	for(int i=0; i<N; i++){
		if(p1[n][i]==1){
			int number = 0;
			for(auto j=0; j<adjList[i].size(); j++){
				int nei = adjList[i][j];
				if(p1[n][nei]==0){
					number++;
				}
			}
			if(number==1){
				s1.push_back(i);
			}
			else{
				s2.push_back(i);
			}
		}
		else{
			notResSet.push_back(i);
		}
	}
}

//***************  12  调用子函数11后，在非解集中，找到进行局部搜索的节点  **************************
int findNode(int n){
	partitionSet(n);
	vector<int> nodeTemp; 
	for(auto i=0; i<notResSet.size(); i++){
		int nodeNum = notResSet[i];
		int number = 0;
		for(auto j=0; j<s1.size(); j++){
			int nodeNum2 = s1[j];
			if(matrix[nodeNum][nodeNum2]==1){
				number++;
			}
			if(number>=2){
				//非解集中的节点nodeNum与s1中的2个或2个以上的节点有连边
				nodeTemp.push_back(nodeNum);
			}
		}
	}
	if(nodeTemp.size()<=0){
		return -1;
	}
	else if(nodeTemp.size()==1){
		return nodeTemp[0];
	}
	else{
		int rNum = (int)(nodeTemp.size()*rand())/(RAND_MAX+1);
		int node = nodeTemp[rNum];
		return node;
	}
} 

//***************  13  进行局部搜索的节点的邻居进行博弈  **************************
void neiGame(int n){
	int i,j;
	int t = 1;
	int size = nodeNei.size();
	int strategy_temp[size];
	double payoff1[size];
	double payoff0[size];
	for(i=0; i<size; i++){
		strategy_temp[i] = 0;
		payoff1[i] = 0;
		payoff0[i] = 0;
	}
	while(t){
		t = 0;
		for(i=0; i<size; i++){
			payoff1[i] = 0;
			payoff0[i] = 0;
		}
		for(i=0; i<size; i++){
			int v = nodeNei[i];
			strategy_temp[i] = p1[n][v];
			for(j=0; j<adjList[v].size(); j++){
				int neibor = adjList[v][j];
				if(p1[n][neibor]==1){
					payoff1[i] = payoff1[i]+1;
					payoff0[i] = payoff0[i]+1+r; 
				}
				else if(p1[n][neibor]==0){
					payoff1[i] = payoff1[i]+1-r;
					payoff0[i] = payoff0[i]+0;
				}
			}
			if(payoff1[i]>=payoff0[i]){
				p1[n][v] = 1;
			}
			else if(payoff1[i]<payoff0[i]){
				p1[n][v] = 0;
			}
			if(t==0){
				if(strategy_temp[i]!=p1[n][v]){
					t = 1;
				}
				else if(strategy_temp[i]==p1[n][v]){
					t = 0;
				}
			}
		}
	}
} 

//***************  14  找到进行局部搜索的节点的邻居集合，调用子函数13  **************************
void gameAgain(int n, int node){
	nodeNei.clear(); 
	for(auto i=0; i<adjList[node].size(); i++){
		int nei = adjList[node][i];
		nodeNei.push_back(nei);
	} 
	neiGame(n);//p1[n]中节点node的邻居节点进行博弈，达到纳什均衡 
}

//***************  15  p1[n]进行局部搜索  **************************
void localSearch(int n){
	//p1中的第n个个体进行局部搜索 
	int localSearchNode = findNode(n);
	int lsTimes = 0;
	if(localSearchNode!=-1){
		p1[n][localSearchNode] = 1;//将进行局部搜索的非解集节点纳入解集中
		gameAgain(n, localSearchNode); //这之后，又得到的新的个体p1[n] 
		lsTimes++;
		if(lsTimes>mlsTime){
			localSearchNode = -1;
		}
		else{
			localSearchNode = findNode(n);
		}
	}
}

//***************  16  p1中的后两个个体进行交叉  **************************
void cross(){
	/*int r1 = (int)(N*rand())/(RAND_MAX+1);
	int r2 = r1 + length;
	for(int i=r1; i<r2; i++){
		int temp = p1[2][i];
		p1[2][i] = p1[3][i];
		p1[3][i] = temp;
	}*/
	
	int r1 = (int)(N*rand())/(RAND_MAX+1);
	int r2 = (int)(N*rand())/(RAND_MAX+1);
	if(r1>=r2){
		int temp = r1;
		r1 = r2;
		r2 = temp;
	}
	for(int i=r1; i<r2; i++){
		int temp = p1[2][i];
		p1[2][i] = p1[3][i];
		p1[3][i] = temp;
	}
}

//***************  17  p1中的后两个个体进行变异  **************************
void mutation(int n){
	int pro;
	for(int i=0; i<N; i++){
		pro = rand()%N;
		if(pro<Proba){
			if(p1[n][i]==1){
				p1[n][i] = 0;
			}
			else{
				p1[n][i] = 1;
			}
		}
	}
	/*int r1 = (int)(N*rand())/(RAND_MAX+1);
	p1[2][r1] = 1-p1[2][r1];
	p1[3][r1] = 1-p1[3][r1];*/
}

//***************  18  计算每个个体代表的解集的大小  **************************
int computeSize(array<int, N> vec){
	int size = 0;
	for(int i=0; i<N; i++){
		if(vec[i]==1){
			size++;
		}
	}
	return size;
}

int computeFitness(array<int, N> vec){
	int fitValue = 0;
	for(int i=0; i<N; i++){
		if(vec[i]==0){
			auto j = adjList[i].begin();
			while(j!=adjList[i].end()){
				if(vec[*j]==0){
					fitValue+=N;
				}
				++j;
			}
		}
	}
	for(int i=0; i<N; i++){
		if(vec[i]==1){
			fitValue++;
		}
	} 
	return fitValue;
}

//***************  19  在父代个体和子代个体中选择最优的两个个体  **************************
void choose(int t){
	vector<int> size;
	int s[4];
	s[0] = computeFitness(p1[0]);
	s[1] = computeFitness(p1[1]);
	s[2] = computeFitness(p1[2]);
	s[3] = computeFitness(p1[3]);
	/*s[0] = computeSize(p1[0]);
	s[1] = computeSize(p1[1]);
	s[2] = computeSize(p1[2]);
	s[3] = computeSize(p1[3]);*/
	size.push_back(s[0]);
	size.push_back(s[1]);
	size.push_back(s[2]);
	size.push_back(s[3]);
	sort(size.begin(),size.end());
	int max1 = size[0];
	int max2 = size[1];
	int n = 0;
	for(int i=0; i<4; i++){
		if(n<=2&&s[i]==max1){
			newPopu[t] = p1[i];
			n++;
		}
		if(n<=2&&s[i]==max2){
			newPopu[t+1] = p1[i];
			//newPopu[t+M/2] = p1[i];
			n++;
		}
	}
} 

//***************  20  新种群替代旧种群  **************************
void newPopu_to_popu(){
	for(int m=0; m<M; m++){
		popu[m] = newPopu[m];
	}
}

//***************  21  随机打乱种群  **************************
void individualSwap(){
	int x = rand()%M;
	int y = rand()%M;
	array<int, N> temp;
	temp = popu[x];
	popu[x] = popu[y];
	popu[y] = temp;
} 

//*******************************  20  计算RunTime次运行的解集大小平均值  *****************************
double computeAverage()
{
	int i;
	int j;
	double average = 0;
	int sum = 0;
	for( i=0; i<RunTime; i++ )
	{
		sum = sum + totalCandidate[i];
	}
	average = (double)( (double)sum/(double)RunTime );
	
	return average;
}

//******************************  21  计算RunTime次运行的解集大小的标准差  ******************************* 
double computeStandardDeviation()
{
	int i;
	double average = computeAverage();
	double standardDeviation = 0;
	double sum = 0;
	for( i=0; i<RunTime; i++ )
	{
		sum = sum + ( totalCandidate[i]-average )*( totalCandidate[i]-average );
	}
	standardDeviation = (double)( sqrt( (double)sum/(double)RunTime ) );
	
	return standardDeviation;
}

//*****************************  22 计算RunTime次运行，得到的最大解  *************************************
int computeMaxCandidate()
{
	int i;
	int max = totalCandidate[0];
	for( i=0; i<RunTime; i++ )
	{
		if( max<totalCandidate[i] )
		{
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
	for( i=0; i<RunTime; i++ )
	{
		if( min>totalCandidate[i] )
		{
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
	for( i=0; i<RunTime; i++ )
	{
		if( totalCandidate[i]==max )
		{
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
	for( i=0; i<RunTime; i++ )
	{
		sum = sum + totalTime[i];
	}
	averageTime = (double)( (double)sum/(double)RunTime );
	
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
	
	ofstream outfile_tm(filename+"_newGEA_Result.txt");
	if (!outfile_tm)
	{
		cerr << "open error!" << endl;
		exit(1);
	}
	outfile_tm << "GEA+LS results : "<<endl<<endl;
	
	outfile_tm << "Network vertices = " << N <<", ";
	outfile_tm << "the total run times is : "<< RunTime << endl<< endl;
	
	outfile_tm << "total candidate of RunTime: "<< " ";
	for( i=0; i<RunTime; i++ )
	{
		outfile_tm << totalCandidate[i] << " ";
	}
	outfile_tm << endl;
	
	outfile_tm<<"averageCandidate: "<<averageCandidate<<endl<< endl;
	outfile_tm<<"stdCandidate: "<<stdCandidate<<endl<< endl;
	outfile_tm<<"maxCandidate: "<<maxCandidate<<endl<< endl;
	outfile_tm<<"minCandidate: "<<minCandidate<<endl<< endl;
	outfile_tm<<"maxCandidateNumber: "<<maxCandidateNumber<<endl<< endl;
	outfile_tm<<"averageTime: "<<averageTime<<endl<< endl;
	
	outfile_tm << "edge number =  : "<< edgenumber << endl<< endl;
	outfile_tm << "the maximum degree of the network is : "<< maxdegree << endl<< endl;
	outfile_tm << "the minimum degree of the network is : "<< mindegree << endl<< endl;
	outfile_tm << "the average degree of the network is : "<< avedegree << endl<< endl;
	
	outfile_tm << "the solution of RunTime: "<<endl;
	for( i=0; i<RunTime; i++ )
	{
		for( j=0; j<N; j++ )
		{
			if( solu[i][j]==1)
			{
				outfile_tm << j <<" ";
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
	
	printf("*********************  共运行 %d 次  *************************：\n\n",RunTime);
	string strName = filename+"_result";
	cout<<strName<<endl;
	printf("Network vertices = %d, run times = %d\n\n",N, RunTime);
	/*sprintf("1、得到的解集中的元素分别为：\n");
	for( i=0; i<RunTime; i++ )
	{
		for( j=0; j<N; j++ )
		{
			if( solu[i][j]==1 )
			{
				printf("%d ",j);
			}			
		}
		printf("\n");
	}
	printf("\n");*/ 
	
	printf("2、得到的解集大小分别为： ");
	for( i=0; i<RunTime; i++ )
	{
		printf("%d ",totalCandidate[i]);
	}
	printf("\n\n");
	
	printf("3、解集大小的平均值为：%f，标准差为：%f \n\n", averageCandidate, stdCandidate);
	
	printf("4、解集大小最大值为：%d，最小值为：%d, 取得最大值的次数为：%d\n\n", maxCandidate, minCandidate, maxCandidateNumber);
	
	printf("5、平均运行时间为：%f \n\n",averageTime); 
	
	
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


