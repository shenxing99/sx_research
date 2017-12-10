#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<fstream>
#include<iostream>
#include<vector>
#include<array>
#include<string>

using namespace std;

#define N 2000  //顶点数
string filename = "BA2000(8).txt";
#define swapTimes 10000
#define b 1.5 // PDG参数 
#define RunTime 30 // 统计RunTime次的运行结果，求平均 
#define T 1000 //最大迭代次数，ILS的终止条件
#define K 4
//********************************************  变量定义  ********************************************* 
vector< vector<int> > adjList; //网络邻接表，存储节点的邻居节点的标号 
vector<int> adjList_line;

vector<int> candidate; //存储所有可能存在2-improvement的节点，存储节点的标号 

int matrix[N][N]; //  1  网络的邻接矩阵
array< array<int, N>, N > matrix_temp;
int deg[N]; //  3  网络中每个节点的度
int maxdegree;
int mindegree;
double avedegree;
int edgeNumber;

array<int, N> solution; // 解集 
array<int, N> solutionBefore; // 解集 
array<int, N> solutionAfter; // 解集 
int perSolution[N]; //  5  用三块结构表示的解，这个结构存储节点的标号，从0开始到N-1 

int tightness[N]; //  7  记录每一个节点的tightness 

int brock1Num; //  8  第一块空间中的节点个数，即独立集的元素个数
int brock2Num; //  9  第二块空间中的节点个数，即自由节点的节点个数 
int brock3Num; //  10  第三块空间中的节点个数，即既不是独立集也不是自由节点的节点个数 

int strategy[N]; // N各节点的策略，0代表背叛，1代表合作
double payoff1[N]; 
double payoff0[N];
int strategy_temp[N];
int set[N]; 

array< array<int, N>, RunTime > solu; // RunTime次运行，每一次得到的解集
int totalcandidate[RunTime]; // RunTime次运行，每一次得到的解集的元素个数，即得到的独立集的元素个数
double totalTime[RunTime]; // RunTime次运行，每一次的运行时间
double averageCandidate;// RunTime次运行，独立集大小的平均值 
double stdCandidate;// RunTime次运行，独立集大小的标准差
int maxCandidate;// RunTime次运行中candidate的最大值 
int minCandidate;// RunTime次运行中candidate的最小值 
int maxCandidateNumber;// RunTime次运行中得到最大candidate的次数
double averageTime;// RunTime次运行的平均运行时间

int BeforeSize;
int AfterSize;

int S;

array< array<int, N>, T > TIterSolution; // T次迭代，每一次得到得到的解集
array< int, T > TiterSolutionCandidate; // T次迭代，记录每一次迭代时，解集的大小
array< int, N > vertexTimes;


//******************************************  子函数声明  *********************************************
void initial(); // 0
void read_Without_e(); 
void initial_RunTime();
void read(); // 1
void complement_graph(); 
void readMatrix();
void obtainAdjList(); // 2
void obtainVerDegree(); // 3
void randomAlgorithm(); // 4
void greedyAlgorithm(); // 5
void computeTightness(); // 8
void permutationSolution(); // 9
void initialCandidate(); // 10
void simpleCandidate(); // 11
void removeTestedVertexFromCandidate( int testedVertexNum ); // 12
void updateCandidateInsert( int insertNum ); // 13  insertNum为插入节点的标号
void updateCandidateRemove( int removeNum ); // 14  removeNum为移出节点的标号 
void removeVerSolution( int removeNum ); // 15
void insertVerSolution( int insertNum ); // 16
void twoImprovement(); // 17  

void initialPDGStrategy();
void PDG();

int computeSolutionCandidate();
double computeAverage(); 
double computeStandardDeviation();
double computeAverageTime();
int computeMaxCandidate();
int computeMinCandidate();
int computeMaxTimes();
void saveResults();
void printResults();

void IndependentSetOrNot();
void maximalIndependentSetOrNot();

void initialVariable_Iteration();
void perturbation(); 
void recordTimes();

//********************************************  main  **************************************************
int main(int argc, char* argv[])
{
	
	int i;
	int j;
	int run;
	int isNumBefore; // 2-improvement LS 前，解集的元素个数 
	int isNumAfter; // 2-improvement LS 后，解集的元素个数 
	int solutionCandidate = 0; // 解集的元素个数 
	
	clock_t start1,end1,start;
	clock_t start2,end2,end; 
	double time1; //记录产生初始解得算法的运行时间 
	double time2; //记录2-improvement的运行时间 
	double runningTime;
	double perturb_time;
	srand( (unsigned int)time(NULL) );
	
	initial(); //变量初始化 
	initial_RunTime();
	
	//read();	//读取网络信息 	
	//complement_graph();
	//read_Without_e();
	readMatrix();
	
	obtainAdjList(); //得到网络的邻接表 	
	obtainVerDegree(); //得到每个节点的度 
			 
	//**************************  2-improvement local search  **************************
	for( run=0; run<RunTime; run++ ) //RunTime次独立运行 
	{
		S = 0;
		initialVariable_Iteration();
		
		//randomAlgorithm(); //用random algorithm产生一个初始解 
		greedyAlgorithm(); //用greedy algorithm产生一个初始解
		//initialPDGStrategy();
		//PDG();
		recordTimes();
		BeforeSize = computeSolutionCandidate();
		
		start = clock();
		perturb_time = 0;
		int itertimes;
		int Stimes;
		for( int iter=0; iter<T; iter++ ) 
		{
			itertimes = iter;
			isNumBefore = 0;
			isNumAfter = 0;
			isNumBefore = computeSolutionCandidate();
			solutionBefore = solution;
			
			start1 = clock();
			perturbation(); 
			end1 = clock();
			perturb_time = perturb_time+(double) (end1-start1)/CLOCKS_PER_SEC;
			recordTimes();
			
			computeTightness(); //计算节点的tightness，只有非独立集的节点才有tightness 		
			permutationSolution(); // 用三块的结构表示独立集、自由节点和非独立集且非自由节点的节点 	
			initialCandidate(); // 根据上一次2-improvement 的结果solution，初始这一次的candidate，重新进入下一轮2-improvement 
			simpleCandidate();
			
			twoImprovement(); // 2-improvement后，得到新的solution，旧solution存储在TIterSolution[t][]中 
			recordTimes();
				
			isNumAfter = computeSolutionCandidate();
			solutionAfter = solution;
			if( isNumAfter>=isNumBefore )
			{
				solution = solutionAfter;
			}
			else if( isNumAfter<isNumBefore )
			{
				int pro = (int)( (100 * rand()) / (RAND_MAX+1) );
				if( pro<5 )
				{
					solution = solutionAfter;
				}
				else
				{
					solution = solutionBefore;
				}
			}
			if(S==swapTimes){
				Stimes = S;
				break;
			}
			
		}	
		//printf("iter = %d, swap = %d,   ", itertimes, S);
		end = clock();
		runningTime = (double) (end-start)/CLOCKS_PER_SEC;
		runningTime = runningTime-perturb_time;
		//time2 = (double) (end2-start2)/CLOCKS_PER_SEC;
		//time2 = time1 + time2; //初始算法+2-improvement的总运行时间 
		
		
		//printf("perturb_time = %f\n", perturb_time);
		solu[run] = solution;
		AfterSize = computeSolutionCandidate();
		totalcandidate[run] = computeSolutionCandidate();
		totalTime[run] = runningTime;
		printf("run = %d, BeforeSize = %d, AfterSize = %d\n", run+1, BeforeSize, AfterSize); 
		//IndependentSetOrNot();
		//maximalIndependentSetOrNot();
	}
	saveResults();
	printResults();
	
	return 0;
}

//***************  0  读取文件  **************************
void read() // 得到了网络邻接矩阵 ， 并记录了边的矩阵edgeMatrix[M][2] 
{
	int i;
	int j;
	FILE* fp = fopen(filename.c_str(),"r"); 
	//FILE* fp = fopen("brock200_2 (12).txt","r");
	
	if( !fp )
	{
		printf("文件打开错误");
	}

	char e;
	int startnode;
	int endnode;
	
	i=0;
	while( !feof(fp) )
	{
		fscanf(fp, "%c %d %d \n", &e, &startnode, &endnode );
		matrix[startnode-1][endnode-1] = 1;
		matrix[endnode-1][startnode-1] = 1;
		//edgeMatrix[i][0] = startnode;
		//edgeMatrix[i][1] = endnode;
		i = i+1;
	}

	fclose( fp );
}

//***************  1  读取文件（网络的边的表示  不带e）  **************************
void read_Without_e()
{
	int i;
	int j;
	FILE* fp = fopen(filename.c_str(),"r"); 
	//FILE* fp = fopen("1dc.1024_24063.txt","r");
	//FILE* fp = fopen("vertex6.txt","r");
	//FILE* fp = fopen("frb30-15-1_450_17827.txt","r");
	
	if( !fp )
	{
		printf("文件打开错误");
	}

	char e;
	int startnode;
	int endnode;
	
	i=0;
	while( !feof(fp) )
	{
		fscanf(fp, "%d %d \n",&startnode, &endnode );
		matrix[startnode-1][endnode-1] = 1;
		matrix[endnode-1][startnode-1] = 1;
		i = i+1;
	}

	fclose( fp );
}

//***************  3  读取邻接矩阵表示的网络  **************************
void readMatrix()
{
	int i;
	int j;
	
	FILE* fp;
	fp = fopen(filename.c_str(),"r"); 
	//fp = fopen( "BA500(10).txt" , "r" );
	if ( fp==NULL )
	{
		printf("文件打开失败");
		exit( 1 );
	}
	char C;
	int row=0;
	int col=0;
	int i1 = 0;//由于matrix的下标需要每次递增1，所以这里另外设置两个迭代器i1和j1
	int j1 = 0;
	while( ( C=fgetc(fp) ) != EOF )
	{
		if ( C=='0' )
		{
			matrix[i1][j1]=0;
			col++;
			j1++;
		}
		else if ( C=='1' )
		{
			matrix[i1][j1]=1;
			col++;
			j1++;
		}
		else if ( C== ' ' )
		{
			col++;
		}
		else if( C == '\n')
		{
			col++;
			row++;
			i1++;
			j1 = 0;
		}
	}
	fclose(fp);
} 

//***************  6  求网络的补图  **************************
void complement_graph()
{
	int i;
	int j;
	
	//int matrix_temp[N][N];
	
	for( i=0; i<N; i++ )
	{
		for( j=0; j<N; j++ )
		{
			matrix_temp[i][j] = matrix[i][j];
		}
	}
	
	for( i=0; i<N; i++ )
	{
		for( j=0; j<N; j++ )
		{
			if( matrix_temp[i][j]==1 )
			{
				matrix[i][j] = 0;
			}
			else if( matrix_temp[i][j]==0 )
			{
				matrix[i][j] = 1;
			}
		}
	}
	
	for( i=0; i<N; i++ )
	{
		matrix[i][i] = 0; 
	}
	
} 

//****************  1  变量初始化  ****************************
void initial()
{
	int i;
	int j;
	int k;
	
	for( i=0; i<N; i++ )
	{
		for( j=0; j<N; j++ )
		{
			matrix[i][j] = 0;
		}
	}

	for( i=0; i<N; i++ )
	{
		deg[i] = 0;
		solution[i] = 3; //对结果表示进行初始化 
		perSolution[i] = 0; //这个结构存储节点的标号 
		tightness[i] = 0; //每个节点的tightness初始化为0 
	}

	for( i=0; i<N; i++ )
	{
		adjList.push_back( adjList_line );
	}
		
	brock1Num = 0;
	brock2Num = 0;
	brock3Num = 0; 
	
	for( i=0; i<N; i++ )
	{
		strategy[i] = 0;
		payoff1[i] = 0;
		payoff0[i] = 0;
		strategy_temp[i] = 0;
		set[i] = 0;
	}
	
	maxdegree = 0;
	mindegree = 0;
	avedegree = 0;
	edgeNumber = 0;
}

//****************  2  变量初始化  ****************************
void initial_RunTime()
{
	int i;
	int j;
	for( i=0; i<RunTime; i++ )
	{
		totalcandidate[i] = 0;
		totalTime[i] = 0;
	}
	for( i=0; i<RunTime; i++ )
	{
		for( j=0; j<N; j++ )
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

//****************  3  得到网络的邻接表(只记录邻居  网络节点编号从0开始)  *************
void obtainAdjList()
{
	int i;
	int j;
	

	for( i=0; i<N; i++ )
	{
		for( j=0; j<N; j++ )
		{
			if( matrix[i][j]==1 )
			{
				adjList[i].push_back(j); //  得到的邻接表中，每一个节点的邻居节点都以标号的升序排列 
			}
		}
	}

}

//**************  4  得到每个节点的度  ****************
void obtainVerDegree()
{
	int i;

	for( i=0; i<N; i++ )
	{
		deg[i] = adjList[i].size();
	}
	
	int max = deg[0];
	int min = deg[0];
	double ave;
	double sum = 0;
	
	for( i=0; i<N; i++ )
	{
		if( max<deg[i] )
		{
			max = deg[i];
		}
	}
	for( i=0; i<N; i++ )
	{
		if( min>deg[i] )
		{
			min = deg[i];
		}
	}
	
	for( i=0; i<N; i++ )
	{
		sum = sum+deg[i];
	}
	ave = (double)( (double)(sum) / (double)(N) );
	
	maxdegree = max;
	mindegree = min;
	avedegree = ave;
	edgeNumber = sum / 2;

}

//**************  5  random algorithm 产生初始解  *************************
void randomAlgorithm()
{
	int i;
	int j;
	int randomNum1;
	int ver;
	vector<int> temp;
		
	for( i=0; i<N; i++ )
	{
		solution[i] = 0;
	}
	
	for( i=0; i<N; i++ )
	{
		tightness[i] = 0;
	}
	
	for( i=0; i<N; i++ )
	{
		if( tightness[i]==0 )
		{
			temp.push_back( i ); // 把自由节点的标号放入temp 
		}
	}
	
	while( (temp.size()) != 0 )
	{
		/*temp.clear();
		
		for( i=0; i<N; i++ )
		{
			if( tightness[i]==0 )
			{
				temp.push_back( i ); // 把自由节点的标号放入temp 
			}
		}*/
		
		//srand( (unsigned int)time(NULL) );
		randomNum1 = (int)( (temp.size() * rand()) / (RAND_MAX+1) );
		ver = temp[randomNum1];
		solution[ver] = 1;
	
		computeTightness(); 
		temp.clear();
		for( i=0; i<N; i++ )
		{
			if( tightness[i]==0 )
			{
				temp.push_back( i ); // 把自由节点的标号放入temp 
			}
		}
	}
		
} 

//**************  6  greedy algorithm 产生初始解  *************************
void greedyAlgorithm()
{
	int i;
	int j;
	int freeVer;
	int costS;
	int minCost;
	int minVer;
	
	vector< vector<int> > freeVertexCost;
	vector<int> freeVertexCost_line;
	
	vector<int> freeVertex; 
	
	for( i=0; i<2; i++ )
	{
		freeVertexCost.push_back( freeVertexCost_line );
	}
	
	for( i=0; i<N; i++ )
	{
		solution[i] = 0;
	}
	
	for( i=0; i<N; i++ )
	{
		tightness[i] = 0;
	}
		
	int r = (int)( N * rand()) / (RAND_MAX+1);
	solution[r] = 1;
	computeTightness();
	for( i=0; i<N; i++ )
	{
		if( tightness[i]==0 )  //记录自由节点的标号 
		{ 
			freeVertex.push_back(i);
		}
	}
	
	while( (freeVertex.size())!=0 )
	{
		freeVertex.clear();
		//freeVertexCost.clear();
		freeVertexCost[0].clear(); 
		freeVertexCost[1].clear();
		
		
		
		for( i=0; i<N; i++ )
		{
			if( tightness[i]==0 )  // 需要找到这个自由节点的邻居中，自由节点的个数 
			{
				freeVertexCost[0].push_back(i); 
				freeVertex.push_back(i);
			}
		}
	
		for( j=0; j<freeVertexCost[0].size(); j++ )
		{
			freeVer = freeVertexCost[0][j];
			costS = 0;
			for( i=0; i<N; i++ )
			{			
				if( (matrix[freeVer][i]==1)&&(tightness[i]==0) )
				{
					costS = costS + 1;
				} 
			}
			freeVertexCost[1].push_back(costS); 
		}
		
		//The greedy algorithm assigns a cost to each free vertex equal to the number of free neighbors it  has,
		//and in each iteration picks a free vertex with the lowest cost.
		minCost = freeVertexCost[1][0];
		minVer = freeVertexCost[0][0]; 
		
		for( i=0; i<freeVertexCost[0].size(); i++ )
		{
			if( minCost > freeVertexCost[1][i] )
			{
				minCost = freeVertexCost[1][i];
				minVer = freeVertexCost[0][i]; //minVer就是要进入解集的节点 
			}		 
		}
	
		solution[minVer] = 1;
	
		computeTightness(); 
	}
	
}

//**************  7  计算节点的tightness  **********************
void computeTightness()
{
	int i;
	int j;
	int k;
	int s;
	for( i=0; i<N; i++ )
	{
		if( solution[i]==0 ) //只有独立集以外的节点才有tightness
		{
			s = 0;
			for( j=0; j<deg[i]; j++ )
			{
				if(solution[ adjList[i][j] ]==1)
				{
					s = s+1;
				}
			}
			tightness[i] = s;
		} 
		else
		{
			tightness[i] = N; //独立集中的节点没有tightness ，此时，将独立集中的节点的tightness设置成一个非常大的数N 
		}
	}
}

//**************  8 将solution分为三块，并计算每一块中的元素个数  ***********************
void permutationSolution()
{
	int i;
	int j;
	int brock1_i;
	int brock2_i;
	int brock3_i;
	
	computeTightness();
	
	brock1_i = 0;
	for( i=0; i<N; i++ )
	{
		if( solution[i]==1 )
		{
			perSolution[brock1_i] = i;
			brock1_i = brock1_i + 1; //brock1_i 是第一块空间中的节点个数，即独立集的元素个数 
		}
	}
	brock1Num = brock1_i;
	
	brock2_i = brock1_i;
	for( i=0; i<N; i++ )
	{
		if( tightness[i]==0 ) //是自由节点
		{
			perSolution[brock2_i] = i;
			brock2_i = brock2_i + 1; //brock2_i-brock1_i 是第二块空间中的节点个数，即自由节点的节点个数 
		} 
	}
	brock2Num = brock2_i-brock1_i;
	
	brock3_i = brock2_i;
	for( i=0; i<N; i++ )
	{
		if( (tightness[i]>0)&&(tightness[i]!=N) )
		{
			perSolution[brock3_i] = i;
			brock3_i = brock3_i + 1; //brock3_i-brock2_i 是第三块空间中的节点个数，即不是独立集也不是自由节点 
		}
	}
	brock3Num = brock3_i-brock2_i;
	
} 

//**************  9  初始化candidate，初始化所有极大独立集中的节点都在candidate中  *********************
void initialCandidate()
{
	int i;
	
	for( i=0; i<N; i++ )
	{
		if( solution[i]==1 )
		{
			candidate.push_back( i );
		}

	}
	
} 

//**************  10  简化candidate，从candidate中删除不可能存在2-improvement的节点  *********************
void simpleCandidate()
{
	int i;
	int j;
	int k;
	int t; 
	int ver;
	int indexS;
	
	vector<int> tight1Neighbor; 
	vector<int> record; 
	vector<int> tempCandidate;

	computeTightness();

	record.clear();
	//solution中没有两个相互独立的1-tight邻居节点的节点，不存在2-improvement，可以从candidate中删除 
	for( i=0; i<candidate.size(); i++ ) //检测candidate中的每一个节点，是否可以被移出candidate
	{
		ver = candidate[i];
		
		tight1Neighbor.clear();
		//record.clear();
		
		for( j=0; j<deg[ ver ]; j++ ) //遍历candidate[i]的所有邻居 
		{
			if( tightness[ adjList[ver][j] ]==1 ) //ver的1-tight邻居节点 
			{
				tight1Neighbor.push_back( adjList[ver][j] ); //找到candidate[i]的所有1-tight邻居 
			} 
		} 
		
		indexS = 0;
		if( tight1Neighbor.size()>=2 )
		{			
			for( k=0; k<tight1Neighbor.size()-1; k++ )
			{
				for( t=k+1; t<tight1Neighbor.size(); t++ )
				{
					if( matrix[ tight1Neighbor[k] ][ tight1Neighbor[t] ]==1 ) //tight1Neighbor[k]
					{
						indexS = indexS+1;  // indexS是指示器，1-tight邻居中，两个邻居如果相邻，则加1
									    // 如果所有这些邻居都相邻，则indexS就等于这些节点全连通时的边数 
					}
				}
			}
	
			if( indexS==( ( (tight1Neighbor.size())*(tight1Neighbor.size()-1) )/2 )  ) //总比较次数   全连通的边数 
			{
			//比较多少次，indexS就是多少，则说明他们相互有连边，即candidate[i]没有相互独立的1-tight邻居，可以移出candidate
			//record.push_back(i); //记录要移出candidate的节点在candidate中的位置
				record.push_back(ver); 			 
			} 
		}
		else
		{
			record.push_back(ver);
		}
	
	} 
	
	tempCandidate.clear();
	
	for( i=0; i<record.size(); i++ )
	{
		for( j=0; j<candidate.size(); j++ )
		{
			if( record[i]==candidate[j] ) //candidate[j]要被删除
			{
				candidate[j] = N; 
			} 
		}
	}
	
	for( i=0; i<candidate.size(); i++ )
	{
		if( candidate[i]!=N ) //tempCandidate存储要保留下来的candidate元素 
		{
			tempCandidate.push_back(candidate[i]);
		}
	}
	
	candidate.clear();
	for( i=0; i<tempCandidate.size(); i++ )
	{
		candidate.push_back( tempCandidate[i] );
	}
	 
} 

//**************  11  遍历某个节点时，将该节点从candidate中移出，不重复遍历节点  ********************
void removeTestedVertexFromCandidate( int testedVertexNum )
{
	int i;
	vector<int> temp;
	
	temp.clear(); 
	for( i=0; i<candidate.size(); i++ )
	{
		if( candidate[i]!=testedVertexNum )
		{
			temp.push_back( candidate[i] );
		}
	}
	
	candidate.clear(); 
	if( temp.size()>0 )
	{
		for( i=0; i<temp.size(); i++ )
		{
			candidate.push_back( temp[i] );
		}
	}
	
} 

//**************  12  有节点插入solution时，更新candidates  *********************
void updateCandidateInsert( int insertNum ) //insertNum为插入节点的标号 
{
	int i;
	int j;
	
	candidate.push_back( insertNum );
	
}

//**************  13  有节点移出solution时，更新candidates  *********************
void updateCandidateRemove( int removeNum ) //removeNum为移出节点的标号 
{
	int i;
	int j;
	
	solution[removeNum] = 0; //removeNum插入到独立集中 
	permutationSolution();
	
	//寻找因为removeNum的移出，removeNum的邻居节点中变为1-tight的邻居节点，
	//这些邻居节点在solution中的邻接节点要进入candidate 
	computeTightness(); //移出一个节点后，更新tightness
	
	vector<int> reTight1;
	for( i=0; i<deg[removeNum]; i++ )
	{
		if( tightness[ adjList[removeNum][i] ]==1 )
		{
			reTight1.push_back( adjList[removeNum][i] ); //找到因为removeNum的移出，removeNum的邻居节点中变为1-tight的邻居节点
		}
	}
	
	for( i=0; i<reTight1.size(); i++ )
	{
		for( j=0; j<brock1Num; j++ )
		{
			if( matrix[ reTight1[i] ][ perSolution[j] ]==1 )
			{
				candidate.push_back( perSolution[j] );
			}
		}
	} 
	
}

//**************  14  从解集中移出一个节点  *********************
void removeVerSolution( int removeNum )
{
	solution[ removeNum ] = 0; 
	
}

//**************  15  解集中插入一个新节点  *********************
void insertVerSolution( int insertNum )
{	
	solution[ insertNum ] = 1; 
	
}

//**************  16  2-improvement  *************************
void twoImprovement()
{
	int i;
	int j;
	int k;
	
	int removeVertex;
	int tempVertex1;
	int tempVertex2;
	
	//for( i=0; i<candidate.size(); i++ )  //对于candidate中的所有节点，都执行一次2-improvement 
	//if( candidate.size()>0 )
	while( candidate.size()>0 )
	{
		S = S+1;
		removeVertex = candidate[0]; //每次检测candidate的第一个元素，直到candidate没有元素，遍历结束 
		//removeVertex = candidate[i];
		
		removeTestedVertexFromCandidate( removeVertex ); // removeVertex被遍历了，把removeVertex从candidate中移出 
		
		removeVerSolution( removeVertex ); // solution改变了，只要solution改变，就要更新candidate 
										   // 从solution中移出removeVertex，solution由S变为S' 
		// removeVertex被临时移走，但是不意味着removeVertex立马变为非独立集中的节点？？？
		 
		permutationSolution(); //solution改变，将solution表示为三块的结构也要相应改变
		computeTightness();
		
		updateCandidateRemove( removeVertex ); // 根据现在的解S'，更新candidate，更新后，再进行简化 
		simpleCandidate(); 
		 
		int index1 = 0; //指示器 
		int index2 = 0;
		
		if( brock2Num>=2 ) //从S中移出一个节点x后，得到 S',如果S'的自由节点数大于等于2，则进行以下操作，否则，说明x的移出不能产生2-improvement 
		{
			for( j=brock1Num; j<brock1Num+brock2Num; j++ ) //从S'的自由节点中找到x=removeVertex的邻居节点 
			{
				tempVertex1 = perSolution[j]; // j要取到perSolution的中间的自由节点那一块
											  // 这一块的起始下标是brock1Num,终止下标是brock1Num+brock2Num-1 
				if( matrix[removeVertex][tempVertex1]==1 )
				{
					index1 = index1 + 1; //记录找到的removeVertex的是S’的自由节点的邻居的个数
					 
					insertVerSolution( tempVertex1 ); //将S'的自由节点中x的邻居节点插入S'，得到S'' 
					permutationSolution();
					computeTightness();
					
					updateCandidateInsert( tempVertex1 ); //根据现在的解S”，更新candidate，更新后，再进行简化 
					simpleCandidate();
					 
					if( brock2Num>=1 ) //如果S''有一个自由节点，则将该自由节点插入S''，完成一个2-improvement 
					{
						tempVertex2 = perSolution[brock1Num]; //取到第二块的第一个元素，也是最后一个元素 
						
						insertVerSolution( tempVertex2 );
						permutationSolution(); 
						computeTightness();
						
						updateCandidateInsert( tempVertex2 );
						simpleCandidate();
						
						j = N; //x完成一个2-improvement，跳出找邻居节点的循环 
					}
					else
					{
						removeVerSolution( tempVertex1 ); //如果S''自由节点数不是1（是0），则
														  //将x从S''移出，重新得到S'，并检查x在S'的自由节点中的另一个邻居节点 
						permutationSolution(); 
						computeTightness();
						
						updateCandidateRemove( tempVertex1 );
						simpleCandidate();	
						
						index2 = index2 + 1; // 当前removeVertex用于2-improvement时，
											 // 遍历其是S'的自由节点的邻居，但是该邻居不能得到2-improvement时，index2加1 					
					}
				}
			}
			if( index1==index2 ) // removeVertex的是S’的自由节点的邻居的个数 = 对removeVertex执行2-improvement失败的次数数 
			{
				insertVerSolution( removeVertex ); 
				permutationSolution();
				computeTightness(); 
						
				updateCandidateInsert( removeVertex );
				simpleCandidate();
			} 
		}  
		else // 把removeVertex重新放回solution中 
		{
			insertVerSolution( removeVertex ); 
			permutationSolution();
			computeTightness(); 
						
			updateCandidateInsert( removeVertex );
			simpleCandidate();
		}
	}
	
	 // 如果还有自由节点，则把自由节点放入解集
	 permutationSolution();
	 while( brock2Num!=0 )
	 {
	 	solution[ perSolution[brock1Num] ] = 1;
	 	permutationSolution();
	 } 
	  
}

//**************  17  PDG game 随机初始化所有节点的策略  ********************
void initialPDGStrategy()
{
	int i;
	int j;
	
	//随机初始化所有节点的策略
	for( i=0; i<N; i++ )
	{
		strategy[i] = 0;
	} 
	//srand( (unsigned int)time(NULL) );
	for( i=0; i<N; i++ )
	{
		strategy[i] = (int)( 100.0*rand() ) / ( RAND_MAX+1 ); //产生0到99的随机数 
	}
	for( i=0; i<N; i++ )
	{
		if(strategy[i]<50)
		{
			strategy[i] = 0;
		}
		else
		{
			strategy[i] = 1;
		}
	} 	
}

//**************  18  PDG game 博弈过程  根据已有的strategy进行博弈  *************
void PDG()
{
	int i;
	int j;
	int t = 1;
	int count = 0;
	
	for( i=0; i<N; i++ )
	{
		strategy_temp[i] = 0;
	} 
	
	for( i=0; i<N; i++ )
	{
		payoff1[i] = 0;
		payoff0[i] = 0;
	}
	
	while(t)
	{
		t = 0;
		count = count+1;
		
		for( i=0; i<N; i++ )
		{
			payoff1[i] = 0;
			payoff0[i] = 0;
		}
		
		for( i=0; i<N; i++ )
		{
			strategy_temp[i] = strategy[i];
			for( j=0; j<N; j++ )
			{
				if( matrix[i][j]==1 )
				{
					if( strategy[j]==1 )
					{
						payoff1[i] = payoff1[i] + 1;
						payoff0[i] = payoff0[i] + b;
					}
					else if( strategy[j]==0 )
					{
						payoff1[i] = payoff1[i] + 0;
						payoff0[i] = payoff0[i] + 0;
					}
				}
			}
			if( payoff1[i]>=payoff0[i] )
			{
				strategy[i] = 1;
			}
			else
			{
				strategy[i] = 0;
			}
			
			if( t==0 )
			{
				if( strategy_temp[i]!=strategy[i] )
				{
					t = 1;
				}
				else if( strategy_temp[i]==strategy[i] )
				{
					t = 0;
				}
			}
		}
	}
	
	for( i=0; i<N; i++ )
	{
		solution[i] = strategy[i];
	}
	
} 

//**************  19  计算solution中元素个数  *************************
int computeSolutionCandidate()
{
	int i;
	int j;
	int solutionCandidate = 0;
	for( i=0; i<N; i++ )
	{
		if( solution[i]==1 )
		{
			solutionCandidate = solutionCandidate + 1;
		}
	}
	
	return solutionCandidate;
	
}

//***************  20  计算RunTime次运行的解集大小平均值**************************
double computeAverage()
{
	int i;
	int j;
	double average = 0;
	int sum = 0;
	for( i=0; i<RunTime; i++ )
	{
		sum = sum + totalcandidate[i];
	}
	average = (double)( (double)sum/(double)RunTime );
	
	return average;
}

//***************  21 计算RunTime次运行的解集大小的标准差**************************
double computeStandardDeviation()
{
	int i;
	double average = computeAverage();
	double standardDeviation = 0;
	double sum = 0;
	for( i=0; i<RunTime; i++ )
	{
		sum = sum + ( totalcandidate[i]-average )*( totalcandidate[i]-average );
	}
	standardDeviation = (double)( sqrt( (double)sum/(double)RunTime ) );
	
	return standardDeviation;
}

//***************  22  计算RunTime次运行的平均运行时间**************************
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

//***************  23  计算RunTime次运行，得到的最大解  **************************
int computeMaxCandidate()
{
	int i;
	int max = totalcandidate[0];
	for( i=0; i<RunTime; i++ )
	{
		if( max<totalcandidate[i] )
		{
			max = totalcandidate[i];
		}
	}
	
	return max; 
}

//***************  24  计算RunTime次运行，得到的最小解  **************************
int computeMinCandidate()
{
	int i;
	int min = totalcandidate[0];
	for( i=0; i<RunTime; i++ )
	{
		if( min>totalcandidate[i] )
		{
			min = totalcandidate[i];
		}
	}
	
	return min;
}

//***************  25  计算RunTime次运行，最优解出现的次数  **************************
int computeMaxTimes()
{
	int i;
	int max = computeMaxCandidate(); 
	int times = 0;
	for( i=0; i<RunTime; i++ )
	{
		if( totalcandidate[i]==max )
		{
			times = times + 1;
		}
	}
	
	return times;
}

//**************  26  将实验结果以文本格式进行存储  *******************
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
	if (!outfile_tm)
	{
		cerr << "open error!" << endl;
		exit(1);
	}
	outfile_tm << "the results of 2-improvement after the greedy algorithm:  "<<endl<<endl;
	outfile_tm << "Network name: " <<filename<<", ";	
	outfile_tm << "Network vertices = " << N <<endl;
	//outfile_tm << "Iteration times = " << T << endl<<endl;
	outfile_tm << "the total run times is : "<< RunTime << endl << endl;
	//outfile_tm << "the number of iteration is : "<< T << endl;
	
	outfile_tm<<"averageCandidate: "<<averageCandidate<<endl;
	outfile_tm<<"stdCandidate: "<<stdCandidate<<endl;
	outfile_tm<<"maxCandidate: "<<maxCandidate<<endl;
	outfile_tm<<"minCandidate: "<<minCandidate<<endl;
	outfile_tm<<"maxCandidateNumber: "<<maxCandidateNumber<<endl;
	outfile_tm<<"averageTime: "<<averageTime<<endl<<endl;
	
	outfile_tm << "total candidate of RunTime: "<< " ";
	for( i=0; i<RunTime; i++ )
	{
		outfile_tm << totalcandidate[i] << " ";
	}
	outfile_tm << endl;
	
	outfile_tm << "the maximum degree of the network is : "<< maxdegree << endl;
	outfile_tm << "the minimum degree of the network is : "<< mindegree << endl;
	outfile_tm << "the average degree of the network is : "<< avedegree << endl;
	outfile_tm << "the number of edges of the network is : "<< edgeNumber << endl<<endl;
	
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

//***************  27  结果输出  **************************
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
	
	cout<<"Network name = "<<filename<<endl;
	printf("Network vertices = %d\n\n",N);
	
	printf("1、得到的解集大小分别为： ");
	for( i=0; i<RunTime; i++ )
	{
		printf("%d ",totalcandidate[i]);
	}
	printf("\n\n");
	
	printf("2、解集大小的平均值为：%f，标准差为：%f \n\n", averageCandidate, stdCandidate);
	
	printf("3、解集大小最大值为：%d，最小值为：%d, 取得最大值的次数为：%d\n\n", maxCandidate, minCandidate, maxCandidateNumber);
	
	printf("4、平均运行时间为：%f \n\n",averageTime); 
	
	
	/*printf("5、得到的解集中的元素分别为：\n");
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
	printf("\n");
	*/
	

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

//**************  28  判断得到的最优解是否是独立集  *******************
void IndependentSetOrNot()
{
	int i;
	int j;
	int index = 1; // //如果index一直等于 1，则说明解集中的节点之间没有连边，即所得结果是独立集 
	
	vector<int> IS;
	IS.clear(); 
	for( i=0; i<N; i++ )
	{
		if( solution[i]==1 )
		{
			IS.push_back(i); // IS存储所有独立集中的节点 
		}
	}
	
	for( i=0; i<IS.size()-1; i++ )
	{
		for( j=i+1; j<IS.size(); j++ )
		{
			if( matrix[ IS[i] ][ IS[j] ]==1 ) // 独立集中的节点IS[i] 与IS[j] 有连边，则说明所得结果不是独立集
			{
				index = 0; //解集不是独立集，index=0 
			} 
		}
	}
	
	//return index; //index=1，是独立集；index=0，不是独立集 
	if( index==1 )
	{
		printf("是独立集,   ");
	}
	else if( index==0 )
	{
		printf("不是独立集,   ");
	}	
	
}

//**************  29  判断是否是极大独立集，直接求出perSolution，看自由节点个数是否为0  ****************
void maximalIndependentSetOrNot()
{
	int i;
	int j;
	int index = 3; //index=1，是极大独立集；index=0，不是极大独立集 
	
	computeTightness();
	permutationSolution(); 
	
	if( brock2Num==0 ) //是极大独立集
	{
		index = 1; 
	} 
	else
	{
		index = 0;
	}
		
	//return index; //index=1，是极大独立集；index=0，不是极大独立集
	if( index==1 )
	{
		printf("是极大独立集\n\n");
	}
	else if( index==0 )
	{
		printf("不是极大独立集\n\n");
	}
}

//*********************  4  初始化迭代过程有关的变量 *****************************
void initialVariable_Iteration()
{
	int i;
	int j;
	
	for( i=0; i<T; i++ )
	{
		for( j=0; j<N; j++ )
		{
			TIterSolution[i][j] = 0;
		}
	}
	for( i=0; i<T; i++ )
	{
		TiterSolutionCandidate[i] = 0;
	}
	
	for( i=0; i<N; i++ )
	{
		solution[i] = 0;
		perSolution[i] = 0;
		tightness[i] = 0;
		strategy_temp[i] = 0;
		payoff1[i] = 0;
		payoff0[i] = 0;
		vertexTimes[i] = 0;
	}
	
} 

//**************  32  perturbations  *************************
void perturbation()
{
	int i;
	int j; 
	
	// 首先确定要插入解集的节点的个数 k 
	int k;
	int solutionSize;
	int alfa;
	double pro;
	solutionSize = computeSolutionCandidate();	
	//srand( (unsigned)time(NULL) );
	alfa = (int)( (2*solutionSize)*rand() ) / ( RAND_MAX+1 ) + 1; //产生1到2*solutionSize之间的的随机数
	if( alfa!=1 )
	{
		k = 1; 
	} 
	else if( alfa==1 ) // k=i+1
	{
		k = 3;	
	}
	
	// 根据 k 的不同取值，选出k个将要被插入的节点,存在pickVerNum中  
	permutationSolution(); 
	
	vector<int> pickVerNum; //被选择的非解节点的标号 
	int r1;
	int r2[K];
	
	int Kver[K]; 
	vector<int> outOfKver;
	int verLongTime;
	int p;
	if( k==1 ) //如果k等于1，则随机地从非解节点中选择一个节点 
	{
		r1 = (int)( ((N-brock1Num) * rand()) / (RAND_MAX+1) ); // N-brock1Num=非解节点总数
		pickVerNum.push_back( perSolution[brock1Num-1+r1] ); // 从perSolution的第二块第三块中选择节点 
	}
	else if( k>1 ) 
	{
		for( i=0; i<K; i++ )
		{
			r2[i] = (int)( ((N-brock1Num) * rand()) / (RAND_MAX+1) ); // 随机选出K个非解集节点 
		}
		for( i=0; i<K; i++ )
		{
			Kver[i] = perSolution [brock1Num-1+r2[i] ];
		}
		// 得到K个节点以外的非解集非自由节点 ，储存在outOfKver中 
		for( i=brock1Num; i<N; i++ )
		{
			p = 0;
			for( j=0; j<K; j++ )
			{
				if( perSolution[i]!=Kver[j] )
				{
					p = p+1;
				}
			}
			if( p==K )
			{
				outOfKver.push_back( perSolution[i] );
			}
		}
		
		// 从outOfKver中选择一个不在解集中时间最久的节点 x=longestTimeVer
		int longestTime = vertexTimes[ outOfKver[0] ];
		int longestTimeVer = outOfKver[0];
		for( i=0; i<outOfKver.size(); i++ )
		{
			if( longestTime<vertexTimes[ outOfKver[i] ] )
			{
				longestTime = vertexTimes[ outOfKver[i] ];
				longestTimeVer = outOfKver[i];
			}
		}
		
		//找到距离longestTimeVer两步的所有节点 
		vector<int> N2x;
		for( i=0; i<adjList[longestTimeVer].size(); i++ )
		{
			int xNeighbor = adjList[longestTimeVer][i];
			for( j=0; j<adjList[xNeighbor].size(); j++ )
			{
				if( solution[ adjList[xNeighbor][j] ]!=1 )
				{
					N2x.push_back( adjList[xNeighbor][j] );
				}
			}
		}
		
		// 从 距离longestTimeVer两步的所有节点 和 longestTimeVer中选择k个节点，插入解集
		vector<int> xN2x;
		for( i=0; i<N2x.size(); i++ )
		{
			xN2x.push_back( N2x[i] );
		} 
		xN2x.push_back(longestTimeVer);
		 
		for( i=0; i<k; i++ )
		{
			int r3 = (int)( (xN2x.size() * rand()) / (RAND_MAX+1) );
			pickVerNum.push_back( xN2x[r3] );			
		}
	}
	
	// 将选出的k个节点插入solution中，并从solution中删除其邻居节点
	int ver;
	if( k==1 )
	{
		ver = pickVerNum[0];
		solution[ver] = 1;
		for( j=0; j<adjList[ver].size(); j++ )
		{
			if( solution[ adjList[ver][j] ]==1 )
			{
				solution[ adjList[ver][j] ] = 0;
			}
		}
	}
	else if( k>1 )
	{
		for( j=0; j<pickVerNum.size(); j++ )
		{
			ver = pickVerNum[j];
			solution[ver] = 1;
			for( j=0; j<adjList[ver].size(); j++ )
			{
				if( solution[ adjList[ver][j] ]==1 )
				{
					solution[ adjList[ver][j] ] = 0;
				}
			}
		}
	} 
	
	// 进行完k个节点的插入以及其邻居的移出后，将自由节点插入solution直至得到极大解
	permutationSolution(); 
	while( brock2Num!=0 )
	{
		i = brock1Num;
		solution[ perSolution[i] ] = 1;
		permutationSolution();
	}
}

//**************  34  记录每个节点在解集中出现过的次数  *************
void recordTimes()
{
	int i;
	int j;
	
	for( i=0; i<N; i++ )
	{
		if( solution[i]==1 )
		{
			vertexTimes[i] = vertexTimes[i]+1;
		}
	}
	
} 
 
