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

#define N 2000  //������
string filename = "BA2000(8).txt";
#define swapTimes 10000
#define b 1.5 // PDG���� 
#define RunTime 30 // ͳ��RunTime�ε����н������ƽ�� 
#define T 1000 //������������ILS����ֹ����
#define K 4
//********************************************  ��������  ********************************************* 
vector< vector<int> > adjList; //�����ڽӱ��洢�ڵ���ھӽڵ�ı�� 
vector<int> adjList_line;

vector<int> candidate; //�洢���п��ܴ���2-improvement�Ľڵ㣬�洢�ڵ�ı�� 

int matrix[N][N]; //  1  ������ڽӾ���
array< array<int, N>, N > matrix_temp;
int deg[N]; //  3  ������ÿ���ڵ�Ķ�
int maxdegree;
int mindegree;
double avedegree;
int edgeNumber;

array<int, N> solution; // �⼯ 
array<int, N> solutionBefore; // �⼯ 
array<int, N> solutionAfter; // �⼯ 
int perSolution[N]; //  5  ������ṹ��ʾ�Ľ⣬����ṹ�洢�ڵ�ı�ţ���0��ʼ��N-1 

int tightness[N]; //  7  ��¼ÿһ���ڵ��tightness 

int brock1Num; //  8  ��һ��ռ��еĽڵ����������������Ԫ�ظ���
int brock2Num; //  9  �ڶ���ռ��еĽڵ�����������ɽڵ�Ľڵ���� 
int brock3Num; //  10  ������ռ��еĽڵ���������Ȳ��Ƕ�����Ҳ�������ɽڵ�Ľڵ���� 

int strategy[N]; // N���ڵ�Ĳ��ԣ�0�����ѣ�1�������
double payoff1[N]; 
double payoff0[N];
int strategy_temp[N];
int set[N]; 

array< array<int, N>, RunTime > solu; // RunTime�����У�ÿһ�εõ��Ľ⼯
int totalcandidate[RunTime]; // RunTime�����У�ÿһ�εõ��Ľ⼯��Ԫ�ظ��������õ��Ķ�������Ԫ�ظ���
double totalTime[RunTime]; // RunTime�����У�ÿһ�ε�����ʱ��
double averageCandidate;// RunTime�����У���������С��ƽ��ֵ 
double stdCandidate;// RunTime�����У���������С�ı�׼��
int maxCandidate;// RunTime��������candidate�����ֵ 
int minCandidate;// RunTime��������candidate����Сֵ 
int maxCandidateNumber;// RunTime�������еõ����candidate�Ĵ���
double averageTime;// RunTime�����е�ƽ������ʱ��

int BeforeSize;
int AfterSize;

int S;

array< array<int, N>, T > TIterSolution; // T�ε�����ÿһ�εõ��õ��Ľ⼯
array< int, T > TiterSolutionCandidate; // T�ε�������¼ÿһ�ε���ʱ���⼯�Ĵ�С
array< int, N > vertexTimes;


//******************************************  �Ӻ�������  *********************************************
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
void updateCandidateInsert( int insertNum ); // 13  insertNumΪ����ڵ�ı��
void updateCandidateRemove( int removeNum ); // 14  removeNumΪ�Ƴ��ڵ�ı�� 
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
	int isNumBefore; // 2-improvement LS ǰ���⼯��Ԫ�ظ��� 
	int isNumAfter; // 2-improvement LS �󣬽⼯��Ԫ�ظ��� 
	int solutionCandidate = 0; // �⼯��Ԫ�ظ��� 
	
	clock_t start1,end1,start;
	clock_t start2,end2,end; 
	double time1; //��¼������ʼ����㷨������ʱ�� 
	double time2; //��¼2-improvement������ʱ�� 
	double runningTime;
	double perturb_time;
	srand( (unsigned int)time(NULL) );
	
	initial(); //������ʼ�� 
	initial_RunTime();
	
	//read();	//��ȡ������Ϣ 	
	//complement_graph();
	//read_Without_e();
	readMatrix();
	
	obtainAdjList(); //�õ�������ڽӱ� 	
	obtainVerDegree(); //�õ�ÿ���ڵ�Ķ� 
			 
	//**************************  2-improvement local search  **************************
	for( run=0; run<RunTime; run++ ) //RunTime�ζ������� 
	{
		S = 0;
		initialVariable_Iteration();
		
		//randomAlgorithm(); //��random algorithm����һ����ʼ�� 
		greedyAlgorithm(); //��greedy algorithm����һ����ʼ��
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
			
			computeTightness(); //����ڵ��tightness��ֻ�зǶ������Ľڵ����tightness 		
			permutationSolution(); // ������Ľṹ��ʾ�����������ɽڵ�ͷǶ������ҷ����ɽڵ�Ľڵ� 	
			initialCandidate(); // ������һ��2-improvement �Ľ��solution����ʼ��һ�ε�candidate�����½�����һ��2-improvement 
			simpleCandidate();
			
			twoImprovement(); // 2-improvement�󣬵õ��µ�solution����solution�洢��TIterSolution[t][]�� 
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
		//time2 = time1 + time2; //��ʼ�㷨+2-improvement��������ʱ�� 
		
		
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

//***************  0  ��ȡ�ļ�  **************************
void read() // �õ��������ڽӾ��� �� ����¼�˱ߵľ���edgeMatrix[M][2] 
{
	int i;
	int j;
	FILE* fp = fopen(filename.c_str(),"r"); 
	//FILE* fp = fopen("brock200_2 (12).txt","r");
	
	if( !fp )
	{
		printf("�ļ��򿪴���");
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

//***************  1  ��ȡ�ļ�������ıߵı�ʾ  ����e��  **************************
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
		printf("�ļ��򿪴���");
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

//***************  3  ��ȡ�ڽӾ����ʾ������  **************************
void readMatrix()
{
	int i;
	int j;
	
	FILE* fp;
	fp = fopen(filename.c_str(),"r"); 
	//fp = fopen( "BA500(10).txt" , "r" );
	if ( fp==NULL )
	{
		printf("�ļ���ʧ��");
		exit( 1 );
	}
	char C;
	int row=0;
	int col=0;
	int i1 = 0;//����matrix���±���Ҫÿ�ε���1����������������������������i1��j1
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

//***************  6  ������Ĳ�ͼ  **************************
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

//****************  1  ������ʼ��  ****************************
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
		solution[i] = 3; //�Խ����ʾ���г�ʼ�� 
		perSolution[i] = 0; //����ṹ�洢�ڵ�ı�� 
		tightness[i] = 0; //ÿ���ڵ��tightness��ʼ��Ϊ0 
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

//****************  2  ������ʼ��  ****************************
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

//****************  3  �õ�������ڽӱ�(ֻ��¼�ھ�  ����ڵ��Ŵ�0��ʼ)  *************
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
				adjList[i].push_back(j); //  �õ����ڽӱ��У�ÿһ���ڵ���ھӽڵ㶼�Ա�ŵ��������� 
			}
		}
	}

}

//**************  4  �õ�ÿ���ڵ�Ķ�  ****************
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

//**************  5  random algorithm ������ʼ��  *************************
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
			temp.push_back( i ); // �����ɽڵ�ı�ŷ���temp 
		}
	}
	
	while( (temp.size()) != 0 )
	{
		/*temp.clear();
		
		for( i=0; i<N; i++ )
		{
			if( tightness[i]==0 )
			{
				temp.push_back( i ); // �����ɽڵ�ı�ŷ���temp 
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
				temp.push_back( i ); // �����ɽڵ�ı�ŷ���temp 
			}
		}
	}
		
} 

//**************  6  greedy algorithm ������ʼ��  *************************
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
		if( tightness[i]==0 )  //��¼���ɽڵ�ı�� 
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
			if( tightness[i]==0 )  // ��Ҫ�ҵ�������ɽڵ���ھ��У����ɽڵ�ĸ��� 
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
				minVer = freeVertexCost[0][i]; //minVer����Ҫ����⼯�Ľڵ� 
			}		 
		}
	
		solution[minVer] = 1;
	
		computeTightness(); 
	}
	
}

//**************  7  ����ڵ��tightness  **********************
void computeTightness()
{
	int i;
	int j;
	int k;
	int s;
	for( i=0; i<N; i++ )
	{
		if( solution[i]==0 ) //ֻ�ж���������Ľڵ����tightness
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
			tightness[i] = N; //�������еĽڵ�û��tightness ����ʱ�����������еĽڵ��tightness���ó�һ���ǳ������N 
		}
	}
}

//**************  8 ��solution��Ϊ���飬������ÿһ���е�Ԫ�ظ���  ***********************
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
			brock1_i = brock1_i + 1; //brock1_i �ǵ�һ��ռ��еĽڵ����������������Ԫ�ظ��� 
		}
	}
	brock1Num = brock1_i;
	
	brock2_i = brock1_i;
	for( i=0; i<N; i++ )
	{
		if( tightness[i]==0 ) //�����ɽڵ�
		{
			perSolution[brock2_i] = i;
			brock2_i = brock2_i + 1; //brock2_i-brock1_i �ǵڶ���ռ��еĽڵ�����������ɽڵ�Ľڵ���� 
		} 
	}
	brock2Num = brock2_i-brock1_i;
	
	brock3_i = brock2_i;
	for( i=0; i<N; i++ )
	{
		if( (tightness[i]>0)&&(tightness[i]!=N) )
		{
			perSolution[brock3_i] = i;
			brock3_i = brock3_i + 1; //brock3_i-brock2_i �ǵ�����ռ��еĽڵ�����������Ƕ�����Ҳ�������ɽڵ� 
		}
	}
	brock3Num = brock3_i-brock2_i;
	
} 

//**************  9  ��ʼ��candidate����ʼ�����м���������еĽڵ㶼��candidate��  *********************
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

//**************  10  ��candidate����candidate��ɾ�������ܴ���2-improvement�Ľڵ�  *********************
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
	//solution��û�������໥������1-tight�ھӽڵ�Ľڵ㣬������2-improvement�����Դ�candidate��ɾ�� 
	for( i=0; i<candidate.size(); i++ ) //���candidate�е�ÿһ���ڵ㣬�Ƿ���Ա��Ƴ�candidate
	{
		ver = candidate[i];
		
		tight1Neighbor.clear();
		//record.clear();
		
		for( j=0; j<deg[ ver ]; j++ ) //����candidate[i]�������ھ� 
		{
			if( tightness[ adjList[ver][j] ]==1 ) //ver��1-tight�ھӽڵ� 
			{
				tight1Neighbor.push_back( adjList[ver][j] ); //�ҵ�candidate[i]������1-tight�ھ� 
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
						indexS = indexS+1;  // indexS��ָʾ����1-tight�ھ��У������ھ�������ڣ����1
									    // ���������Щ�ھӶ����ڣ���indexS�͵�����Щ�ڵ�ȫ��ͨʱ�ı��� 
					}
				}
			}
	
			if( indexS==( ( (tight1Neighbor.size())*(tight1Neighbor.size()-1) )/2 )  ) //�ܱȽϴ���   ȫ��ͨ�ı��� 
			{
			//�Ƚ϶��ٴΣ�indexS���Ƕ��٣���˵�������໥�����ߣ���candidate[i]û���໥������1-tight�ھӣ������Ƴ�candidate
			//record.push_back(i); //��¼Ҫ�Ƴ�candidate�Ľڵ���candidate�е�λ��
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
			if( record[i]==candidate[j] ) //candidate[j]Ҫ��ɾ��
			{
				candidate[j] = N; 
			} 
		}
	}
	
	for( i=0; i<candidate.size(); i++ )
	{
		if( candidate[i]!=N ) //tempCandidate�洢Ҫ����������candidateԪ�� 
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

//**************  11  ����ĳ���ڵ�ʱ�����ýڵ��candidate���Ƴ������ظ������ڵ�  ********************
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

//**************  12  �нڵ����solutionʱ������candidates  *********************
void updateCandidateInsert( int insertNum ) //insertNumΪ����ڵ�ı�� 
{
	int i;
	int j;
	
	candidate.push_back( insertNum );
	
}

//**************  13  �нڵ��Ƴ�solutionʱ������candidates  *********************
void updateCandidateRemove( int removeNum ) //removeNumΪ�Ƴ��ڵ�ı�� 
{
	int i;
	int j;
	
	solution[removeNum] = 0; //removeNum���뵽�������� 
	permutationSolution();
	
	//Ѱ����ΪremoveNum���Ƴ���removeNum���ھӽڵ��б�Ϊ1-tight���ھӽڵ㣬
	//��Щ�ھӽڵ���solution�е��ڽӽڵ�Ҫ����candidate 
	computeTightness(); //�Ƴ�һ���ڵ�󣬸���tightness
	
	vector<int> reTight1;
	for( i=0; i<deg[removeNum]; i++ )
	{
		if( tightness[ adjList[removeNum][i] ]==1 )
		{
			reTight1.push_back( adjList[removeNum][i] ); //�ҵ���ΪremoveNum���Ƴ���removeNum���ھӽڵ��б�Ϊ1-tight���ھӽڵ�
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

//**************  14  �ӽ⼯���Ƴ�һ���ڵ�  *********************
void removeVerSolution( int removeNum )
{
	solution[ removeNum ] = 0; 
	
}

//**************  15  �⼯�в���һ���½ڵ�  *********************
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
	
	//for( i=0; i<candidate.size(); i++ )  //����candidate�е����нڵ㣬��ִ��һ��2-improvement 
	//if( candidate.size()>0 )
	while( candidate.size()>0 )
	{
		S = S+1;
		removeVertex = candidate[0]; //ÿ�μ��candidate�ĵ�һ��Ԫ�أ�ֱ��candidateû��Ԫ�أ��������� 
		//removeVertex = candidate[i];
		
		removeTestedVertexFromCandidate( removeVertex ); // removeVertex�������ˣ���removeVertex��candidate���Ƴ� 
		
		removeVerSolution( removeVertex ); // solution�ı��ˣ�ֻҪsolution�ı䣬��Ҫ����candidate 
										   // ��solution���Ƴ�removeVertex��solution��S��ΪS' 
		// removeVertex����ʱ���ߣ����ǲ���ζ��removeVertex�����Ϊ�Ƕ������еĽڵ㣿����
		 
		permutationSolution(); //solution�ı䣬��solution��ʾΪ����ĽṹҲҪ��Ӧ�ı�
		computeTightness();
		
		updateCandidateRemove( removeVertex ); // �������ڵĽ�S'������candidate�����º��ٽ��м� 
		simpleCandidate(); 
		 
		int index1 = 0; //ָʾ�� 
		int index2 = 0;
		
		if( brock2Num>=2 ) //��S���Ƴ�һ���ڵ�x�󣬵õ� S',���S'�����ɽڵ������ڵ���2����������²���������˵��x���Ƴ����ܲ���2-improvement 
		{
			for( j=brock1Num; j<brock1Num+brock2Num; j++ ) //��S'�����ɽڵ����ҵ�x=removeVertex���ھӽڵ� 
			{
				tempVertex1 = perSolution[j]; // jҪȡ��perSolution���м�����ɽڵ���һ��
											  // ��һ�����ʼ�±���brock1Num,��ֹ�±���brock1Num+brock2Num-1 
				if( matrix[removeVertex][tempVertex1]==1 )
				{
					index1 = index1 + 1; //��¼�ҵ���removeVertex����S�������ɽڵ���ھӵĸ���
					 
					insertVerSolution( tempVertex1 ); //��S'�����ɽڵ���x���ھӽڵ����S'���õ�S'' 
					permutationSolution();
					computeTightness();
					
					updateCandidateInsert( tempVertex1 ); //�������ڵĽ�S��������candidate�����º��ٽ��м� 
					simpleCandidate();
					 
					if( brock2Num>=1 ) //���S''��һ�����ɽڵ㣬�򽫸����ɽڵ����S''�����һ��2-improvement 
					{
						tempVertex2 = perSolution[brock1Num]; //ȡ���ڶ���ĵ�һ��Ԫ�أ�Ҳ�����һ��Ԫ�� 
						
						insertVerSolution( tempVertex2 );
						permutationSolution(); 
						computeTightness();
						
						updateCandidateInsert( tempVertex2 );
						simpleCandidate();
						
						j = N; //x���һ��2-improvement���������ھӽڵ��ѭ�� 
					}
					else
					{
						removeVerSolution( tempVertex1 ); //���S''���ɽڵ�������1����0������
														  //��x��S''�Ƴ������µõ�S'�������x��S'�����ɽڵ��е���һ���ھӽڵ� 
						permutationSolution(); 
						computeTightness();
						
						updateCandidateRemove( tempVertex1 );
						simpleCandidate();	
						
						index2 = index2 + 1; // ��ǰremoveVertex����2-improvementʱ��
											 // ��������S'�����ɽڵ���ھӣ����Ǹ��ھӲ��ܵõ�2-improvementʱ��index2��1 					
					}
				}
			}
			if( index1==index2 ) // removeVertex����S�������ɽڵ���ھӵĸ��� = ��removeVertexִ��2-improvementʧ�ܵĴ����� 
			{
				insertVerSolution( removeVertex ); 
				permutationSolution();
				computeTightness(); 
						
				updateCandidateInsert( removeVertex );
				simpleCandidate();
			} 
		}  
		else // ��removeVertex���·Ż�solution�� 
		{
			insertVerSolution( removeVertex ); 
			permutationSolution();
			computeTightness(); 
						
			updateCandidateInsert( removeVertex );
			simpleCandidate();
		}
	}
	
	 // ����������ɽڵ㣬������ɽڵ����⼯
	 permutationSolution();
	 while( brock2Num!=0 )
	 {
	 	solution[ perSolution[brock1Num] ] = 1;
	 	permutationSolution();
	 } 
	  
}

//**************  17  PDG game �����ʼ�����нڵ�Ĳ���  ********************
void initialPDGStrategy()
{
	int i;
	int j;
	
	//�����ʼ�����нڵ�Ĳ���
	for( i=0; i<N; i++ )
	{
		strategy[i] = 0;
	} 
	//srand( (unsigned int)time(NULL) );
	for( i=0; i<N; i++ )
	{
		strategy[i] = (int)( 100.0*rand() ) / ( RAND_MAX+1 ); //����0��99������� 
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

//**************  18  PDG game ���Ĺ���  �������е�strategy���в���  *************
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

//**************  19  ����solution��Ԫ�ظ���  *************************
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

//***************  20  ����RunTime�����еĽ⼯��Сƽ��ֵ**************************
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

//***************  21 ����RunTime�����еĽ⼯��С�ı�׼��**************************
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

//***************  22  ����RunTime�����е�ƽ������ʱ��**************************
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

//***************  23  ����RunTime�����У��õ�������  **************************
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

//***************  24  ����RunTime�����У��õ�����С��  **************************
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

//***************  25  ����RunTime�����У����Ž���ֵĴ���  **************************
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

//**************  26  ��ʵ�������ı���ʽ���д洢  *******************
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

//***************  27  ������  **************************
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
	
	printf("*********************  ������ %d ��  *************************��\n\n",RunTime);
	
	cout<<"Network name = "<<filename<<endl;
	printf("Network vertices = %d\n\n",N);
	
	printf("1���õ��Ľ⼯��С�ֱ�Ϊ�� ");
	for( i=0; i<RunTime; i++ )
	{
		printf("%d ",totalcandidate[i]);
	}
	printf("\n\n");
	
	printf("2���⼯��С��ƽ��ֵΪ��%f����׼��Ϊ��%f \n\n", averageCandidate, stdCandidate);
	
	printf("3���⼯��С���ֵΪ��%d����СֵΪ��%d, ȡ�����ֵ�Ĵ���Ϊ��%d\n\n", maxCandidate, minCandidate, maxCandidateNumber);
	
	printf("4��ƽ������ʱ��Ϊ��%f \n\n",averageTime); 
	
	
	/*printf("5���õ��Ľ⼯�е�Ԫ�طֱ�Ϊ��\n");
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

//**************  28  �жϵõ������Ž��Ƿ��Ƕ�����  *******************
void IndependentSetOrNot()
{
	int i;
	int j;
	int index = 1; // //���indexһֱ���� 1����˵���⼯�еĽڵ�֮��û�����ߣ������ý���Ƕ����� 
	
	vector<int> IS;
	IS.clear(); 
	for( i=0; i<N; i++ )
	{
		if( solution[i]==1 )
		{
			IS.push_back(i); // IS�洢���ж������еĽڵ� 
		}
	}
	
	for( i=0; i<IS.size()-1; i++ )
	{
		for( j=i+1; j<IS.size(); j++ )
		{
			if( matrix[ IS[i] ][ IS[j] ]==1 ) // �������еĽڵ�IS[i] ��IS[j] �����ߣ���˵�����ý�����Ƕ�����
			{
				index = 0; //�⼯���Ƕ�������index=0 
			} 
		}
	}
	
	//return index; //index=1���Ƕ�������index=0�����Ƕ����� 
	if( index==1 )
	{
		printf("�Ƕ�����,   ");
	}
	else if( index==0 )
	{
		printf("���Ƕ�����,   ");
	}	
	
}

//**************  29  �ж��Ƿ��Ǽ����������ֱ�����perSolution�������ɽڵ�����Ƿ�Ϊ0  ****************
void maximalIndependentSetOrNot()
{
	int i;
	int j;
	int index = 3; //index=1���Ǽ����������index=0�����Ǽ�������� 
	
	computeTightness();
	permutationSolution(); 
	
	if( brock2Num==0 ) //�Ǽ��������
	{
		index = 1; 
	} 
	else
	{
		index = 0;
	}
		
	//return index; //index=1���Ǽ����������index=0�����Ǽ��������
	if( index==1 )
	{
		printf("�Ǽ��������\n\n");
	}
	else if( index==0 )
	{
		printf("���Ǽ��������\n\n");
	}
}

//*********************  4  ��ʼ�����������йصı��� *****************************
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
	
	// ����ȷ��Ҫ����⼯�Ľڵ�ĸ��� k 
	int k;
	int solutionSize;
	int alfa;
	double pro;
	solutionSize = computeSolutionCandidate();	
	//srand( (unsigned)time(NULL) );
	alfa = (int)( (2*solutionSize)*rand() ) / ( RAND_MAX+1 ) + 1; //����1��2*solutionSize֮��ĵ������
	if( alfa!=1 )
	{
		k = 1; 
	} 
	else if( alfa==1 ) // k=i+1
	{
		k = 3;	
	}
	
	// ���� k �Ĳ�ͬȡֵ��ѡ��k����Ҫ������Ľڵ�,����pickVerNum��  
	permutationSolution(); 
	
	vector<int> pickVerNum; //��ѡ��ķǽ�ڵ�ı�� 
	int r1;
	int r2[K];
	
	int Kver[K]; 
	vector<int> outOfKver;
	int verLongTime;
	int p;
	if( k==1 ) //���k����1��������شӷǽ�ڵ���ѡ��һ���ڵ� 
	{
		r1 = (int)( ((N-brock1Num) * rand()) / (RAND_MAX+1) ); // N-brock1Num=�ǽ�ڵ�����
		pickVerNum.push_back( perSolution[brock1Num-1+r1] ); // ��perSolution�ĵڶ����������ѡ��ڵ� 
	}
	else if( k>1 ) 
	{
		for( i=0; i<K; i++ )
		{
			r2[i] = (int)( ((N-brock1Num) * rand()) / (RAND_MAX+1) ); // ���ѡ��K���ǽ⼯�ڵ� 
		}
		for( i=0; i<K; i++ )
		{
			Kver[i] = perSolution [brock1Num-1+r2[i] ];
		}
		// �õ�K���ڵ�����ķǽ⼯�����ɽڵ� ��������outOfKver�� 
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
		
		// ��outOfKver��ѡ��һ�����ڽ⼯��ʱ����õĽڵ� x=longestTimeVer
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
		
		//�ҵ�����longestTimeVer���������нڵ� 
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
		
		// �� ����longestTimeVer���������нڵ� �� longestTimeVer��ѡ��k���ڵ㣬����⼯
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
	
	// ��ѡ����k���ڵ����solution�У�����solution��ɾ�����ھӽڵ�
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
	
	// ������k���ڵ�Ĳ����Լ����ھӵ��Ƴ��󣬽����ɽڵ����solutionֱ���õ������
	permutationSolution(); 
	while( brock2Num!=0 )
	{
		i = brock1Num;
		solution[ perSolution[i] ] = 1;
		permutationSolution();
	}
}

//**************  34  ��¼ÿ���ڵ��ڽ⼯�г��ֹ��Ĵ���  *************
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
 
