#include<iostream>
#include<vector>
#include<array> 
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string>

using namespace std;

#define N 6
#define E 10
#define cutofftime 3600
#define RunTime 10

string filename = "p_hat1500-2 (65).txt";
  
array< array<int,N>, N > matrix;//邻接矩阵 
array< array<int, N>, N > matrix_temp; 
vector< vector<int> > adjList; 
vector<int> adjList_line;
int deg[N]; 
int maxdegree; 
int mindegree; 
double avedegree;  
int edgenumber; 

array< int, N > global_solution;
array<int,N> solution;
int weight[E];
int dscore[N];


array<int,RunTime> runSolutionSize;
array< array<int, N>, RunTime > runSolution;
int maxSize;
int maxNumber;
int minSize;
int minNumber;
double aveSize;
double aveTime;
double std;

//********************************************  函数声明  **********************************************
void read();
void readMatrix();
void complement_graph();
void obtainAdjList();
void obtainVerDegree();
 

//********************************************  main  **********************************************
int main(int argc, char* argv[]){
	int i;
	int j;
	
	globalInitialVariable();
	read();
	//readMatrix();
	//complement_graph();
	obtainAdjList();
	obtainVerDegree(); 
	
	return 0;
}

void globalInitialVariable(){
	int i;
	int j;
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			matrix[i][j] = 0;
			matrix_temp[i][j] = 0;
		}
	}
	for( i=0; i<N; i++ ){
		adjList.push_back(adjList_line);
		deg[i] = 0;	
	}
	maxdegree = 0;
	mindegree = 0;
	avedegree = 0;
	edgenumber = 0;
	for(i=0; i<RunTime; i++){
		runSolutionSize[i] = 0;
		for(j=0; j<N; j++){
			runSolution[i][j] = 0;
		}
	}
	maxSize = 0;
	maxNumber = 0;
	minSize = 0;
	minNumber = 0;
	aveSize = 0;
	aveTime = 0;
	std = 0;
}

//********************************************  子函数定义  **********************************************
//***************  1  读取文件（网络的边的表示）  **************************
void read() // 得到了网络邻接矩阵 ， 并记录了边的矩阵edgeMatrix[M][2] 
{
	int i;
	int j;
	FILE* fp = fopen(filename.c_str(),"r"); 	
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
		i = i+1;
	}

	fclose( fp );
}

//***************  1  读取文件（网络的边的表示）  **************************
void read_without_e() // 得到了网络邻接矩阵 ， 并记录了边的矩阵edgeMatrix[M][2] 
{
	int i;
	int j;
	FILE* fp = fopen(filename.c_str(),"r"); 
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
		fscanf(fp, "%d %d \n", &startnode, &endnode );
		matrix[startnode-1][endnode-1] = 1;
		matrix[endnode-1][startnode-1] = 1;
		i = i+1;
	}

	fclose( fp );
}

//***************  2  读取邻接矩阵表示的网络  **************************
void readMatrix()
{
	int i;
	int j;
	FILE* fp = fopen(filename.c_str(),"r"); 
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

//***************************  5  求网络的补图  *****************************
void complement_graph()
{
	int i;
	int j;	
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

//************************************  6  得到网络邻接表  ************************************
void obtainAdjList()
{
	int i;
	int j;
	
	for( i=0; i<N; i++ )
	{
		adjList[i].clear();
	}
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

//********************************  7  得到网络的度信息  *****************************
void obtainVerDegree()
{
	int i;

	for( i=0; i<N; i++ )
	{
		deg[i] = adjList[i].size();
	}
	
	int max = deg[0];
	for( i=0; i<N; i++ )
	{
		if( max<deg[i] )
		{
			max = deg[i];
		}
	} 
	maxdegree = max; // 最大度 
	
	int min = deg[0];
	for( i=0; i<N; i++ )
	{
		if( min>deg[i] )
		{
			min = deg[i];
		}
	}
	mindegree = min; // 最小度  
	
	double sumdegree = 0;
	for( i=0; i<N; i++ )
	{
		sumdegree = sumdegree + deg[i];
	}
	avedegree = (double)( (double)sumdegree/ (double)N); // 平均度 
	
	edgenumber = (int)( sumdegree/ 2); // 总边数 
}

