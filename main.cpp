#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <queue>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <array>

//using namespace std;

#define NUMBER_OF_TESTS 1

typedef std::array<short int, 16> TBoard;

const TBoard SOLUTION = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0 };

struct CAlgorithm {
	std::string(*Method)(const  TBoard &startPos);
	std::string Name;
	// Процент ненайденных решений
	double FailureRate;
	// Среднее время работыы
	double MeanTime;
	// На сколько шагов в среднем найденное решение хуже наилучшего
	double MeaтLag;
};

std::vector<CAlgorithm> Algorithms;

struct CNode{
	TBoard Board;
	short int GapPos;
	std::string Path;
	void FindGapPos();
	void FirstInit(const TBoard &startPos);
};

void CNode::FindGapPos()
{
	for( int i = 0; i < 16; i++ ) {
		if (Board[i] == 0) {
			GapPos = i;
			return;
		}
	}
	// throw exception there
	throw 1;
}

void CNode::FirstInit(const TBoard &startPos)
{
	Board = startPos;
	FindGapPos();
	Path = "";
}

#define NO_SOLUTIONS "No solutions"

TBoard start_pos;

int Distances[16][16];

void CalcDistances()
{
    for (int truePazzle = 1; truePazzle < 16; truePazzle++){
        for (int pos = 0; pos < 16; pos++){
            int i1, j1, i2, j2;
            i1 = (truePazzle-1) / 4;
            j1 = (truePazzle-1) % 4;

            i2 = pos / 4;
            j2 = pos % 4;
            int dist = abs(i1 - i2) + abs(j1 - j2);
            Distances[pos][truePazzle] = dist;
        }
    }
}

int GetManhattanDistance( const TBoard& board )
{
    int res = 0;
    for( int i=0; i < 16; i++ )
        res += Distances[i][board[i]];
    return res;
}

bool IsSolveExist( const TBoard& startPos ) 
{
	int inv = 0;
	for( int i = 0; i<16; i++ )
		if( startPos[i] )
			for( int j = 0; j<i; ++j )
				if( startPos[j] > startPos[i] )
					inv++;
	for( int i = 0; i<16; i++ )
		if( startPos[i] == 0 )
			inv += 1 + i / 4;
	if( inv % 2 ) {
		std::cout << "Solution not exist!!!";
	}
	return (inv % 2 == 0);
}

// generate next playing positions
void GetSuccessors( const CNode& curPos, std::vector<CNode>& successors )
{
    int gp = curPos.GapPos;
    if( gp > 3 ) { // can move up
        CNode newPos = curPos;
        std::swap(newPos.Board[gp], newPos.Board[gp - 4]);
        newPos.GapPos = gp - 4;
        newPos.Path += 'U';
        successors.push_back(newPos);
    }
    if( gp < 4 * 3 ) { // can move down
        CNode newPos = curPos;
        std::swap(newPos.Board[gp], newPos.Board[gp + 4]);
        newPos.GapPos = gp + 4;
        newPos.Path += 'D';
        successors.push_back(newPos);
    }
    if( gp % 4 < 3 ) { // can move right
        CNode newPos = curPos;
        std::swap(newPos.Board[gp], newPos.Board[gp + 1]);
        newPos.GapPos = gp + 1;
        newPos.Path += 'R';
        successors.push_back(newPos);
    }
    if( gp % 4 > 0 ) { //can move left
        CNode newPos = curPos;
        std::swap(newPos.Board[gp], newPos.Board[gp - 1]);
        newPos.GapPos = gp - 1;
        newPos.Path += 'L';
        successors.push_back(newPos);
    }
}

bool have_solution;

//----------------------IDA* Algorithm-----------------

std::string idaRec( CNode curNode, int curDepth, int maxDepth )
{
    static bool haveSolution;
    static std::string solve;

    if( curDepth == 0 ) {
        haveSolution = false;
        solve = NO_SOLUTIONS;
    }

    int dist = GetManhattanDistance(curNode.Board);

    if( dist == 0 ) {
        haveSolution = true;
        solve = curNode.Path;
        return solve;
    }

    if( curDepth + dist > maxDepth )
        return solve;

    std::vector<CNode> nextNodes;
    GetSuccessors(curNode, nextNodes);
    for( int i=0, maxi = nextNodes.size(); i < maxi && !haveSolution; i++ ) {
        idaRec(nextNodes[i], curDepth + 1, maxDepth);
    }

    nextNodes.clear();

    return solve;
}


std::string IdaStar(const TBoard& startPos ) 
{
	const int MAX_DEPTH = 30;
    std::string res;
	CNode startNode;
	startNode.FirstInit(startPos);
    for( int i = 0; i <= MAX_DEPTH; i++ ) {
        res = idaRec(startNode, 0, i);
        if (res != NO_SOLUTIONS)
            return res;
    }
    return NO_SOLUTIONS;
}

//----------------------A star algorithm---------------------------
// he will be here:)

TBoard GenerateStartPos() 
{
	TBoard res;
	//res = { 1, 2, 3, 4 , 5, 6, 7, 8, 9, 10, 15, 11, 13, 14, 0, 12 };
	res = {00, 01, 04, 8, 05, 15, 02, 03, 9, 13, 07, 12, 11, 10, 06, 14};
	return res;
}


void InitAlgorithms()
{
	CAlgorithm alg1;
	alg1.Name = "IDA*";
	alg1.Method = &IdaStar;

	Algorithms.push_back(alg1);
}


int main(){
	CalcDistances();
	InitAlgorithms();

	for( int testNumber = 0; testNumber < NUMBER_OF_TESTS; testNumber++ ) {
		TBoard startPos;
		std::string path;
		startPos = GenerateStartPos();
		for (auto pAlg = Algorithms.begin(); pAlg != Algorithms.end(); pAlg++) {
			path = (*pAlg).Method(startPos);
			std::cout << path << std::endl;
		}
	}
	system("PAUSE");
    return 0;
};