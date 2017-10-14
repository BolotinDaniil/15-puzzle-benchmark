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
#include <chrono>
#include <climits>
#include <functional>
#include <algorithm>
#include <utility>


//using namespace std;

#define NUMBER_OF_TESTS 10

//store playing board in a linear array;
typedef std::array<short int, 16> TBoard;

const TBoard SOLUTION = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0 };


// struct for description and statictics of using algorithms
struct CAlgorithm {
	std::string(*SolverAlgorithm)(const  TBoard &startPos);
	std::string Name;
	// share of not founds solutions
	double FailureShare;
	// mean work time in milliseconds
	double MeanTime;
	// mean distance between best solution and solution of this method
	double MeanLag;

	CAlgorithm( char *name, std::string(*solverAlgorithm)(const  TBoard &startPos) )
	{
		Name = name;
	 	SolverAlgorithm = solverAlgorithm;
		MeanLag = MeanTime = FailureShare = 0;
	}
};

// Add here you algorithm in InitAlgorithms function
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

void CNode::FirstInit( const TBoard &startPos )
{
	Board = startPos;
	FindGapPos();
	Path = "";
}

#define NO_SOLUTIONS "No solutions"

int Distances[16][16];

void CalcDistances()
{
    for( int truePazzle = 1; truePazzle < 16; truePazzle++ ) {
        for( int pos = 0; pos < 16; pos++ ){
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
	for( int i = 0; i < 16; i++ ) {
		res += Distances[i][board[i]];
	}
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
		//std::cout << "Solution not exist!!!";
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

//----------------------IDA* algorithm-----------------

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

	if( curDepth + dist > maxDepth ) {
		return solve;
	}

	std::vector<CNode> nextNodes;
	GetSuccessors(curNode, nextNodes);
	for( int i=0, maxi = nextNodes.size(); i < maxi && !haveSolution; i++ ) {
		idaRec(nextNodes[i], curDepth + 1, maxDepth);
	}

	nextNodes.clear();

	return solve;
}


std::string IdaStar( const TBoard& startPos ) 
{
	const int MAX_DEPTH = 36;
	std::string res;
	CNode startNode;
	startNode.FirstInit(startPos);
	for( int i = 0; i <= MAX_DEPTH; i++ ) {
		res = idaRec(startNode, 0, i);
		if( res != NO_SOLUTIONS ) {
			return res;
		}
	}
	return NO_SOLUTIONS;
}

//----------------------A* algorithm---------------------------
std::string AStar( const TBoard& startPos ) 
{
	int dist;
	int numVisited = 0;
	int maxFrontier = 0;
	CNode startNode;
	startNode.FirstInit(startPos);

	double magicConstant = 1; // warning: value more than 1 will make unacceptable estimation function, set it only when debug

	std::multimap<int, CNode> pqueue; //list for next nodes to expand
	dist = GetManhattanDistance(startPos);
	pqueue.insert(std::pair <int, CNode>(int(magicConstant * dist), startNode));
	std::set<TBoard> visited;
	while( !pqueue.empty() ) {
		CNode curNode = pqueue.begin()->second;
		pqueue.erase(pqueue.begin());
		numVisited += 1;
		//        if (numVisited % 1000000 == 0)
		//            std::cout << numVisited <<" : " << pqueue.size()<< std::endl;
		visited.insert(curNode.Board);

		if( curNode.Board == SOLUTION ) {
			//std::cout << "Visited: " << numVisited << std::endl;
			//std::cout << "Max Frontier: " << maxFrontier << std::endl;
			return curNode.Path;
		}

		std::vector<CNode> nextNodes;
		nextNodes.clear();
		GetSuccessors(curNode, nextNodes);
		for( int i = 0, maxi = nextNodes.size(); i < maxi; i++ ) {
			if( visited.find(nextNodes[i].Board) != visited.end() ) {
				continue;
			}
			dist = GetManhattanDistance(nextNodes[i].Board);
			int cost = int(dist * magicConstant) + int(nextNodes[i].Path.length());
			pqueue.insert(std::pair<int, CNode>(cost, nextNodes[i]));
			if( pqueue.size() > maxFrontier ) {
				maxFrontier = pqueue.size();
			}
		}
		// freeze protection
		if( maxFrontier > 1e5 ) {
			break;
		}
	}
	return NO_SOLUTIONS;
}

//----------------------Beam search algorithm START ---------------------------

std::vector<CNode> getNewFrontFromPossibleFront( std::vector<CNode>& possibleFrontNodes, const int frontCapacity )
{
	std::vector<std::pair<int, int> > possibleFrontDistances;
	std::vector<CNode> newFront;
	for (int possibleFrontIndex = 0; possibleFrontIndex < possibleFrontNodes.size(); ++possibleFrontIndex) {
		possibleFrontDistances.emplace_back(GetManhattanDistance(possibleFrontNodes[possibleFrontIndex].Board), possibleFrontIndex);
	}

	std::sort(possibleFrontDistances.begin(), possibleFrontDistances.end());

	for (int possibleFrontIndex = 0; possibleFrontIndex < std::min((int)possibleFrontNodes.size(), frontCapacity); ++possibleFrontIndex) {
		newFront.push_back(possibleFrontNodes[possibleFrontDistances[possibleFrontIndex].second]);
	}

	return newFront;
}

std::string BeamSearch( const TBoard& startPos, const int frontCapacity ) 
{
	int stepCount = 0;
	const int STEP_LIMIT = 1e2;

	CNode startNode;
	startNode.FirstInit(startPos);

	std::vector<CNode> frontNodes;
	frontNodes.push_back(startNode);

	while( stepCount < STEP_LIMIT && frontNodes.size() != 0) {
		++stepCount;

		for (int index = 0; index < frontNodes.size(); ++index) {
			if (GetManhattanDistance(frontNodes[index].Board) == 0) {
				return frontNodes[index].Path;
			}
		}

		std::vector<CNode> possibleFrontNodes;
		for (int frontIndex = 0; frontIndex < frontNodes.size(); ++frontIndex) {
			std::vector<CNode> children;
			GetSuccessors(frontNodes[frontIndex], children);
			for (int childIndex = 0; childIndex < children.size(); ++childIndex) {
				possibleFrontNodes.push_back(children[childIndex]);
			}
		}

		if (possibleFrontNodes.size() == 0) {
			return NO_SOLUTIONS;
		}
		
		frontNodes = getNewFrontFromPossibleFront(possibleFrontNodes, frontCapacity);
	}

	return NO_SOLUTIONS;
}

std::string BeamSearchFrontCap1( const TBoard& startPos ) {
	return BeamSearch(startPos, 1);
}

std::string BeamSearchFrontCap128( const TBoard& startPos ) {
	return BeamSearch(startPos, 128);
}

std::string BeamSearchFrontCap1024( const TBoard& startPos ) {
	return BeamSearch(startPos, 1024);
}

std::string BeamSearchFrontCap262144( const TBoard& startPos ) {
	return BeamSearch(startPos, 262144);
}

//----------------------Beam search algorithm END ---------------------------

TBoard GenerateStartPos() 
{
	// generate new random start position
	// this function have to generate only correct start positions (use func IsSolveExist( const TBoard& startPos ) )
	TBoard res;
	bool ok = false;
	while( !ok ) {
		std::array<bool, 16> was;
		for( int i = 0; i < 16; i++ ) {
			was[i] = false;
		}
		int k;
		for( int i = 0; i < 16; i++ ) {
			int k;
			while( was[k = rand() % 16] );
			was[k] = true;
			res[i] = k;
		}
		if( IsSolveExist(res) ) {
			ok = true;
		}
	}

	//res = { 1, 2, 3, 4 , 5, 6, 7, 8, 9, 10, 15, 11, 13, 14, 0, 12 };
	//res = {00, 01, 04, 8, 05, 15, 02, 03, 9, 13, 07, 12, 11, 10, 06, 14};
	return res;
}

void InitAlgorithms()
{
	CAlgorithm alg1("IDA*", &IdaStar);
	Algorithms.push_back(alg1);

	CAlgorithm alg2("A*", &AStar);
	Algorithms.push_back(alg2);

	// Не находит ничего
	//CAlgorithm beamSearch1024("Beam search front capacity: 1024", &BeamSearchFrontCap1024);
	//Algorithms.push_back(beamSearch1024);

	// Работает 70 сек на одном тесте
	CAlgorithm beamSearch262144("Beam search front capacity: 262144", &BeamSearchFrontCap262144);
	Algorithms.push_back(beamSearch262144);
}

void RunBenchmark() 
{
	for( int testNumber = 0; testNumber < NUMBER_OF_TESTS; testNumber++ ) {
		TBoard startPos;
		std::string path;
		startPos = GenerateStartPos();
		int minPathLen = INT_MAX;
		std::vector<std::string> paths;
		paths.clear();
		for( auto pAlg = Algorithms.begin(); pAlg != Algorithms.end(); pAlg++ ) {
			std::chrono::high_resolution_clock::time_point startTime, finTime;
			startTime = std::chrono::high_resolution_clock::now();
			path = pAlg->SolverAlgorithm(startPos);
			finTime = std::chrono::high_resolution_clock::now();

			pAlg->MeanTime += (std::chrono::duration_cast<std::chrono::nanoseconds> (finTime - startTime)).count();

			//std::cout << path << std::endl;
			paths.push_back(path);
			if( path == NO_SOLUTIONS ) {
				pAlg->FailureShare += 1;
			}
			else if( minPathLen > path.length() ) {
				minPathLen = path.length();
			}
		}
		for( int i = 0; i < Algorithms.size(); i++ ) {
			if( paths[i] != NO_SOLUTIONS) {
				Algorithms[i].MeanLag += paths[i].length() - minPathLen;
			}
		}
	}
	for (auto pAlg = Algorithms.begin(); pAlg != Algorithms.end(); pAlg++) {
		// calc mean and convert nanoseconds in milliseconds
		pAlg->MeanTime /= NUMBER_OF_TESTS * 1e6;
		if( pAlg->FailureShare < NUMBER_OF_TESTS ) {
			pAlg->MeanLag /= (NUMBER_OF_TESTS - pAlg->FailureShare);
		}
		pAlg->FailureShare /= NUMBER_OF_TESTS;
	}
}

void OutAlgorithmsStatistics()
{
	for( auto pAlg = Algorithms.begin(); pAlg != Algorithms.end(); pAlg++ ) {
		std::cout << "---------------------------------------------" << std::endl;
		std::cout << "Name: " << pAlg->Name << std::endl;
		std::cout << "Mean Time: " << pAlg->MeanTime << " ms" << std::endl;
		std::cout << "Failures Rate: " << pAlg->FailureShare << std::endl;
		std::cout << "Mean distance to best solve: " << pAlg->MeanLag << std::endl;
	}
}

int main()
{
	// automatic randomization
	// srand(time(0));
	// Results should be reproducible
	srand(42);

	CalcDistances();
	InitAlgorithms();
	RunBenchmark();
	OutAlgorithmsStatistics();

	system("PAUSE");
	return 0;
}
