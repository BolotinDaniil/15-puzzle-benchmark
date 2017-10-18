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

#define NUMBER_OF_TESTS 5

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
	const int MAX_DEPTH = 40;
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
		if( maxFrontier > 1e6 ) {
			break;
		}
	}
	return NO_SOLUTIONS;
}

//----------------------A* with memory limit algorithm---------------------------

CNode popNode( const int MAX_DEPTH, std::array <std::multimap<int, CNode>, 100> &front )
{
	int minCost = INT_MAX, bestLen = -1;
	for( int i = 0; i < MAX_DEPTH; i++ ) {
		if( front[i].size() > 0 && front[i].begin()->first < minCost ) {
			minCost = front[i].begin()->first;
			bestLen = i;
		}
	}
	CNode res = front[bestLen].begin()->second;
	front[bestLen].erase(front[bestLen].begin());
	return res;
}

void removeNode( const int MAX_DEPTH, std::array <std::multimap<int, CNode>, 100> &front )
{
	int i = 0;
	while( front[i++].size() == 0 ) {}
	i--;
	front[i].erase(--front[i].end());
	return ;
}

std::string AStarMemLimit( const TBoard& startPos )
{
	int dist;
	int numVisited = 0;
	int maxFrontier = 0;
	CNode startNode, curNode;
	startNode.FirstInit(startPos);

	double magicConstant = 1; // warning: value more than 1 will make unacceptable estimation function, set it only when debug
	const int MAX_FRONT_SIZE = 1e5;
	const int MAX_DEPTH = 100;
	std::array <std::multimap<int, CNode>, MAX_DEPTH> front;
	for (int i = 0; i < MAX_DEPTH; i++) {
		front[i].clear();
		//std::cout << front[i].size() << std::endl;
	}

	int frontLen = 1;
	dist = GetManhattanDistance(startPos);
	front[startNode.Path.length()].insert(std::pair <int, CNode>(int(magicConstant * dist), startNode));

	std::set<TBoard> visited;
	while( frontLen > 0 && numVisited < 1e7 ) {
		curNode = popNode(MAX_DEPTH, front);
		frontLen--;
		numVisited++;

		visited.insert(curNode.Board);
		while (frontLen > MAX_FRONT_SIZE) {
			removeNode(MAX_DEPTH, front);
			frontLen--;
		}

		if (curNode.Board == SOLUTION) {
			//std::cout << "Visited: " << numVisited << std::endl;
			//std::cout << "Max Frontier: " << maxFrontier << std::endl;
			return curNode.Path;
		}

		std::vector<CNode> nextNodes;
		nextNodes.clear();
		GetSuccessors(curNode, nextNodes);
		for (int i = 0, maxi = nextNodes.size(); i < maxi; i++) {
			if (visited.find(nextNodes[i].Board) != visited.end() || nextNodes[i].Path.length() >= MAX_DEPTH) {
				continue;
			}
			dist = GetManhattanDistance(nextNodes[i].Board);
			int cost = int(dist * magicConstant) + int(nextNodes[i].Path.length());

			front[nextNodes[i].Path.length()].insert(std::pair <int, CNode>(cost, nextNodes[i]));
			frontLen++;

			if( frontLen > maxFrontier ) {
				maxFrontier = frontLen;
			}
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

//----------------------RBFS START-------------------------------------------
std::map<TBoard, int> elevatedDistances;
std::set<TBoard> currentPath;
int min_distance;

bool distancesComparer(std::pair<int, CNode> a, std::pair<int, CNode> b) {
	return a.first < b.first;
}

int GetDistance(const TBoard& pos) {
	auto it = elevatedDistances.find(pos);
	if (it != elevatedDistances.end()) {
		return it->second;
	}
	else {
		return GetManhattanDistance(pos);
	}
}

struct CReturnResult
{
	CNode solution;
	bool isNull;
	int bound;
};

void SetDistance(const TBoard& pos, int value) {
	elevatedDistances.insert(std::make_pair(pos, value));
}

CReturnResult RBFS(CNode currentNode, int bound, CReturnResult result) {
	CReturnResult resTemp;
	if (GetManhattanDistance(currentNode.Board) == 0) {
		resTemp.bound = bound;
		resTemp.isNull = false;
		resTemp.solution = currentNode;
		return resTemp;
	}

	std::vector<CNode> successors;
	GetSuccessors(currentNode, successors);

	std::vector<CNode> filteredSuccessors;
	for (auto it = successors.begin(); it != successors.end(); it++) {
		if (currentPath.find(it->Board) == currentPath.end()) {
			filteredSuccessors.push_back(*it);
		}
	}

	if (filteredSuccessors.size() == 0) {
		resTemp.isNull = true;
		resTemp.bound = GetDistance(currentNode.Board);
		return resTemp;
	}
	
	TBoard previousBest;
	for (;;) {
		
		int alternative = 0;

		std::vector<std::pair<int, CNode>> distances;
		for (int i = 0; i < filteredSuccessors.size(); i++) {
			distances.push_back(std::make_pair(GetDistance(filteredSuccessors[i].Board), filteredSuccessors[i]));
		}
		std::sort(distances.begin(), distances.end(), distancesComparer);

		std::pair<int, CNode> best = distances[0];
		
		int bestDistance = best.first;
		if (previousBest != best.second.Board) {
			previousBest = best.second.Board;
		}
		else {
			resTemp.isNull = true;
			resTemp.bound = bestDistance;
			return resTemp;
		}
		if (bestDistance < min_distance) {
			min_distance = bestDistance;
			//std::cout << min_distance << std::endl;
		}
		if (bestDistance > bound) {
			resTemp.isNull = true;
			resTemp.bound = bestDistance;
			return resTemp;
		}

		if (distances.size() > 1) {
			alternative = distances[1].first;
		}
		else {
			alternative = INT_MAX;
		}

		bound = bound < alternative ? bound : alternative;

		currentPath.insert(best.second.Board);
		resTemp = RBFS(best.second, bound, result);
		currentPath.erase(best.second.Board);

		SetDistance(best.second.Board, resTemp.bound);
		bound = resTemp.bound;

		if (!resTemp.isNull) {
			return resTemp;
		}
	}
}

std::string RecursiveBestFirstSearch(const TBoard& startPos) {
	CNode startNode;
	startNode.FirstInit(startPos);
	elevatedDistances.clear();
	currentPath.clear();
	CReturnResult result;
	result.isNull = false;
	result.solution = startNode;
	result.bound = 0;
	currentPath.insert(startNode.Board);
	min_distance = INT_MAX;
	CReturnResult resTemp = RBFS(startNode, INT_MAX, result);
	
	if (resTemp.isNull) {
		return NO_SOLUTIONS;
	}
	else {
		return resTemp.solution.Path;
	}
}

//----------------------------------RBFS END---------------------------------

//------------------------BRANCH AND BOUNDARY--------------------------------

std::string BranchAndBoundaryBFS( const TBoard& startPos, int problemUpperBound, int queueStopSize )
{
	int dist;
	int numVisited = 0;
	int maxFrontier = 0;
	CNode startNode;
	startNode.FirstInit( startPos );

	std::multimap<int, CNode> pqueue;
	dist = GetManhattanDistance( startPos );
	pqueue.insert( std::pair <int, CNode>( dist, startNode ) );
	std::set<TBoard> visited;
	std::vector<CNode> nextNodes;

	while( !pqueue.empty() ) {
		CNode curNode = pqueue.begin()->second;
		pqueue.erase( pqueue.begin() );
		numVisited += 1;
		visited.insert( curNode.Board );

		if( curNode.Board == SOLUTION ) {
			return curNode.Path;
		}

		nextNodes.clear();
		GetSuccessors( curNode, nextNodes );
		for( int i = 0; i < nextNodes.size(); i++ ) {
			if( visited.find( nextNodes[i].Board ) != visited.end() ) {
				continue;
			}
			dist = GetManhattanDistance( nextNodes[i].Board );
			int solutionLowerBound = dist + int( nextNodes[i].Path.length() );
			
			if( solutionLowerBound > problemUpperBound ) {
				continue;
			}

			pqueue.insert( std::pair<int, CNode>( solutionLowerBound, nextNodes[i] ) );
		}

		if( pqueue.size() > queueStopSize ) {
			break;
		}
	}
	return NO_SOLUTIONS;
}

std::string BranchAndBoundaryFromRBFSHeuristic( const TBoard& startPos )
{
	int problemUpperBound = RecursiveBestFirstSearch( startPos ).length();
	return BranchAndBoundaryBFS( startPos, problemUpperBound, 1000000000 );
}

//------------------------BRANCH AND BOUNDARY--------------------------------

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
	//CAlgorithm alg1("IDA*", &IdaStar);
	//Algorithms.push_back(alg1);

	//CAlgorithm alg2("A*", &AStar);
	//Algorithms.push_back(alg2);

	// Не находит ничего
	//CAlgorithm beamSearch1024("Beam search front capacity: 1024", &BeamSearchFrontCap1024);
	//Algorithms.push_back(beamSearch1024);

	// Работает 70 сек на одном тесте
	//CAlgorithm beamSearch262144("Beam search front capacity: 262144", &BeamSearchFrontCap262144);
	//Algorithms.push_back(beamSearch262144);

	CAlgorithm aStarMemLimit("A* with memory limit", AStarMemLimit);
	Algorithms.push_back(aStarMemLimit);
	
	// Работает 16мс на одном прогоне, правда результат отстает по длине на 806.
	CAlgorithm rbfs("RBFS", &RecursiveBestFirstSearch);
	Algorithms.push_back(rbfs);

	// С этой эвристикой 20с на тест, гарантирует оптимальность решения
	CAlgorithm branchAndBoundary( "B&B", &BranchAndBoundaryFromRBFSHeuristic );
	Algorithms.push_back( branchAndBoundary );
}

void RunBenchmark() 
{
	for( int testNumber = 0; testNumber < NUMBER_OF_TESTS; testNumber++ ) {

		std::cout << "Start test " << testNumber << std::endl;

		TBoard startPos;
		std::string path;
		startPos = GenerateStartPos();
		int minPathLen = INT_MAX;
		std::vector<std::string> paths;
		paths.clear();
		for( auto pAlg = Algorithms.begin(); pAlg != Algorithms.end(); pAlg++ ) {

			std::cout << "Testing " << pAlg->Name << std::endl;

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
	std::cout << NUMBER_OF_TESTS << " tests\n";
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
