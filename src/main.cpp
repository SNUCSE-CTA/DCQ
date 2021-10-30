#include <cassert>

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <random>

#if defined(_METHOD_C)
#include <map>
#include <set>
#endif

#include "graph.h"
#include "algorithm.h"
#include "timer.h"

using namespace std;

enum Mode { INDEX, QUERY, NEIGHBOR, STAT, ALL, FILTER };

struct IndexUnit {
	int32_t id;
	Graph* g;
	Coloring* coloring = NULL;
	pair<int32_t, int32_t>* cldist = NULL; //color-label distribution
	int32_t cllen; //length of cldist
    #ifdef _METHOD_C
    vector<int> degseq;
    #endif
};

//INPUT
void printHelp();
int32_t parseMode(char* aString);

//INDEX and QUERY
void initReadGraphs();
bool readGraphs(string aFileName, Graph* aGraph, int& aId);
void copyColoring(Coloring* aCol1, Coloring* aCol2);
void printGraph(ofstream& out, int32_t aId, Graph* aGraph);
void printColoring(ofstream& aFile, int32_t aId, Coloring* aColoring, Graph* aGraph);
void printCLDist(ofstream& aFile, int32_t aId, pair<int32_t,int32_t>* aDist, int32_t aLen);
#if defined(_METHOD_C)
double readIndex(string aInputFile, string aColorFile, string aDistFile, vector<IndexUnit>& index, vector<pair<vector<int32_t>, int32_t>>& spec);
#else
double readIndex(string aInputFile, string aColorFile, string aDistFile, vector<IndexUnit>& index);
#endif
pair<int, int> binarySearch(vector<IndexUnit>& index, pair<int32_t,int32_t>* aDist, int32_t aLen, int32_t lo=-1, int32_t hi=-1);

//NEIGHBOR and STAT
void extractNeighborGraph(int32_t aVertex, int32_t aD, Graph* aG, int32_t& aGraphId, string aOutputFile);
void printStat(string aInputFile);


//used in readGraphs
streampos currPos;
bool hasCurrPos = false;

#ifdef _METHOD_C
bool operator<(const vector<int32_t>& a, const vector<int32_t>& b){
    if (a.size() != b.size()) {
        return a.size() < b.size();
    }

    const size_t size = a.size();
    for (size_t i = 0; i < size; ++i) {
        if (a[i] != b[i]) {
            return a[i] < b[i];
        }
    }

    return false;
}

int compare(const vector<int32_t>& a, const vector<int32_t>& b) {
    if (a.size() != b.size()) {
        return static_cast<int>(a.size()) - static_cast<int>(b.size());
    }

    const size_t size = a.size();
    for (size_t i = 0; i < size; ++i) {
        if (a[i] != b[i]) {
            return static_cast<int>(a[i]) - static_cast<int>(b[i]);
        }
    }

    return 0;
}

int binary_search(vector<pair<vector<int32_t>, int32_t>>& spec, vector<int32_t>& degseq) {
    int lo = 0;
    int hi = spec.size() - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        int c = compare(spec[mid].first, degseq);
        if (c == 0) {
            return mid;
        } else if (c < 0) {
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }

    return -1;
}
#endif

int main(int argc, char* argv[])
{
	Timer timer;
	int32_t mode = -1; //type of job to do
	int32_t d = -1;
	int32_t limit = -1;
	string inputFile;
	string outputFile;
	string queryFile;
	string distFile;
	string colorFile;

	//parse input arguments
	for(int32_t i = 1; i < argc; ++i) {
		if(argv[i][0] == '-') {
			switch(argv[i][1]) {
				case 'm': //mode
					mode = parseMode(argv[i+1]);
					break;
				case 'd': //d-neighborhoods
					d = atoi(argv[i+1]);
					break;
				case 'f':
					inputFile = argv[i+1];
					break;
				case 'o':
					outputFile = argv[i+1];
					break;
				case 'l':
					limit = atoi(argv[i+1]);
					break;
				case 'i':
					distFile = argv[i+1];
					break;
				case 'c':
					colorFile = argv[i+1];
					break;
				case 'q':
					queryFile = argv[i+1];
					break;

			}
		} // if
	} // for

	if( mode == -1 ) {
		printHelp();
		return -1;
	}

	switch(mode) {
        case FILTER: {
            Refinement R;
            vector<pair<vector<int>, vector<pair<int, int>>>> S;
            vector<pair<Graph*, Coloring*>> I;

            // Read data graphs.
            bool hasGraph = false;
			do {
                int id;
                Graph* G = new Graph;
				hasGraph = readGraphs(inputFile, G, id);
				if (hasGraph) {
                    // Generate degree sequence.
                    vector<int> degseq;
                    degseq.reserve(G->numNode);
                    for (int i = 0; i < G->numNode; ++i) {
                        degseq.push_back(G->d[i]);
                    }
                    sort(degseq.rbegin(), degseq.rend());

                    // Generate color-label distribution.
					initGlobal(G->numNode);
					R.runSingle(G);

					Coloring* col = new Coloring;
					copyColoring(R.getStableColoring(), col);

                    vector<pair<int, int>> cldist;
                    for (int i = 0; i < G->numNode; ++i) {
                        cldist.push_back(make_pair(col->color[col->inv[i]], G->l[i]));
                    }
                    sort(cldist.begin(), cldist.end());

                    S.push_back(make_pair(degseq, cldist));
                    I.push_back(make_pair(G, col));
				}
			} while (hasGraph);

            int ACC = 0;
            int BCC = 0;
            int CCC = 0;
            int DDD = 0;

            double rA = 0;
            double rB = 0;
            double rC = 0;

            // Read query graphs.
			initReadGraphs();
			do {
                int id;
                Graph* G = new Graph;
				hasGraph = readGraphs(queryFile, G, id);
                if (!hasGraph) {
                    break;
                }

                vector<int> degseq;
                degseq.reserve(G->numNode);
                for (int i = 0; i < G->numNode; ++i) {
                    degseq.push_back(G->d[i]);
                }
                sort(degseq.rbegin(), degseq.rend());

                initGlobal(G->numNode);
                R.runSingle(G);
                Coloring* col = R.getStableColoring();

                vector<pair<int, int>> cldist;
                for (int i = 0; i < G->numNode; ++i) {
                    cldist.push_back(make_pair(col->color[col->inv[i]], G->l[i]));
                }
                sort(cldist.begin(), cldist.end());

                int acc = 0;
                int bcc = 0;
                int ccc = 0;
                int ddd = 0;
                Algorithm alg;
                for (size_t i = 0; i < S.size(); ++i) {
                    auto& s = S[i];
                    bool x = (s.first == degseq);
                    bool y = (s.second == cldist);

                    if (x) ++acc;
                    if (y) ++bcc;
                    if (x && y) {
                        ++ccc;
					    if (alg.run(G, I[i].first, col, I[i].second, R.getNumTreeNodeSingle())) {
                            ++ddd;
                        }
                    }
                }
                ACC += acc;
                BCC += bcc;
                CCC += ccc;
                DDD += ddd;

                if (acc != 0) rA += static_cast<double>(acc - ddd) / static_cast<double>(acc);
                if (bcc != 0) rB += static_cast<double>(bcc - ddd) / static_cast<double>(bcc);
                if (ccc != 0) rC += static_cast<double>(ccc - ddd) / static_cast<double>(ccc);

            } while (hasGraph);

            cout << "degree sequence FPR: " << rA/200 << endl;
            cout << "color-label distribution FPR: " << rB/200 << endl;
            cout << "combination FPR: " << rC/200 << endl;

            break;
        }
        case ALL: {
			Graph g;
			Refinement cr;
			vector<IndexUnit> index;
			double refineTime = 0;
			double sumRefineTime = 0;
			double sumSortTime = 0;

			//READ graphs and COLOR REFINEMENT
			initReadGraphs();
			bool hasGraph = false;
			do {
				IndexUnit unit;
				unit.g = new Graph;

				hasGraph = readGraphs(inputFile, unit.g, unit.id);
				if (hasGraph) {
timer.start();
					//REFINEMENT
                    #ifdef _METHOD_C
                    unit.degseq.reserve(unit.g->numNode);
                    for (int i = 0; i < unit.g->numNode; ++i) {
                        unit.degseq.push_back(unit.g->d[i]);
                    }
                    // sort(unit.degseq.begin(), unit.degseq.end());
                    sort(unit.degseq.rbegin(), unit.degseq.rend());
                    #endif

					initGlobal(unit.g->numNode);
					cr.runSingle(unit.g);

					unit.coloring = new Coloring;
					copyColoring(cr.getStableColoring(), unit.coloring);
					
					//color-label distribution
					unit.cllen = unit.g->numNode;
					unit.cldist = new pair<int32_t, int32_t>[unit.cllen];
					for(int i = 0; i < unit.g->numNode; ++i) {
						unit.cldist[i].first = unit.coloring->color[ unit.coloring->inv[i] ]; //color
						unit.cldist[i].second = unit.g->l[i];//label
					}
					sort(unit.cldist, unit.cldist + unit.cllen, 
					[](const pair<int32_t, int32_t>& a, const pair<int32_t, int32_t>& b) -> bool {
						if(a.first != b.first)
							return a.first < b.first;
						else
							return a.second < b.second;
					});
					
					index.push_back(unit);
refineTime = timer.end();
sumRefineTime += refineTime;
				}
				else {
					delete unit.g;
				}
			} while(hasGraph);

timer.start();
			clearGlobal(); 
			cr.clearSingle();
sumRefineTime += timer.end();
			
			//SORT Index by color-label distribution
timer.start();
			sort(index.begin(), index.end(),
			[](const IndexUnit& a, const IndexUnit& b) -> bool {
                #if defined(_METHOD_C)
                int c = compare(a.degseq, b.degseq);
                if (c == 0) {
                #endif
                    if (a.cllen != b.cllen) {
                        return a.cllen < b.cllen;
                    } else {
                        int32_t i = 0;
                        while (i < a.cllen - 1) {
                            if (a.cldist[i].first != b.cldist[i].first)
                                break;
                            if (a.cldist[i].second != b.cldist[i].second)
                                break;
                            ++i;
                        }
                        
                        if (a.cldist[i].first != b.cldist[i].first) {
                            return a.cldist[i].first < b.cldist[i].first;
                        } else {
                            return a.cldist[i].second < b.cldist[i].second;
                        }
                    }
                #if defined(_METHOD_C)
                } else {
                    return c < 0;
                }
                #endif
			});

            #if defined(_METHOD_C)
            vector<pair<vector<int32_t>, int32_t>> spec;

            const int size = index.size();
            for (int i = 0; i < size; ++i) {
                const auto& unit = index[i];
                bool is_new = false;
                if (spec.empty()) {
                    is_new = true;
                } else {
                    const vector<int32_t>& last_degseq = spec.back().first;
                    if (unit.degseq.size() != last_degseq.size()) {
                        is_new = true;
                    } else {
                        bool is_identical = true;
                        for (int32_t i = 0; i < unit.degseq.size(); ++i) {
                            if (unit.degseq[i] != last_degseq[i]) {
                                is_identical = false;
                                break;
                            }
                        }

                        if (!is_identical) {
                            is_new = true;
                        }
                    }
                }

                if (is_new) {
                    spec.emplace_back(unit.degseq, i);
                }
            }
            #endif
sumSortTime += timer.end();

			Graph* query = new Graph;
			pair<int32_t, int32_t>* cldist = NULL;
			int32_t cllen = 0;
			int32_t id;
			double candTime = 0;
			double checkTime = 0;
			double sumCandTime = 0;
			double sumCheckTime = 0;

            #ifdef BREAKDOWN
            Timer tt;
            double sum_degseq_time = 0;
            double sum_degseq_search_time = 0;
            double sum_collab_time = 0;
            double sum_collab_search_time = 0;
            double sum_backtrack_time = 0;

            int32_t bkd_candidates_degseq = 0;
            int32_t bkd_candidates_collab = 0;
            int32_t bkd_matches = 0;
            #endif

			int num_query_graphs_not_found = 0;
            int not_found_degseq = 0;
            int no_candidates = 0;
            int num_matches = 0;
            int num_queries = 0;
            int num_candidates = 0;

			vector<vector<int> > answer;

			initReadGraphs();
			while( readGraphs(queryFile, query, id) ) {
                ++num_queries;
                vector<int> ans;

                Refinement cr;
                Coloring* coloring = nullptr;
                pair<int, int> range = make_pair(0, -1);
timer.start();
                #if defined(_METHOD_C)
                #ifdef BREAKDOWN
                tt.start();
                #endif
                vector<int> degseq;
                degseq.reserve(query->numNode);
                for (int i = 0; i < query->numNode; ++i) {
                    degseq.push_back(query->d[i]);
                }
                // sort(degseq.begin(), degseq.end());
                sort(degseq.rbegin(), degseq.rend());
                #ifdef BREAKDOWN
                sum_degseq_time += tt.end();
                tt.start();
                #endif
                int x = binary_search(spec, degseq);
                #ifdef BREAKDOWN
                sum_degseq_search_time += tt.end();
                #endif

                if (x == -1) {
                    #ifdef BREAKDOWN
                    cerr << 0 << endl;
                    #endif
                    ++not_found_degseq;
                } else {
                    int lo = spec[x].second;
                    int hi = (x+1 != spec.size()) ? spec[x+1].second : index.size();
                    #ifdef BREAKDOWN
                    cerr << hi - lo << ' ';
                    bkd_candidates_degseq += (hi - lo);
                    #endif
                #endif
                    //REFINE
                    #ifdef BREAKDOWN
                    tt.start();
                    #endif
                    initGlobal(query->numNode); //init markCell and markNode;
                    cr.runSingle(query);
                    coloring = cr.getStableColoring();

                    //color-label distribution of query
                    if( cllen < query->numNode ) { //cllen denotes the length of cldist array (for recycling array)
                        if( cldist != NULL )
                            delete[] cldist;
                        cldist = new pair<int32_t, int32_t>[query->numNode];
                        cllen = query->numNode;
                    }
                    for(int i = 0; i < query->numNode; ++i) {
                        cldist[i].first = coloring->color[ coloring->inv[i] ];
                        cldist[i].second = query->l[i];
                    }
                    sort(cldist, cldist + query->numNode, 
                    [](const pair<int32_t, int32_t>& a, const pair<int32_t, int32_t>& b) -> bool {
                        if(a.first != b.first)
                            return a.first < b.first;
                        else
                            return a.second < b.second;
                    });
                    #ifdef BREAKDOWN
                    sum_collab_time += tt.end();
                    tt.start();
                    #endif

                    #if defined(_METHOD_C)
                    range = binarySearch(index, cldist, query->numNode, lo, hi - 1);
                    #else
                    range = binarySearch(index, cldist, query->numNode);
                    #endif

                    #ifdef BREAKDOWN
                    sum_collab_search_time += tt.end();
                    #endif

                    #ifdef BREAKDOWN
                    if (range.first >= range.second) {
                        cerr << 0 << endl;
                    } else {
                        cerr << (range.second - range.first) << ' ';
                        bkd_candidates_collab += (range.second - range.first);
                    }
                    #endif

                #if defined(_METHOD_C)
                }
                #endif
candTime = timer.end();
sumCandTime += candTime;

                if (range.first >= range.second) { 
                    checkTime = 0; 
                    ++num_query_graphs_not_found;
                    ++no_candidates;
                    answer.push_back(ans);
                    continue;
                } else {
                    num_candidates += (range.second - range.first);
                }
                    
				//VERIFICATION (Backtrack)
                int _num_matches = 0;
timer.start();
					Algorithm alg;
                #ifdef BREAKDOWN
                tt.start();
                #endif

				for (int32_t i = range.first; i < range.second; ++i) {
					bool same = alg.run(query, index[i].g, coloring, index[i].coloring, cr.getNumTreeNodeSingle());
					if (same) {
                        ++_num_matches;
                        ans.push_back(index[i].id);
					}
				}

                #ifdef BREAKDOWN
                sum_backtrack_time += tt.end();
                #endif
checkTime = timer.end();
sumCheckTime += checkTime;

                answer.push_back(ans);

				clearGlobal();

                if (_num_matches == 0) {
					++num_query_graphs_not_found;
                    #ifdef BREAKDOWN
                    cerr << 0 << endl;
                    #endif
				} else {
                    num_matches += _num_matches;
                    #ifdef BREAKDOWN
                    cerr << _num_matches << endl;
                    bkd_matches += _num_matches;
                    #endif
                }

			} // while
			delete query;
			delete[] cldist;
			
			cout << "File: " << queryFile << endl;
			cout << "Refine Build BinarySearch Verification" << endl;
			cout << setprecision(4) << fixed << sumSortTime << ' ' << sumRefineTime << ' ' << sumCandTime << " " << sumCheckTime << endl;
			cout << "# query graphs not found: " << num_query_graphs_not_found << endl;
			cout << "# not found degseq: " << not_found_degseq << endl;
            cout << "# not found collab: " << no_candidates - not_found_degseq << endl;
            cout << "# not found btrack: " << num_query_graphs_not_found - no_candidates << endl;
            cout << "# candidates: " << num_candidates << endl;
            // cout << "# check colorings: " << num_check_colorings << endl;
            // cout << "# backtracks: " << num_backtracks << endl;
            cout << "avg. # matches: " << static_cast<double>(num_matches) / static_cast<double>(num_queries - num_query_graphs_not_found) << endl;

			ofstream out("answer.txt");
			for(int i = 0; i < answer.size(); ++i) {
				out << i << ": ";
				sort(answer[i].begin(), answer[i].end());
				for(int j = 0; j < answer[i].size(); ++j) {
					out << answer[i][j];
					if( j < answer[i].size() - 1 )
						out << " ";
				}
				out << endl;
			}
			out.close();
            
            #ifdef BREAKDOWN
            cout << "breakdown" << endl;
            cout << "  total degseq time: " << sum_degseq_time << ' ' << sum_degseq_search_time << endl;
            cout << "  total collab time: " << sum_collab_time << ' ' << sum_collab_search_time << endl;
            cout << "  total btrack time: " << sum_backtrack_time << endl;
            // cout << not_found_degseq << ' ' << no_candidates - not_found_degseq << ' ' << num_query_graphs_not_found - no_candidates << endl;
            cout << bkd_candidates_degseq << ' ' << bkd_candidates_collab << ' ' << bkd_matches << endl;
            cout << sum_degseq_time << ' ' << sum_degseq_search_time << ' ' << sum_collab_time << ' ' << sum_collab_search_time << ' ' << sum_backtrack_time << endl;
            #endif

            // cout << endl << t1 << ' ' << t2 << ' ' << t3 << ' ' << t4 << endl;

			//FREE
			for(int32_t i = 0; i < index.size(); ++i) {
				delete index[i].g;
				delete index[i].coloring;
				delete[] index[i].cldist;
			}
            break;
        }
		case INDEX:
		{
			Graph g;
			Refinement cr;
			vector<IndexUnit> index;
			double refineTime = 0;
			double sumRefineTime = 0;
			double sumSortTime = 0;

			//READ graphs and COLOR REFINEMENT
			initReadGraphs();
			bool hasGraph = false;
			do {
				IndexUnit unit;
				unit.g = new Graph;

				hasGraph = readGraphs(inputFile, unit.g, unit.id);
				if(hasGraph) {
timer.start();
					//REFINEMENT
                    #ifdef _METHOD_C
                    unit.degseq.reserve(unit.g->numNode);
                    for (int i = 0; i < unit.g->numNode; ++i) {
                        unit.degseq.push_back(unit.g->d[i]);
                    }
                    sort(unit.degseq.begin(), unit.degseq.end());
                    #endif

					initGlobal(unit.g->numNode);
					cr.runSingle(unit.g);
// cout << setprecision(4) << fixed << "Graph " << unit.id << ": " << refineTime << " sec" << endl;

					unit.coloring = new Coloring;
					copyColoring(cr.getStableColoring(), unit.coloring);
					
					//color-label distribution
					unit.cllen = unit.g->numNode;
					unit.cldist = new pair<int32_t, int32_t>[unit.cllen];
					for(int i = 0; i < unit.g->numNode; ++i) {
						unit.cldist[i].first = unit.coloring->color[ unit.coloring->inv[i] ]; //color
						unit.cldist[i].second = unit.g->l[i];//label
					}
					sort(unit.cldist, unit.cldist + unit.cllen, 
					[](const pair<int32_t, int32_t>& a, const pair<int32_t, int32_t>& b) -> bool {
						if(a.first != b.first)
							return a.first < b.first;
						else
							return a.second < b.second;
					});
					
					index.push_back(unit);
refineTime = timer.end();
sumRefineTime += refineTime;
				}
				else {
					delete unit.g;
				}
			} while(hasGraph);

timer.start();
			clearGlobal(); 
			cr.clearSingle();
sumRefineTime += timer.end();
			
			//SORT Index by color-label distribution
timer.start();
			sort(index.begin(), index.end(),
			[](const IndexUnit& a, const IndexUnit& b) -> bool {
                #if defined(_METHOD_C)
                int c = compare(a.degseq, b.degseq);
                if (c == 0) {
                #endif
                    if (a.cllen != b.cllen) {
                        return a.cllen < b.cllen;
                    } else {
                        int32_t i = 0;
                        while (i < a.cllen - 1) {
                            if (a.cldist[i].first != b.cldist[i].first)
                                break;
                            if (a.cldist[i].second != b.cldist[i].second)
                                break;
                            ++i;
                        }
                        
                        if (a.cldist[i].first != b.cldist[i].first) {
                            return a.cldist[i].first < b.cldist[i].first;
                        } else {
                            return a.cldist[i].second < b.cldist[i].second;
                        }
                    }
                #if defined(_METHOD_C)
                } else {
                    return c < 0;
                }
                #endif
			});
sumSortTime += timer.end();

cout << "File: " << inputFile << endl;
cout << "Refinement Sort Total" << endl;
cout << setprecision(4) << fixed << sumRefineTime << " " << sumSortTime << " " << sumRefineTime + sumSortTime << endl;

			//SAVE to files
			ofstream graph_out(outputFile);
			ofstream color_out(colorFile);
			ofstream dist_out(distFile);
			for(int32_t i = 0; i < index.size(); ++i) {
				printGraph(graph_out, index[i].id, index[i].g); //print original graph
				printColoring(color_out, index[i].id, index[i].coloring, index[i].g); //print coloring with modified degree
				printCLDist(dist_out, index[i].id, index[i].cldist, index[i].cllen);
			}
			graph_out.close();
			color_out.close();
			dist_out.close();

			//FREE
			for(int32_t i = 0; i < index.size(); ++i) {
				delete index[i].g;
				delete index[i].coloring;
				delete[] index[i].cldist;
			}
		} //INDEX
			break;

		case QUERY:
		{
			Graph* query = new Graph;
			pair<int32_t, int32_t>* cldist = NULL;
			int32_t cllen = 0;
			vector<IndexUnit> index;
			int32_t id;
			double candTime = 0;
			double checkTime = 0;
			double sumCandTime = 0;
			double sumCheckTime = 0;

            #if defined(_METHOD_C)
            // Binary search using deg. seq.
            vector<pair<vector<int32_t>, int32_t>> spec;
            double index_time = readIndex(inputFile, colorFile, distFile, index, spec);
            #else
			double index_time = readIndex(inputFile, colorFile, distFile, index);
            #endif

			int num_query_graphs_not_found = 0;
            int num_query_graphs_not_found_only_VE = 0;
            int no_candidates = 0;
            int num_matches = 0;
            int num_queries = 0;
            int num_candidates = 0;

			initReadGraphs();
			while( readGraphs(queryFile, query, id) ) {
                ++num_queries;
                Refinement cr;
                Coloring* coloring = nullptr;
                pair<int, int> range = make_pair(0, -1);
timer.start();
                #if defined(_METHOD_C)
                vector<int> degseq;
                degseq.reserve(query->numNode);
                for (int i = 0; i < query->numNode; ++i) {
                    degseq.push_back(query->d[i]);
                }
                sort(degseq.begin(), degseq.end());

                int x = binary_search(spec, degseq);
                if (x == -1) {
                    ++num_query_graphs_not_found_only_VE;
                } else {
                    int lo = spec[x].second;
                    int hi = (x+1 != spec.size()) ? spec[x+1].second : index.size()-1;
                #endif
                    //REFINE
                    initGlobal(query->numNode); //init markCell and markNode;
                    cr.runSingle(query);
                    coloring = cr.getStableColoring();

                    //color-label distribution of query
                    if( cllen < query->numNode ) { //cllen denotes the length of cldist array (for recycling array)
                        if( cldist != NULL )
                            delete[] cldist;
                        cldist = new pair<int32_t, int32_t>[query->numNode];
                        cllen = query->numNode;
                    }
                    for(int i = 0; i < query->numNode; ++i) {
                        cldist[i].first = coloring->color[ coloring->inv[i] ];
                        cldist[i].second = query->l[i];
                    }
                    sort(cldist, cldist + query->numNode, 
                    [](const pair<int32_t, int32_t>& a, const pair<int32_t, int32_t>& b) -> bool {
                        if(a.first != b.first)
                            return a.first < b.first;
                        else
                            return a.second < b.second;
                    });

                    #if defined(_METHOD_C)
                    range = binarySearch(index, cldist, query->numNode, lo, hi);
                    #else
                    range = binarySearch(index, cldist, query->numNode);
                    #endif
                #if defined(_METHOD_C)
                }
                #endif
candTime = timer.end();
sumCandTime += candTime;

                if( range.first >= range.second) { 
                    checkTime = 0; 
                    //cout << setprecision(4) << fixed << "Query " << id << ": " << candTime + checkTime << " sec" << endl; 
                    ++num_query_graphs_not_found;
                    ++no_candidates;
                    continue;
                } else {
                    num_candidates += (range.second - range.first);
                }
                    
timer.start();
				//VERIFICATION (Backtrack)
				// vector<int32_t> identical;
                int _num_matches = 0;
				for(int32_t i = range.first; i < range.second; ++i) {
					Algorithm alg;
					bool same = alg.run(query, index[i].g, coloring, index[i].coloring, cr.getNumTreeNodeSingle());
					if (same) {
						// identical.push_back(index[i].id);
                        ++_num_matches;
					}
				}
checkTime = timer.end();
sumCheckTime += checkTime;
				clearGlobal();

				// if (identical.empty()) {
                if (_num_matches == 0) {
					++num_query_graphs_not_found;
				} else {
                    // num_matches += static_cast<int>(identical.size());
                    num_matches += _num_matches;
                }
/*
cout << setprecision(4) << fixed << "Query " << id << ": " << candTime + checkTime << " sec: " << endl;
cout << id << ": ";
sort(identical.begin(), identical.end());
for(int32_t i = 0; i < identical.size(); ++i)
	cout << identical[i] << " ";
cout << endl;
*/

			} // while
			delete query;
			delete[] cldist;
			
			cout << "File: " << queryFile << endl;
			cout << "ReadIndex BinarySearch Verification Total" << endl;
			cout << setprecision(4) << fixed << index_time << ' ' << sumCandTime << " " << sumCheckTime << " " << sumCandTime + sumCheckTime << endl;
			cout << "# query graphs not found: " << num_query_graphs_not_found << endl;
			cout << "# query graphs not found (only VE): " << num_query_graphs_not_found_only_VE << endl;
            cout << "No candidates: " << no_candidates << endl;
            cout << "# candidates: " << num_candidates << endl;
            // cout << "# check colorings: " << num_check_colorings << endl;
            // cout << "# backtracks: " << num_backtracks << endl;
            // cout << "avg. # candidates: " << static_cast<double>(num_candidates) / static_cast<double>(num_queries) << endl;
            cout << "avg. # matches: " << static_cast<double>(num_matches) / static_cast<double>(num_queries) << endl;

			//FREE
			for(int32_t i = 0; i < index.size(); ++i) {
				delete index[i].g;
				delete index[i].coloring;
				delete[] index[i].cldist;
			}
		}//QUERY
			break;

		case NEIGHBOR: 
		{
			Graph* g = new Graph(inputFile);
			if(g->fail()) {
				cout << "There is no file " << inputFile << endl;
				delete g;
				return -1;
			}

			vector<int32_t> vt;
			for(int32_t i = 0; i < g->numNode; ++i) {
				vt.push_back(i);
			}
			if( limit > 0 ) { //shuffle when we randomly select vertices
				random_device rd;
				mt19937 g(rd());
				shuffle(vt.begin(), vt.end(), g);
			}
			if( limit > 0 )
				limit = (limit < g->numNode) ? limit : g->numNode;
			else
				limit = g->numNode;

			int32_t graphId = 0;
			for(int32_t i = 0; i < limit; ++i) {
				int32_t v = vt[i];
				extractNeighborGraph(v, d, g, graphId, outputFile);
			}
		
			delete g;
		}
			break;

		case STAT:
		{
			printStat(inputFile);
		}
			break;

	} // switch

	return 0;
}

void printHelp()
{
	cout << "**********Usage**********" << endl;
	cout << "-m: select a mode from the following options" << endl;
	cout << "    " << "neighbor: extract d-neighborhood graph" << endl;
	cout << "    " << "          (required parameters: -d, -f, -o)" << endl;
	cout << "    " << "          (optional parameter: -l)" << endl;
	cout << "    " << "stat: statistics of the graphs" << endl;
	cout << "    " << "      (required parameter: -f)" << endl;
	cout << "    " << "index: build index of the input set of graphs" << endl;
	cout << "    " << "       (required parameters: -f, -i, -c, -o)" << endl;
	cout << "    " << "query: isomorphism query" << endl;
	cout << "    " << "       (required parameters: -q, -i, -c)" << endl;
	cout << "-d: set d (distance) of d-neighborhood" << endl;
	cout << "-f: input file name" << endl;
	cout << "-i: invariant file name" << endl;
	cout << "-c: coloring file name" << endl;
	cout << "-o: output file name" << endl;
	cout << "-q: query file name" << endl;
	cout << "-l: limit" << endl;
}

int32_t parseMode(char* aString)
{
	if( strcmp(aString, "neighbor") == 0 )
		return NEIGHBOR;
	if( strcmp(aString, "stat") == 0 )
		return STAT;
	if( strcmp(aString, "index") == 0 )
		return INDEX;
	if( strcmp(aString, "query") == 0 )
		return QUERY;
	if( strcmp(aString, "all") == 0 )
		return ALL;
	if( strcmp(aString, "filter") == 0 )
		return FILTER;
	
	return -1;
}

//PRINT d-neighborhood of aVertex
void extractNeighborGraph(int32_t aVertex, int32_t aD, Graph* aG, int32_t& aGraphId, string aOutputFile)
{
	int32_t n = aG->numNode;
	int32_t* dneighbor = global_memory.getLLArray(n);
	int32_t neighCount = 0;
	int32_t* visited = global_memory.getLLArray(n);
	memset(visited, 0, sizeof(int32_t) * n);
	int32_t* queue = global_memory.getLLArray(n);
	int32_t queueStart = 0;
	int32_t queueEnd = 0;
	int32_t nextQueueEnd = 0;
	int32_t depth = 0;

	//Step A. find d-hop neighbors
	visited[aVertex] = 1;
	dneighbor[neighCount] = aVertex; //insert itself to d-neighborhood
	++neighCount;
	queue[queueEnd] = aVertex;
	++queueEnd;
	nextQueueEnd = queueEnd;

	while(true) {
		while( queueStart < queueEnd ) {
			int32_t curr = queue[queueStart];
			++queueStart;

			int32_t* neigh = aG->e[curr];
			for(int32_t i = 0; i < aG->d[curr]; ++i) {
				int32_t child = neigh[i];
				if( visited[child] == 0 ) {
					visited[child] = 1;
					queue[nextQueueEnd] = child;
					++nextQueueEnd;

					dneighbor[neighCount] = child;
					++neighCount;
				}
			} //for(i)
		} //while(queueStart < queueEnd)

		if( queueEnd == nextQueueEnd )
			break;

		queueStart = queueEnd;
		queueEnd = nextQueueEnd;

		++depth;
		if(depth == aD)
			break;
	} //while(true)
	global_memory.returnLLArray(queue, n); //we do not need a queue after BFS

	//Step B. Print d-neighborhood
	ofstream out(aOutputFile, fstream::out | fstream::app);
	int32_t* newname = global_memory.getLLArray(n);
		
	//vertex ID should be in [0, #vertices)
	for(int32_t i = 0; i < neighCount; ++i) {
		int32_t node = dneighbor[i];
		newname[node] = i;
	}
	
	//print components to the output file
	out << "t " << aGraphId << " " << neighCount << endl;
	++aGraphId;
	for(int32_t i = 0; i < neighCount; ++i)
	out << "v " << i << " " << aG->l[ dneighbor[i] ] << endl;
	for(int32_t i = 0; i < neighCount; ++i) {
		int32_t node = dneighbor[i];
		for(int32_t j = 0; j < aG->d[node]; ++j) {
			int32_t adj = aG->e[node][j];
			if(visited[adj] == 1) {
				if( newname[node] < newname[adj] )
					out << "e " << newname[node] << " " << newname[adj] << " 0" << endl;
			}
		} //for(j)
	} //for(i)
	out.close();
	global_memory.returnLLArray(dneighbor, n);
	global_memory.returnLLArray(visited, n);
	global_memory.returnLLArray(newname, n);
}

void printStat(string aInputFile)
{
	ifstream in(aInputFile);
	if(in.fail()) {
		cout << "There is no file " << aInputFile << endl;
		return;
	}

	int64_t sumNumNode = 0;
	int64_t sumNumEdge = 0;
	int32_t numGraph = -1;
	int32_t maxNumNode = 0;
	int32_t minNumNode = 2100000000;
	int32_t numNode = 0;
	int32_t numEdge = 0;
	int32_t atleastthree = 0;
	int32_t al100 = 0;
	int32_t al1000 = 0;
	int32_t al10000 = 0;

	char tag;
	int32_t vid, vlabel, left, right, elabel, gid;
	string line;
	
	while(in >> tag) {
		if(tag == 't') {
			++numGraph;
			sumNumNode += numNode;
			sumNumEdge += numEdge;
			if( maxNumNode < numNode )
				maxNumNode = numNode;
			if( numNode != 0 && minNumNode > numNode )
				minNumNode = numNode;
			if( numNode >= 3 )
				++atleastthree;
			if( numNode >= 100 )
				++al100;
			if( numNode >= 1000 )
				++al1000;
			if( numNode >= 10000 )
				++al10000;

			numNode = 0;
			numEdge = 0;
			getline(in, line);
		}
		else if(tag == 'v') {
			++numNode;
			getline(in, line);
		}
		else if(tag == 'e') {
			++numEdge;
			getline(in, line);
		}
	}
	++numGraph;
	sumNumNode += numNode;
	sumNumEdge += numEdge;
	if( maxNumNode < numNode )
		maxNumNode = numNode;
	if( numNode != 0 && minNumNode > numNode )
		minNumNode = numNode;
	if( numNode >= 3 )
		++atleastthree;

	in.close();

	cout << "#graphs: " << numGraph << endl;
	cout << "avg. #vertices: " << sumNumNode / (double)numGraph << endl;
	cout << "avg. #edges: " << sumNumEdge / (double)numGraph << endl;
	cout << "min #vertices: " << minNumNode << endl;
	cout << "max #vertices: " << maxNumNode << endl;
	cout << "at least 3 vertices: " << atleastthree << endl;
	cout << "at least 100, 1000, 10000: " << endl;
	cout << al100 << "\t" << al1000 << "\t" << al10000 << endl;
}


void initReadGraphs()
{
	hasCurrPos = false;
}

bool readGraphs(string aFileName, Graph* aGraph, int& aId)
{
	if(aGraph == NULL) {
		cout << "ERROR in readGraphs: aGraph == NULL" << endl;
		return false;
	}

	ifstream infile(aFileName);
	if(infile.fail()) {
		cout << "ERROR in readGraphs: " << aFileName << " not exists" << endl;
		return false;
	}

	if( hasCurrPos ) {
		infile.clear();
		infile.seekg(currPos);
	}

	char tag;
	int32_t vid, vlabel, left, right, elabel, gid;
	int32_t numNode = 0;
	int32_t numEdge = 0;
	bool hasGraph = false;
	string line;
	streampos edgePos = infile.tellg();

	//first read: read vertices and count degree for each vertex
	while( infile >> tag ) {
		if( tag == 't' ) {
			if(hasGraph)
				break;
			else
				hasGraph = true;

			infile >> gid >> numNode;
			aId = gid;
			//allocate memory for vertices
			if(aGraph->nLen < numNode) {
				aGraph->clear();
				aGraph->od = new int32_t[numNode];
				aGraph->d = new int32_t[numNode];
				aGraph->l = new int32_t[numNode];
				aGraph->one = new char[numNode];
				aGraph->adjPos = new int32_t[numNode];
				aGraph->e = new int32_t*[numNode];
				aGraph->eLen = new int32_t[numNode];
				for(int32_t i = 0; i < numNode; ++i) {
					aGraph->e[i] = NULL;
					aGraph->eLen[i] = 0;
				}
				aGraph->nLen = numNode;
			}
			aGraph->numNode = numNode;
		}
		else if( tag == 'v' ) {
			infile >> vid >> vlabel;
			if( numNode <= vid )
				cout << "ERROR in " << __FUNCTION__ << "(): vid >= #nodes" << endl;
			aGraph->d[vid] = 0;
			aGraph->l[vid] = vlabel;
			aGraph->one[vid] = 0;

			edgePos = infile.tellg(); //record the last line that tag 'v' appears.
		}
		else if( tag == 'e' ) {
			infile >> left >> right >> elabel;
			if( numNode <= left )
				cout << "ERROR in " << __FUNCTION__ << "(): left >= #nodes" << endl;
			if( numNode <= right )
				cout << "ERROR in " << __FUNCTION__ << "(): right >= #nodes" << endl;
			++(aGraph->d[left]);
			++(aGraph->d[right]);
			++numEdge;
		}

		currPos = infile.tellg(); //in case there are no edges
		hasCurrPos = true;
	}

	if(hasGraph == false)
		return false;

	//allocate memory for edges
	for(int32_t i = 0; i < numNode; ++i) {
		if( aGraph->eLen[i] < aGraph->d[i] ) {
			if( aGraph->e[i] != NULL ) {
				delete[] aGraph->e[i];
				aGraph->e[i] = NULL;
			}
			aGraph->e[i] = new int32_t[ aGraph->d[i] ];
			aGraph->eLen[i] = aGraph->d[i];
		}
		aGraph->d[i] = 0;
	}
	aGraph->numEdge = numEdge;

	//second read: read edges
	infile.clear();
	infile.seekg(edgePos);

	while( infile >> tag ) {
		if( tag == 't' ) {
			break;
		}
		else if( tag == 'e' ) {
			infile >> left >> right >> elabel;
			aGraph->e[left][aGraph->d[left]] = right;
			aGraph->e[right][aGraph->d[right]] = left;
			++(aGraph->d[left]);
			++(aGraph->d[right]);
		}
		else {
			cout << "ERROR in " << __FUNCTION__ << "(): invalid file format" << endl;
			break;
		}

		currPos = infile.tellg();
		hasCurrPos = true;
	}
	//store original degree
	for(vid = 0; vid < numNode; ++vid)
		aGraph->od[vid] = aGraph->d[vid];

	//calculate the start position of edge list for each vertex
	aGraph->adjPos[0] = 0;
	for(vid = 1; vid < numNode; ++vid)
		aGraph->adjPos[vid] = aGraph->adjPos[vid - 1] + aGraph->d[vid - 1];

	infile.close();
	
	return true;
}

void copyColoring(Coloring* aCol1, Coloring* aCol2)
{
	if(aCol1 == NULL || aCol2 == NULL)
		return;

	if(aCol2 != NULL)
		aCol2->clear();
	
	aCol2->init(aCol1->numNode);
	aCol2->numNode = aCol1->numNode;
	aCol2->numCell = aCol1->numCell;
	for(int32_t i = 0; i < aCol1->numNode; ++i) {
		aCol2->color[i] = aCol1->color[i];
		aCol2->perm[i] = aCol1->perm[i];
		aCol2->inv[i] = aCol1->inv[i];
		aCol2->cellSize[i] = aCol1->cellSize[i];
	}
}

void printGraph(ofstream& out, int32_t aId, Graph* aGraph)
{
	out << "t " << aId << " " << aGraph->numNode << endl;
	for(int32_t i = 0; i < aGraph->numNode; ++i) {
		out << "v " << i << " " << aGraph->l[i] << endl;
	}
	for(int32_t i = 0; i < aGraph->numNode; ++i) {
		for(int32_t j = 0; j < aGraph->od[i]; ++j) { //note that we use original degree
			int32_t adj = aGraph->e[i][j];
			if(i < adj) {
				out << "e " << i << " " << adj << " 0" << endl;
			}
		}
	}
}

void printColoring(ofstream& out, int32_t aId, Coloring* aColoring, Graph* aGraph)
{
	out << "t " << aId << " " << aColoring->numNode << " ";
	for(int32_t i = 0; i < aColoring->numNode; ++i) {
		out << aColoring->color[ aColoring->inv[i] ] << " ";
	}

	for(int32_t i = 0; i < aGraph->numNode; ++i) {
		out << (int)(aGraph->one[i]);
		if( i < aGraph->numNode - 1)
			out << " ";
		else
			out << endl;
	}
}

void printCLDist(ofstream& out, int32_t id, pair<int32_t, int32_t>* dist, int32_t len)
{
	out << "t " << id << " " << len << " ";
	for(int i = 0; i < len; ++i) {
		out << dist[i].first << " " << dist[i].second;
		if( i < len - 1)
			out << " ";
		else
			out << endl;
	}
}

#if defined(_METHOD_C)
double readIndex(string aInputFile, string aColorFile, string aDistFile, vector<IndexUnit>& index, vector<pair<vector<int32_t>, int32_t>>& spec)
#else
double readIndex(string aInputFile, string aColorFile, string aDistFile, vector<IndexUnit>& index)
#endif
{
	int32_t id;
	char tag;
	int32_t ind = 0;
	Graph g;
	//READ GRAPHS
	initReadGraphs();
	bool hasGraph = false;

    Timer t;
    double index_time = 0.0;
	do {
		IndexUnit unit;
		unit.g = new Graph;
		hasGraph = readGraphs(aInputFile, unit.g, unit.id);

        if (hasGraph) {
            #if defined(_METHOD_C)
            unit.degseq.reserve(unit.g->numNode);
            for (int i = 0; i < unit.g->numNode; ++i) {
                unit.degseq.push_back(unit.g->d[i]);
            }
            sort(unit.degseq.begin(), unit.degseq.end());

    t.start();
            bool is_new = false;
            if (spec.empty()) {
                is_new = true;
            } else {
                const vector<int32_t>& last_degseq = spec.back().first;
                if (unit.degseq.size() != last_degseq.size()) {
                    is_new = true;
                } else {
                    bool is_identical = true;
                    for (int32_t i = 0; i < unit.degseq.size(); ++i) {
                        if (unit.degseq[i] != last_degseq[i]) {
                            is_identical = false;
                            break;
                        }
                    }

                    if (!is_identical) {
                        is_new = true;
                    }
                }
            }

            if (is_new) {
                spec.emplace_back(unit.degseq, static_cast<int32_t>(index.size()));
            }
    index_time += t.end();
            #endif

			index.push_back(unit);
        }

	} while(hasGraph);
	
	//READ color-label distribution
	//we assume that the order of graphs and color-label distribution are the same.
	ifstream in(aDistFile);
	while( in >> tag ) {
		if( tag == 't' ) {
			in >> id >> index[ind].cllen;
			index[ind].cldist = new pair<int32_t, int32_t>[index[ind].cllen];
			for(int32_t i = 0; i < index[ind].cllen; ++i) {
				in >> index[ind].cldist[i].first >> index[ind].cldist[i].second;
			}
			++ind;
		}
		else {
			cout << "INVALID FILE FORMAT " << aDistFile << endl;
			break;
		}
	}
	in.close();

	//READ COLORING & ONE
	in.open(aColorFile);
	int32_t* color = NULL;
	int32_t colorLen = 0;
	ind = 0;
	while( in >> tag ) {
		if( tag == 't' ) {
			int32_t n;
			in >> id >> n;

			if( colorLen < n ) {
				if(color != NULL) {
					delete[] color;
					color = NULL;
				}
				color = new int32_t[n];
				colorLen = n;
			}
			//READ COLORING
			for(int32_t i = 0; i < n; ++i)
				in >> color[i]; //color of ui

			//Since color[i] != coloring->color, we need to do the following job.
			index[ind].coloring = new Coloring(n);
			Coloring* coloring = index[ind].coloring;
			memset(coloring->cellSize, 0, sizeof(int32_t) * n);

			for(int32_t i = 0; i < n; ++i) {
				int32_t c = color[i];
				coloring->perm[c + coloring->cellSize[c]] = i;
				coloring->inv[i] = c + coloring->cellSize[c];
				coloring->color[c + coloring->cellSize[c]] = c;
				++(coloring->cellSize[c]);
			}
			coloring->numCell = 0;
			for(int32_t i = 0; i < n; i += coloring->cellSize[i]) {
				++(coloring->numCell);
			}

			//READ ONE
			int32_t ione = 0;
			for(int32_t i = 0; i < n; ++i) {
				in >> ione;
				index[ind].g->one[i] = (char)ione;
			}
			
			//Process one vertices
			for(int32_t i = 0; i < index[ind].g->numNode; ++i) {
				int32_t l = 0;
				int32_t r = index[ind].g->od[i] - 1;
				while( l <= r ) {
					int32_t adj = index[ind].g->e[i][l];
					if( index[ind].g->one[adj] ) {
						index[ind].g->e[i][l] = index[ind].g->e[i][r];
						index[ind].g->e[i][r] = adj;
						--r;
					}
					else {
						++l;
					}
				}
				if( index[ind].g->one[i] )
					index[ind].g->d[i] = -1;
				else
					index[ind].g->d[i] = r + 1;
			}

			++ind;
		}
		else {
			cout << "INVALID FILE FORMAT " << aColorFile << endl;
			break;
		}
	}
	in.close();

	if( color != NULL ) {
		delete[] color;
	}

    return index_time;
}

pair<int32_t, int32_t> binarySearch(vector<IndexUnit>& index, pair<int32_t,int32_t>* aDist, int32_t aLen, int32_t lo, int32_t hi)
{
	int32_t start = -1;
	int32_t end = -1;
	int32_t l = (lo != -1) ? lo : 0;
	int32_t h = (hi != -1) ? hi : index.size() - 1;

	// get the start index
	while( l <= h ) {
		int32_t mid = (h - l)/2 + l;
		if( index[mid].cllen > aLen ) {
			h = mid - 1;
		}
		else if( index[mid].cllen < aLen ) {
			l = mid + 1;
		}
		else { //cllen == aLen
			int i = 0;
			for(i = 0; i < aLen - 1; ++i) {
				if(index[mid].cldist[i].first != aDist[i].first || 
				   index[mid].cldist[i].second != aDist[i].second)
					break;
			}
			if(index[mid].cldist[i].first != aDist[i].first) {
				if(index[mid].cldist[i].first > aDist[i].first)
					h = mid - 1;
				else
					l = mid + 1;
			}
			else { //index.first == aDist.first
				if(index[mid].cldist[i].second > aDist[i].second)
					h = mid - 1;
				else if(index[mid].cldist[i].second < aDist[i].second)
					l = mid + 1;
				else { //index.first & second == aDist.first & second
					start = mid;
					h = mid - 1;
				}
			}
		}
	} //while

	//get the end index
	l = (lo != -1) ? lo : 0;
	h = (hi != -1) ? hi : index.size() - 1;
	while( l <= h ) {
		int32_t mid = (h - l)/2 + l;
		if( index[mid].cllen > aLen ) {
			h = mid - 1;
		}
		else if( index[mid].cllen < aLen ) {
			l = mid + 1;
		}
		else { //cllen == aLen
			int i = 0;
			for(i = 0; i < aLen - 1; ++i) {
				if(index[mid].cldist[i].first != aDist[i].first ||
				   index[mid].cldist[i].second != aDist[i].second)
				   break;
			}
			if(index[mid].cldist[i].first != aDist[i].first) {
				if(index[mid].cldist[i].first > aDist[i].first)
					h = mid - 1;
				else
					l = mid + 1;
			}
			else { //index.first == aDist.first
				if(index[mid].cldist[i].second > aDist[i].second)
					h = mid - 1;
				else if(index[mid].cldist[i].second < aDist[i].second)
					l = mid + 1;
				else {
					end = mid;
					l = mid + 1;
				}
			}
		}
	} //while

	if( start == -1 || end == -1 )
		return make_pair(0, -1);
	else
		return make_pair(start, end + 1);
}
