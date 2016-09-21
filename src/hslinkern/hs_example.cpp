#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <fstream>
#include "hslinkern.hpp"

using namespace std;

void usage(const char *th) {
  cout << hslinkern_NAME <<", (C) 2011 Rene van Bevern <rene.vanbevern@tu-berlin.de>" << endl
       << "-------------------------------------------------------------------------" << endl;
  cout << "Usage: " << th << " [options] hypergraph_file.txt" << endl;
  cout << endl << "Valid options are" << endl;
  cout << "-k K: kernelize the instance asking for a hitting set of size K." << endl;
  cout << "      Per default, an upper bound on the minimum hitting set size" << endl;
  cout << "      is used as K." << endl;
  cout << "-o OUTPUT: write kernelized hypergraph into the file OUTPUT." << endl;
  cout << "           By default, no output is generated, only statistics." << endl;
  cout << "-s: find the optimum hitting set in the kernelized instance using" << endl;
  cout << "    an exponential-time exact search-tree algorithm, " << endl;
  cout << "      time O(d^k + |V| + |E|)" << endl;
  exit(1);
}    

int main(int argc, char *argv[]) {
  int c;
  bool searchtree = false;
  int k = -1;
  ofstream of;
  
  while((c=getopt(argc, argv, "k:o:s")) != -1)
    switch (c) {
    case 'k': {
      k = atoi(optarg);
      assert(k > 0);
      break;
    }
    case 'o': {
      of.open(optarg);
      break;
    }
    case 's': {
      searchtree = true;
      break;
    }
    default: usage(argv[0]);
    }
  
  if(optind >= argc) {
    cout << "Error: no filename given" << endl;
    usage(argv[0]);
  }

  cout << "Reading hypergraph " << argv[optind] << endl;
  Hypergraph G(read_hypergraph(argv[optind]));
  cout << "Graph statistics" <<endl;
  Graphstats stats(get_stats(G));
  print_stats(stats);

  edges_by_size(G, stats);
  Hypergraph apx(approx_hs(G, stats));
  

  int apxbound = get_stats(apx).edgesums;
  int greedybound = greedy_hs(G, stats);
  if(k < 0) {
    k = apxbound;
    if(greedybound < k) k = greedybound;
  }

  cout << "Upper bound HS     : " << k << endl;
  cout << "Lower bound HS     : " << apx.size() << endl;

    //    return 0;
    //  }
  Hypergraph H(kernelize(G, stats, k));
  Graphstats hstats(get_stats(H));
  cout << endl << "Kernel statistics" << endl;
  print_stats(hstats);

  if(searchtree) {
    vector<bool> hit(hstats.vertices+1);
    cout << endl << "... solving ..." << endl;
    cout << endl << "Optimal solution size: " << searchtree_hs(H,hstats,hit,k) << endl;
  }

  if(of.is_open()) {
    for(Hypergraph::iterator i = H.begin(); i != H.end(); ++i) {
      copy(i->begin(), i->end(), ostream_iterator<Vertex>(of, " "));
      of << endl;
    }
    of.close();
  }
  return 0;
}
