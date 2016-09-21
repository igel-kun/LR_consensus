#include <iterator>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <cassert>
#include "Hypergraph.hpp"
#include "db_map.hpp"

typedef map<Edge, MLeaf> DB;
DB db;

bool MLeaf::used(Vertex n) {
  set<Vertex>::iterator i = usedb.find(n);
  if(i == usedb.end()) return false;
  else return true;
}

void MLeaf::setused(Vertex n) {
  usedb.insert(n);
}

Leaf* lookup(const Edge &e) {
  DB::iterator i = db.find(e);
  if(i == db.end())
    i = db.insert(pair<Edge, MLeaf>(e, MLeaf())).first;
  return &i->second;
}

void initialize_db(const Hypergraph &G, const Graphstats &s) {
}
