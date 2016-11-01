

#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/labeled_graph.hpp>
#include <iostream>
#include <utility>

using namespace std;
using namespace boost;

#define LEFT false
#define RIGHT true


typedef pair<int, bool> Label;

struct VertexProperty
{
  MyNode* corresponding_node;
  unsigned incident_non_singleton = 0;
};

struct EdgeProperty
{
  int white_weight = 0
  int red_weight = 0;
  int green_weight = 0;

//  EdgeProperty(const int w_weight, const int r_weight, const int g_weight):
//    white_weight(w_weight), red_weight(r_weight), green_weight(g_weight) {}
};


class MultiGraph
{	
private:
  typedef boost::adjacency_list< listS , vecS, undirectedS, VertexProperty, EdgeProperty > GraphType;
  typedef GraphType::vertex_descriptor Vertex;
	
	GraphType Graph;
  
  // map StIds to Vertices
  vector<Vertex> StId_to_Vertex;
  // a vector of vertices on the right side in order of addition to Graph
  list<Vertex> Rx;
public:
  const CentroidPath& centroid_path;
	
	MultiGraph(const CentroidPath& cp, const unsigned max_StId = 0):
    StId_to_Vertex(max_StId), centroid_path(cp),
  {
    for(const MyNode* node: cp.get_nodes_top_down())
      addVertex(node, RIGHT);
  }

	const GraphType& getGraph() const{
		return Graph;
	}
	
  const VertexProperty& operator[](const StId id) const
  {
    return Graph[StId_to_Vertex[id]];
  }
  const VertexProperty& operator[](const Vertex v) const
  {
    return Graph[v];
  }

	//! add a vertex to the graph on side LR
	void addVertex(const MyNode* v, const bool side){
    const StId v_id = stid(v);
    if(vertex >= StId_to_Vertex.size()) StId_to_Vertex.resize(vertex + 1);
    Vertex v_in_g = boost::add_vertex(VertexProperty{v, 0}, Graph);
		StId_to_Vertex[v_id] = v_in_g;
    if(side == RIGHT) Rx.push_back(v_in_g);
	}
	
  //! add an edge of weight "weight" between v1 (must be on the left!) and v2 (must be on the right!) 
  void add_edge(const StId v1, const StId v2, const int weight_w, const int weight_g, const int weight_r)
  {
    const Vertex& l(StId_to_Vertex[v1]);
    const Vertex& r(StId_to_Vertex[v2]);

    bool success = boost::add_edge(l, r, EdgeProperty{ 0, 0, 0 }, Graph).second;

    // update non-singleton-incident vertices
    if(success && (boost::degree(l, Graph) == 2))
      Graph[r].incident_non_singleton++;
  }

  void add_edge(const StId v1, const StId v2, const int weight)
  {
    add_edge(v1, v2, weight, weight, weight);
  }


  //! return the weight of an agreement matching containing only edges below (not strictly) x
  unsigned agreement_matching_below(const Vertex& u, const bool side)
  {
    if(boost::degree(u, Graph) == 1){
      // step 1: get adjacent vertex and find it in the search tree
      const Edge& uv = *boost::out_edges(u, Graph).first;
      const Vertex& v = *boost::target(uv, Graph);
      list<AM_ST_Node*> path;
      search_tree.lookup(v, path);

      // step 2: process white edge
      // step 2.1: get agreement matching with e as topmost edge
      const unsigned am_uv_topmost = agreement_matching_below(uv, rv_path);
      // step 2.2: update vertices on the root-v path in the search tree
      for(AM_ST_Node* z: path) z->m = am_uv_topmost;

      // step 3: process red edge
      unsigned max_g = 0;
      for(AM_ST_Node* z: path) {
        // update y
        update(max_g, z->g);
        update(z->y, max_g + z->r);
        // update g
        if(z->get_left()) update(z->get_left()->g, max_g);
        if(z->get_right()) update(z->get_right()->g, max_g);
        // clear g
        z->g = 0;
        // update r
        update(z->r, Graph[uv].red_weight);
      }

      //step 4: process green edge

    } else {
    }
#error continue here
  }

  //! return the weight of an agreement matching with e as topmost edge
  unsigned agreement_matching_below(const Edge& uv, const list<AM_ST_Node*>& rv_path)
  {
    // case1: the best matching contains another white edge: compute 1 + max{m(z) | z in lfringe(v)}
    // case2: the best matching contains only uv and a proper crossing
    //NOTE: slightly differing from the paper, we pull the max_{z in lfringe} out of the other max'es
    unsigned max_m = 0;
    unsigned max_crossing = 0;
    unsigned max_g_on_path = 0;
    for(const AM_ST_Node* w: rv_path){
      const AM_ST_Node* left = w.get_left();
      if(left) {
        update(max_m, left->m);
        // get max { g(z') | z' is a proper ancestor of left }
        update(max_g_on_path, w->g); //NOTE: this is w and not left!
        update(max_crossing, std::max({left->x, max_g_on_path + left->r, left->y});
      }
    }
    // TODO: research whether it's correct to add the white weight of uv to both of them
    // and what in the hell is the +1 about? The paper does not explain any of this...
    return Graph[uv].white_weight + std::max(max_m + 1, max_crossing);
  }

};	

// int main()
// {
//   using namespace boost;
//   typedef adjacency_list<multisetS , vecS, undirectedS, Vertex, Edge> Graph;
//                                          //^^^^^^^^^^^ <-- fix
//   
//   Graph g;
//   auto a= add_vertex(Vertex{ "A" }, g);
//   auto b= add_vertex(Vertex{ "B" }, g);
//   auto c= add_vertex(Vertex{ "C" }, g);
//   add_edge(a, b, Edge{ 10, "E1" }, g);
//   add_edge(a, b, Edge{ 10, "E2" }, g);
//   add_edge(a, b, Edge{ 10, "E3" }, g);
//   add_edge(a, c, Edge{ 10, "E4" }, g);
// 
//   // checking number of edges
//   std::cout<< num_edges(g)<< std::endl;
// 
//   // printing edges branching from A
//   auto erange= out_edges(a, g);
//   for(auto i= erange.first; i!= erange.second; ++ i)
//     std::cout<< g[*i].code<< std::endl;
// 
//   // now we want to iterate over edges that connect A and B
//   auto wtf= boost::edge_range(a, b, g);
//   for(auto i= wtf.first; i!= wtf.second; ++ i)
//     std::cout<< g[*i].code<< std::endl;
// }
