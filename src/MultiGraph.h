#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/labeled_graph.hpp>
#include <iostream>
#include <utility>

using namespace std;
using namespace boost;


struct VertexMultiGraph
{
  pair<int, bool> infos ;
//  int id;
//  bool rightSide;
};

struct EdgeMultiGraph
{
  double weight;
  std::string color;
};

typedef boost::labeled_graph<
    boost::adjacency_list< multisetS , vecS, undirectedS, VertexMultiGraph, EdgeMultiGraph >,
    pair<int, bool> > MultiGraphType;

typedef MultiGraphType::vertex_descriptor vertexType;

	//typedef boost::adjacency_list<multisetS , vecS, undirectedS, VertexMultiGraph, EdgeMultiGraph> MultiGraphType;
	vector <int> Rx;
	vector <int> Lx;


class MultiGraph
{	
	private:
	
	MultiGraphType * Graph;
	
	public:
	
	 MultiGraph(){new MultiGraphType(); }


	MultiGraphType * getGraph(){
		return Graph;
	}
	
	vector <int> getRxVertices(){
		return Rx;
	}
	
	vector <int> getLxVertices(){
		return Lx;
	}
	
	void setRxVertices (vector <int> vertices){
		Rx = vertices;
	}
	
	void addRxVertex (int vertex){
		Rx.push_back(vertex);
	}
	
	void setLxVertices (vector <int> vertices){
		Lx = vertices;
	}
	
	void addLxVertex (int vertex){
		Lx.push_back(vertex);
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