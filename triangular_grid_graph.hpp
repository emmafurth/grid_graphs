#ifndef TRIANGULAR_GRID_GRAPH_HPP
#define TRIANGULAR_GRID_GRAPH_HPP

#include <utility> 
#include <vector>
#include <unistd.h>
#include "boost/config.hpp" // put this first to suppress some VC++ warnings

#include <boost/utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/property_maps/null_property_map.hpp>
#include <boost/property_map/shared_array_property_map.hpp>
#include <boost/graph/property_maps/constant_property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>

enum edge_boundary_t { edge_boundary };
enum edge_in_two_factor_t { edge_in_two_factor };
enum vertex_in_two_factor_t { vertex_in_two_factor };

namespace boost {
  BOOST_INSTALL_PROPERTY(edge, boundary);
  BOOST_INSTALL_PROPERTY(edge, in_two_factor);
  BOOST_INSTALL_PROPERTY(vertex, in_two_factor);
}

using namespace std;
using namespace boost;

const double epsilon = 0.001;

/*template<typename edge_t>
struct Triple{
	edge_t first;
	edge_t second;
	edge_t third;
	Triple() {}
	Triple(edge_t e1, edge_t e2, edge_t e3) : first(e1), second(e2), third(e3) {}
};*/

/*template <typename Vertex, class Graph>
struct neighbor_of {
  neighbor_of() { }
  neighbor_of(Vertex _v, Graph _G) : v(_v), G(_G) { }
  template <typename vertex_t>
  bool operator()(const vertex_t& u) const {
    return edge(u, v, G).first;
  }
  Vertex v;
  Graph G;
};*/

/*template<class graph_t, class CoordMap>
void get_coord_grid(graph_t &G, CoordMap &coord_map){
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	typedef typename graph_traits<graph_t>::adjacency_iterator vertex_iter_t;
	
	map<vertex_t, bool> recorded;
	vertex_t v;
	vertex_iter_t vi, vi_end;
	v = *(vertices(G).first);
	for (tie(vi,vi_end) = adjacent_vertices;


}*/
template <typename bool_map_t>
class bool_writer {
public:
	bool_writer(bool_map_t _bmap) : bmap(_bmap) {}
	template <typename Edge>
	void operator()(std::ostream& out, const Edge& e) const {
		out << "[style=" << (bmap[e] ? "solid":"dashed") << "]";
		/*bool is_true;
		try {
			is_true = bmap.at(e);
		}
		catch (const std::out_of_range& oor) {
			is_true = false;
		}
		out << "[style=" << (is_true ? "solid":"dashed") << "]";*/
	}
private:
	bool_map_t bmap;
};
template <typename bool_map_t>
inline bool_writer<bool_map_t>
make_bool_writer(bool_map_t t) {
	return bool_writer<bool_map_t>(t);
}
template<typename edge_t, typename graph_t>
string edge_str(edge_t e, graph_t &g){
	stringstream s;
	s << "(" << get(vertex_name, g, source(e, g)) << "," << get(vertex_name, g, target(e, g)) << ")";
	return s.str();
}

template <typename two_factor_map_t>
struct two_factor_filter {
  two_factor_filter() { }
  two_factor_filter(two_factor_map_t tfmap) : m_tf(tfmap) { }
  template <typename Edge>
  bool operator()(const Edge& e) const {
  	//if( m_tf[e] >1)
  	//	cout << "m_tf[" << e << "] = " << m_tf[e] << endl;
    return m_tf[e];
  }
  two_factor_map_t m_tf;
};

// Returns true if e1 and e2 contain the same vertex
template<typename edge_t, typename graph_t>
bool share_vertex(edge_t e1, edge_t e2, graph_t &g){
	//typename graph_traits<graph_t>::vertex_descriptor v = graph_traits<graph_t>::null_vertex();	
	return (source(e1, g) == source(e2, g) || source(e1, g) == target(e2, g) || target(e1, g) == source(e2, g) || target(e1, g) == target(e2, g));
}

// Returns the vertex shared by e1 and e2. Returns graph_traits<graph_t>::null_vertex(); if they do not share a vertex.
template<typename edge_t, typename graph_t>
typename graph_traits<graph_t>::vertex_descriptor shared_vertex(edge_t e1, edge_t e2, graph_t &g){
	//typename graph_traits<graph_t>::vertex_descriptor v = graph_traits<graph_t>::null_vertex();	
	if (source(e1, g) == source(e2, g) || source(e1, g) == target(e2, g))
		return source(e1, g);
	if (target(e1, g) == source(e2, g) || target(e1, g) == target(e2, g))
		return target(e1, g);
	
	return graph_traits<graph_t>::null_vertex();	
}

// If e1 and e2 are part of a triangle (i.e.: a cycle of length 3), returns the third edge in that triangle.
// That is, if e1 = (u, v) and e2 = (v, w), and (u, w) is an edge in g, then this returns (u,w).
template<typename edge_t, typename graph_t>
optional<edge_t> triangle(edge_t e1, edge_t e2, graph_t &g){
	typename graph_traits<graph_t>::vertex_descriptor v = shared_vertex(e1, e2, g);
	if (v == graph_traits<graph_t>::null_vertex())
		return optional<edge_t>();
	typename graph_traits<graph_t>::vertex_descriptor v1, v2;
	v1 = (source(e1,g) == v ? target(e1,g) : source(e1, g));
	v2 = (source(e2,g) == v ? target(e2,g) : source(e2, g));
	
	edge_t e3;
	bool check;
	tie(e3, check) = edge(v1, v2, g);
	if (check)
		return optional<edge_t>(e3);
	
	return optional<edge_t>();
}

template<typename edge_t, typename graph_t>
bool is_parallelogram(edge_t e1, edge_t e2, graph_t &g){
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	/*pair<vertex_t, vertex_t> shared1, shared2;
	shared1 = shared_neighbors(e1, g);
	shared2 = shared_neighbors(e2, g);*/
	
	vertex_t v1 = source(e1, g);
	vertex_t v2 = target(e1, g);
	vertex_t v3 = source(e2, g);
	vertex_t v4 = target(e2, g);

	if (v1 == v3 || v1 == v4 || v2 == v3 || v2 == v4)
		return false;
	
	return (edge(v1, v3, g).second && edge(v1, v4, g).second && (edge(v2, v3, g).second || edge(v2, v4, g).second)) ||
								(edge(v2, v3, g).second && edge(v2, v4, g).second && (edge(v1, v3, g).second || edge(v1, v4, g).second));
}
/*template<typename edge_t, typename graph_t>
optional<Triple<edge_t> > parallelogram(edge_t e1, edge_t e2, graph_t &g){
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	
	vertex_t v1 = source(e1, g);
	vertex_t v2 = target(e1, g);
	vertex_t v3 = source(e2, g);
	vertex_t v4 = target(e2, g);

	if (v1 == v3 || v1 == v4 || v2 == v3 || v2 == v4)
		return optional<Triple<edge_t> >();
	Triple<edge_t> tmp;
	if (edge(v1, v3, g).second && edge(v1, v4, g).second){
		if (edge(v2, v3, g).second){
			tmp = Triple<edge_t>(edge(v1, v3, g).first, edge(v1, v4, g).first, edge(v2, v3, g).first);
			return optional<Triple<edge_t> >(tmp);
		}
		else if (edge(v2, v4, g).second)
			return optional<Triple<edge_t> >(Triple<edge_t>(edge(v1, v3, g).first, edge(v1, v4, g).first, edge(v2, v4, g).first));
		else
			return optional<Triple<edge_t> >();
	}
	else if (edge(v2, v3, g).second && edge(v2, v4, g).second){
		if (edge(v1, v3, g).second)
			return optional<Triple<edge_t> >(Triple<edge_t>(edge(v2, v3, g).first, edge(v2, v4, g).first, edge(v1, v3, g).first));
		else if (edge(v1, v4, g).second)
			return optional<Triple<edge_t> >(Triple<edge_t>(edge(v2, v3, g).first, edge(v2, v4, g).first, edge(v1, v4, g).first));
		else
			return optional<Triple<edge_t> >();
	}
	else
		return optional<Triple<edge_t> >();
}
*/

template<typename edge_t, typename graph_t>
vector<edge_t> parallelogram(edge_t e1, edge_t e2, graph_t &g){
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	
	vector<edge_t> result;
	
	vertex_t v1 = source(e1, g);
	vertex_t v2 = target(e1, g);
	vertex_t v3 = source(e2, g);
	vertex_t v4 = target(e2, g);

	if (v1 == v3 || v1 == v4 || v2 == v3 || v2 == v4)
		return result;
	
	
	if (edge(v1, v3, g).second && edge(v2, v4, g).second){
		result.push_back(e1);
		result.push_back(e2);
		result.push_back(edge(v1, v3, g).first);
		result.push_back(edge(v2, v4, g).first);
	}
	else if (edge(v1, v4, g).second && edge(v2, v3, g).second){
		result.push_back(e1);
		result.push_back(e2);
		result.push_back(edge(v1, v4, g).first);
		result.push_back(edge(v2, v3, g).first);
	}
	
	return result;
}
// Let e = (u,v). This returns the vertices that are in the intersection of N(u) and N(v). note that in 
// triangular grid graphs there are at most two of these. If there is only one of these, then the second element
// the pair will be graph_traits<graph_t>::null_vertex();	If there are none, the first will also be graph_traits<graph_t>::null_vertex();	
template<typename Edge, typename graph_t>
pair<typename graph_traits<graph_t>::vertex_descriptor, typename graph_traits<graph_t>::vertex_descriptor  > 
shared_neighbors(Edge e, graph_t &G) {
	typedef typename graph_traits<graph_t>::adjacency_iterator adj_iter_t;
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	vertex_t v1 = source(e, G);
	vertex_t v2 = target(e, G);
	vertex_t null_vertex = graph_traits<graph_t>::null_vertex();
	std::pair<vertex_t, vertex_t> result(null_vertex,null_vertex);
	adj_iter_t vi, vj, vi_end, vj_end;
	for (tie(vi, vi_end) = adjacent_vertices(v1, G); vi != vi_end; ++vi){
		if (*vi == v2)
			continue;		
		for (tie(vj, vj_end) = adjacent_vertices(v2, G); vj != vj_end; ++vj){  			
			if (*vj == v1)
				continue;
			if (*vi == *vj){
				if (result.first == null_vertex){
					result.first = *vi;
				}
				else if (*vi != result.first){
					result.second = *vi;
					// Since there are at most two intermediate vertices, you can return
					// as soon you as you find the second.
					return result;
				}
			}
		}
	}
	return result;
}
// Returns all edges (u, w) such that u and w are adjacent to v.
template<typename vertex_t, typename graph_t>
vector<typename graph_traits<graph_t>::edge_descriptor > edge_neighbors(vertex_t v, graph_t &g){
	typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
	vector<edge_t> e_neighbors;
	map<vertex_t, bool> in_vector;
	pair<vertex_t, vertex_t> shared;
	
	typename graph_traits<graph_t>::adjacency_iterator vi, vi_end;
	for (tie(vi, vi_end) = adjacent_vertices(v, g); vi != vi_end; ++vi){
		shared = shared_neighbors(edge(v, *vi, g).first, g);
		if (!in_vector[shared.first]){
			e_neighbors.push_back(edge(shared.first, *vi, g).first);
			in_vector[shared.first] = true;
		}
		if (shared.second != graph_traits<graph_t>::null_vertex() && !in_vector[shared.second]){
			e_neighbors.push_back(edge(shared.second, *vi, g).first);
			in_vector[shared.second] = true;
		}
	}
	return e_neighbors;
}

// Returns true if the subgraph induced by the neighbors of v is a locally connected graph.
template<typename vertex_t, typename graph_t> 
bool locally_connected(vertex_t v, graph_t &G){
	int n = num_vertices(G);
	if (n <= 1)
		return true;
		
	typedef typename graph_traits<graph_t>::out_edge_iterator edge_iter_t;
	typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
	vertex_t null_vertex = graph_traits<graph_t>::null_vertex();
	edge_iter_t ei, ei_end; 
	pair<vertex_t, vertex_t> shared;
	int i = 0;
	for (tie(ei, ei_end) = out_edges(v, G); ei != ei_end; ++ei){
		shared = shared_neighbors<edge_t, graph_t>(*ei, G);
		if (shared.first == null_vertex)
			return false;
		else if (shared.second == null_vertex)
			i++;
	}
	return (i <= 2);
}

// Returns true if g is polygonal. Does so by checking that all vertices are locally connected.
template<typename graph_t> 
bool polygonal(graph_t &g){

	std::vector<int> component(num_vertices(g));
  int num = connected_components(g, &component[0]);
  if (num > 1)
  	return false;
	typedef typename graph_traits<graph_t>::vertex_iterator vertex_iter_t;
	vertex_iter_t vi, vi_end;
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi){
		if (!locally_connected(*vi, g))
			return false;	
	}
	return true;
}

// Returns true if the graph is the Star of David graph.
template<typename graph_t> 
bool star_of_david(graph_t &g){
	if (num_vertices(g) != 13)
		return false;
		
	typename graph_traits<graph_t>::vertex_descriptor v;
	typename graph_traits<graph_t>::vertex_iterator vi, vi_end;
	
	int deg_2_count = 0;
	int deg_5_count = 0;
	int deg_6_count = 0;
	
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi){
		if (degree(*vi, g) == 6){
			if (deg_6_count >= 1)
				return false;
			v = *vi;
			deg_6_count++;
		}
		else if (degree(*vi, g) == 5){
			if (deg_5_count >= 6)
				return false;
			deg_5_count++;
		}
		else if (degree(*vi, g) == 2){
			if (deg_2_count >= 6)
				return false;
			deg_2_count++;
		}
		else
			return false;
	}
	if (deg_6_count != 1 || deg_5_count != 6 || deg_2_count != 6)
		return false;
	
	vector<typename graph_traits<graph_t>::edge_descriptor > c_neighbors = edge_neighbors(v, g);
	if (c_neighbors.size() != 6)
		return false;
	
	typename graph_traits<graph_t>::adjacency_iterator vj, vj_end;
	for(tie(vj, vj_end) = adjacent_vertices(v, g); vj != vj_end; ++vj){
		if (degree(*vj, g) != 5)
			return false;
	}
	return true;
}

// Returns a list of all vertices u such that 1) u is not in the cycle C, and 2) there 
// exist vertices v and w in C such that u, v, and w form a triangle. 
template<typename vertex_cycle_map_t, typename graph_t>
list<typename graph_traits<graph_t>::vertex_descriptor > cycle_extendors(vertex_cycle_map_t &vc_map, graph_t &g){
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	typedef typename graph_traits<graph_t>::vertex_iterator vertex_iter_t;
	//typedef typename graph_traits<graph_t>::out_edge_iterator edge_iter_t;
	typedef typename graph_traits<graph_t>::adjacency_iterator adj_iter_t;
	
	list<typename graph_traits<graph_t>::vertex_descriptor > extendors;
	
	vertex_iter_t vi, vi_end;
	adj_iter_t vj, vj_end;
	pair<vertex_t, vertex_t> shared;
	
	
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi){
		if (vc_map[*vi])
			continue;
		for (tie(vj, vj_end) = adjacent_vertices(*vi, g); vj != vj_end; ++vj){
			if (vc_map[*vj]){
				shared = shared_neighbors(edge(*vi, *vj, g).first, g);
				if (vc_map[shared.first] || vc_map[shared.second]){
					extendors.push_back(*vi);
					break;
				}	 
			}
		}
	}		
		
	return extendors;
}

// Takes v, and checks to see if there are any edges (u,w) in C such that u and w are both
// adjacent to v. If so, this removes (u,w) from C, adds (u,v) and (v, w), and returns true.
// Otherwise returns false.
template<typename vertex_t, typename vertex_cycle_map_t, typename edge_cycle_map_t, typename graph_t>
bool triangulate(vertex_t v, vertex_cycle_map_t &vc_map, edge_cycle_map_t &ec_map, graph_t &g){
	if (vc_map[v])
		return false;
	
	typedef typename graph_traits<graph_t>::adjacency_iterator adj_iter_t;
	
	adj_iter_t vi, vi_end;
	pair<vertex_t, vertex_t> shared;
	for (tie(vi, vi_end) = adjacent_vertices(v, g); vi != vi_end; ++vi){
		shared = shared_neighbors(edge(v, *vi, g).first, g);
		if (vc_map[shared.first] && ec_map[edge(*vi, shared.first, g).first]){
				vc_map[v] = true;
				ec_map[edge(*vi, shared.first, g).first] = false;
				ec_map[edge(*vi, v, g).first] = true;
				ec_map[edge(v, shared.first, g).first] = true;
				return true;
		}
		else if (shared.second != graph_traits<graph_t>::null_vertex() && vc_map[shared.second] && ec_map[edge(*vi, shared.second, g).first]){
				vc_map[v] = true;
				ec_map[edge(*vi, shared.second, g).first] = false;
				ec_map[edge(*vi, v, g).first] = true;
				ec_map[edge(v, shared.second, g).first] = true;
				return true;
		}
	}
	return false;
}

// Returns 60, 120, or 180 if (v_minus, v) and (v_plus, v) form a 60, 120, or 180 degree
// angle respectively. If none of the above are true (e.g.: if (v_minus, v) and/or 
// (v, v_plus) are not edges), returns -1. 
template<typename vertex_t, typename graph_t>
int angle(vertex_t &v_minus, vertex_t &v, vertex_t &v_plus, graph_t &g){
	typename graph_traits<graph_t>::edge_descriptor e1, e2;

	if (!edge(v_minus, v, g).second)
		return -1;
	if (!edge(v, v_plus, g).second)
		return -1;
	
	if (edge(v_minus, v_plus, g).second)
		return 60;
	
	
	e1 = edge(v_minus, v, g).first;
	e2 = edge(v, v_plus, g).first;
	
	pair<vertex_t, vertex_t> shared1, shared2;
	
	shared1 = shared_neighbors(e1, g);
	shared2 = shared_neighbors(e2, g);
	
	if (shared1.first == shared2.first || shared1.first == shared2.second || shared1.second == shared2.first || shared1.second == shared2.second)
		return 120;
	
	else
		return 180; 
}

// If v is in C, then this finds its incident edges in C and returns the degree of the
// angle they form (60, 120, or 180). Otherwise returns -1. 
template<typename vertex_t, typename edge_cycle_map_t, typename graph_t>
int angle(vertex_t &v, edge_cycle_map_t &ec_map, graph_t &g){
	vertex_t v_minus, v_plus;
	typename graph_traits<graph_t>::adjacency_iterator vi, vi_end;
	int i = 0;
	for (tie(vi, vi_end) = adjacent_vertices(v, g); vi != vi_end; ++vi){
		if (ec_map[edge(v, *vi, g).first]){
			if (i == 0){
				v_minus = *vi;
				i++;
			}
			else {
				v_plus = *vi; 
				return angle(v_minus, v, v_plus, g);
			}
		}
	}
	return -1;
}

// TODO: add better check for when cycle does not represent cycle. Also, ordering?
// If v is in C, then finds and returns the edges incident to v in C.
template<typename vertex_t, typename vertex_cycle_map_t, typename edge_cycle_map_t, typename graph_t>
optional<pair<typename graph_traits<graph_t>::edge_descriptor, typename graph_traits<graph_t>::edge_descriptor > > 
cycle_out_edges(vertex_t v, vertex_cycle_map_t &vc_map, edge_cycle_map_t &ec_map, graph_t &g){
	typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
	if (!vc_map[v])
		return optional<pair<edge_t, edge_t> >();
	
	typename graph_traits<graph_t>::out_edge_iterator ei, ei_end;
	typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
	edge_t e;
	int i = 0;
	for(tie(ei,ei_end) = out_edges(v, g); ei != ei_end; ++ei){
		if (ec_map[*ei]){
			if (i == 0){
				e = *ei;
				i++;
			}
			else 
				return optional<pair<edge_t, edge_t> >(make_pair(e, *ei));
		}
	}
	return optional<pair<edge_t, edge_t> >();
}

// Assuming e1 and e2 form a triangle, like this: 
//   o-----o
//   \\   //
// e1 \\ // e2
//      o    
// Checks to see if e1 and e2 are part of a bolt/hourglass; if so, returns vector of edges 
// in the bolt/hourglass.
template<typename edge_t, typename edge_cycle_map_t, typename graph_t>
vector<edge_t> bolt_or_hourglass_from_triangle(edge_t e1, edge_t e2, edge_cycle_map_t &ec_map, graph_t &g){
	vector<edge_t> b_or_h;
	if (!share_vertex(e1, e2, g) || !triangle(e1, e2, g).is_initialized())
		return b_or_h;
	
	typename graph_traits<graph_t>::vertex_descriptor v = shared_vertex(e1, e2, g);
	edge_t e = triangle(e1, e2, g).get();
	
	if (!(ec_map[e1] && ec_map[e2] && !ec_map[e]))
		return b_or_h;
	
	vector<edge_t> e_neighbors = edge_neighbors(v, g);
	for(typename vector<edge_t>::iterator ei = e_neighbors.begin(); ei != e_neighbors.end(); ++ei){
		if (ec_map[*ei]){
			b_or_h.push_back(e);
			b_or_h.push_back(*ei);
			if (share_vertex(e, *ei, g) && share_vertex(*ei, e1, g)){
				b_or_h.push_back(e2);
				b_or_h.push_back(triangle(e1, *ei, g).get());
				return b_or_h;
			}
			else if (share_vertex(e, *ei, g) && share_vertex(*ei, e2, g)){
				b_or_h.push_back(e1);
				b_or_h.push_back(triangle(e2, *ei, g).get());
				return b_or_h;
			}
			else {
				b_or_h.push_back(e1);
				b_or_h.push_back(e2);
				b_or_h.push_back(edge(source(*ei, g), v, g).first);
				b_or_h.push_back(edge(target(*ei, g), v, g).first);
				return b_or_h;
			}
		}
	}
	// Will be an empty vector if you get to this point.
	return b_or_h;
}

// Checks to see if e is part of a bolt/hourglass; if so, returns vector of edges in the 
// bolt/hourglass. For now, this only works if e is not in C.
template<typename edge_t, typename vertex_cycle_map_t, typename edge_cycle_map_t, typename graph_t>
vector<edge_t > bolt_or_hourglass(edge_t e, vertex_cycle_map_t &vc_map, edge_cycle_map_t &ec_map, graph_t &g){
	// For now, this only works for edges not in the cycle.
	if (ec_map[e])
		return vector<edge_t>();
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	vector<edge_t> b_or_h;
	vertex_t v1, v2, v;
	v1 = source(e, g);
	v2 = target(e, g);
	if (!vc_map[v1] || !vc_map[v2])
		return b_or_h;
	
	edge_t e1, e2;
	pair<vertex_t, vertex_t> shared = shared_neighbors(e, g);
	vector<edge_t> e_neighbors;
	if (vc_map[shared.first])
		v = shared.first;
	else if (vc_map[shared.second])
		v = shared.second;
	else
		return b_or_h;
	
	e1 = edge(v, v1, g).first;
	e2 = edge(v, v2, g).first;
	if (ec_map[e1] && ec_map[e2]){
		b_or_h = bolt_or_hourglass_from_triangle(e1, e2, ec_map, g);
		if (!b_or_h.empty())
			return b_or_h;
	}
	if (angle(v1, ec_map, g) == 60){
		tie(e1, e2) = cycle_out_edges(v2, vc_map, ec_map, g).get();
		optional<edge_t> e3 = triangle(e, e1, g);
		if (e3.is_initialized()){
			b_or_h.push_back(e);
			b_or_h.push_back(e1);
			b_or_h.push_back(e3.get());
			pair<edge_t, edge_t> c_out = cycle_out_edges(v1, vc_map, ec_map, g).get();
			b_or_h.push_back(c_out.first);
			b_or_h.push_back(c_out.second);
			b_or_h.push_back(triangle(c_out.first, c_out.second, g).get());
			return b_or_h;
		}
		e3 = triangle(e, e2, g);
		if (e3.is_initialized()){
			b_or_h.push_back(e);
			b_or_h.push_back(e2);
			b_or_h.push_back(e3.get());
			pair<edge_t, edge_t> c_out = cycle_out_edges(v1, vc_map, ec_map, g).get();
			b_or_h.push_back(c_out.first);
			b_or_h.push_back(c_out.second);
			b_or_h.push_back(triangle(c_out.first, c_out.second, g).get());
			return b_or_h;
		}
	}
	if (angle(v2, ec_map, g) == 60){
		tie(e1, e2) = cycle_out_edges(v1, vc_map, ec_map, g).get();
		optional<edge_t> e3 = triangle(e, e1, g);
		if (e3.is_initialized()){
			b_or_h.push_back(e);
			b_or_h.push_back(e1);
			b_or_h.push_back(e3.get());
			pair<edge_t, edge_t> c_out = cycle_out_edges(v2, vc_map, ec_map, g).get();
			b_or_h.push_back(c_out.first);
			b_or_h.push_back(c_out.second);
			b_or_h.push_back(triangle(c_out.first, c_out.second, g).get());
			return b_or_h;
		}
		e3 = triangle(e, e2, g);
		if (e3.is_initialized()){
			b_or_h.push_back(e);
			b_or_h.push_back(e2);
			b_or_h.push_back(e3.get());
			pair<edge_t, edge_t> c_out = cycle_out_edges(v2, vc_map, ec_map, g).get();
			b_or_h.push_back(c_out.first);
			b_or_h.push_back(c_out.second);
			b_or_h.push_back(triangle(c_out.first, c_out.second, g).get());
			return b_or_h;
		}
	}
	return b_or_h;
	
}

// Note that this function assumes the following configuration around e:
//         e
// o====o-----o====o
//  \  //\   /\\  /
//   \//  \ /  \\/
//    o    o    o
// Checks to see if e is part of an alternating kissing cycle containing at most two
// bolts/hourglasses. If so, returns vector of edges in the alternating kissing cycle. 
// (edges shared by the two bolts/hourglasses will be included twice; this will not affect
// performance of flip). For now, this only works if e is not in C.
template<typename edge_t, typename vertex_cycle_map_t, typename edge_cycle_map_t, typename graph_t>
vector<edge_t > double_bolt_or_hourglass(edge_t e, vertex_cycle_map_t &vc_map, edge_cycle_map_t &ec_map, graph_t &g){
	// For now, this only works for edges not in the cycle.
	if (ec_map[e])
		return vector<edge_t>();
	
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	vector<edge_t> b_or_h1, b_or_h2, b_or_h;
	vertex_t v1, v2, v;
	v1 = source(e, g);
	v2 = target(e, g);
	if (!vc_map[v1] || !vc_map[v2])
		return b_or_h;
	
	if (!(angle(v1, ec_map, g) == 60 && angle(v2, ec_map, g) == 60))
		return b_or_h;
		
	edge_t e1, e2, e3, e4;
	pair<vertex_t, vertex_t> shared = shared_neighbors(e, g);
	tie(e1, e2) = cycle_out_edges(v1, vc_map, ec_map, g).get();
	tie(e3, e4) = cycle_out_edges(v2, vc_map, ec_map, g).get();
	
	// Returns true if it looks something like this (let shared.first = v3) (e1 and e2 may swap places)
	//         v1   e    v2 
	//         	o-------o
	//         //\     /\\
	//     e1 //  \   /  \\ e2
	//       //    \ /    \\
	//      o       o       o
	//              v3
	bool check1 = ((triangle(e1, edge(v1, shared.first, g).first, g).is_initialized() || 
									triangle(e2, edge(v1, shared.first, g).first, g).is_initialized()) && 
									(triangle(e3, edge(v2, shared.first, g).first, g).is_initialized() || 
									triangle(e4, edge(v2, shared.first, g).first, g).is_initialized()));
	// Like check1, but v3 = shared.second			
	bool check2 = ((triangle(e1, edge(v1, shared.second, g).first, g).is_initialized() || 
									triangle(e2, edge(v1, shared.second, g).first, g).is_initialized()) && 
									(triangle(e3, edge(v2, shared.second, g).first, g).is_initialized() || 
									triangle(e4, edge(v2, shared.second, g).first, g).is_initialized()));
	if (vc_map[shared.first] && check1)
		v = shared.first;
	else if (vc_map[shared.second] && check2)
		v = shared.second;
	else
		return b_or_h;
		
	b_or_h1 = bolt_or_hourglass(edge(v, v1, g).first, vc_map, ec_map, g);
	if (!b_or_h1.empty()){
		flip(b_or_h1, ec_map, g);
		b_or_h2 = bolt_or_hourglass(e, vc_map, ec_map, g);
		if (!b_or_h2.empty()){
			flip(b_or_h1, ec_map, g);
			b_or_h.reserve(b_or_h1.size()+b_or_h2.size());
			b_or_h.insert( b_or_h.end(), b_or_h1.begin(), b_or_h1.end() );
			b_or_h.insert( b_or_h.end(), b_or_h2.begin(), b_or_h2.end() );
			return b_or_h;
		}
		flip(b_or_h1, ec_map, g);
		b_or_h1.clear();
	}
	b_or_h1 = bolt_or_hourglass(edge(v, v2, g).first, vc_map, ec_map, g);
	if (!b_or_h1.empty()){
		flip(b_or_h1, ec_map, g);
		b_or_h2 = bolt_or_hourglass(e, vc_map, ec_map, g);
		if (!b_or_h2.empty()){
			flip(b_or_h1, ec_map, g);
			b_or_h.reserve(b_or_h1.size()+b_or_h2.size());
			b_or_h.insert( b_or_h.end(), b_or_h1.begin(), b_or_h1.end() );
			b_or_h.insert( b_or_h.end(), b_or_h2.begin(), b_or_h2.end() );
			return b_or_h;
		}
		flip(b_or_h1, ec_map, g);
		b_or_h1.clear();
	}
	return b_or_h;
}

// Note that this function assumes the following configuration around e:
//    e
// o-----o
// \\   //
//  \\ // 
//    o    
// Checks to see if e is part of an alternating kissing cycle containing at most two
// bolts/hourglasses. If so, returns vector of edges in the alternating kissing cycle. 
// (edges shared by the two bolts/hourglasses will be included twice; this will not affect
// performance of flip). For now, this only works if e is not in C.
//
// TODO: Come up with better function name.
template<typename edge_t, typename vertex_cycle_map_t, typename edge_cycle_map_t, typename graph_t>
vector<edge_t > double_bolt_or_hourglass2(edge_t e, vertex_cycle_map_t &vc_map, edge_cycle_map_t &ec_map, graph_t &g){
	// For now, this only works for edges not in the cycle.
	if (ec_map[e])
		return vector<edge_t>();
		
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	vector<edge_t> b_or_h1, b_or_h2, b_or_h;
	
	// Check to make sure there isn't a single bolt or hourglass nearby.
	b_or_h = bolt_or_hourglass(e, vc_map, ec_map, g);
	if (!b_or_h.empty())
		return b_or_h;
	
	vertex_t v1, v2, v;
	v1 = source(e, g);
	v2 = target(e, g);
	if (!vc_map[v1] || !vc_map[v2])
		return b_or_h;
	
	edge_t e1, e2;
	pair<vertex_t, vertex_t> shared = shared_neighbors(e, g);
	if (ec_map[edge(shared.first, v1, g).first] && ec_map[edge(shared.first, v2, g).first]){
		v = shared.first;
		e1 = edge(shared.first, v1, g).first;
		e2 = edge(shared.first, v2, g).first;
	}
	else if (ec_map[edge(shared.second, v1, g).first] && ec_map[edge(shared.second, v2, g).first]){
		v = shared.second;
		e1 = edge(shared.second, v1, g).first;
		e2 = edge(shared.second, v2, g).first;
	}
	else
		return b_or_h;
		
	vector<edge_t> e_neighbors = edge_neighbors(v, g);
	for(typename vector<edge_t>::iterator ei = e_neighbors.begin(); ei != e_neighbors.end(); ++ei){
		b_or_h1 = bolt_or_hourglass(*ei, vc_map, ec_map, g);
		if (b_or_h1.empty())
			continue;
		
		flip(b_or_h1, ec_map, g);
		b_or_h2 = bolt_or_hourglass_from_triangle(e1, e2, ec_map, g);
		flip(b_or_h1, ec_map, g);
		if (!b_or_h2.empty()){
			b_or_h.reserve(b_or_h1.size()+b_or_h2.size());
			b_or_h.insert( b_or_h.end(), b_or_h1.begin(), b_or_h1.end() );
			b_or_h.insert( b_or_h.end(), b_or_h2.begin(), b_or_h2.end() );
			return b_or_h;
		}
	}
	return b_or_h;
}


// Note that this function assumes the following configuration around e:
//    e
// o-----o
// \\   //
//  \\ // 
//    o    
// Checks to see if e is part of an alternating kissing cycle containing an arbitrary
// number of bolts/hourglasses (see thesis). If so, returns vector of edges in the 
// alternating kissing cycle (edges shared by two bolts/hourglasses will be included twice;
// this will not affect performance of flip). For now, this only works if e is not in C.
template<typename edge_t, typename vertex_cycle_map_t, typename edge_cycle_map_t, typename graph_t>
bool find_linear_alt_kiss_cycle(edge_t e, vertex_cycle_map_t &vc_map, edge_cycle_map_t &ec_map, graph_t &g){
	// For now, this only works for edges not in the cycle.
	if (ec_map[e])
		return false;
		
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	vector<edge_t> double_b_or_h;
	vertex_t v1, v2, v;
	v1 = source(e, g);
	v2 = target(e, g);
	if (!vc_map[v1] || !vc_map[v2])
		return false;
	
	edge_t e1, e2;
	pair<vertex_t, vertex_t> shared = shared_neighbors(e, g);
	if (ec_map[edge(shared.first, v1, g).first] && ec_map[edge(shared.first, v2, g).first]){
		v = shared.first;
		e1 = edge(shared.first, v1, g).first;
		e2 = edge(shared.first, v2, g).first;
	}
	else if (ec_map[edge(shared.second, v1, g).first] && ec_map[edge(shared.second, v2, g).first]){
		v = shared.second;
		e1 = edge(shared.second, v1, g).first;
		e2 = edge(shared.second, v2, g).first;
	}
	else
		return false;
		
	vector<vertex_t> vc_vector;
	vector<edge_t> ec_vector;
	cycle_as_vector(vc_vector, ec_vector, vc_map, ec_map, g, v);
	
	optional<pair<edge_t, edge_t> > es;
	optional<edge_t> opt_edge;
	int i = -2;
	
	do{
		i = i+2;
		es = cycle_out_edges(vc_vector[i], vc_map, ec_map, g);
		assert(es.is_initialized());
		
		opt_edge = triangle(es.get().first, es.get().second, g);
		assert(opt_edge.is_initialized());
		e = opt_edge.get();
		
		double_b_or_h = double_bolt_or_hourglass2(e, vc_map, ec_map, g);
	}	while(i < vc_vector.size() && double_b_or_h.empty());
	if (double_b_or_h.empty())
		return false;
	v = vc_vector[i];
	
	do {
		es = cycle_out_edges(v, vc_map, ec_map, g);
		assert(es.is_initialized());
	
		opt_edge = triangle(es.get().first, es.get().second, g);
		assert(opt_edge.is_initialized());
		e = opt_edge.get();
		
		double_b_or_h = double_bolt_or_hourglass2(e, vc_map, ec_map, g);
		
		flip(double_b_or_h, ec_map, g);
		shared = shared_neighbors(e, g);
		v = (shared.first == v ? shared.second : shared.first);
	} while (vc_map[v]);
	bool test = triangulate(v, vc_map, ec_map, g);
	return test;
}

// Returns true if the edges in alt_kiss_cycle are an alternating kissing cycle.
template<typename edge_t, typename edge_cycle_map_t, typename graph_t>
bool is_alt_kiss_cycle(vector<edge_t> &alt_kiss_cycle, edge_cycle_map_t &ec_map, graph_t &g){
	typedef adjacency_list<vecS, vecS, undirectedS, 
												property<vertex_in_two_factor_t, bool>, 
												property<edge_in_two_factor_t, bool> 
												> Graph;
	
	map<typename graph_traits<graph_t>::vertex_descriptor, bool> in_new_graph;
	map<typename graph_traits<graph_t>::vertex_descriptor, typename graph_traits<Graph>::vertex_descriptor> counterpart;
	Graph G;
	typename graph_traits<Graph>::vertex_descriptor v1, v2;
	for(typename vector<edge_t>::iterator ei = alt_kiss_cycle.begin(); ei != alt_kiss_cycle.end(); ++ei){
		if (!in_new_graph[source(*ei, g)]){
			v1 = add_vertex(ec_map[*ei], G);
			in_new_graph[source(*ei, g)] = true;
			counterpart[source(*ei, g)] = v1;
		}
		else {
			v1 = counterpart[source(*ei, g)];
			if (!get(vertex_in_two_factor, G, v1) && ec_map[*ei])
				put(vertex_in_two_factor, G, v1, true);
		}
		if (!in_new_graph[target(*ei, g)]){
			v2 = add_vertex(ec_map[*ei], G);
			in_new_graph[target(*ei, g)] = true;
			counterpart[target(*ei, g)] = v2;
		}
		else {
			v2 = counterpart[target(*ei, g)];
			if (!get(vertex_in_two_factor, G, v2) && ec_map[*ei])
				put(vertex_in_two_factor, G, v2, true);
		}
		add_edge(v1, v2, ec_map[*ei], G);
	}
	
	typename graph_traits<Graph>::vertex_iterator vi, vi_end;
	for(tie(vi, vi_end) = vertices(G); vi != vi_end; ++vi){
		if (!get(vertex_in_two_factor, G, *vi))
			return false;
		if (degree(*vi, G) == 2 || degree(*vi, G) == 4){
			int in_tf_count = 0;
			int not_in_tf_count = 0;
			typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
			for (tie(ei, ei_end) = out_edges(*vi, G); ei != ei_end; ++ei){
				if (get(edge_in_two_factor, G, *ei))
					in_tf_count++;
				else
					not_in_tf_count++;
			}
			if (!(in_tf_count == not_in_tf_count && in_tf_count == (int)(degree(*vi, G)/2)))
				return false;
		}
	}
	return true;
}

// Flips alt_kiss_cycle. That is, if some edge e is in C and alt_kiss_cycle, then it will 
// be removed from C. If e is not in C but is in alt_kiss_cycle, then it will be added to
// C. If e is in alt_kiss_cycle an even number of times, then it will be unchanged.
// If e is in alt_kiss_cycle an odd number of times, then it will be changed as if it were 
// only included once. 
template<typename edge_t, typename edge_cycle_map_t, typename graph_t>
void flip(vector<edge_t> &alt_kiss_cycle, edge_cycle_map_t &ec_map, graph_t &g){
	for(typename vector<edge_t>::iterator ei = alt_kiss_cycle.begin(); ei != alt_kiss_cycle.end(); ++ei)
		ec_map[*ei] = !ec_map[*ei];
}

// Takes the property maps indicating whether a particular vertex/edge is in C, and 
// populates vc_vector and ec_vector with the vertices/edges in C in order.
template<typename vertex_t, typename edge_t, typename vertex_cycle_map_t, typename edge_cycle_map_t, typename graph_t>
void cycle_as_vector(vector<vertex_t>& vc_vector, vector<edge_t>& ec_vector, vertex_cycle_map_t &vc_map, edge_cycle_map_t &ec_map, graph_t &g, vertex_t v_first = typename graph_traits<graph_t>::null_vertex()){
	vc_vector.clear();
	ec_vector.clear();
	
	if (v_first == graph_traits<graph_t>::null_vertex()){
		typename graph_traits<graph_t>::vertex_iterator vi, vi_end;
		tie(vi, vi_end) = vertices(g);
		while (!vc_map[*vi] && vi != vi_end)
			++vi;
	
		if (!vc_map[*vi])
			return;
	
		v_first = *vi;
	}
	else if(!vc_map[v_first])
		return;
		
	edge_t e;
	
	vc_vector.push_back(v_first); 
	
	optional<pair<edge_t, edge_t> > es = cycle_out_edges(v_first, vc_map, ec_map, g); 
	assert(es.is_initialized());
	e = es.get().first;
	
	ec_vector.push_back(e);
	vertex_t v = (source(e, g) == v_first ? target(e, g) : source(e, g));
	
	while (v != v_first){
				
		vc_vector.push_back(v); 
		es = cycle_out_edges(v, vc_map, ec_map, g);
		assert(es.is_initialized()); 

		e = (e == es.get().first ? es.get().second : es.get().first);
		v = (source(e, g) == v ? target(e, g) : source(e, g));
		ec_vector.push_back(e);
	}
}

// If a cycle is non-Hamiltonian, extends it to a Hamiltonian cycle
// (only works if g is polygonal and not the Star of David)
template<typename vertex_cycle_map_t, typename edge_cycle_map_t, typename graph_t>
void extend_cycle(vertex_cycle_map_t &vc_map, edge_cycle_map_t &ec_map, graph_t &g){
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
	typename graph_traits<graph_t>::adjacency_iterator vj, vj_end;
	
	list<vertex_t> extendors = cycle_extendors(vc_map, g);
	pair<vertex_t, vertex_t> shared;
	typename list<vertex_t>::iterator vi;
	vector<edge_t> e_neighbors, b_or_h;
	bool test;
	while (!extendors.empty()){
		vi = extendors.begin();
		// This while loop iterates through the extendors and extends to all the vertices that
		// can be "simply triangulated." That is, they could not
		// 		1) be triangulated immediately
		//		2) be triangulated with one bolt/hourglass flip
		// Extendors requiring more flips are skipped for later.
		while(vi != extendors.end()){
			// Extendors may have duplicates, so this check is necessary.  
			if (vc_map[*vi]){
				vi = extendors.erase(vi);
				continue;
			}
			test = triangulate(*vi, vc_map, ec_map, g);
			if (!test){
				e_neighbors = edge_neighbors(*vi, g);
				for (typename vector<edge_t>::iterator ei = e_neighbors.begin(); ei != e_neighbors.end(); ++ei){
					b_or_h = bolt_or_hourglass(*ei, vc_map, ec_map, g);
					if (!b_or_h.empty()){
						flip(b_or_h, ec_map, g);
						b_or_h.clear();
						test = triangulate(*vi, vc_map, ec_map, g);
						break;
					}
				}
				if (!test){
					++vi;
					continue;
				}
			}
			// Adds neighbors of vi to extendors as necessary. May add extendors that are already in there.
			// TODO: Move this boiler plate to seperate function
			for (tie(vj, vj_end) = adjacent_vertices(*vi, g); vj != vj_end; ++vj){
				if (!vc_map[*vj]){
					shared = shared_neighbors(edge(*vi, *vj, g).first, g);
					if (vc_map[shared.first] || vc_map[shared.second]){
						extendors.push_back(*vj);
					}
				
				}
			}
			vi = extendors.erase(vi);
		}
		test = false;
		vi = extendors.begin();
		// This while loop iterates through all the extendors that could not be simply triangulated,
		// and looks for either a "double" bolt/hourglass, or a linear alternating kissing cycle.
		while(vi != extendors.end()){ 
			// Extendors may have duplicates, so this check is necessary.  
			if (vc_map[*vi]){
				vi = extendors.erase(vi);
				continue;
			}
			e_neighbors = edge_neighbors(*vi, g);
			for (typename vector<edge_t>::iterator ei = e_neighbors.begin(); ei != e_neighbors.end(); ++ei){
				b_or_h = double_bolt_or_hourglass(*ei, vc_map, ec_map, g);
				if (!b_or_h.empty()){
					flip(b_or_h, ec_map, g);
					b_or_h.clear();
					test = triangulate(*vi, vc_map, ec_map, g);
					break;
				}
				test = find_linear_alt_kiss_cycle(*ei, vc_map, ec_map, g);
				if (test){
					break;
				}
			}
			if (test){
				// Adds neighbors of vi to extendors as necessary. May add extendors that are already in there.
				for (tie(vj, vj_end) = adjacent_vertices(*vi, g); vj != vj_end; ++vj){
					if (!vc_map[*vj]){
						shared = shared_neighbors(edge(*vi, *vj, g).first, g);
						if (vc_map[shared.first] || vc_map[shared.second]){
							extendors.push_back(*vj);
						}
			
					}
				}
				vi = extendors.erase(vi);
				break;
			}
			else
				++vi;
		}
	}
}


// (x1,y1) = x*(1,0) + y*(0.5, sqrt(3)/2)
pair<float, float> convert_coord(int x, int y){
	float x1 = x+((float)y)/2.0;
	float y1 = y*sqrt(3.0)/2.0;
	return make_pair(x1, y1);
}

pair<float, float> convert_coord(pair<int, int> coord){
	return convert_coord(coord.first, coord.second);
}

float coord_distance(float x1, float y1, float x2, float y2){
	float result = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	// If the distance is less then epsilon away from an integer, then round it. 
	if (abs(floor(result+0.5)-result) < epsilon)
		return floor(result+0.5);
	else 
		return result;
}
float coord_distance(pair<float, float> coord1, pair<float, float> coord2){
	return coord_distance(coord1.first, coord1.second, coord2.first, coord2.second);
}

float coord_distance(pair<int, int> coord1, pair<int, int> coord2){
	pair<float, float> converted_coord1 = convert_coord(coord1);
	pair<float, float> converted_coord2 = convert_coord(coord2);
	return coord_distance(converted_coord1, converted_coord2);
}

vector<pair<int, int> > get_random_coords(int x, int y){
	int num_vertices = rand()%(x*y/3)+(x*y/3);
	//cout << "Num vertices: " << num_vertices << endl;
	vector<pair<int, int> > coords;
	
	int already_picked[x*y];
	int a, b, c;
	for(int i=0; i < num_vertices; i++){
		do{
			a = rand()%(x*y);
		}while(already_picked[a] == 1);
		b = a % x;
		c = a / x;
		//cout << "Picked (" << b << "," << c << ")\n";
		already_picked[a] = 1;
		coords.push_back(make_pair(b, c));
	}
	return coords;
}

template<typename graph_t>
graph_t get_triangular_grid_graph_from_coords(vector<pair<int, int> > coords){
	graph_t g;
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	map<vertex_t, pair<int, int> > vertex_to_coord;
	vertex_t v;
	string name;
	for (int i=0; i< coords.size(); i++){
	//for (vector<pair<int, int> >::iterator coord = coords.begin() ; coord != coords.end(); ++coord){
		name = "(" + to_string(coords[i].first) + "," + to_string(coords[i].second) + ")";
		v = add_vertex(name, g);
		//cout << "Adding " << name << endl;
		//put(vertex_name, g, v, name);
		vertex_to_coord[v] = coords[i];
	}
	typename graph_traits<graph_t>::vertex_iterator vi, vj, vi_end, vj_end;
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi){
		for(tie(vj, vj_end) = vertices(g); vj != vj_end; ++vj){
			//cout << "Distance between (" << vertex_to_coord[*vi].first << "," << vertex_to_coord[*vi].second << ") and (" << vertex_to_coord[*vj].first << "," << vertex_to_coord[*vj].second << ") is " << coord_distance(vertex_to_coord[*vi], vertex_to_coord[*vj]) << endl;
			if (coord_distance(vertex_to_coord[*vi], vertex_to_coord[*vj]) == 1.0 && !edge(*vi, *vj, g).second){
				//cout << "Adding (" << vertex_to_coord[*vi].first << "," << vertex_to_coord[*vi].second << ") and (" << vertex_to_coord[*vj].first << "," << vertex_to_coord[*vj].second << ")\n";
				add_edge(*vi, *vj, g);
			}
		}
	}
	return g;
}

template<typename graph_t>
graph_t random_triangular_grid_graph(int x, int y){
	vector<pair<int, int> > coords;
	graph_t g;
	int i = 0;
	do {
		cout << "Attempt " << i << endl;
		coords = get_random_coords(x, y);
		cout << "Got coords\n";
		g = get_triangular_grid_graph_from_coords<graph_t>(coords);
		i++;
	}while(!polygonal(g));
	return g;
}

// Sets a property map indicating whether a given edge is a boundary edge or not.
template<typename boundary_vertex_map_t, typename boundary_edge_map_t, typename graph_t>
void boundary_graph(boundary_vertex_map_t &vb_map, boundary_edge_map_t &eb_map, graph_t &g){
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	vertex_t null_vertex = graph_traits<graph_t>::null_vertex();
	typedef typename graph_traits<graph_t>::edge_iterator edge_iter_t;
	edge_iter_t ei, ei_end; 
	pair<vertex_t, vertex_t> shared;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
		shared = shared_neighbors(*ei, g);
		eb_map[*ei] = (shared.first == null_vertex || shared.second == null_vertex);
		if (eb_map[*ei]){
			vb_map[source(*ei, g)] = true;
			vb_map[target(*ei, g)] = true;
		}
	}			
}

template<typename boundary_vertex_map_t, typename boundary_edge_map_t, typename graph_t>
void boundary_based_two_factor(boundary_vertex_map_t &vb_map, boundary_edge_map_t &eb_map, graph_t &g){
	boundary_graph(vb_map, eb_map, g);
	extend_cycle(vb_map, eb_map, g);	
}

// Determines the number of components in the twofactor

template<typename tf_edge_map_t, typename graph_t>
int two_factor_components(tf_edge_map_t &etf_map, graph_t &g){
	//two_factor_map_t tfmap = get(edge_in_two_factor, G);
	two_factor_filter<tf_edge_map_t> filter(etf_map);
	filtered_graph<graph_t, two_factor_filter<tf_edge_map_t> >
		fg(g, filter);
	
	std::vector<int> component(num_vertices(fg));
	return connected_components(fg, &component[0]);
}

template<typename tf_edge_map_t, typename graph_t>
int two_factor_components(tf_edge_map_t &etf_map, graph_t &g, vector<int> &component){
	//two_factor_map_t tfmap = get(edge_in_two_factor, G);
	two_factor_filter<tf_edge_map_t> filter(etf_map);
	filtered_graph<graph_t, two_factor_filter<tf_edge_map_t> >
		fg(g, filter);
	
	component.clear();
	component.resize(num_vertices(fg));
	return connected_components(fg, &component[0]);
}

template<typename tf_vertex_map_t, typename tf_edge_map_t, typename graph_t>
list<pair<typename graph_traits<graph_t>::edge_descriptor, typename graph_traits<graph_t>::edge_descriptor> > 
get_boundary_parallelograms(tf_vertex_map_t &vtf_map, tf_edge_map_t &etf_map, graph_t &g){
	typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
	vector<int> component;
	int n = two_factor_components(etf_map, g, component);
	list<pair<edge_t, edge_t> > parallelograms;
	if (n <= 1)
		return parallelograms;
	
	typedef typename graph_traits<graph_t>::edge_iterator edge_iter_t;
	edge_iter_t ei, ei_end, ej, ej_end; 
	
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
		for (tie(ej, ej_end) = edges(g); ej != ei; ++ej){
			if (is_parallelogram(*ei, *ej, g) && 
					component[source(*ei, g)] == component[target(*ei, g)] &&
					component[source(*ej, g)] == component[target(*ej, g)] &&
					component[source(*ei, g)] != component[source(*ej, g)])
				parallelograms.push_back(make_pair(*ei, *ej));
		}
	}
	return parallelograms;
}

template<typename component_graph_t, typename tf_vertex_map_t, typename tf_edge_map_t, typename graph_t>
component_graph_t get_component_graph(tf_vertex_map_t &vtf_map, tf_edge_map_t &etf_map, graph_t &g){
	//typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_name_t, int> > component_graph_t;
	typedef typename graph_traits<component_graph_t>::vertex_descriptor component_vertex_t;
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
	
	vector<int> component;
	int k = two_factor_components(etf_map, g, component);
	cout << "k = " << k << endl;
	component_graph_t cg;
	vector<component_vertex_t> component_id_to_vertex(k);
	
	component_vertex_t v;
	for(int i=0; i<k; i++){
		v = add_vertex(i, cg);
		cout << "Added " << v << " corresponding to " << i << endl;
		component_id_to_vertex[i] = v;
	}
	cout << "Done!\n";
	list<pair<edge_t, edge_t> > pars = get_boundary_parallelograms(vtf_map, etf_map, g);
	cout << "Got boundary parallelograms.\n";
	int x, y;
	int i = 0;
	for(typename list<pair<edge_t, edge_t> >::iterator pi = pars.begin(); pi != pars.end(); ++pi){
		x = component[source((*pi).first, g)];
		y = component[source((*pi).second, g)];
		if (!edge(component_id_to_vertex[x], component_id_to_vertex[y], cg).second){
			cout << "Adding edge (" << x << "," << y << ")\n";
			add_edge(component_id_to_vertex[x], component_id_to_vertex[y], i, cg);
			i++;
		}
	}
	return cg;
}

template<typename edge_t, typename tf_vertex_map_t, typename tf_edge_map_t, typename graph_t>
bool merge_at_boundary_parallelogram(edge_t &e1, edge_t &e2, tf_vertex_map_t &vtf_map, tf_edge_map_t &etf_map, graph_t &g){
	if (!is_parallelogram(e1, e2, g))
		return false;
	
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
		
	vector<edge_t> alt_kiss_cycle;
	
	if (!etf_map[e1]){
		alt_kiss_cycle = bolt_or_hourglass(e1, vtf_map, etf_map, g);
		if (alt_kiss_cycle.empty()){
			optional<pair<edge_t, edge_t> > es = cycle_out_edges(source(e1, g), vtf_map, etf_map, g);
			assert(es.is_initialized());
			if (is_parallelogram(e2, es.get().first, g))
				e1 = es.get().first;
			else if (is_parallelogram(e2, es.get().second, g))
				e1 = es.get().second;
			else {
				es = cycle_out_edges(target(e1, g), vtf_map, etf_map, g);
				assert(es.is_initialized());
				if (is_parallelogram(e2, es.get().first, g))
					e1 = es.get().first;
				else if (is_parallelogram(e2, es.get().second, g))
					e1 = es.get().second;
				else
					return false;
			}
		}
		else
			flip(alt_kiss_cycle, etf_map, g);
	}
	
	if (!etf_map[e2]){
		alt_kiss_cycle = bolt_or_hourglass(e2, vtf_map, etf_map, g);
		if (alt_kiss_cycle.empty()){
			optional<pair<edge_t, edge_t> > es = cycle_out_edges(source(e2, g), vtf_map, etf_map, g);
			assert(es.is_initialized());
			if (is_parallelogram(e1, es.get().first, g))
				e2 = es.get().first;
			else if (is_parallelogram(e1, es.get().second, g))
				e2 = es.get().second;
			else {
				es = cycle_out_edges(target(e2, g), vtf_map, etf_map, g);
				assert(es.is_initialized());
				if (is_parallelogram(e1, es.get().first, g))
					e2 = es.get().first;
				else if (is_parallelogram(e1, es.get().second, g))
					e2 = es.get().second;
				else
					return false;
			}
		}
		else
			flip(alt_kiss_cycle, etf_map, g);
	}
	
	if (!(etf_map[e1] && etf_map[e2]))
		return false;
	
	alt_kiss_cycle = parallelogram(e1, e2, g);
	assert(!alt_kiss_cycle.empty());
	flip(alt_kiss_cycle, etf_map, g);
	return true;
}

template<typename tf_vertex_map_t, typename tf_edge_map_t, typename graph_t>
bool merge_all_components(tf_vertex_map_t &vtf_map, tf_edge_map_t &etf_map, graph_t &g){
	
	vector<int> component;
	int k = two_factor_components(etf_map, g, component);
	if (k == 1)
		return true;
	if (k == 0)
		return false;
	
	typedef typename graph_traits<graph_t>::edge_descriptor edge_t;
	typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, int, 
																									property<vertex_name_t, int> >,
																									property<edge_index_t, size_t> 
																									> c_graph_t;
																									
	c_graph_t cg;
	vector<graph_traits<c_graph_t>::vertex_descriptor> component_id_to_vertex(k);
	
	graph_traits<c_graph_t>::vertex_descriptor v;
	for(int i=0; i<k; i++){
		v = add_vertex(i, cg);
		put(vertex_name, cg, v, i);
		component_id_to_vertex[i] = v;
	}
	
	list<pair<edge_t, edge_t> > pars = get_boundary_parallelograms(vtf_map, etf_map, g);

	int x, y;
	int i = 0;
	
	//shared_array_property_map<bool, property_map<c_graph_t, edge_index_t>::const_type> is_set(num_vertices(cg), get(edge_index, cg));
	
	
	//shared_array_property_map<string, property_map<c_graph_t, edge_index_t> > par_map(num_edges(cg), get(edge_index, cg));
	vector_property_map<pair<edge_t, edge_t> > par_map;
	
	graph_traits<c_graph_t>::edge_descriptor e;
	pair<edge_t, edge_t> p;
	for(typename list<pair<edge_t, edge_t> >::iterator pi = pars.begin(); pi != pars.end(); ++pi){
		x = component[source((*pi).first, g)];
		y = component[source((*pi).second, g)];
		if (!edge(component_id_to_vertex[x], component_id_to_vertex[y], cg).second){
			e = add_edge(component_id_to_vertex[x], component_id_to_vertex[y], i, cg).first;
			par_map[get(edge_index, cg, e)] = *pi;
			i++;
		}
		else{
			e = edge(component_id_to_vertex[x], component_id_to_vertex[y], cg).first;
			p = par_map[get(edge_index, cg, e)];
			if ((!etf_map[p.first] && !etf_map[p.second]) && (etf_map[(*pi).first] || etf_map[(*pi).second]))
				par_map[get(edge_index, cg, e)] = *pi;
			else if ((!etf_map[p.first] || !etf_map[p.second]) && (etf_map[(*pi).first] && etf_map[(*pi).second]))
				par_map[get(edge_index, cg, e)] = *pi;
		}
	}
	
	shared_array_property_map<graph_traits<c_graph_t>::vertex_descriptor, property_map<c_graph_t, vertex_index_t>::const_type> pred(num_vertices(cg), get(vertex_index, cg));
	mt19937 random_generator;		
	random_spanning_tree(cg, random_generator, predecessor_map(pred));
	
	graph_traits<c_graph_t>::edge_iterator ei, ei_end;
	for(tie(ei,ei_end)=edges(cg); ei!=ei_end; ++ei){
		if (!(get(pred, source(*ei, cg)) == target(*ei, cg) || get(pred, target(*ei, cg)) == source(*ei, cg)))
			continue;
		
		p = par_map[get(edge_index, cg, *ei)];
		cout << "Merging at " << edge_str(p.first, g) << " and " << edge_str(p.second, g) << endl;
		if (!merge_at_boundary_parallelogram(p.first, p.second, vtf_map, etf_map, g)){
			cout << "Failed at " << edge_str(p.first, g) << " and " << edge_str(p.second, g) << endl;
			return false;
		}
	}
	
	return true;
}
#endif
