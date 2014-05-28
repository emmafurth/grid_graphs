#ifndef TRIANGULAR_GRID_GRAPH_HPP
#define TRIANGULAR_GRID_GRAPH_HPP

#include <utility> 
#include <vector>
#include <unistd.h>
#include "boost/config.hpp" // put this first to suppress some VC++ warnings

#include <boost/utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/property_maps/null_property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>

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
bool share_vertex(edge_t e1, edge_t e2, graph_t &g){
	//typename graph_traits<graph_t>::vertex_descriptor v = graph_traits<graph_t>::null_vertex();	
	return (source(e1, g) == source(e2, g) || source(e1, g) == target(e2, g) || target(e1, g) == source(e2, g) || target(e1, g) == target(e2, g));
}

template<typename edge_t, typename graph_t>
typename graph_traits<graph_t>::vertex_descriptor shared_vertex(edge_t e1, edge_t e2, graph_t &g){
	//typename graph_traits<graph_t>::vertex_descriptor v = graph_traits<graph_t>::null_vertex();	
	if (source(e1, g) == source(e2, g) || source(e1, g) == target(e2, g))
		return source(e1, g);
	if (target(e1, g) == source(e2, g) || target(e1, g) == target(e2, g))
		return target(e1, g);
	
	return graph_traits<graph_t>::null_vertex();	
}

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

template<typename Edge, typename graph_t>
pair<typename graph_traits<graph_t>::vertex_descriptor, typename graph_traits<graph_t>::vertex_descriptor  > shared_neighbors(Edge e, graph_t &G) {
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
					// as you as you find the second.
					return result;
				}
			}
		}
	}
	return result;
}

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

template<typename graph_t> 
bool polygonal(graph_t &g){
	typedef typename graph_traits<graph_t>::vertex_iterator vertex_iter_t;
	vertex_iter_t vi, vi_end;
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi){
		if (!locally_connected(*vi, g))
			return false;	
	}
	return true;
}

template<typename boundary_map_t, typename graph_t>
void boundary_graph(boundary_map_t &b_map, graph_t &g){
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	vertex_t null_vertex = graph_traits<graph_t>::null_vertex();
	typedef typename graph_traits<graph_t>::edge_iterator edge_iter_t;
	edge_iter_t ei, ei_end; 
	pair<vertex_t, vertex_t> shared;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
		shared = shared_neighbors(*ei, g);
		b_map[*ei] = (shared.first == null_vertex || shared.second == null_vertex);
	}			
}

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

template<typename edge_t, typename edge_cycle_map_t, typename graph_t>
void flip(vector<edge_t> &alt_kiss_cycle, edge_cycle_map_t &ec_map, graph_t &g){
	for(typename vector<edge_t>::iterator ei = alt_kiss_cycle.begin(); ei != alt_kiss_cycle.end(); ++ei)
		ec_map[*ei] = !ec_map[*ei];
}

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

#endif