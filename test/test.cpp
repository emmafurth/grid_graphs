#include "boost/config.hpp" // put this first to suppress some VC++ warnings

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <time.h>
#include <sys/time.h>

#include <boost/utility.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/random_spanning_tree.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/property_map/shared_array_property_map.hpp>
#include <boost/graph/property_maps/constant_property_map.hpp>
#include <boost/lexical_cast.hpp>
/*#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>*/
#include "../triangular_grid_graph.hpp"

using namespace std;
using namespace boost;

//typedef adjacency_list<vecS, vecS, undirectedS> Graph;
typedef adjacency_list<vecS, vecS, undirectedS, 
												property<vertex_name_t, int, property<vertex_in_two_factor_t, bool> >, 
												property<edge_in_two_factor_t, bool> 
												> Graph;
typedef adjacency_list<vecS, vecS, undirectedS, 
												property<vertex_name_t, string, property<vertex_in_two_factor_t, bool> >, 
												property<edge_in_two_factor_t, bool> 
												> graph_t;												
typedef graph_traits<Graph>::vertex_descriptor vertex_t;
typedef graph_traits<Graph>::edge_descriptor edge_t;
typedef graph_traits<Graph>::vertex_iterator vertex_iter_t;
typedef graph_traits<Graph>::edge_iterator edge_iter_t;
typedef graph_traits<Graph>::adjacency_iterator adj_iter_t;
typedef graph_traits<Graph>::out_edge_iterator out_edge_iter_t;

vertex_t null_vertex = graph_traits<Graph>::null_vertex();
boost::mt19937 random_generator;

template<typename graph_t, typename edge_cycle_map_t>
void print(string filename, graph_t &G, edge_cycle_map_t &ec_map){
			//std::map<std::string,std::string> graph_attr, vertex_attr, edge_attr;
			//vertex_attr["shape"] = "point";
			std::map<std::string,std::string> graph_attr, vertex_attr, edge_attr;
			//vertex_attr["shape"] = "point";
			std::ofstream outfile(filename.c_str());
			write_graphviz(outfile, G, make_label_writer(get(vertex_name, G)), 
																 make_bool_writer(ec_map),
																 make_graph_attributes_writer(graph_attr, vertex_attr,
                                                     edge_attr));
		}

template<typename graph_t, typename edge_cycle_map_t, typename vertex_name_map_t>
void print(string filename, graph_t &G, edge_cycle_map_t &ec_map, vertex_name_map_t& name){
			//std::map<std::string,std::string> graph_attr, vertex_attr, edge_attr;
			//vertex_attr["shape"] = "point";
			std::map<std::string,std::string> graph_attr, vertex_attr, edge_attr;
			//vertex_attr["shape"] = "point";
			std::ofstream outfile(filename.c_str());
			write_graphviz(outfile, G, make_label_writer(name), 
																 make_bool_writer(ec_map),
																 make_graph_attributes_writer(graph_attr, vertex_attr,
                                                     edge_attr));
		}
		
template <typename Graph_T, typename PredMap>
void write_spanning_tree(const Graph_T& g, PredMap pred, string filename) {
  shared_array_property_map<string, typename property_map<Graph_T, edge_index_t>::const_type> edge_style(num_edges(g), get(edge_index, g));
  cout << num_edges(g) << endl;
  //shared_array_property_map<string, typename property_map<Graph_T, vertex_index_t>::const_type> vertex_pos(num_vertices(g), get(vertex_index, g));
  
  BGL_FORALL_EDGES_T(e, g, Graph_T) {
  	cout << get(edge_style, e) << endl;
  }
  
  cout << "--------\n";
  
  BGL_FORALL_EDGES_T(e, g, Graph_T) {
  	cout << "Edge " << e << ": " << (get(pred, target(e, g)) == source(e, g) || get(pred, source(e, g)) == target(e, g)) << endl;
    put(edge_style, e, (get(pred, target(e, g)) == source(e, g) || get(pred, source(e, g)) == target(e, g)) ? "bold" : "dotted");
  	//edge_style[e] = (get(pred, target(e, g)) == source(e, g) || get(pred, source(e, g)) == target(e, g)) ? "bold" : "dotted";
  	
  	cout << get(edge_style, e) << endl;
  }
  
  cout << "--------\n";
  
  BGL_FORALL_EDGES_T(e, g, Graph_T) {
  	cout << e << ": " << edge_style[e] << endl;
  }
  /*BGL_FORALL_VERTICES_T(v, g, Graph_T) {
    put(vertex_pos, v, lexical_cast<string>(v[0]) + "," + lexical_cast<string>(v[1]));
  }*/
  dynamic_properties dp;
  dp.property("style", edge_style);
  //dp.property("pos", vertex_pos);
  //dp.property("len", weight);
  dp.property("node_id", get(vertex_index, g));
  ofstream out(filename.c_str());
  write_graphviz_dp(out, g, dp);
}

// Taken from here: http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c
std::string exec(char* cmd) {
		strcat(cmd, " 2>&1");
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 128, pipe) != NULL)
    		result += buffer;
    }
    pclose(pipe);
    return result;
}
//template<typename graph_t, typename vertex_cycle_map_t, typename edge_cycle_map_t>
//void read_two_factor(string filename, graph_t &g, vertex_cycle_map_t &vc_map, edge_cycle_map_t &ec_map){
template<typename graph_t>
void read_two_factor(string filename, graph_t &g){
	//UnlabeledGraph G_prime;
	dynamic_properties dp;
	ifstream in(filename.c_str());
	dp.property("id", get(vertex_name, g)); 
	dp.property("v_in_cycle", get(vertex_in_two_factor, g));
	dp.property("e_in_cycle", get(edge_in_two_factor, g));
	//cout << "Test1\n";
	read_graphviz(in, g, dp, "id");
	
	//copy_to_labeled_graph(G_prime, G);
	//cout << "Test2\n";
}

vector<pair<int, int> > read_coord_file(string filename){
	vector<pair<int, int> > coords;
	
	if (filename.length() <= 6)
		return coords;
	if (0 != filename.compare (filename.length() - 6, 6, ".coord"))
		return coords;
	ifstream infile(filename.c_str());
	int x, y;
	char dummy;
	
	while (infile >> x >> dummy >> y)
	{
			cout << x << "," << y <<endl;
			coords.push_back(make_pair(x,y));
	}
	return coords;
}

/*void test_b_or_h(){
	enum { A, B, C, D, E, F, G , N};
	const char* name = "ABCDEFG";
	Graph g(N);
	add_edge(A, B, g);
	add_edge(A, C, g);
	add_edge(A, D, g);
	add_edge(B, D, g);
	add_edge(B, E, g);
	add_edge(C, D, g);
	add_edge(C, F, g);
	add_edge(D, E, g);
	add_edge(D, F, g);
	add_edge(D, G, g);
	add_edge(E, G, g);
	add_edge(F, G, g);
	
	map<vertex_t, bool> vc_map;
	vc_map[A] = false;
	vc_map[B] = false;
	vc_map[C] = false;
	vc_map[D] = true;
	vc_map[E] = false;
	vc_map[F] = true;
	vc_map[G] = true;
	
	
	list<vertex_t> extendors;
	extendors = cycle_extendors(vc_map, g);
	cout << "Extendors:\n";
	for(list<vertex_t>::iterator vi = extendors.begin(); vi != extendors.end(); ++vi)
		cout << *vi << endl; 
	
	map<edge_t, bool> ec_map;
	ec_map[edge(A, B, g).first] = false;
	ec_map[edge(A, C, g).first] = false;
	ec_map[edge(A, D, g).first] = false;
	ec_map[edge(B, D, g).first] = false;
	ec_map[edge(B, E, g).first] = false;
	ec_map[edge(C, D, g).first] = false;
	ec_map[edge(C, F, g).first] = false;
	ec_map[edge(D, E, g).first] = false;
	ec_map[edge(D, F, g).first] = true;
	ec_map[edge(D, G, g).first] = true;
	ec_map[edge(E, G, g).first] = false;
	ec_map[edge(F, G, g).first] = true;
	
	extend_cycle(vc_map, ec_map, g);
	
	vector<edge_t> b_or_h;
	edge_iter_t ei, ei_end;
	cout << "The following are in the cycle:\n";
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
		if (ec_map[*ei]) 
				cout << *ei << endl;
	}
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
		b_or_h = bolt_or_hourglass(*ei, vc_map, ec_map, g);
		if (!b_or_h.empty()){
			cout << "AKC for " << *ei << endl;
			for(vector<edge_t>::iterator ej = b_or_h.begin(); ej != b_or_h.end(); ++ej)
				cout << "   " << *ej << " is " << (ec_map[*ej] ? " BLACK." : " white.") << endl;
				
			cout << "Is " << (is_alt_kiss_cycle(b_or_h, ec_map, g) ? "" : "NOT ") << "an AKC\n";
		}
	}
	
	b_or_h.push_back(edge(D, F, g).first);
	
	cout << "Is " << (is_alt_kiss_cycle(b_or_h, ec_map, g) ? "" : "NOT ") << "an AKC\n";
	
	
}*/

void test_b_or_h_case_1(string filename){
	Graph g;
	edge_iter_t ei, ei_end;
	/*map<vertex_t, bool> vc_map;
	map<edge_t, bool> ec_map;*/
	
	
	//property_map<Graph, vertex_in_two_factor_t>::type vc_map;// = get(vertex_in_two_factor, g);
	//property_map<Graph, edge_in_two_factor_t>::type ec_map;// = get(edge_in_two_factor, g);
	
	
	/*map<vertex_t, bool> v_in_c;
	associative_property_map<map<vertex_t, bool> >
    vc_map(v_in_c);

	map<edge_t, bool> e_in_c;
	associative_property_map<map<edge_t, bool> >
    ec_map(e_in_c);*/
	read_two_factor(filename, g);//, vc_map, ec_map);
	
	property_map<Graph, vertex_in_two_factor_t>::type vc_map = get(vertex_in_two_factor, g);
	property_map<Graph, edge_in_two_factor_t>::type ec_map = get(edge_in_two_factor, g);
	
	vertex_iter_t vi, vi_end;
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
		cout << get(vertex_name, g, *vi) << " is in cycle: " << (vc_map[*vi] ? "YES" : "NO") << endl;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
		cout << "(" << get(vertex_name, g, source(*ei, g)) << "," << get(vertex_name, g, target(*ei, g)) << ") is in cycle: " << (ec_map[*ei] ? "YES" : "NO") << endl;
	
	
	list<vertex_t> extendors;
	extendors = cycle_extendors(vc_map, g);
	cout << "Extendors:\n";
	for(list<vertex_t>::iterator vi = extendors.begin(); vi != extendors.end(); ++vi)
		cout << *vi << endl;

	property_map<Graph, vertex_name_t>::type name = get(vertex_name, g);
	
	print("twofactor.dot", g, ec_map, name);

	if (star_of_david(g)){
		cout << "Is Star of David graph, not Hamiltonian.\n";
		return;
	}
	cout << "Test\n";
	extend_cycle(vc_map, ec_map, g);
	cout << "Test\n";
	vector<edge_t> b_or_h;
	cout << "The following are in the cycle:\n";
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
		if (ec_map[*ei]) 
				cout << *ei << endl;
	}
	
	/*vector<edge_t> e_neighbors = edge_neighbors(extendors.front(), g);
	
	find_linear_alt_kiss_cycle(e_neighbors[0], vc_map, ec_map, g);*/
	
	print("extended.dot", g, ec_map);
	
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi){
		if (get(vertex_name, g, *vi) == 8)
			break;
	}
	cout << "Start at: " << name[*vi] << endl;
	vector<vertex_t> vc_vector;
	vector<edge_t> ec_vector;
	cycle_as_vector(vc_vector, ec_vector, vc_map, ec_map, g, *vi);	
	cout << "Vertices in the cycle in order:\n";
	for(vector<vertex_t>::iterator vi = vc_vector.begin(); vi != vc_vector.end(); ++vi)
		cout << get(vertex_name, g, *vi) << endl;
	
	/*cout << "Edges in the cycle in order:\n";
	for(vector<edge_t>::iterator ei = ec_vector.begin(); ei != ec_vector.end(); ++ei)
		cout << "(" << get(vertex_name, g, source(*ei,g)) << "," << get(vertex_name, g, target(*ei,g)) << ")\n";*/
	
}


void test_extend(string filename){
	Graph g;
	edge_iter_t ei, ei_end;
	dynamic_properties dp;
	ifstream in(filename.c_str());
	
	dp.property("id", get(vertex_name, g)); 
	//dp.property("v_in_cycle", tmp1);
	dp.property("v_in_cycle", get(vertex_in_two_factor, g));
	//dp.property("e_in_cycle", tmp2);
	dp.property("e_in_cycle", get(edge_in_two_factor, g));

	cout << "Test1\n";
	read_graphviz(in, g, dp, "id");
	
	vertex_iter_t vi, vi_end;
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
		put(vertex_in_two_factor, g, *vi, false);
	
	int m = num_edges(g);
	cout << "Test\n";
	int e_id = rand()%m;
	int e_count = 0;
	edge_t e;
	
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
		if (e_count == e_id){
			e = *ei;
		}
		put(edge_in_two_factor, g, *ei, false);
		e_count++;
	}
	
	pair<vertex_t, vertex_t> shared = shared_neighbors(e, g);
	put(vertex_in_two_factor, g, source(e, g), true);
	put(vertex_in_two_factor, g, target(e, g), true);
	put(vertex_in_two_factor, g, shared.first, true);
	
	put(edge_in_two_factor, g, e, true);
	put(edge_in_two_factor, g, edge(shared.first, source(e, g), g).first, true);
	put(edge_in_two_factor, g, edge(shared.first, target(e, g), g).first, true);
	
	property_map<Graph, vertex_in_two_factor_t>::type vc_map = get(vertex_in_two_factor, g);
	property_map<Graph, edge_in_two_factor_t>::type ec_map = get(edge_in_two_factor, g);
	
	print("twofactor.dot", g, ec_map);
	// Should comment out these next two lines if dot isn't installed
	char cmd1[] = "dot -Kneato -Tpng -o twofactor.png twofactor.dot";
	exec(cmd1);
	
	cout << "Starting...\n";
	extend_cycle(vc_map, ec_map, g);
	
	print("extended.dot", g, ec_map);
	// Should comment out these next two lines if dot isn't installed
	char cmd2[] = "dot -Kneato -Tpng -o extended.png extended.dot";
	exec(cmd2);
}

void random_polygonal_graph(){
	graph_t g = random_triangular_grid_graph<graph_t>(5, 5);
	std::ofstream outfile("random.dot");
	write_graphviz(outfile, g, make_label_writer(get(vertex_name, g)));
	char cmd[] = "dot -Kneato -Tpng -o random.png random.dot";
	exec(cmd);
}

void test_boundary_two_factor(string filename){
	edge_iter_t ei, ei_end, ej, ej_end;
	/*Graph g;
	dynamic_properties dp;
	ifstream in(filename.c_str());
	
	dp.property("id", get(vertex_name, g)); 
	//dp.property("v_in_cycle", tmp1);
	dp.property("v_in_cycle", get(vertex_in_two_factor, g));
	//dp.property("e_in_cycle", tmp2);
	dp.property("e_in_cycle", get(edge_in_two_factor, g));

	read_graphviz(in, g, dp, "id");*/
	
	vector<pair<int, int> > coords = read_coord_file(filename);
	graph_t g = get_triangular_grid_graph_from_coords<graph_t>(coords);
	
	vertex_iter_t vi, vi_end;
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
		put(vertex_in_two_factor, g, *vi, false);
	
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
		put(edge_in_two_factor, g, *ei, false);
	
	property_map<graph_t, vertex_in_two_factor_t>::type vc_map = get(vertex_in_two_factor, g);
	property_map<graph_t, edge_in_two_factor_t>::type ec_map = get(edge_in_two_factor, g);
	boundary_based_two_factor(vc_map, ec_map, g);
	
	print("twofactor.dot", g, ec_map);
	char cmd[] = "dot -Kneato -Tpng -o twofactor.png twofactor.dot";
	exec(cmd);
	
	vector<int> component;
	int k = two_factor_components(ec_map, g, component);
	
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
		cout << "component[" << get(vertex_name, g, *vi) << "] = " << component[*vi] << endl;
		
	/*cout << "Parallelograms:\n";
	optional<Triple<edge_t> > result;
	for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
		for (tie(ej, ej_end) = edges(g); ej != ei; ++ej){
			if (ec_map[*ei] && ec_map[*ej] && is_parallelogram(*ei, *ej, g)){
				result = parallelogram(*ei, *ej, g);
				assert(result.is_initialized());
				cout << edge_str(*ei, g) << " and " << edge_str(*ej, g) << " along with " << edge_str(result.get().first, g) << ", " << edge_str(result.get().second, g) << ", and " << edge_str(result.get().third, g) << endl;
			}
		}
	}
	*/
	list<pair<edge_t, edge_t> > pars;
	pars = get_boundary_parallelograms(vc_map, ec_map, g);
	cout << "Boundary parallelograms:\n";
	for(list<pair<edge_t, edge_t> >::iterator pi = pars.begin(); pi != pars.end(); ++pi)
		cout << edge_str((*pi).first, g) << " and " << edge_str((*pi).second, g) << endl;
	
	typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, int, 
																									property<vertex_name_t, int> >,
																									property<edge_index_t, size_t> 
																									> c_graph_t;
	c_graph_t cg = get_component_graph<c_graph_t>(vc_map, ec_map, g);
	
	/*enum { A, B, C, N};
	const char* name = "ABC";
	c_graph_t cg(N);
	
	add_edge(A, B, 0, cg);
	add_edge(A, C, 1, cg);
	add_edge(B, C, 2, cg);*/
	
	
	std::ofstream outfile("component.dot");
	write_graphviz(outfile, cg, make_label_writer(get(vertex_name, cg))); 
	char cmd2[] = "dot -Kneato -Tpng -o component.png component.dot";
	exec(cmd2);
	
	//shared_array_property_map<gt::vertex_descriptor, property_map<graph_type, vertex_index_t>::const_type> pred(num_vertices(g), get(vertex_index, g));
	/*shared_array_property_map<graph_traits<c_graph_t>::vertex_descriptor, property_map<c_graph_t, vertex_index_t>::const_type> pred(num_vertices(cg), get(vertex_index, cg));
	
	random_spanning_tree(cg, random_generator, predecessor_map(pred));
	write_spanning_tree(cg, pred, "unweight_random_st.dot");

	graph_traits<c_graph_t>::vertex_iterator vj, vj_end;
	for(tie(vj, vj_end) = vertices(cg); vj != vj_end; ++vj)
		cout << "Predecessor of " << *vj << " is  " << pred[*vj] << endl;
		*/
		
	cout << merge_all_components(vc_map, ec_map, g) << endl;
	
	print("extended.dot", g, ec_map);
	char cmd3[] = "dot -Kneato -Tpng -o extended.png extended.dot";
	exec(cmd3);
	
}

int main(int argc,char* argv[]){
	
	srand(time(NULL));
	string filename = argv[1];
	
	//random_polygonal_graph();
	test_boundary_two_factor(filename);
	
	/*pair<int, int> coord1(2,0);
	pair<int, int> coord2(1,1);
	pair<float, float> float_coord1 = convert_coord(coord1);
	pair<float, float> float_coord2 = convert_coord(coord2);
	
	cout << "Convert: (" << float_coord1.first << "," << float_coord1.second << ") and ("  << float_coord2.first << "," << float_coord2.second << ")\n";
	cout << "Distance: " << coord_distance(coord1, coord2) << endl;*/
	
	//read_coord_file(filename);
	/*int i = 0;
	enum { A, B, C, D, E, F, G , N};
	const char* name = "ABCDEFG";
	Graph g(N);
	if (i == 0) {
		add_edge(A, B, g);
		add_edge(A, C, g);
		add_edge(A, D, g);
		add_edge(B, D, g);
		add_edge(B, E, g);
		add_edge(C, D, g);
		add_edge(C, F, g);
		add_edge(D, E, g);
		add_edge(D, F, g);
		add_edge(D, G, g);
		add_edge(E, G, g);
		add_edge(F, G, g);
	}
	else {
	
		add_edge(A, B, g);
		add_edge(A, C, g);
		add_edge(A, D, g);
		add_edge(B, D, g);
		add_edge(B, E, g);
		add_edge(C, D, g);
		add_edge(D, E, g);
		add_edge(D, F, g);
		add_edge(D, G, g);
		add_edge(F, G, g);
	}
	
	out_edge_iter_t ei, ei_end;
	for (tie(ei, ei_end) = out_edges(D, g); ei != ei_end; ++ei){
		pair<vertex_t, vertex_t> result = shared_neighbors(*ei, g);
		cout << "Shared neighbors of " << *ei << " are " << (result.first == null_vertex ? -1 : result.first)  << ", " << (result.second == null_vertex ? -1 : result.second) << endl;
	}
	
	vertex_iter_t vi, vi_end;
	for(tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
		cout << *vi << " is " << (locally_connected(*vi, g) ? "" : " NOT ") << "locally connected.\n";
	
	cout << "The graph is " << (polygonal(g) ? "" : "NOT ") << "polygonal\n";
	
	map<edge_t, bool> b_map;

	boundary_graph(b_map, g);
	edge_iter_t ej, ej_end;
	cout << "The following are in the boundary map:\n";
	for (tie(ej, ej_end) = edges(g); ej != ej_end; ++ej){
		if (b_map[*ej]) 
			cout << *ej << endl;
	}
	
	map<vertex_t, bool> vc_map;
	vc_map[A] = false;
	vc_map[B] = false;
	vc_map[C] = false;
	vc_map[D] = true;
	vc_map[E] = false;
	vc_map[F] = true;
	vc_map[G] = true;
	
	
	list<vertex_t> extendors;
	extendors = cycle_extendors(vc_map, g);
	cout << "Extendors:\n";
	for(list<vertex_t>::iterator vi = extendors.begin(); vi != extendors.end(); ++vi)
		cout << *vi << endl; 
	
	map<edge_t, bool> ec_map;
	if (i == 0){
		ec_map[edge(A, B, g).first] = false;
		ec_map[edge(A, C, g).first] = false;
		ec_map[edge(A, D, g).first] = false;
		ec_map[edge(B, D, g).first] = false;
		ec_map[edge(B, E, g).first] = false;
		ec_map[edge(C, D, g).first] = false;
		ec_map[edge(C, F, g).first] = false;
		ec_map[edge(D, E, g).first] = false;
		ec_map[edge(D, F, g).first] = true;
		ec_map[edge(D, G, g).first] = true;
		ec_map[edge(E, G, g).first] = false;
		ec_map[edge(F, G, g).first] = true;
		
		extend_cycle(vc_map, ec_map, g);
		
		cout << "The following are in the cycle:\n";
		for (tie(ej, ej_end) = edges(g); ej != ej_end; ++ej){
			if (ec_map[*ej]) 
				cout << *ej << endl;
		}
		
	}

	
	*/
	
	
}