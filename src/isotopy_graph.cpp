#include "isotopy_graph.h"
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace Isotopy {

Graph::Graph(int delta, const std::vector<bool>&sign_vector, const std::set<std::pair<int, int>>& edges)
: delta(delta) {
  std::size_t nverts = (delta + 1) * (delta + 2) / 2;
  assert(sign_vector.size() == nverts && "sign vector length does not match number of vertices");

  size_t nverts_triang = (delta + 1) * (delta + 2) / 2;
  size_t nvert_formula = 4 * nverts_triang - 4 * (delta + 1) +1;

  std::map<int,int> vertex_map;

  sign.resize(nvert_formula);
  side_points.resize(4 * delta);
  sides.resize(4 * delta, std::vector<bool>(4, false)); 

  int side_point_count = 0;
  int global_count = 0;
  int count = 0;
  //First Quadrant
  for (int limit = delta; limit >= 0; --limit) {
    for (int i = 0; i <= limit; ++i) {
      if (i == limit) {
        side_points[side_point_count] = count;
        if (limit == 0) {
          sides[side_point_count][0] = true;
          sides[side_point_count][1] = true;
        } else if (limit == delta) {
          sides[side_point_count][0] = true;
          sides[side_point_count][3] = true;
        } else {
          sides[side_point_count][0] = true;
        }
        side_point_count++;

      }

      vertex_map[count] = count;
      sign[count] =sign_vector[count];
      count++;
      global_count++;
    }
  }
  //Second Quadrant
  count = 0;
  for (int limit = delta; limit >= 0; --limit) {
    for (int i = 0; i <= limit; ++i) {
      if (i == 0) {
        vertex_map[count+nverts] = count;
      } else {
        if (i == limit) {
          side_points[side_point_count +  limit -1 ] = global_count;
          if (limit == delta) {
            sides[side_point_count + limit -1][1] = true;
            sides[side_point_count + limit -1][2] = true;
          } else {
            sides[side_point_count + limit -1][1] = true;
          }
        }
        vertex_map[count+nverts] = global_count;
        if (i % 2 == 0) {
          sign[global_count] = sign_vector[count];
        } else {
          sign[global_count] = !sign_vector[count];
        }
        global_count++;
      }
      count++;
    }
  }
  side_point_count += delta;
  
  //Third and Fourth Quadrants
  count = 0;
  for (int limit = delta; limit >= 0; --limit) {
    for (int i = 0; i <= limit; ++i) {
      if (limit == delta) {
        vertex_map[count+2*nverts] = vertex_map[count+nverts];
        vertex_map[count+3*nverts] = vertex_map[count];
      } else if (i == 0) {
        vertex_map[count+2*nverts] = global_count;
        vertex_map[count+3*nverts] = global_count;
        if (limit % 2 == 0) {
          sign[global_count] = sign_vector[count];
        } else {
          sign[global_count] = !sign_vector[count];
        }
        if (limit == 0) {
          side_points[side_point_count] = global_count;
          sides[side_point_count][2] = true;
          sides[side_point_count][3] = true;
          side_point_count++;
        }
        global_count++;
      } else {
        vertex_map[count+2*nverts] = global_count;
        if (i % 2 == 0) {
          sign[global_count] = sign_vector[count];
        } else {
          sign[global_count] = !sign_vector[count];
        }
        if (limit % 2 == 0) {
          sign[global_count] = sign[global_count];
        } else {
          sign[global_count] = !sign[global_count];
        }
        //Save new side points 
        if (i == limit) {
          side_points[side_point_count] = global_count;
          sides[side_point_count][2] = true;
          side_points[side_point_count + 2*limit] = global_count+1;
          sides[side_point_count + 2*limit][3] = true;
          side_point_count++;
        }



        global_count++;
        vertex_map[count+3*nverts] = global_count;
        if (limit % 2 == 0) {
          sign[global_count] = sign_vector[count];
        } else {
          sign[global_count] = !sign_vector[count];
        }
        global_count++;
      }
      count++;
    }
  }
  /*
  std::cout << "Vertex Map:\n";
  for (const auto& [key, value] : vertex_map) {
    std::cout << "Original: " << key << " -> New: " << value << "\n";
  }
  */

  adjacency.resize(global_count);
  for (const auto& edge : edges) {
    assert(edge.first < static_cast<int>(nverts) && edge.second < static_cast<int>(nverts) && "edge vertex index out of range");
    for (size_t k = 0; k < 4; ++k) {
      int v1 = vertex_map[edge.first + k * nverts];
      int v2 = vertex_map[edge.second + k * nverts];

      adjacency[v1].insert(v2);
      adjacency[v2].insert(v1);
    }
  }
}

void Graph::connected_components() {
  std::set<int> visited;
  int component_count = 0;

  component.resize(adjacency.size(), -1); // Ensure enough space for all vertex indices

  for (size_t vertex = 0; vertex < adjacency.size(); ++vertex) {
    if (visited.find(vertex) == visited.end()) {
      std::set<int> stack;
      stack.insert(vertex);
      std::set<int> temp_adjacency =  std::set<int>();

      while (!stack.empty()) {
        int current = *stack.begin();
        stack.erase(stack.begin());
        visited.insert(current);
        component[current] = component_count;

        for (int neighbor : adjacency[current]) {
          if (visited.find(neighbor) == visited.end()) {
            if (sign[neighbor] != sign[current]) continue; // Only traverse edges with same sign
              stack.insert(neighbor);
          } else {
            if (sign[neighbor] != sign[current]){
              auto c2 = component[neighbor];
              temp_adjacency.insert(c2);
              component_adjacency[c2].insert(component_count);

            }
          }
        }
      }
      component_adjacency.push_back(temp_adjacency);
      component_count++;
    }
  }
  std::cout << "Total components found: " << component_count << "\n";
  for (size_t c = 0; c < component_adjacency.size(); ++c) {
    std::cout << "Component " << c << " is adjacent to components: ";
    for (int adj : component_adjacency[c]) {
      std::cout << adj << " ";
    }
    std::cout << "\n";
  }

}

void Graph::pre_isotopy_root() {
  // Placeholder for pre-isotopy root logic
  /*
  We need to compute: 
  nbs: A vector of sets of other sidepoints a point is connected to via its component
  */

  std::vector<std::set<int>> nbs(side_points.size());

  for (size_t i = 0; i < side_points.size(); ++i) {
    int vertex = side_points[i];
    int comp = component[vertex];
  }



  /*
  auto component2side_point = std::vector<std::set<int>>(component_adjacency.size());
  for (size_t i = 0; i < side_points.size(); ++i) {
    int vertex = side_points[i];
    //int antipode = side_points[(i + 2 * delta) % (4 * delta)];

    int comp = component[vertex];
    component2side_point[comp].insert(vertex);
  }
  */

}


void Graph::isotopy_root() {

}




}
