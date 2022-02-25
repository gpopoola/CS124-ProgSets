#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>

using namespace std;

struct vertex {
    float coord[4];
    vertex* parent;
    int rank;
};

struct edge {
    float weight;
    vertex* src;
    vertex* dst;
};

default_random_engine generator;
uniform_real_distribution<double> weight_gen(0.0,1.0);

// Create vertices for the graph
vector<vertex*> vertex_constructor(int n, int dim){
    vector<vertex*> vertices;

    for (int i = 0; i < n; i++){
        vertex* v = new vertex;
        if (dim >= 2){
            v->coord[0] = weight_gen(generator);
            v->coord[1] = weight_gen(generator);
        }
        if (dim >= 3){
            v->coord[2] = weight_gen(generator);
        }
        if (dim == 4){
            v->coord[3] = weight_gen(generator);
        }
        v->parent = v;
        v->rank = 0;
        vertices.push_back(v);
    }

    return vertices;
}

bool compareEdges(edge* a, edge* b) {
    return b->weight > a->weight;
}

float euclidean_dist(edge* curr_e, int dim) {
    float diff = 0;

    for (int i = 0; i < dim; i++){
        diff += pow(curr_e->dst->coord[i] - curr_e->src->coord[i], 2.0);
    }

    return sqrt(diff);
}

vertex* find(vertex* v){
    if (v->parent == v){
        return v;
    }

    v->parent = find(v->parent);
    return v->parent;
}

void link(vertex* v, vertex* w){
    if (v->rank > w->rank){
        swap(v,w);

    } else if (v->rank == w->rank) {
    
        w->rank++;
    }
    v->parent = w;
}

void unite(vertex* v, vertex* w){
    link(find(v), find(w));
}

// Generate n(n-1)/2 edge weights for n vertices according to the dimension passed in
vector<edge*> G_builder(float n, int dim, vector<vertex*> vertices){
    
    vector<edge*> edge_arr;

    /*Thresholds prunes for edge weights at or below values set by the
    below functions with a given input n */
    float threshold_0 = (100)/(n-1);
    float threshold_2 = 2.25/pow(n, 1.0/2.0);
    float threshold_3 = 2.25/pow(n, 1.0/3.0);
    float threshold_4 = 2.25/pow(n, 1.0/4.0);

    /* Fill each spot in the array with an edge weight */
    /* Use nested for loops to facilitate setting the endpoints of edges */

        for (int i = 0; i < n-1; i++){
            for (int j = i; j < n-1; j++){
                edge* curr_e = new edge;
                curr_e->src = vertices[i];
                curr_e->dst = vertices[j+1];
                if (dim > 0){
                    curr_e->weight = euclidean_dist(curr_e, dim);
                }
                else {
                    curr_e->weight = weight_gen(generator);
                }

                if ((dim == 0 && curr_e->weight <= threshold_0) ||
                    (dim == 2 && curr_e->weight <= threshold_2) ||
                    (dim == 3 && curr_e->weight <= threshold_3) ||
                    (dim == 4 && curr_e->weight <= threshold_4)) {
                    edge_arr.push_back(curr_e);
                } else {
                    delete curr_e;
                }

            }
        }

    return edge_arr;
}

void G_destructor(vector<edge*> edge_arr, vector<vertex*> vertices){

    // Loop through every edge pointer in G and free its memory
    for (int i = 0; i < edge_arr.size(); i++){
        delete edge_arr[i];  
    }

    // Loop through every vertex pointer in G and free its memory
    for (int i = 0; i < vertices.size(); i++){
        delete vertices[i];
    }
}

float kruskal(vector<edge*> edge_arr, int n){
    int num_e = edge_arr.size();

    vector<edge*> X;

    float X_sum = 0.0;

    sort(edge_arr.begin(), edge_arr.end(), compareEdges);

    for (int i = 0; i < num_e; i++){
        if (find(edge_arr[i]->src) != find(edge_arr[i]->dst)) {
            X.push_back(edge_arr[i]);
            X_sum += edge_arr[i]->weight;

            unite(edge_arr[i]->src, edge_arr[i]->dst);
        }
    }
    return X_sum;
        
}

int main(int argc, char** argv){
    int n = atoi(argv[2]);
    int trials = atoi(argv[3]);
    int dim = atoi(argv[4]);

    if (n == 0 || n == 1){
        cout << "AVG: ";
        cout << 0 << endl;  
        return 0;
    }
    
    vector<edge*> edge_arr;
    vector<vertex*> vertices;
    float result = 0;

    for (int i = 0; i < trials; i++){
        vertices = vertex_constructor(n, dim);
        edge_arr = G_builder(n, dim, vertices);

        result += kruskal(edge_arr, n); 
        G_destructor(edge_arr, vertices);
    }

    cout << n << endl << trials << endl << dim << endl;
    cout << "SUM: ";
    cout << result << endl;

    result /= trials;

    cout << "AVG: ";
    cout << result << endl;
}