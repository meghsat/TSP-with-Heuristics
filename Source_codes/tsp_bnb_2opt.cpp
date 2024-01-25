#include <bits/stdc++.h>
#include <ctime>
#include <chrono>
using namespace std;
#define INF INT_MAX
map<vector<int>, double> ans; 
double mindist = 1e9; 
double upperbound=0;
double lowerbound=0;

int get_minkey(double key[], bool mstarr[], int V) {
    double min = INT_MAX;
    int min_index;
    for (int v = 0; v < V; v++) {
        if (mstarr[v] == false && key[v] < min)
            min = key[v], min_index = v;
    }
    return min_index;
}
double Weight_MST(int parent[], vector<vector<double>>& graph,int V)
{
    double weight=0;
    for (int i = 1; i < V; i++)
    {
        weight+=graph[i][parent[i]];
    }
        return weight;
}
int* MST(vector<vector<double>>& graph, int V) {
    int* parent = new int[V];
    double* key = new double[V];
    bool* mstarr = new bool[V];
    for (int i = 0; i < V; i++) {
        key[i] = INT_MAX, mstarr[i] = false;
    }
    key[0] = 0;
    parent[0] = -1;
    for (int count = 0; count < V - 1; count++) {
        int u = get_minkey(key, mstarr, V);
        mstarr[u] = true;
        for (int v = 0; v < V; v++) {
            if (graph[u][v] && mstarr[v] == false && graph[u][v] < key[v])
                parent[v] = u, key[v] = graph[u][v];
        }
    }
   
    return parent;
}
void lowerbound_MST(vector<vector<double>>& graph,int n)
{
    for(int i=0;i<n;i++)
    {
        vector<vector<double>> graph2(n, vector<double>(n, 0));
        int a=0,b=0;
        for (int j = 0; j < n; j++) {
        if(j!=i){
            b=0;
        for (int k = 0; k < n; k++) {
            if(k!=i){
            graph2[a][b]=graph[j][k];
            b++;
            }
        }
        a++;
        }
    }
    int first_min,second_min=0;
    vector<double> graph_sub={};
    for(int p=0;p<n;p++)
    {
        if(graph[i][p]!=0)
        {
            graph_sub.push_back(graph[i][p]);
        }
    }
    sort(graph_sub.begin(),graph_sub.end());
    first_min=graph_sub[0],second_min=graph_sub[1];
    int* parent=MST(graph2,n-1);
lowerbound= max(Weight_MST(parent, graph2,n-1)+first_min+second_min,lowerbound);
}}

double calculateTourDistance(const vector<int>& tour, const vector<vector<double>>& distanceMatrix) {
    double distance = 0;
    for (size_t i = 0; i < tour.size() - 1; ++i) {
        distance += distanceMatrix[tour[i]][tour[i + 1]];
    }
    distance += distanceMatrix[tour.back()][tour.front()]; 
    return distance;
}

vector<int> twoOptSwap(const vector<int>& tour, int i, int k) {
    vector<int> newTour = tour;
    while (i < k) {
        swap(newTour[i], newTour[k]);
        i++;
        k--;
    }
    return newTour;
}
vector<int> twoOptHeuristic(const vector<int>& tour, const vector<vector<double>>& distanceMatrix) {
    int n = tour.size();
    vector<int> currentTour = tour;
    for (int i = 0; i < n - 1; ++i) {
        for (int k = i + 1; k < n ; ++k) {
            vector<int> newTour = twoOptSwap(currentTour, i , k);
            // cout<<i+1<<" "<<k<<endl;
            // for(int m=0;m<newTour.size();m++)
            // {
            //     cout<<newTour[m]<<" ";
            // }      
            double currentDistance = calculateTourDistance(currentTour, distanceMatrix);
            double newDistance = calculateTourDistance(newTour, distanceMatrix);
          //  cout<<currentDistance<<" "<<newDistance<<endl;
            if (newDistance < currentDistance) {
                currentTour = newTour;
            }
    //         auto end = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end-start;
    // if(elapsed_seconds.count()/60>=9)
    // {
    //     cout<<elapsed_seconds.count()/60<<" ";
    //     return currentTour;
    // }
        }
    }
    cout<<"done";
    return currentTour;
}

void bnbdfs(int u, bool visited[], vector<vector<double>>& graph, double curdist, int cntr, vector<int> paths, int n, int source) 
{
    visited[u] = 1;
    cntr++;
    if (cntr == n) {
      // cout<<"distance "<<curdist + graph[u][source]<<endl;
        if (curdist + graph[u][source] < mindist &&  curdist + graph[u][source] <= upperbound*2 && curdist + graph[u][source] >= lowerbound) {
            mindist = curdist + graph[u][source];
            vector<int> optimizedTour = twoOptHeuristic(paths, graph);
            double currentDistance = calculateTourDistance(optimizedTour, graph);
            cout<<" mintdist: "<<mindist<<" "<<currentDistance<<endl;
            if(mindist<=currentDistance)
            {
                ans[paths] = curdist + graph[u][source];
                for(int i=0;i<paths.size();i++)
            {
                cout<<paths[i]<<" ";
            }
           cout<<" mintdist: "<<mindist<<endl;
            }
            else{
                mindist=currentDistance;
                ans[optimizedTour]=mindist;
                 cout<<"optimized tour";
            for(int i=0;i<optimizedTour.size();i++)
            {
                cout<<optimizedTour[i]<<" ";
            }
             cout<<" mintdist: "<<mindist<<endl;
            }
    //         auto end = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end-start;
    // if(elapsed_seconds.count()/60>=9)
    // {
    //     cout<<elapsed_seconds.count()/60<<" ";
    //     return;
    // }
        }
    }
    for (int i = 0; i < n; i++) {
        if (visited[i] == 0) {
            paths.push_back(i);
          //  cout<<"distance "<<curdist + graph[u][i]<<endl;
            if (curdist + graph[u][i] < mindist && curdist + graph[u][i] <= upperbound*2 ) {
                bnbdfs(i, visited, graph, curdist + graph[u][i], cntr, paths, n, source);
    //              auto end = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end-start;
    // if(elapsed_seconds.count()/60>=9)
    // {
    //     cout<<elapsed_seconds.count()/60<<" ";
    //     return;
    // }
            }
            paths.pop_back();
        }
    }
    visited[u] = 0;
}

int main() {
    string filename;
    cout << "Enter the input file name with extension: ";
    cin >> filename;

    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return 1;
    }

    int n;
    inputFile >> n;

    vector<vector<double>> graph(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inputFile >> graph[i][j];
        }
    }
    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < n; j++) {
    //         cout << graph[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    // Heuristic: MST
    int* parent = MST(graph, n);
    upperbound= 2*Weight_MST(parent, graph,n);
    lowerbound_MST(graph, n);
    cout<<upperbound<<" "<<lowerbound<<endl;
    // for(int i=0;i<n;i++)
    // {
    //     cout<<parent[i]<<" ";
    // }
    bool* visited = new bool[n]{0};
   
    for (int i = 0; i < n; i++) {
        int s = i;
        visited[s] = 1;
        bnbdfs(s, visited, graph, 0, 0, {s}, n, s);
    //     auto end = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end-start;
    // if(elapsed_seconds.count()/60>=9)
    // {
    //     cout<<elapsed_seconds.count()/60<<" ";
    //     break;
    // }
    }

    cout << mindist << endl;
    map<vector<int>, double>::iterator it = ans.begin();
    while (it != ans.end()) {
        cout << "Key: ";
        for (const int& value : it->first) {
            cout << value << " ";
        }
        cout << "Value: " << it->second << endl;
        ++it;
    }
    return 0;
}
