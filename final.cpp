#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

struct edge{
    int truss;
    int support;
    int settle;
};

bool order(int a, int b, unordered_map<int, int>& degrees){
    bool var = false;
    if(degrees[a] < degrees[b]){
        var = true;
    } 
    else if(degrees[a] == degrees[b] && a < b){
        var = true;
    }
    return var;
}

void update_support(unordered_map<int, unordered_map<int, int>>& graph, unordered_map<int, int>& degrees, 
unordered_map<int, int>& vtop, int rank, int size){
    unordered_map<int, set<vector<int>>> pinfo; 
    for(auto x : graph){
        for(auto y : x.second){
            if(graph.find(y.first) != graph.end()){
                for(auto z : graph[x.first]){
                    if(order(y.first, z.first, degrees)){
                        if(graph[y.first].find(z.first) != graph[y.first].end()){
                            graph[x.first][y.first]++;          
                            graph[x.first][z.first]++;
                            graph[y.first][z.first]++;
                        }    
                    }
                    else if(graph.find(z.first) == graph.end()){
                        vector<int> v;
                        v.push_back(x.first);
                        v.push_back(z.first);
                        v.push_back(y.first); 
                        pinfo[vtop[z.first]].insert(v);                                               
                    }
                }
            }
            else{
                for(auto z : graph[x.first]){
                    if(order(y.first, z.first, degrees)){
                        vector<int> v;
                        v.push_back(x.first);
                        v.push_back(y.first);
                        v.push_back(z.first); 
                        pinfo[vtop[y.first]].insert(v);
                    }
                }                
            }
        }
    }

    int sendbuf[size];
    int recvbuf[size];
    memset(recvbuf, 0, size);
    memset(sendbuf, 0, size);

    int sendcount = 0;
    for(int i = 0; i < size; i++){
        sendbuf[i] = pinfo[i].size();
        sendcount += pinfo[i].size();
    }
    
    MPI_Alltoall(sendbuf, 1, MPI_INT, recvbuf, 1, MPI_INT, MPI_COMM_WORLD);

    int recvcount = 0;
    for(int i = 0; i < size; i++){
        recvcount += recvbuf[i];
    }

    int sendreq[sendcount*3], recvreq[recvcount*3];
    
    int send_counts[size], recv_counts[size], send_disp[size], recv_disp[size];

    int k = 0;
    for(int i = 0; i < size; i++){
        send_disp[i] = k;
        send_counts[i] = 3*pinfo[i].size();
        for(auto a : pinfo[i]){
            sendreq[k++] = a[0];
            sendreq[k++] = a[1];
            sendreq[k++] = a[2];
        }
    }

    k = 0;
    for(int i = 0; i < size; i++){
        recv_disp[i] = k;
        recv_counts[i] = 3*recvbuf[i];
        k += recv_counts[i];
    }

    MPI_Alltoallv(sendreq, send_counts, send_disp, MPI_INT, recvreq, recv_counts, recv_disp, MPI_INT, MPI_COMM_WORLD);

    unordered_map<int, set<vector<int>>> responses;

    for(int i = 0; i < size; i++){
        int l = recv_disp[i];
        for(int j = 0; j < recv_counts[i]; j += 3){
            int u = recvreq[l++];
            int v = recvreq[l++];
            int w = recvreq[l++];
           
            if(graph[v].find(w) != graph[v].end()){
                responses[i].insert({u, v, w, 1});
                graph[v][w]++;
            }
            else{
                responses[i].insert({u, v, w, 0});
            }
        }
    }

    int response_send[4*recvcount], response_recv[4*sendcount];

    int send_response_disp[size], send_response_counts[size], recv_response_disp[size], recv_response_counts[size];

    k = 0;

    for(int i = 0; i < size; i++){
        send_response_disp[i] = k;
        send_response_counts[i] = 4*responses[i].size();
        for(auto a : responses[i]){
            response_send[k++] = a[0];
            response_send[k++] = a[1];
            response_send[k++] = a[2];
            response_send[k++] = a[3];
        }
    }
    
    k = 0;

    for(int i = 0; i < size; i++){
        recv_response_disp[i] = k;
        recv_response_counts[i] = 4*pinfo[i].size();
        k += recv_response_counts[i];
    }

    MPI_Alltoallv(response_send, send_response_counts, send_response_disp, MPI_INT, response_recv, recv_response_counts, 
    recv_response_disp, MPI_INT, MPI_COMM_WORLD);

    for(int i = 0; i < size; i++){
        int l = recv_response_disp[i];
        for(int j = 0; j < recv_response_counts[i]; j += 4){
            int u = response_recv[l++];
            int v = response_recv[l++];
            int w = response_recv[l++];
            int ans = response_recv[l++];

            if(ans == 1){
                graph[u][v]++;
                graph[u][w]++;
            }
        }
    }
    // cout << rank << " " << "update_support done" << endl;
}

// struct CompareEdge {
//     bool operator()(pair<pair<int, int>, edge*>& a, pair<pair<int, int>, edge*>& b) {
//         return a.second->truss > b.second->truss;
//     }
// };

// int find_min(unordered_map<int, unordered_map<int, edge*>>& graph){
//     int minimum = INT_MAX;
//     priority_queue<pair<pair<int, int>, edge*>, vector<pair<pair<int, int>, edge*>>, CompareEdge> pq;

//     for(auto x : graph){
//         for(auto y : x.second){
//             if(y.second->settle == 0){
//                 pq.push({{x.first, y.first}, y.second});
//             }
//         }
//     }

//     while(!pq.empty()){
//         auto curr = pq.top();
//         pq.pop();
//         if(curr.second->settle == 0){
//             minimum = curr.second->truss;
//             break;
//         }
//     }

//     return minimum;
// }

// int initialise(unordered_map<int, unordered_map<int, edge*>>& graph){
//     int total = 0;
//     for(auto a : graph){
//         for(auto b : a.second){
//             b.second->truss = b.second->support + 2;
//         }
//         total += graph[a.first].size();
//     }
//     return total;
// }

// void min_truss(unordered_map<int, unordered_map<int, edge*>>& graph){
//     int count = initialise(graph);
//     int min;
//     while(true){
//         min = find_min(graph);
//         for(auto x : graph){
//             for(auto y : x.second){
//                 if(y.second->truss == min && y.second->settle == 0){
//                     for(auto z : graph[x.first]){
//                         if(graph[y.first].find(z.first) != graph[y.first].end()){
//                             if(graph[x.first][z.first]->settle == 0 && graph[y.first][z.first]->settle == 0){
//                                 if(graph[x.first][z.first]->truss != min){
//                                     graph[x.first][z.first]->truss--;
//                                 }
//                                 if(graph[y.first][z.first]->truss != min){
//                                     graph[y.first][z.first]->truss--;                                    
//                                 }    
//                             }
//                         }
//                         else if(graph[z.first].find(y.first) != graph[z.first].end()){
//                             if(graph[x.first][z.first]->settle == 0 && graph[z.first][y.first]->settle == 0){
//                                 if(graph[x.first][z.first]->truss != min){
//                                     graph[x.first][z.first]->truss--;
//                                 }
//                                 if(graph[z.first][y.first]->truss != min){
//                                     graph[z.first][y.first]->truss--;                                    
//                                 }    
//                             }
//                         }
//                     }
//                     graph[x.first][y.first]->settle = 1;
//                     count -= 1;
//                 }
//             }
//         }
//         if(count == 0){
//             break;
//         }
//     }
// }

void quick_truss(unordered_map<int, unordered_map<int, int>>& graph, int k){
    queue<pair<int, int>> deletable;
    for(auto a : graph){
        for(auto b : a.second){
            if(b.second < k && a.first < b.first){
                deletable.push(make_pair(a.first, b.first));
            }
        }
    }
    // cout << "step1 " << endl;
    while(!deletable.empty()){
        int x = deletable.front().first;
        int y = deletable.front().second;
        deletable.pop();
        
        if(graph[x].find(y) != graph[x].end()){
            graph[x].erase(y);
            graph[y].erase(x);
            for(auto z : graph[x]){
                // cout << "step2 " << z.first << endl;
                if(graph[y].find(z.first) != graph[y].end()){
                    
                    graph[y][z.first]--;
                    // cout << "step3 " << z.first << " " << y << endl;
                    graph[z.first][y]--;
                    
                    if(graph[y][z.first] < k){
                        deletable.push(make_pair(y, z.first));
                    }
                    
                    graph[x][z.first]--;
                    graph[z.first][x]--;
                    if(graph[x][z.first] < k){
                        deletable.push(make_pair(x, z.first));
                    }
                }
            }
        }
    }

    // cout << "quick truss done" << endl;
}

void quick_filter(unordered_map<int, unordered_map<int, int>>& graph){
    vector<int> todelete;
    for(auto a : graph){
        if(a.second.empty()){
            todelete.push_back(a.first);
        }
    }
    for(auto a : todelete){
        graph.erase(a);
    }
    // cout << "quick filter done" << endl;
}

void dfs(int u, unordered_map<int, unordered_map<int, int>>& graph, unordered_map<int, bool>& visited, vector<int>& component, int k) {
    visited[u] = true;
    component.push_back(u);
    for (auto v : graph[u]) {
        if(!visited[v.first]) {
            dfs(v.first, graph, visited, component, k);
        }
    }
}

vector<vector<int>> find_Connected_Components(unordered_map<int, unordered_map<int, int>>& graph, int k) {
    quick_truss(graph, k);
    quick_filter(graph);
    vector<vector<int>> components;
    unordered_map<int, bool> visited;
    for(auto a : graph){
        visited[a.first] = false;
    }
    for(auto u : graph) {
        if(!visited[u.first]) {
            vector<int> component;
            dfs(u.first, graph, visited, component, k);
            components.push_back(component);
        }
    }
    // cout << "found connected components " << components.size() << endl;
    return components;
}

unordered_map<int, vector<int>> influencer_vertices(unordered_map<int, unordered_map <int, int>>& graph, int k, int p) {
    unordered_map<int, unordered_map<int, int>> graph_temp = graph;
    vector<vector<int>> k_truss_groups = find_Connected_Components(graph, k);

    unordered_map<int, vector<int>> influencers;

    unordered_map<int, int> vtocomp;
    for(int i = 0; i < k_truss_groups.size(); i++){
        for(int j = 0; j < k_truss_groups[i].size(); j++){
            vtocomp[k_truss_groups[i][j]] = i;
        }
    }

    vector<int> vertices;
    for(auto i : graph_temp){
        vertices.push_back(i.first);
    }
    for(int i = 0; i < vertices.size(); i++){
        unordered_set<int> ntocomp;
        for(auto a : graph_temp[vertices[i]]){
            if(vtocomp.find(a.first) != vtocomp.end()){
                ntocomp.insert(vtocomp[a.first]);
            }
        }
        if(ntocomp.size() >= p){
            // cout << ntocomp.size() << endl;
            vector<int> influenced;
            for(auto a : ntocomp){
                for(auto b : k_truss_groups[a]){
                    influenced.push_back(b);
                }
            }
            influencers[vertices[i]] = influenced;
        }
    }
    // cout << "influencers found" << endl;
    return influencers;
}

int main(int argc, char* argv[]) {

    string t = argv[1];
    int task = stoi(t.substr(t.find("=")+1));

    string inputpath = argv[2];
    inputpath = inputpath.substr(inputpath.find("=") + 1);

    string headerpath = argv[3];
    headerpath = headerpath.substr(headerpath.find("=") + 1);

    string outputpath = argv[4];
    outputpath = outputpath.substr(outputpath.find("=")+1);

    string v = argv[5];
    int verbose = stoi(v.substr(v.find("=")+1));

    string start_f = argv[6];
    int start_k = stoi(start_f.substr(start_f.find("=")+1));

    string end_f = argv[7];
    int end_k = stoi(end_f.substr(end_f.find("=")+1));

    string p_num = argv[8];
    int p_ego = stoi(p_num.substr(p_num.find("=")+1));

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 0;
    int m = 0;

    ifstream inputfile(inputpath, ios::binary);
    inputfile.read((char*)&n, 4);
    inputfile.read((char*)&m, 4);

    ifstream headerfile(headerpath, ios::binary);

    int subgraphSize = n / size;
    int vertex_offset, vertex, degree;
    
    unordered_map<int, int> degrees;    
    unordered_map<int, int> vtop;   
    unordered_map<int, unordered_map<int, int>> subgraph; 

    for(int i = 0; i < n; i++){    
        inputfile.read((char*)&vertex, 4);
        inputfile.read((char*)&degree, 4);
        degrees[vertex] = degree;
        inputfile.ignore(4*degree);
    }

    int rank_offset;
    int start = rank*subgraphSize;
    int end;
    if(rank < size - 1){
        end = (rank + 1)*subgraphSize - 1;
    }
    else{    
        end = n - 1;
    }

    for(int i = 0; i < size; i++){
        int start_p = i*subgraphSize;
        int end_p;
        if(i < size - 1){
            end_p = (i + 1)*subgraphSize - 1;
        }
        else{
            end_p = n - 1;
        }
        for(int j = start_p; j <= end_p; j++){
            vtop[j] = i;
        }
    }
    headerfile.seekg(4*start, ios::beg);
    headerfile.read((char*)&rank_offset, 4);
    inputfile.seekg(rank_offset, ios::beg);
    for(int i = start; i <= end; i++){
        int node = 0;
        inputfile.read((char*)&node, 4);
        int deg = 0;
        inputfile.read((char*)&deg, 4);
        for(int j = 0; j < deg; j++){
            int neighbour = 0;
            inputfile.read((char*)&neighbour, 4);
            if(order(node, neighbour, degrees)){
                subgraph[node][neighbour] = 0;
            }
        }
    }
    inputfile.close();
    headerfile.close();

    ofstream outfile(outputpath);

    // cout << "read done" << endl;

    // cout << "help: " << subgraph.size() << endl;

    update_support(subgraph, degrees, vtop, rank, size);  

    int recvbuf[size];
    memset(recvbuf, 0, size);

    int sendcount = 0;
    for(auto a : subgraph){
        sendcount += a.second.size();
        // cout << rank << " debug " << a.first << "->" << a.second.size();
    }
    

    // cout << rank << " debug " << sendcount << endl;

    // string s = "";
    // for(int i = 0; i < size; i++){
    //     s += to_string(recvbuf[i]);
    //     s += " ";
    // }
    
    // cout << rank << " " << sendcount << ", " << s;

    MPI_Gather(&sendcount, 1, MPI_INT, recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // cout << rank << " " << "step 1" << endl;
    // cout << "recvbuf[0]" << recvbuf[0] << endl;
    int recvcount = 0;
    if(rank == 0){    
        for(int i = 0; i < size; i++){
            recvcount += recvbuf[i];
            // cout << recvbuf[i];
        }
        // cout << endl;
    } 

    // cout << rank << "debug hoja mc" << endl;
    // cout << rank << " " << sendcount << " " << recvcount << endl;

    int sendreq[sendcount*3];
    int* recvreq;
    int* recv_counts;
    int* recv_disp;

    if(rank == 0){
        recvreq = new int[recvcount*3];
        recv_counts = new int[size];
        recv_disp = new int[size];
    }

    // cout << rank << " " << "step 2" << endl;
    
    int k = 0;
    for(auto a : subgraph){
        for(auto b : a.second){
            sendreq[k++] = a.first;
            sendreq[k++] = b.first;
            sendreq[k++] = b.second;
        }
    }

    k = 0;
    if(rank == 0){
        for(int i = 0; i < size; i++){
            recv_disp[i] = k;
            recv_counts[i] = 3*recvbuf[i];
            k += recv_counts[i];
        }
    }

    MPI_Gatherv(sendreq, sendcount*3, MPI_INT, recvreq, recv_counts, recv_disp, MPI_INT, 0, MPI_COMM_WORLD);

    // cout << rank << " done till here" << endl;
    
    if(rank == 0){

        for(int i = 1; i < size; i++){
            int l = recv_disp[i];
            for(int j = 0; j < recv_counts[i]; j += 3){
                int u = recvreq[l++];
                int v = recvreq[l++];
                int support = recvreq[l++];
                subgraph[u][v] = support;
            }
        }

        for(auto a : subgraph){
            for(auto b : a.second){
                subgraph[b.first][a.first] = subgraph[a.first][b.first];
            }
        }

        // cout << "total graph ready to be processed" << endl;

        // cout << endl;

        if(task == 1){

            if(verbose == 0){
                for(int s = start_k; s <= end_k; s++){
                    vector<vector<int>> components;
                    components = find_Connected_Components(subgraph, s);
                    if(s < end_k){
                        if(components.empty()){
                            outfile << 0 << " ";
                        }
                        else{
                            outfile << 1 << " ";
                        }
                    }
                    else{
                        if(components.empty()){
                            outfile << 0;
                        }
                        else{
                            outfile << 1;
                        }    
                    }
                }
            }

            else if(verbose == 1){
                for(int s = start_k; s <= end_k; s++){
                    vector<vector<int>> components;
                    components = find_Connected_Components(subgraph, s);
                    if(components.empty()){
                        outfile << 0 << "\n";
                    }
                    else{
                        outfile << 1 << "\n";
                        outfile << components.size() << "\n";
                        for(int i = 0; i < components.size(); i++){
                            for(int j = 0; j < components[i].size(); j++){
                                if(j < components[i].size() - 1){
                                    outfile << components[i][j] << " ";
                                }
                                else{
                                    outfile << components[i][j];
                                }
                            }
                            outfile << "\n";
                        }
                    }   
                }    
            }
        }
        else if(task == 2){
            
            unordered_map<int, vector<int>> influence;
            influence = influencer_vertices(subgraph, end_k, p_ego);

            if(verbose == 0){
                outfile << influence.size() << "\n";
                int isize = 0;
                if(!influence.empty()){
                    for(auto a : influence){
                        if(isize < influence.size() - 1){
                            outfile << a.first << " ";
                        }
                        else{
                            outfile << a.first;
                        }
                        isize++;
                    }
                }
            }
            else if(verbose == 1){
                int isize = influence.size();
                outfile << isize << "\n";
                if(isize != 0){
                    int count = 0;
                    for(auto a : influence){
                        outfile << a.first << "\n";
                        int linesize = 0;
                        for(auto b : a.second){
                            if(linesize < a.second.size() - 1){
                                outfile << b << " ";
                            }
                            else{
                                if(count < n - 1){
                                    outfile << b << "\n";
                                }
                                else{
                                    outfile << b;
                                }
                            }
                            linesize++;
                        }
                        count++;
                    }
                }
            }
        }
    }

    outfile.close();
    
    MPI_Finalize();
    
    return 0;

}