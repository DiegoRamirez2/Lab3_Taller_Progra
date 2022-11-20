#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(){
    vector<vector<int> > v;
    v.push_back({1,2,3});
    v.push_back({4,5,6});
    v.push_back({7,8,9});

    for(int i=0; i<v.size(); i++){
        for(int j=0; j<v[i].size(); j++){
            cout << v[i][j] << " ";
        }
        cout << endl;
    }
}