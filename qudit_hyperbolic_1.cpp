//By Yi Jiang, April 2021
#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<cmath>
#include<random>
#include<string>
#include<string.h>
#include <sys/types.h>
#include <algorithm>
#include<climits>
using namespace std;
    int d=7;//degree
    int dd=5;//degree of dual graph
    string latticelist[]={"5_7_1260.mtx","5_7_6090.mtx"};
    string duallatticelist[]={"7_5_1260.mtx","7_5_6090.mtx"};

int** sb;//bonds list for spins
int** bs;//spins list for bonds
//vector<vector<char> > bc;//coefficient
int ns,nb,nc=0;

bool read_incidence(string readname){
    if(readname.length()){
        ifstream mtxin(readname.c_str());
        if(!mtxin.good()){
            cout<<"error opening file"<<endl;
            return false;
        }
        else{
            while (mtxin.peek() == '%') mtxin.ignore(2048, '\n');
            mtxin>>ns>>nb>>nc;
            sb=new int*[ns];
            bs=new int*[nb];
            for(int i=0;i<ns;++i){
                sb[i]=new int[d];
                for(int j=0;j<d;++j)
                    sb[i][j]=-1;
            }
            for(int i=0;i<nb;++i){
                bs[i]=new int[2]{-1,-1};
            }
            for(int i=0;i<nc;++i){
                int s,b,c;
                mtxin>>s>>b>>c;
                for(int i=0;i<d;++i){
                    if(sb[s-1][i]>=0)
                        continue;
                    sb[s-1][i]=b;
                    break;
                }
                for(int i=0;i<2;++i){
                    if(bs[b-1][i]>=0)
                        continue;
                    bs[b-1][i]=s;
                    break;
                }
            }
            mtxin.close();
        }
        return true;
    }
    else{
        cout<<"error in function \"read_incidence\""<<endl;
        return false;
    }
    /*
    d=4;
    lattsize="square_"+to_string(l);
    sb=new int*[l*l];
    bs=new int*[l*l*2];
    for(int i=0;i<l*l;++i){
        sb[i]=new int[d];
    }
    for(int i=0;i<l*l*2;++i){
        bs[i]=new int[2];
    }
    for(int i=0;i<l;++i){
        for(int j=0;j<l;++j){
            sb[i*l+j][0]=i*l+j;
            sb[i*l+j][1]=i*l+(j+l-1)%l;
            sb[i*l+j][2]=(i+l-1)%l*l+j+l*l;
            sb[i*l+j][3]=i*l+j+l*l;
        }
    }
    for(int i=0;i<l;++i){
        for(int j=0;j<l;++j){
            bs[i*l+j][0]=i*l+j;
            bs[i*l+j][1]=i*l+(j+1)%l;
            bs[i*l+j+l*l][0]=i*l+j;
            bs[i*l+j+l*l][1]=(i+1)%l*l+j;
        }
    }
    ns=l*l;
    nb=l*l*2;
    return true;
    */
}
int** fb;//bonds list for faces
int** bf;//faces list for bonds
int nf;//# of faces
bool read_incidence_d(string readname){
    if(readname.length()){
        ifstream mtxin(readname.c_str());
        if(!mtxin.good()){
            cout<<"error opening file"<<endl;
            return false;
        }
        else{
            while (mtxin.peek() == '%') mtxin.ignore(2048, '\n');
            mtxin>>nf>>nb>>nc;
            fb=new int*[nf];
            bf=new int*[nb];
            for(int i=0;i<nf;++i){
                fb[i]=new int[dd];
                for(int j=0;j<dd;++j)
                    fb[i][j]=-1;
            }
            for(int i=0;i<nb;++i){
                bf[i]=new int[2]{-1,-1};
            }
            for(int i=0;i<nc;++i){
                int s,b,c;
                mtxin>>s>>b>>c;
                for(int i=0;i<dd;++i){
                    if(fb[s-1][i]>=0)
                        continue;
                    fb[s-1][i]=b;
                    break;
                }
                for(int i=0;i<2;++i){
                    if(bf[b-1][i]>=0)
                        continue;
                    bf[b-1][i]=s;
                    break;
                }
            }
            mtxin.close();
        }
        return true;
    }
    /*
    d=4;
    lattsize="square_"+to_string(l);
    sb2=new int*[l*l];
    bs2=new int*[l*l*2];
    for(int i=0;i<l*l;++i){
        sb2[i]=new int[d];
    }
    for(int i=0;i<l*l*2;++i){
        bs2[i]=new int[2];
    }
    for(int i=0;i<l;++i){
        for(int j=0;j<l;++j){
            sb2[i*l+j][0]=i*l+j;
            sb2[i*l+j][1]=i*l+(j+l-1)%l;
            sb2[i*l+j][2]=(i+l-1)%l*l+j+l*l;
            sb2[i*l+j][3]=i*l+j+l*l;
        }
    }
    for(int i=0;i<l;++i){
        for(int j=0;j<l;++j){
            bs2[i*l+j][0]=i*l+j;
            bs2[i*l+j][1]=i*l+(j+1)%l;
            bs2[i*l+j+l*l][0]=i*l+j;
            bs2[i*l+j+l*l][1]=(i+1)%l*l+j;
        }
    }
    ns=l*l;
    nb=l*l*2;
    return true;*/
    else{
        cout<<"error in function \"read_incidence_d\""<<endl;
        return false;
    }
}
int delete_lattice(){
    for(int i=0;i<ns;++i)
        delete[] sb[i];
    for(int i=0;i<nb;++i)
        delete[] bs[i];
    for(int i=0;i<nf;++i)
        delete[] fb[i];
    for(int i=0;i<nb;++i)
        delete[] bf[i];
    delete[] sb,bs,fb,bf;
    return 0;
}
int main(){
    if(sizeof(latticelist)/sizeof(latticelist[0])!=sizeof(duallatticelist)/sizeof(duallatticelist[0])){
        cout<<"Number of lattices and dual lattices don't match"<<endl;
        return -1;
    }
    for(int k=0;k<sizeof(latticelist)/sizeof(latticelist[0]);k++){
        if(!read_incidence(latticelist[k])||!read_incidence_d(duallatticelist[k]))
            return -1;
        for(int i=0;i<nb;i++){
            bs[i][1]*=-1;
        }
        bf[0][1]*=-1;
        vector<pair<int,int> > ef_curr;//current layer of edge-face
        vector<pair<int,int> > ef_next;//next layer of edge-face
        //vector<int> faces_next;//next layer of faces
        
        ef_next.push_back(make_pair(1,bf[0][0]));
        ef_curr.push_back(ef_next[0]);
        
        int layercount=0;
        while(ef_next.size()){
            layercount++;
            if(layercount>=100){
                cout<<"too many layers"<<endl;
                return -1;
            }
            for(int i=0;i<ef_next.size();i++){
                pair<int,int> ef=ef_next[i];
                //current face is ef.second, work on all the other edges, create current layer.
                vector<int> edges;
                for(int j=0;j<dd;j++){
                    edges.push_back(fb[ef.second-1][j]);
                }
                pair<int,int> e1;
                if(abs(bf[ef.first-1][0])==ef.second){
                    if(bf[ef.first-1][0]<0){
                        e1=make_pair(abs(bs[ef.first-1][1]),abs(bs[ef.first-1][0]));
                    }
                    else if(bf[ef.first-1][1]<0){
                        e1=make_pair(abs(bs[ef.first-1][0]),abs(bs[ef.first-1][1]));
                    }
                    else{
                        cout<<"something wrong in part 000";
                        return -1;
                    }
                }
                else if(abs(bf[ef.first-1][1])==ef.second){
                    if(bf[ef.first-1][1]<0){
                        e1=make_pair(abs(bs[ef.first-1][1]),abs(bs[ef.first-1][0]));
                    }
                    else if(bf[ef.first-1][0]<0){
                        e1=make_pair(abs(bs[ef.first-1][0]),abs(bs[ef.first-1][1]));
                    }
                    else{
                        cout<<"something wrong in part 001";
                        return -1;
                    }
                }
                else{
                    cout<<"something wrong in part 002";
                    return -1;
                }
                for(int j=0;j<edges.size();j++){
                    if(abs(bs[edges[j]-1][0])==e1.first&&abs(bs[edges[j]-1][1])==e1.second || abs(bs[edges[j]-1][1])==e1.first&&abs(bs[edges[j]-1][0])==e1.second){
                        edges.erase(edges.begin()+j);
                    }
                }
                /*int loopcount=0;*/
                while(edges.size()){
                    /*loopcount++;
                    cout<<"edges.size() is "<<edges.size()<<endl;
                    int flag=0;*/
                    for(int j=0;j<edges.size();j++){
                        if(abs(bs[edges[j]-1][0])==e1.second || abs(bs[edges[j]-1][1])==e1.second){
                            /*flag=1;*/
                            if(abs(bs[edges[j]-1][0])==e1.second){//right order
                                e1.first=e1.second;
                                e1.second=abs(bs[edges[j]-1][1]);
                                if(abs(bf[edges[j]-1][0])==ef.second){
                                    if(bf[edges[j]-1][0]<0){
                                        cout<<"graph not orientable"<<endl;
                                        return -1;
                                        edges.erase(edges.begin()+j);
                                    }
                                    if(bf[edges[j]-1][1]>0){
                                        bf[edges[j]-1][1]*=-1;
                                        ef_curr.push_back(make_pair(edges[j],ef.second));
                                    }
                                }
                                else{
                                    if(bf[edges[j]-1][1]<0){
                                        cout<<"graph not orientable"<<endl;
                                        return -1;
                                        edges.erase(edges.begin()+j);
                                    }
                                    if(bf[edges[j]-1][0]>0){
                                        bf[edges[j]-1][0]*=-1;
                                        ef_curr.push_back(make_pair(edges[j],ef.second));
                                    }
                                }
                            }
                            else{//wrong order
                                e1.first=e1.second;
                                e1.second=abs(bs[edges[j]-1][0]);
                                if(abs(bf[edges[j]-1][0])==ef.second){
                                    if(bf[edges[j]-1][1]<0){
                                        cout<<"graph not orientable"<<endl;
                                        return -1;
                                        edges.erase(edges.begin()+j);
                                    }
                                    if(bf[edges[j]-1][0]>0){
                                        bf[edges[j]-1][0]*=-1;
                                        ef_curr.push_back(make_pair(edges[j],ef.second));
                                    }
                                }
                                else{
                                    if(bf[edges[j]-1][0]<0){
                                        cout<<"graph not orientable"<<endl;
                                        return -1;
                                        edges.erase(edges.begin()+j);
                                    }
                                    if(bf[edges[j]-1][1]>0){
                                        bf[edges[j]-1][1]*=-1;
                                        ef_curr.push_back(make_pair(edges[j],ef.second));
                                    }
                                }
                            }
                            edges.erase(edges.begin()+j);
                            break;
                        }
                    }
                    /*
                    if(flag==0){
                        cout<<"edge not found?"<<endl;
                    }
                    if(loopcount>dd){
                        cout<<"stuck in loop"<<endl;
                        return -1;
                    }
                    */
                }
            }
            //cout<<"ef_next.size() is "<<ef_next.size()<<endl;
            //cout<<"ef_curr.size() is "<<ef_curr.size()<<endl;
            //create next layer based on current layer
            ef_next.clear();
            for(int j=0;j<ef_curr.size();j++){
                if(abs(bf[ef_curr[j].first-1][0])==ef_curr[j].second){
                    ef_next.push_back(make_pair(ef_curr[j].first,abs(bf[ef_curr[j].first-1][1])));
                }
                else if(abs(bf[ef_curr[j].first-1][1])==ef_curr[j].second){
                    ef_next.push_back(make_pair(ef_curr[j].first,abs(bf[ef_curr[j].first-1][0])));
                }
            }
            ef_curr.clear();
        }
        //output new incidence matrix
        ofstream outfile;
        string outname="q_"+latticelist[k];
        outfile.open(outname);
        outfile<<"%%MatrixMarket matrix coordinate integer general"<<endl<<"%"<<endl<<ns<<" "<<nb<<" "<<2*nb<<endl;
        for(int i=0;i<ns;++i){
            for(int j=0;j<d;j++){
                outfile<<i+1<<" "<<sb[i][j]<<" ";
                if(abs(bs[sb[i][j]-1][0])==i+1){
                    if(bs[sb[i][j]-1][0]>0)
                        outfile<<1<<endl;
                    else
                        outfile<<-1<<endl;
                }
                else if(abs(bs[sb[i][j]-1][1])==i+1){
                    if(bs[sb[i][j]-1][1]>0)
                        outfile<<1<<endl;
                    else
                        outfile<<-1<<endl;
                }
            }
        }
        outfile.close();
        //output new dual incidence matrix
        outname="q_"+duallatticelist[k];
        outfile.open(outname);
        outfile<<"%%MatrixMarket matrix coordinate integer general"<<endl<<"%"<<endl<<nf<<" "<<nb<<" "<<2*nb<<endl;
        for(int i=0;i<nf;++i){
            for(int j=0;j<dd;j++){
                outfile<<i+1<<" "<<fb[i][j]<<" ";
                if(abs(bf[fb[i][j]-1][0])==i+1){
                    if(bf[fb[i][j]-1][0]>0)
                        outfile<<1<<endl;
                    else
                        outfile<<-1<<endl;
                }
                else if(abs(bf[fb[i][j]-1][1])==i+1){
                    if(bf[fb[i][j]-1][1]>0)
                        outfile<<1<<endl;
                    else
                        outfile<<-1<<endl;
                }
            }
        }
        outfile.close();
        delete_lattice();
    }
    return 0;
}
