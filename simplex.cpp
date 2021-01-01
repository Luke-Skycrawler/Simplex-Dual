#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <iomanip>
#include <fstream>
#define NDEBUG   // assertions enabled
using namespace std;
/********************************************
 * input specification
 * n, m #of variables, #of restrictions respectively
 * 'max'/'min' , 'primal'/'dual'
 * c0 ,..., cn-1
 * a01 ,..., a1_n-1, '>='|'<='|'=' , b0
 * ...
 * a_m-1_1 ,..., am-1_n-1, '>='|'<='|'=' , bm-1
*********************************************
TODO: simple Euler elimination, disregard Gauss-Siedel in the beginning 
*********************************************/
static int m,n;
struct simplex_table_collum{
    static int M;   // number of all variables, equal to the vector size
    int var;        // index of the non-base variable
    vector<double> &a;
    double b;
    simplex_table_collum(vector<double> &a=*new vector<double>,double &b=*new double(0.0),int var=-1):a(a),b(b),var(var){}
    void replace_with(int j){
        assert(a[j]!=0);
        var=j;
    }
    void normalize(){
        *this/=a[var];
    }
    
    simplex_table_collum& operator *(double k) const{
        if(tmp_vec)delete tmp_vec;
        if(tmp_double)delete tmp_double;
        if(tmp)delete tmp;
        tmp_vec=new vector<double>(a);
        tmp_double=new double(b);
        tmp=new simplex_table_collum(*tmp_vec,*tmp_double,-1);
        // FIXME: #1 CLOSED (but not beautifully) leakage needed to tend to, don't know how
        for(int i=0;i<M;i++)tmp->a[i]*=k;
        tmp->b*=k;
        return *tmp;
    }
    // simplex_table_collum& operator +(const simplex_table_collum &j) const{
    // }
    simplex_table_collum& operator -=(const simplex_table_collum &j){
        for(int i=0;i<M;i++)a[i]-=j.a[i];
        b-=j.b;
        return *this;
    }
    simplex_table_collum& operator /=(double aij){
        for(int i=0;i<M;i++)a[i]/=aij;
        b/=aij;
        return *this;
    }
    private:
    static simplex_table_collum *tmp;
    static vector<double> *tmp_vec;
    static double *tmp_double;
};
int simplex_table_collum::M;
simplex_table_collum* simplex_table_collum::tmp=NULL;
vector<double>* simplex_table_collum::tmp_vec=NULL;
double* simplex_table_collum::tmp_double=NULL; 
static vector<simplex_table_collum> collums;
static int check(vector<simplex_table_collum> &collums,simplex_table_collum &z);
static int select_base_in(vector<simplex_table_collum> &collums,simplex_table_collum &z);
static int select_base_out(int _in,vector<simplex_table_collum> &collums);
int main(int argc,char **argv){
    ifstream _input;
    #ifdef UNIT_TEST
    _input.open(argv[1]);
    #else
    _input.open("test2.0.txt");
    #endif
    _input>>n>>m;
    vector<double> *a=new vector<double>[m],c;
    // FIXME: #0 target is min
    double *b=new double[m];
	c.reserve(m+n);
    c.resize(0);
    // simplex_table_collum::M=m+n;
    collums.reserve(m);
    for(int i=0;i<m+n;i++){
        double tmp;
        if(i<n){
            _input>>tmp;
            c.push_back(tmp);
        }
        else c.push_back(0.0);
    }
    simplex_table_collum z(c,*new double(0.0),-1);
    for(int j=0;j<m;j++)a[j].resize(m+n);
    for(int j=0;j<m;j++){
        int op;
        // a[j].reserve(m+n);
        for(int i=0;i<n;i++){
            double tmp;
            _input>>tmp;
            a[j][i]=tmp;
        }
        assert(op!=1);
        _input>>b[j]>>op;
        if(op==-1){
            a[j][m+j]=1.0;
            collums.push_back(simplex_table_collum(a[j],b[j],j+m));
        }
        else if(op==0){
            a[j][m+j]=0.0;
            assert(a[j][j]!=0);
            // TODO: select a base variable
            collums.push_back(simplex_table_collum(a[j],b[j],j));
            collums[j].normalize();
        }
    }
    int e;
    for(int i=0;i<n;i++){
        _input>>e;
        assert(e==1);
    }
    _input.close();
    simplex_table_collum::M=m+n;
/********************************************
 * end of input section
*********************************************/
    int ret,iter=0;
    do{
        if(ret=check(collums,z))break;
        int _in=select_base_in(collums,z);
        int _out=select_base_out(_in,collums);
        collums[_out].replace_with(_in);
        for(int j=0;j<m;j++){
            if(j==_out)continue;
            double coeff = collums[j].a[_in]/collums[_out].a[_in];
            collums[j]-=collums[_out]*coeff;
        }
        z-=collums[_out]*(z.a[_in]/collums[_out].a[_in]);
        collums[_out]/=collums[_out].a[_in];

        cout<<"iteration #"<<iter++<<endl;
        #ifdef NDEBUG
        cout<<"swaped in: "<<_in<<"\tswaped out: "<<_out<<endl;
        for(int j=0;j<m;j++){
            for(int i=0;i<m+n;i++)cout<<setw(5)<<collums[j].a[i];
            cout<<"|\t"<<collums[j].b<<endl;
        }
        cout<<"************************************"<<endl;
        #endif
    }
    while(1);
    switch (ret){
        case 1:
            cout<<"max z="<<-z.b<<endl;
            for(int j=0;j<m;j++){       // FIXME: #2 use rank instead of m to be more exact
                cout<<"x"<<collums[j].var<<"="<<setw(5)<<collums[j].b/collums[j].a[collums[j].var]<<endl;
            }                
            break;
        case -1:cout<<"infinate maximum"<<endl;break;
        // TODO: notify which variable is the cause
        case -2:cout<<"no solution"<<endl;break;
        default:break;
    }
    return 0;
}
static int check(vector<simplex_table_collum> &collums,simplex_table_collum &z){
    for(int i=0;i<m+n;i++)if(z.a[i]>0.0){
        // non-base variable
        int j=0;
        for(j=0;j<m;j++){
            if(collums[j].a[i]>0.0)break;
            else if(collums[j].b<0.0)
                return -2;      // no solution
            // FIXME: #3 possible rank reduction, inspection needed
        }
        if(j==m)return -1;      // infinete optimal solution
        else{
            z.var=i;
            return 0;           // continue swapping out
        }
    }
    return 1;
}
inline static int select_base_in(vector<simplex_table_collum> &collums,simplex_table_collum &z){
    return z.var;               // refer to check for information
}
static int select_base_out(int i,vector<simplex_table_collum> &collums){
    static int min,index_min;
    int init=1;
    for(int j=0;j<m;j++){
        if(collums[j].a[i]<=0)continue;
        else {
            double threshold=collums[j].b/collums[j].a[i];
            if(init||threshold<min){
                min=threshold;
                index_min=j;
                if(init)init=0;
            }
        }
    }
    return index_min;
}
