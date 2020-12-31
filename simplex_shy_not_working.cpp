#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <iomanip>
// #define NDEBUG   // assertions enabled
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
#define M (m+n)
struct simplex_table_collum{
    // static int M;   // number of all variables, equal to the vector size
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
        simplex_table_collum *tmp=new simplex_table_collum(*new vector<double>(a),*new double(b),-1);
        // FIXME: leakage needed to tend to, don't know how
        for(int i=0;i<M;i++)tmp->a[i]*=k;
        tmp->b*=k;
        return *tmp;
    }
    // simplex_table_collum& operator +(const simplex_table_collum &j) const{
    // }
    simplex_table_collum& operator -=(const simplex_table_collum &j){
        for(int i=0;i<M;i++)a[i]+=j.a[i];
        b+=j.b;
        return *this;
    }
    simplex_table_collum& operator /=(double aij){
        for(int i=0;i<M;i++)a[i]/=aij;
        b/=aij;
        return *this;
    }
};
static vector<simplex_table_collum> collums;
static int check(vector<simplex_table_collum> &collums,simplex_table_collum &z);
static int select_base_in(vector<simplex_table_collum> &collums,simplex_table_collum &z);
static int select_base_out(int _in,vector<simplex_table_collum> &collums);
int main(void){
    string arg1,arg2;
    cin>>n>>m;
    cin>>arg1>>arg2;
    assert(arg1=="max");
    assert(arg2=="primal");
    vector<double> *a=new vector<double>[m],c;
    double *b=new double[m];
    c.resize(m+n);
    // simplex_table_collum::M=m+n;
    collums.reserve(m);
    for(int i=0;i<m+n;i++){
        double tmp;
        if(i<n){
            cin>>tmp;
            c.push_back(tmp);
        }
        else c.push_back(0.0);
    }
    simplex_table_collum z(c,*new double(0.0),-1);
    for(int j=0;j<m;j++)a[j].reserve(m+n);
    for(int j=0;j<m;j++){
        string op;
        // a[j].reserve(m+n);
        for(int i=0;i<n;i++){
            double tmp;
            cin>>tmp;
            a[j][i]=tmp;
        }
        cin>>op;
        assert(op!=">=");
        cin>>b[j];
        if(op=="<="){
            a[j][m+j]=1.0;
            collums.push_back(simplex_table_collum(a[j],b[j],j+m));
        }
        else if(op=="="){
            a[j][m+j]=0.0;
            assert(a[j][j]!=0);
            // TODO: select a base variable
            collums.push_back(simplex_table_collum(a[j],b[j],j));
            collums[j].normalize();
        }
    }
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
            for(int j=0;j<m;j++){       // FIXME: use rank instead of m to be more exact
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
            // FIXME: possible rank reduction, inspection needed
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
