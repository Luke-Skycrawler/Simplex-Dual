#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <map>
// #define NDEBUG   // assertions enabled
// #define CASE__NDEBUG
#define UNIT_TEST
#define EPS 5e-7
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
TODO: FIXME: #6 fix peculiar cases: nexample
*********************************************/
static int m,n;
struct simplex_table_row{
    static int M;   // number of all variables, equal to the vector size
    int var;        // index of the non-base variable
    vector<double> &a;
    double &b;
    simplex_table_row(vector<double> &a=*new vector<double>,double &b=*new double(0.0),int var=-1):a(a),b(b),var(var){}
    void replace_with(int j){
        assert(a[j]!=0);
        var=j;
    }
    void normalize(){
        *this/=a[var];
    }
    
    simplex_table_row& operator *(double k) const{
        if(tmp_vec)delete tmp_vec;
        if(tmp_double)delete tmp_double;
        if(tmp)delete tmp;
        tmp_vec=new vector<double>(a);
        tmp_double=new double(b);
        tmp=new simplex_table_row(*tmp_vec,*tmp_double,-1);
        // FIXME: #1 CLOSED (but not beautifully) leakage needed to tend to, don't know how
        for(int i=0;i<M;i++)tmp->a[i]*=k;
        tmp->b*=k;
        return *tmp;
    }
    // simplex_table_row& operator +(const simplex_table_row &j) const{
    // }
    simplex_table_row& operator =(const simplex_table_row &j){
        a=j.a;
        b=j.b;
        return *this;
    }
    simplex_table_row& operator -=(const simplex_table_row &j){
        for(int i=0;i<M;i++)a[i]-=j.a[i];
        b-=j.b;
        return *this;
    }
    simplex_table_row& operator /=(double aij){
        for(int i=0;i<M;i++)a[i]/=aij;
        b/=aij;
        return *this;
    }
    private:
    static simplex_table_row *tmp;
    static vector<double> *tmp_vec;
    static double *tmp_double;
};
int simplex_table_row::M;
simplex_table_row* simplex_table_row::tmp=NULL;
vector<double>* simplex_table_row::tmp_vec=NULL;
double* simplex_table_row::tmp_double=NULL; 
static vector<simplex_table_row> rows;
static int check(vector<simplex_table_row> &rows,simplex_table_row &z);
static int select_base_in(vector<simplex_table_row> &rows,simplex_table_row &z);
static int select_base_out(int _in,vector<simplex_table_row> &rows);
int main(int argc,char **argv){
    ifstream _input;
    #ifdef UNIT_TEST
    _input.open(argv[1]);
    #else
    _input.open("adlittle.txt");
    #endif
    _input>>n>>m;
    vector<double> *a=new vector<double>[m],c;
    // FIXME: #0 CLOSED target is min
    // FIXME: #5 CLOSED output format
    double *b=new double[m];
	c.reserve(2*m+n);
    c.resize(0);
    // simplex_table_row::M=m+n;
    rows.reserve(m);
    for(int i=0;i<m+n;i++){
        double tmp;
        if(i<n){
            _input>>tmp;
            c.push_back(tmp);
        }
        else c.push_back(0.0);
    }
    simplex_table_row::M=m+n;
    simplex_table_row z(c,*new double(0.0),-1);
    for(int j=0;j<m;j++){
        a[j].reserve(2*m+n);
        a[j].resize(m+n);
    }
    for(int j=0;j<m;j++){
        int op;
        // a[j].reserve(m+n);
        for(int i=0;i<n;i++){
            double tmp;
            _input>>tmp;
            a[j][i]=tmp;
        }
        _input>>b[j]>>op;
        if(op==-1||op==1){
            a[j][n+j]=-op;
            rows.push_back(simplex_table_row(a[j],b[j],j+n));
            if(op==1)
                rows[j].normalize();
        }
        else if(op==0){
            a[j][m-1+j]=0.0;
            rows.push_back(simplex_table_row(a[j],b[j],-1));
        }
    }
    int e;
    map<int,int> counterpart;
    for(int i=0;i<n;i++){
        _input>>e;
        if(e==-1)
            for(int j=0;j<m;j++){
                rows[j].a.push_back(-rows[j].a[i]);
                rows[j].a[i]=0.0;
                if(rows[j].var==i){
                    rows[j].var=simplex_table_row::M;
                    rows[j].normalize();
                }
                counterpart[simplex_table_row::M++]=i;
                // FIXME: CLOSED output should know this change in sign
            }
        else if(e==0){
            for(int j=0;j<m;j++)
                rows[j].a.push_back(-rows[j].a[i]);
            z.a.push_back(-z.a[i]);
            counterpart[simplex_table_row::M++]=i;
        }
    }
    _input.close();
/********************************************
 * end of input section
/********************************************
 * select the base variable bundle
 * naive, without pivoting
 * numerical error could be overwhelming
*********************************************/
    int base_cnt;
    for(int i=0;i<n;i++){
        base_cnt=0;
        for(int j=0;j<m;j++){
            if(rows[j].var!=-1)base_cnt++;
            if(rows[j].var!=-1||fabs(rows[j].a[i])<EPS)continue;
            rows[j].var=i;
            base_cnt++;
            for(int k=0;k<m;k++){
                if(j==k)continue;
                rows[k]-=rows[j]*(rows[k].a[i]/rows[j].a[i]);
            }
            z-=rows[j]*(z.a[i]/rows[j].a[i]);
            rows[j].normalize();
            break;
        }
        if(base_cnt==m)break;
    }
    base_cnt=0;
    for(int j=0;j<m;j++)if(rows[j].var!=-1)base_cnt++;
    int ret=0;
    if(base_cnt<m){
        // search for zero row
        for(int j=0;j<m;j++){
            int i;
            for(i=0;i<n;i++)if(fabs(rows[j].a[i])>EPS)break;
            if(i==n){
                if(rows[j].b==0.0)
                    rows[j--]=rows[--m];    //elegant
                else ret=-2;
            }
        }
    }
    // TODO: examine if any b[i]<0 with var!=-1
            // for(int j=0;j<m;j++)
            //     if(rows[j].b<0.0)
            //         ret=-2;
            // // FIXME: selecting a feasible base not that easy
    // facing bj<0 scenario
/********************************************
 * feasible solution not built yet (with some negative base variables)
********************************************/

    int iter=0;
    do{
        int _in,_out,j;
        if(ret)break;
        #ifdef CASE__NDEBUG
        try{
            for(j=0;j<m;j++){
                if(b[j]<-EPS){
                    cout<<j<<" "<<b[j]<<" "<<iter<<endl;
                }
                assert(b[j]>=-EPS);
                // assert(rows[j].a[rows[j].var]==1.0);
            }
        }
        catch(...){
            cout<<"error"<<endl;
        }
        #endif
        for(j=0;j<m;j++)if(rows[j].b<-EPS){
            int k;
            for(k=0;k<(simplex_table_row::M);k++)
                if(rows[j].a[k]<-EPS)break;
            if(k==simplex_table_row::M)
                ret=-2;      // no solution if all the variables have positive coeffients
            else{
                _in=k;
                _out=j;
                ret=0;
            } 
            break;
        }
        if(ret)break;
        if(j==m){
            if(ret=check(rows,z))break;
            _in=select_base_in(rows,z);
            _out=select_base_out(_in,rows);
        }
        rows[_out].replace_with(_in);
        for(j=0;j<m;j++){
            if(j==_out)continue;
            double coeff = rows[j].a[_in]/rows[_out].a[_in];
            rows[j]-=rows[_out]*coeff;
        }
        z-=rows[_out]*(z.a[_in]/rows[_out].a[_in]);
        rows[_out]/=rows[_out].a[_in];
        cout<<"iteration #"<<iter++<<"  "<<-z.b<<endl;
        #ifdef NDEBUG
        cout<<"swaped in: "<<_in<<"\tswaped out: "<<_out<<endl;
        for(int j=0;j<m;j++){
            for(int i=0;i<simplex_table_row::M;i++)cout<<setw(5)<<rows[j].a[i];
            cout<<"|\t"<<rows[j].b<<endl;
        }
        cout<<"************************************"<<endl;
        #endif
    }
    while(1);
    switch (ret){
        case 1:
            cout<<"min z="<<(-z.b)<<endl;   // FIXME: #4 CLOSED output -0
            for(int j=0;j<m;j++){       // FIXME: CLOSED #2 use rank instead of m to be more exact
                cout<<"x"<<rows[j].var<<"="<<setw(5)<<rows[j].b/rows[j].a[rows[j].var]<<endl;
            }                
            break;
        case -1:ret++;cout<<"infinate maximum"<<endl;break;
        // TODO: notify which variable is the cause
        case -2:ret++;cout<<"no solution"<<endl;break;
        default:break;
    }
    #ifdef UNIT_TEST
        string _output_filename(argv[1]);
        _output_filename.replace(_output_filename.find("txt"),3,"out1");
        ofstream fout;
        fout.open(_output_filename);
        fout<<ret<<endl;
        if(ret==1){
            if(z.b<EPS&&z.b>=0.0)fout<<0<<endl;
            else 
                fout<<fixed<<setprecision(6)<<-z.b<<endl;
            double* map=new double[n];
            for(int j=0;j<m;j++){
                if(rows[j].var<n)
                    map[rows[j].var]+=rows[j].b/rows[j].a[rows[j].var];
                else if(counterpart.count(rows[j].var))
                    map[counterpart[rows[j].var]]-=rows[j].b/rows[j].a[rows[j].var];
            }
            for(int i=0;i<n;i++){
                if(map[i]<0.0&&map[i]>-EPS)map[i]=0.0;
                fout<<fixed<<setprecision(6)<<map[i]<<" ";
            }
            delete []map;
        }
        fout.close();
    #endif
    return 0;
}
static int check(vector<simplex_table_row> &rows,simplex_table_row &z){
    for(int i=0;i<simplex_table_row::M;i++)if(z.a[i]<-EPS){   // NOTE: #0 modified '<' to change min/max
        // non-base variable
        int j=0;
        for(j=0;j<m;j++){
            if(rows[j].a[i]>EPS)break;
            // FIXME: #3 CLOSED possible rank reduction, inspection needed
        }
        if(j==m)return -1;      // infinete optimal solution
        else{
            z.var=i;
            return 0;           // continue swapping out
        }
    }
    return 1;
}
inline static int select_base_in(vector<simplex_table_row> &rows,simplex_table_row &z){
    return z.var;               // refer to check for information
}
static int select_base_out(int i,vector<simplex_table_row> &rows){
    static int min,index_min;
    int init=1;
    for(int j=0;j<m;j++){
        if(rows[j].a[i]<=EPS)continue;
        else {
            double threshold=rows[j].b/rows[j].a[i];
            if(init||threshold<min){
                min=threshold;
                index_min=j;
                if(init)init=0;
            }
        }
    }
    return index_min;
}
