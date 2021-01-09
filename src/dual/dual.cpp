#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <map>
#include <float.h>
#include "PreProcess_dual.h"
#define EPS 5E-5
#define INF (DBL_MAX)
#define INFINITY_SOLUTION 0
#define NO_SOLUTION -1
#define NORMAL 1
using namespace std;


struct simplex_table_row{
    static int M;   // number of all variables, equal to the vector size
    int var;        // index of the non-base variable
    vector<double> &a;
    double &b;
    
    simplex_table_row(vector<double> &a = *new vector<double>,\
                      double &b = *new double(0.0),\
                      int var = -1)\
                      :a(a) ,b(b), var(var){}
    
    simplex_table_row(const simplex_table_row &s)\
                      :a(*new vector<double>),\
                      b(*new double(s.b)),\
                      var(s.var)
    {
        for(int i=0;i<M;i++)a.push_back(s.a[i]);
    }

    void replace_with(int j)
    {
        assert(a[j]!=0);
        var=j;
    }

    void normalize()
    {
        *this/=a[var];
    }

    simplex_table_row& operator *(double k) const
    {
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
    
    simplex_table_row& operator =(const simplex_table_row &j)
    {
        a=j.a;
        b=j.b;
        return *this;
    }
    
    simplex_table_row& operator -=(const simplex_table_row &j)
    {
        for(int i=0;i<M;i++)a[i]-=j.a[i];
        b-=j.b;
        return *this;
    }
    
    simplex_table_row& operator /=(double aij)
    {
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
simplex_table_row* simplex_table_row::tmp = NULL;
vector<double>* simplex_table_row::tmp_vec = NULL;
double* simplex_table_row::tmp_double = NULL; 

static vector<simplex_table_row> rows_main;
static vector<simplex_table_row> rows_2phase;
simplex_table_row z_main, z_2phase;

vector<double> *a_main, *a_2phase;
vector<double> b_main, b_2phase;
vector<double> c_main, c_2phase;
static int m, n;

void printDebug(simplex_table_row& z, vector<simplex_table_row>& rows, int m);
int dualCore(simplex_table_row& z, vector<simplex_table_row>& rows, int m);
int beginDual(vector<double>* _a, vector<double>& _b, vector<double>& _c);
void printResult(simplex_table_row& z, vector<simplex_table_row>& rows, int n, int m, outputInfo* output_info, int output_num);

inline int lessThanZero(double t)
{
    return t < -EPS;
}
inline int greaterThanZero(double t)
{
    return t > EPS;
}
inline int equalZero(double t)
{
    return t >= -EPS && t <= EPS;
}
inline int lessThan(double a, double b)
{
    return a < (b - EPS);
}
inline int greaterThan(double a, double b)
{
    return a > (b + EPS);
}
inline int equal(double a, double b)
{
    return (!lessThan(a, b)) && (!greaterThan(a, b));
}

int main(int argc,char **argv)
{
    //get input
    struct preProcessDualResult* input; //a pointer to standardized input
    input = preProcessDual(argv[1]);
    n = input->var_num;
    m = input->constraint_num;

    //If not dual-feasible, we do not use dual. (why? it might be dual-feasible after transformation)
    for(int i = 1; i <= n; i++)
    {
        if(greaterThanZero(input->c[i]))
        {
            cout << "Input is not dual-feasible! Please use the **primal** simplex method." << endl;
            return -1;
        }
    }

    //Solve the original problem
    int dual_ret = beginDual(input->a, input->b, input->c);
    switch(dual_ret)
    {
        case NO_SOLUTION:
            cout << "-1" << endl;
            break;
        case NORMAL:
            cout << "1" << endl;
            printResult(z_main, rows_main, n, m, input->output_info, input->output_num);
            break;
        case INFINITY_SOLUTION:
            cout << "0" << endl;
            break;
    }
    return 0;
}

int beginDual(vector<double>* _a, vector<double>& _b, vector<double>& _c)
{
    /** Initialization **/
    simplex_table_row::M = m + n;
    //a, b, c
    //`a_main`: a (m)*(n + m) matrix
    //`b_main`: a (m) vector
    //`c_main`: a (n + m) vector
    a_main = new vector<double>[m];
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < n; j++) //copy
            a_main[i].push_back(_a[i + 1][j + 1]);
        for(int k = 0; k < m; k++) //slack variables
            a_main[i].push_back(i==k);
    }
    for(int i = 0; i < m; i++) //copy
        b_main.push_back(_b[i + 1]);
    for(int j = 0; j < n; j++) //copy
        c_main.push_back(_c[j + 1]);
    for(int k = 0; k < m; k++) //slack variables
        c_main.push_back(0.0);

    for(int i = 0; i < m; i++)
    {
        /** Push each line into the `rows`,
         * and the initial base for row `i`
         * is exactly the added slack variable **/
        rows_main.push_back(simplex_table_row (a_main[i], b_main[i], m + i));
    }
    z_main = simplex_table_row(c_main, *new double(0.0));

    // printDebug(z_main, rows_main, m);
    return dualCore(z_main, rows_main, m);
}



int dualCore(simplex_table_row& z, vector<simplex_table_row>& rows, int m)
{
    while(true)
    {
        int outBase = -1, inBase = -1, outBaseLine = -1;
        double b_max = 0;
        for(int i = 0; i < m; i++)
        {
            if(lessThan(rows[i].b, b_max))
            {
                b_max = rows[i].b;
                outBase = rows[i].var;
                outBaseLine = i;
            }
            else if(equal(rows[i].b, b_max) && rows[i].var < outBase)
            {
                b_max = rows[i].b;
                outBase = rows[i].var;
                outBaseLine = i;
            }
        }
        if(outBase == -1)
        {
            // cout << "dualCore finished." << endl; 
            // printDebug(z, rows, m);
            return NORMAL;
        }
    
        double cdiva_min = INF;
        for(int j = 0; j < simplex_table_row::M; j++)
        {
            if(greaterThanZero(rows[outBaseLine].a[j]) || equalZero(rows[outBaseLine].a[j]))
                continue;
            double div = z.a[j] / rows[outBaseLine].a[j];
            if(lessThan(div, cdiva_min))
            {
                cdiva_min = div;
                inBase = j;
            }
        }
        if(inBase == -1)
        {
            return NO_SOLUTION;
        }

        //Change base
        rows[outBaseLine].replace_with(inBase);
        rows[outBaseLine].normalize();
        //Transform `z`
        for(int j = 0; j < simplex_table_row::M; j++)
        {
            double coeff = z.a[inBase];
            z -= rows[outBaseLine] * coeff;
        }
        //Transform other rows
        for(int i = 0; i < m; i++)
        {
            if(i == outBaseLine) continue; //skip itself
            double coeff = rows[i].a[inBase];
            rows[i] -= rows[outBaseLine] * coeff;
        }
    }
    
    return NORMAL;
}

void printResult(simplex_table_row& z, vector<simplex_table_row>& rows, int n, int m, outputInfo* output_info, int output_num)
{
    cout << fixed << setprecision(6) << z_main.b  << endl;
    vector<double> solution(n, 0);
    for(int i = 0; i < m; i++)
    {
        solution[rows[i].var] = rows[i].b;
    }

    if(output_info[1].type == 0)
    {
        cout << fixed << setprecision(6) << solution[0];
    }
    else if(output_info[1].type == 1)
    {
        cout << fixed << setprecision(6) << -solution[0];
    }
    else if(output_info[1].type == -2)
    {
        cout << fixed << setprecision(6) << solution[output_info[1].first_aux - 1] - solution[output_info[1].second_aux - 1];
    }
    for(int i = 2; i <= output_num; i++)
    {
        if(output_info[i].type == 0)
        {
            cout << fixed << setprecision(6) << " " << solution[i - 1];
        }
        else if(output_info[i].type == 1)
        {
            cout << fixed << setprecision(6) << " " << -solution[i - 1];
        }
        else if(output_info[i].type == -2)
        {
            cout << fixed << setprecision(6) << " " << solution[output_info[i].first_aux - 1] - solution[output_info[i].second_aux - 1];
        }
    }
    cout << endl;
}

void printDebug(simplex_table_row& z, vector<simplex_table_row>& rows, int m)
{
    cout << "/************ DEBUG ************/" << endl;
    cout << "Target: max z = " << endl;
    for(int j = 0; j < simplex_table_row::M; j++)
    {
        cout << z.a[j] << " ";
    }
    cout << "| "  << z.b << endl;
    cout << "------Constraints------" << endl;
    for(int i = 0; i < m; i++)
    {
        for(int j = 0; j < simplex_table_row::M; j++)
            cout << rows[i].a[j] << " ";
        cout << "| " << rows[i].b;
        cout << "\t(Base: x" << rows[i].var << ")" << endl;
    }
    cout << "/*******************************/" << endl << endl;
}
