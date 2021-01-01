#ifndef __PRE_PROCESS__
#define __PRE_PROCESS__

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

#define MAXN 201
#define MAXM 501

/********************************************
 * struct definitions
*********************************************/
struct outputInfo //Output Information
{
    int type;/*
            0: xj is exactly xj
            1: xj is actually -xj
            2: xj is actually xm - xn
            */
    int first_aux, second_aux; // m, n
};

struct preProcessResult
{
    int var_num; //total number of variables after standardize()
    int new_var_num; //number of newly added variables
    vector<double> *a;
    vector<double> b;
    vector<double> c;
    struct outputInfo *output_info;

    preProcessResult(\
        int _var_num,\
        int _new_var_num,\
        vector<double> *_a,\
        vector<double> &_b,\
        vector<double> &_c,\
        struct outputInfo *_output_info):\
            var_num(_var_num), new_var_num(_new_var_num),\
            a(_a), b(_b), c(_c), output_info(_output_info){}
};

/********************************************
 * function declarations
*********************************************/
static void readInput(string file_name);
static void standardize();
struct preProcessResult* preProcess(string file_name); //the function you should call only

#endif

// struct simplex_table_row{
// public:
//     static int M;   // number of all variables, equal to the vector size
//     int var;        // index of the non-base variable
//     vector<double> &a;
//     double b;
    
//     simplex_table_row(
//         vector<double> &a = *new vector<double>,\
//         double &b = *new double(0.0),\
//         int var=-1)\
//         :a(a), b(b), var(var){}
    
//     void replace_with(int j)
//     {
//         assert(a[j] != 0);
//         var = j;
//     }
    
//     void normalize()
//     {
//         *this /= a[var];
//     }
    
//     simplex_table_row& operator *(double k) const
//     {
//         if(tmp_vec)
//             delete tmp_vec;
//         if(tmp_double)
//             delete tmp_double;
//         if(tmp)
//             delete tmp;
//         tmp_vec = new vector<double>(a);
//         tmp_double = new double(b);
//         tmp = new simplex_table_row(*tmp_vec, *tmp_double, -1);
//         // FIXME: #1 CLOSED (but not beautifully) leakage needed to tend to, don't know how
//         for(int i = 0; i < M; i++)
//             tmp->a[i] *= k;
//         tmp->b *= k;
//         return *tmp;
//     }
    
//     // simplex_table_row& operator +(const simplex_table_row &j) const{
//     // }
    
//     simplex_table_row& operator -=(const simplex_table_row &j)
//     {
//         for(int i = 0; i < M; i++)
//             a[i] -= j.a[i];
//         b -= j.b;
//         return *this;
//     }
    
//     simplex_table_row& operator /=(double aij)
//     {
//         for(int i = 0; i < M; i++)
//             a[i] /= aij;
//         b /= aij;
//         return *this;
//     }
    
// private:
//     static simplex_table_row *tmp;
//     static vector<double> *tmp_vec;
//     static double *tmp_double;
// };