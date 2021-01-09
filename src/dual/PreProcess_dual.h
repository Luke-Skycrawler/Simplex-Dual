#ifndef __PRE_PROCESS_DUAL__
#define __PRE_PROCESS_DUAL__

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <assert.h>

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

struct preProcessDualResult
{
    int var_num; //total number of variables after standardize()
    int constraint_num;
    vector<double> *a;
    vector<double> b;
    vector<double> c;
    struct outputInfo *output_info;
    int output_num;

    preProcessDualResult(\
        int _var_num,\
        int _constraint_num,\
        vector<double> *_a,\
        vector<double> &_b,\
        vector<double> &_c,\
        struct outputInfo *_output_info,
        int _output_num):\
            var_num(_var_num), constraint_num(_constraint_num),\
            a(_a), b(_b), c(_c), output_info(_output_info), output_num(_output_num){}
};

/********************************************
 * function declarations
*********************************************/
static void readInput(string file_name);
static void standardize();
struct preProcessDualResult* preProcessDual(string file_name); //the function you should call only

#endif