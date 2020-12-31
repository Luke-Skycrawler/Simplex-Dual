#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#define MAXN 201
#define MAXM 501

using namespace std;

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
} output_info[MAXN];

/********************************************
 * global variables
*********************************************/
//input
int n, m;
int in_c[MAXN] = {0};
int in_a[MAXM][MAXN] = {0};
int in_b[MAXM] = {0}, in_d[MAXM] = {0}, in_e[MAXN] = {0};
//others
int var_num;
vector<double> *a;
vector<double> b;
vector<double> c;

/********************************************
 * function declarations
*********************************************/
void readInput(string file_name);
void standardize();

/********************************************
 * main()
*********************************************/
int main(int argc, char* argv[])
{
    
    if(argc != 2)
    {
        cout << "Usage: simplex [input-file]" << endl;
        return -1;
    }
    readInput(argv[1]);
    if(n < m)
    {
        cout << "[Error] Row rank is not full! Developing..." << endl;
        return -2;
    }
    standardize();
    return 0;
}


/********************************************
 * function definitions
*********************************************/
void readInput(string file_name)
{
    int i, j;
    ifstream in_file;
    in_file.open(file_name, ios::in);
    in_file >> n >> m;
    for(i = 1; i <= n; i++)
    {
        in_file >> in_c[i];
    }
    for(i = 1; i <= m; i++)
    {
        for(j = 1; j <= n; j++)
        {
            in_file >> in_a[i][j];
        }
        in_file >> in_b[i] >> in_d[i];
    }
    for(i = 1; i <= n; i++)
    {
        in_file >> in_e[i];
    }
    in_file.close();

    /******************** debug ********************/
    cout << endl << "Load problem from " << file_name << " successfully!\n**Original LP**" << endl;
    cout << "min z = ";
    cout << in_c[1] << "*[x_1]";
    for(i = 2; i <= n; i++)
    {
        cout << " + " << in_c[i] << "*[x_" << i << "]";
    }
    cout << endl;
    for(i = 1; i <= m; i++)
    {
        cout << in_a[i][1] << "*[x_1]";
        for(j = 2; j <= n; j++)
        {
            cout << " + " << in_a[i][j] << "*[x_" << j << "]";
        }
        if(in_d[i] < 0) 
            cout << " ≤ " << in_b[i] << endl;
        else if(in_d[i] == 0) 
            cout << " = " << in_b[i] << endl;
        else
            cout << " ≥ " << in_b[i] << endl;
    }
    for(j = 1; j <= n; j++)
    {
        if(in_e[j] == -1)
            cout << "[x_" << j << "]" << " ≤ 0" << endl;
        else if(in_e[j] == 1)
            cout << "[x_" << j << "]" << " ≥ 0" << endl;
        else
            cout << "[x_" << j << "]" << " <> 0" << endl;
    }
}

void standardize()
{
    int i, j;
    int new_var_num = 0; //record of newly added variables

    //initialize global variables
    a = new vector<double>[m + 1]; //malloc for matrix `a`
    b.push_back(0.0);
    c.push_back(0.0);

    /*  copy `in_a` -> `a`
        copy `in_b` -> `b`
        copy -`in_c` -> `c` (min -> max)
    */
    for(i = 1; i <= m; i++)
    {
        a[i].push_back(0.0);
        for(j = 1; j <= n; j++)
        {
            a[i].push_back(in_a[i][j]);
        }
        
        b.push_back(in_b[i]);
    }
    for(j = 1; j <= n; j++)
    {
        c.push_back(-in_c[j]);
    }

    //Check e_j
    for(j = 1; j <= n; j++)
    {
        if(in_e[j] == -1) //flip over c_j & a_ij, i = 1, 2, ..., m
        {
            c[j] *= -1;
            for(i = 1; i <= m; i++)
            {
                a[i][j] *= -1;
            }
            //mark on output info
            output_info[j].type = 1; 
        }
        else if(in_e[j] == 0) //transform x_j into 2 nonnegative auxiliary variables -- x_j = x_m - x_n
        {
            c.push_back(-c[j]);
            for(i = 1; i <= m; i++)
            {
                a[i].push_back(-a[i][j]);
            }
            new_var_num += 1; //record the newly added variable
            //mark on output info
            output_info[j].type = 2;
            output_info[j].first_aux = j;
            output_info[j].second_aux = n + new_var_num;
        }
        else
        {
            output_info[j].type = 0;
        }
    }

    //Check d_i
    for(i = 1; i <= m; i++)
    {
        if(in_d[i] == -1) // <=; add a nonnegative auxiliary to the left 
        {
            c.push_back(0.0);
            for(int k = 1; k <= m; k++)
            {
                a[k].push_back((k==i)); //for row i, new coefficient `1` is pushed
                                       //otherwise, `0` is pushed
            }
            new_var_num += 1; //record the newly added variable
        }
        else if(in_d[i] == 1) // >=; substract a nonnegative auxiliary to the left 
        {
            c.push_back(0.0);
            for(int k = 1; k <= m; k++)
            {
                a[k].push_back(-(k==i)); //for row i, new coefficient `-1` is pushed
                                        //otherwise, `0` is pushed
            }
            new_var_num += 1; //record the newly added variable
        }
    }

    var_num = n + new_var_num; //update global record for variables' number

    /******************** debug ********************/
    cout << endl << "**Standard LP**" << endl;
    cout << "max -z = ";
    cout << c[1] << "*[x_1]";
    for(i = 2; i <= var_num; i++)
    {
        cout << " + " << c[i] << "*[x_" << i << "]";
    }
    cout << endl;
    for(i = 1; i <= m; i++)
    {
        cout << a[i][1] << "*[x_1]";
        for(j = 2; j <= var_num; j++)
        {
            cout << " + " << a[i][j] << "*[x_" << j << "]";
        }
        cout << " = " << b[i] << endl;
    }
}