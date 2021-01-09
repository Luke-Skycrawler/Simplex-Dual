#include "PreProcess_dual.h"

/********************************************
 * global variables definitions
*********************************************/
//Input variables -- static
static int n, m;
static double in_c[MAXN] = {0};
static double in_a[MAXM][MAXN] = {0};
static double in_b[MAXM] = {0};
static int in_d[MAXM] = {0}, in_e[MAXN] = {0};

//Variables for **export** -- global
int var_num; //total number of variables after standardize()
int new_var_num = 0; //number of newly added variables
int constraint_num = 0;
int new_constraint_num = 0;
vector<double> *a;
vector<double> b;
vector<double> c;
struct outputInfo output_info[MAXN];


/********************************************
 * function definitions
*********************************************/
struct preProcessDualResult* preProcessDual(string file_name)
{
    readInput(file_name);
    standardize();
    return new preProcessDualResult(var_num, constraint_num, a, b, c, output_info, n);
}

void readInput(string file_name)
{
    int i, j;
    ifstream in_file;
    in_file.open(file_name, ios::in);
    in_file >> n >> m;
    assert(n < MAXN);
    assert(m < MAXM);
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
    // cout << endl << "Load problem from " << file_name << " successfully!\n**Original LP**" << endl;
    // cout << "min z = ";
    // cout << in_c[1] << "*[x_1]";
    // for(i = 2; i <= n; i++)
    // {
    //     cout << " + " << in_c[i] << "*[x_" << i << "]";
    // }
    // cout << endl;
    // for(i = 1; i <= m; i++)
    // {
    //     cout << in_a[i][1] << "*[x_1]";
    //     for(j = 2; j <= n; j++)
    //     {
    //         cout << " + " << in_a[i][j] << "*[x_" << j << "]";
    //     }
    //     if(in_d[i] < 0) 
    //         cout << " <= " << in_b[i] << endl;
    //     else if(in_d[i] == 0) 
    //         cout << " = " << in_b[i] << endl;
    //     else
    //         cout << " >= " << in_b[i] << endl;
    // }
    // for(j = 1; j <= n; j++)
    // {
    //     if(in_e[j] == -1)
    //         cout << "[x_" << j << "]" << " <= 0" << endl;
    //     else if(in_e[j] == 1)
    //         cout << "[x_" << j << "]" << " >= 0" << endl;
    //     else
    //         cout << "[x_" << j << "]" << " <> 0" << endl;
    // }
}

void standardize()
{
    int i, j;

    //initialize global variables
    a = new vector<double>[MAXM + MAXN]; //malloc for matrix `a`
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
    //Difference from PRIMAL:
    //we prefer **less equal**!
    for(i = 1; i <= m; i++)
    {
        if(in_d[i] == 1) // >=; reverse both sides
        {
            for(j = 1; j <= n + new_var_num; j++)
            {
                a[i][j] *= -1;
            }
            b[i] *= -1;
        }
        else if(in_d[i] == 0) //split = into <= and >=
        {
            new_constraint_num++;
            for(j = 0; j <= n + new_var_num; j++)
            {
                a[m + new_constraint_num].push_back(-a[i][j]);
            }
            b.push_back(-b[i]);
        }
    }


    var_num = n + new_var_num; //update global record for variables' number
    constraint_num = m + new_constraint_num;

    /******************** debug ********************/
    // cout << endl << "**Standard LP**" << endl;
    // cout << "max -z = ";
    // cout << c[1] << "*[x_1]";
    // for(i = 2; i <= var_num; i++)
    // {
    //     cout << " + " << c[i] << "*[x_" << i << "]";
    // }
    // cout << endl;
    // for(i = 1; i <= constraint_num; i++)
    // {
    //     cout << a[i][1] << "*[x_1]";
    //     for(j = 2; j <= var_num; j++)
    //     {
    //         cout << " + " << a[i][j] << "*[x_" << j << "]";
    //     }
    //     cout << " <= " << b[i] << endl;
    // }
}