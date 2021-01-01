#include "PreProcess.h"

int main()
{
    struct preProcessResult* input; //a pointer as a result acceptor

    input = preProcess("./test.txt"); /** the one and only one call you should make **/

    cout << endl << "-------" << endl;
    cout << "`var_num`: " << input->var_num << endl; //total variable number
    cout << "`new_var_num`: " << input->new_var_num << endl; //newly added variable number
    cout << "number of elements in `a[1]`: " << input->a[1].size() - 1 << endl; //the column number of `a`
    cout << "number of elements in `b`: " << input->b.size() - 1 << endl; //also the row number of `a`
    cout << "number of elements in `c`: " << input->c.size() - 1 << endl;
    //output information explanation
    for(int i = 1; i <= input->var_num - input->new_var_num; i++)
    {
        cout << "`output_info[" << i << "]`: ";
        cout << "\ttype = " << input->output_info[i].type << endl;
        if(input->output_info[i].type == 1)
        {
            cout << "\t[x_" << i << "] = -[x_" << i << "]" << endl;
        }
        if(input->output_info[i].type == 2)
        {
            cout << "\t[x_" << i << "] = [x_" << input->output_info[i].first_aux << "] - [x_" << input->output_info[i].second_aux << "]" << endl;
        }
    }
    
    return 0;
}