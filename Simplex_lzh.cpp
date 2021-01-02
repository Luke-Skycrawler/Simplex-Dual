#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <cfloat>

#define MAX_ROW 500
#define MAX_COL 200
#define ZERO 1e-4
using namespace std;

class Simplex
{
private:
    int n, m;
    int *e;
    int *d;
    double *b;
    double *c;
    double **a;
    double *x;
    void prepare()
    {
        int i, j;
        for(i=0; i<m; i++)
        {
            if(b[i]<0)
            {
                b[i] = -b[i];
                d[i] = -d[i];
                for(j=0; j<n; j++)
                    a[i][j] = -a[i][j];
            }
        }
        vector<double> m1[MAX_ROW];
        vector<double> m2[MAX_ROW];
        for(i=0; i<m; i++)
        {
            for(j=0; j<n; j++)
                m1[i].push_back(a[i][j]);
        }
        for(i=0; i<n; i++)
        {
            if(e[i]<0)
            {
                for(j=0; j<m; j++)
                    m1[j][i] = -m1[j][i];
                c[i] = -c[i];
            }
            else if(e[i]==0)
            {
                for(j=0; j<m; j++)
                    m1[j].push_back(-1);
            }
        }
        for(i=0; i<m; i++)
        {
            if(d[i]>0)
            {
                for(j=0; j<m; j++)
                {
                    if(j==i)
                        m1[j].push_back(-1);
                    else
                        m1[j].push_back(0);
                }
            }
            if(d[i]<0)
            {
                for(j=0; j<m; j++)
                {
                    if(j==i)
                        m1[j].push_back(1);
                    else
                        m1[j].push_back(0);
                }
            }
        }
        int newCol = m1[0].size();
        
        for(i=0; i<m; i++)
        {
            delete[] this->a[i];
            this->a[i] = new double[newCol];
            for(j=0; j<newCol; j++)
                a[i][j] = m1[i][j];
        }

        double *tmpc = new double[newCol];
        for(j=0; j<newCol; j++)
        {
            if(j<this->n)   
                tmpc[j] = this->c[j];
            else
                tmpc[j] = 0;
        }
        delete[] this->c;
        this->c = tmpc;

        x = new double[newCol];
        for(i=0; i<newCol; i++)
            x[i] = 0;
        
        this->n = newCol;
    }
public:
    Simplex(): n(0), m(0), z(0), e(NULL), d(NULL), b(NULL), c(NULL), a(NULL), x(NULL){};
    void readFile(const char *path)
    {
        int i, j;
        ifstream file(path);
        if(!file.is_open())
        {
            cout << "Open file error" << endl;
        }
        file >> this->n >> this->m;
        this->a = new double*[this->m];
        for(i=0; i<this->m; i++)
            this->a[i] = new double[this->n];
        this->b = new double[this->m];
	    this->c = new double[this->n];
	    this->d = new int[this->m];
	    this->e = new int[this->n];
        for(i=0; i<this->n; i++)
            file >> this->c[i];
        for(i=0; i<this->m; i++)
        {
            for(j=0; j<this->n; j++)
                file >> this->a[i][j];
            file >> this->b[i];
            file >> this->d[i];
        }
        for(i=0; i<this->n; i++)
            file >> this->e[i];
        prepare();
    }
    void print()
    {
        int i, j;
        cout << "m = " << m << endl;
        cout << "n = " << n << endl;
        cout << "c = ";
        for(i=0; i<n; i++)
            cout << c[i] << "\t";
        cout << endl;
        cout << "m = " << endl;
        for(i=0; i<m; i++)
        {
            for(j=0; j<n; j++)
                cout << a[i][j] << "\t";
            cout << b[i] << endl;
        }
    }
    void print_root()
    {
        if(x==NULL)
        {
            cout << "No root" << endl;
            return;
        }
        int i;
        for(i=0; i<this->n; i++)
            cout << x[i] << " ";
        cout << endl;
    }
    Simplex transform()
    {
        int i, j;
        Simplex M;
        M.n = this->m+this->n;
        M.m = this->m;
        M.b = new double[M.m];
        for(i=0; i<M.m; i++)
            M.b[i] = this->b[i];
        M.c = new double[M.n];
        for(i=0; i<M.n; i++)
        {
            if(i<this->n)
                M.c[i] = 0;
            else
                M.c[i] = 1;
        }
        M.a = new double*[M.m];
        for(i=0; i<M.m; i++)
        {
            M.a[i] = new double[M.n];
            for(j=0; j<M.n; j++)
            {
                if(j<this->n)
                    M.a[i][j] = this->a[i][j];
                else if(i==j-this->n)
                    M.a[i][j] = 1;
                else
                    M.a[i][j] = 0;
            }
        }
        M.x = new double[M.n];
        for(i=0; i<M.n; i++)
            M.x[i] = 0;
        for(i=this->n; i<M.n; i++)
        {
            for(j=0; j<M.n; j++)
                M.c[j] -= M.a[i-this->n][j];
            M.z -= b[i-this->n];
        }
        return M;
    }
    Simplex transport_back_ab()
    {
        int i, j;
        Simplex m;
        m.m = this->m;
        m.n = this->n-this->m;

        m.b = new double[m.m];
        for(i=0; i<m.m; i++)
            m.b[i] = this->b[i];
        m.c = new double[m.n];
        for(i=0; i<m.n; i++)
            m.c[i] = 0;
        m.x = new double[m.n];
        for(i=0; i<m.n; i++)
            m.x[i] = 0;
        m.a = new double*[m.m];
        for(i=0; i<m.m; i++)
        {
            m.a[i] = new double[m.n];
            for (j=0; j<m.n; j++)
            {
                m.a[i][j] = a[i][j];
            }
        }
        m.z = 0;
        return m;
    }
    void copy_c(Simplex orig)
    {
        int i;
        for(i = 0; i < orig.n; i++)
            this->c[i] = orig.c[i];
    }
    void make_normal(Simplex orig)
    {
        int i, j, k;
        int index;
        double factor;
        for(i=0; i<this->n; i++)
        {
            if(fabs(orig.x[i])>ZERO)
            {
                factor = c[i];
                for(j=0; j<this->m; j++)
                {
                    if(this->a[j][i] == 1)
                    {
                        for(k=0; k<this->n; k++)
                            this->c[k] -= this->a[j][k]*factor;
                        this->z += factor*this->b[j];
                        break;
                    }
                }
            }
        }
    }
    int do_counting()
    {
        int i, j, k;
        double **table;
        double **tmptable;
        table = new double*[m+1];
        tmptable = new double*[m+1];
        for(i=0; i<m+1; i++)
        {
            table[i] = new double[n];
            tmptable[i] = new double[n];
        }
        for(i=0; i<m; i++)
        {
            for(j=0; j<n; j++)
                table[i][j] = a[i][j];
        }
        for(i=0; i<m; i++)
            table[i][n] = b[i];
        for(j=0; j<n; j++)
            table[m][j] = c[j];
        table[m][n] = -z;

        int in_base, out_base;
        int argmin_outbase;
        double min_outbase = 0;
        double temp;
        double cur_min;
        bool is_base;
        while(1)
        {
            in_base = -1;
            cur_min = 0;
            for(i=0; i<n; i++)
            {
                if(table[m][i]<0)
                {
                    in_base = i;
                    break;
                }
            }
            if(in_base<0)
            {
                for(i=0; i<m; i++)
                {
                    for(j=0; j<n; j++)
                    {
                        if(table[i][j]==1)
                        {
                            is_base = true;
                            for(k=0; k<m; k++)
                            {
                                if(k!=i && table[k][j]!=0)
                                {
                                    is_base = false;
                                    break;
                                }
                            }
                            if(is_base == true)
                                x[j] = table[i][n];   
                        }
                    }
                }
                z = -table[m][n];

                for(i=0; i<m; i++)
                {
                    for(j=0; j<n; j++)
                        a[i][j] = table[i][j];
                }
                for(i=0; i<m; i++)
                    b[i] = table[i][n];
                return 1;
            }
            else
            {
                argmin_outbase = -1;
                min_outbase = DBL_MAX;
                for(i=0; i<m; i++)
                {
                    if(table[i][in_base]>0)
                    {
                        temp = table[i][n]/table[i][in_base];
                        if(temp<min_outbase)
                        {
                            min_outbase = temp;
                            argmin_outbase = i;
                        }
                    }
                }
                if(argmin_outbase>=0)
                    out_base = argmin_outbase;
                else
                    return 0;
                for(i=0; i<m+1; i++)
                {
                    for(j=0; j<n+1; j++)
                    {
                        if(i!=out_base)
                            tmptable[i][j] = table[i][j]-(table[i][in_base]/table[out_base][in_base])*table[out_base][j];
                        else
                            tmptable[i][j] = table[i][j]/table[out_base][in_base];
                    }
                }
                for(i=0; i<m+1; i++)
                {
                    for(j=0; j<n+1; j++)
                        table[i][j] = tmptable[i][j];
                }
            }
        }
    }
    double z;
};

int main()
{
	Simplex M;
    M.readFile("ÁÙ½ç.txt");
	Simplex M1;
	M1 = M.transform();
	M1.do_counting();
	if(M1.z>ZERO)
	{
		cout << -1 << endl;
	}
	else
	{
		Simplex M2;
		M2 = M1.transport_back_ab();
		M2.copy_c(M);
		M2.make_normal(M1);
		if(M2.do_counting()==0)
		{
			cout << 0 << endl;
			return 0;
		}
		cout << 1 << endl;
		cout << M2.z << endl;
		M2.print_root();
	}
}