/*
This program is implemented for numeric genotype data that is coded as 0, 1, and 2.
*/


#include <iostream>
#include <fstream>
#include <string>

#include <vector>
#include <map>
#include <set>

#include <math.h>
#include <algorithm>

#include <omp.h>

#include "gurobi_c++.h"
using namespace std;

typedef vector<int> intv;
typedef vector<string> strv;
typedef vector<double> douv;
typedef set<int> ints;
typedef map<int, ints> intsm;
typedef map<int, int> intm;

const double TOL = 1e-5;
const int MAX_ITR = 1000;
const double LEAST_THRE = 0.0;
const bool INMEM = true;
const string INPUT = "f1";

const int TIMELIMIT = 3600;
const string THRE[] = {"0.7","0.5","0.3","0.2"};	// has to be sorted from large to small!!!
const string K[] = {"1","2","3"};


struct wedge
{
	int b;
	double w;
};
typedef vector<wedge> wedgev;
typedef map<int, wedgev> wedgevm;

struct domsol
{
	int totaln;
	int isopt;
	int pron;	// non-isolated nodes
	double soln;	// solution on unioslated graph only

	int ison;
	intv allsol;	// allsol consists of solution unisolated graph and isolated nodes. 
					// 1: in the dom, 2: isolated, 3: initialization, 0: otherwise 
};

// split function
void splitstr(strv& splits, const string& str, const string& deli = " ")
{
	size_t start, end;
	start = str.find_first_not_of(deli);
	splits.clear();

	while (start != string::npos)
	{
		end = str.find_first_of(deli, start + 1);
		if (end != string::npos)
		{
			splits.push_back(str.substr(start, end - start));
			start = str.find_first_not_of(deli, end + 1);
		}
		else
		{
			splits.push_back(str.substr(start));
			break;
		}
	}
}
void get_input(strv& fnames)
{
	string str, ful_file;
	ful_file = "input/" + INPUT + ".txt";

	ifstream fs(ful_file.c_str());

	if (!fs)
	{
		std::cout << "Cannot open filelist!" << endl;
		return;
	}

	while (!fs.eof())
	{
		getline(fs, str);
		if (str.size()>0)
			fnames.push_back(str);
	}

	fs.close();
}

//-----------------------------------------------------------------------------------------//
// Function to calculate LD (r^2) 
/*
* This part borrowed the idea of linkage disequilibrium calculation method from
* Rogers, Alan R.  and Huff, Chad. 2009. Linkage Disequilibrium in Loci with Unknown Phase. Genetics. 182(3):839-844.

* For any two biallelic markers, there are 4 types of gametes, which are AB, Ab, aB, and ab.
* genotypes AA (BB), Aa (Bb), and aa (bb) are numbered as 0, 1, and 2

* In order to calculate LD with missing data, the missing value is coded as 5.

* Input:
*
* hh is a vector of size 4, where the estimated haplotype by EM will be put
*
* h is a vector of 4 haplotype frequencies.
*
* x is a 3X3 matrix of genotype counts.  The phenotypes are coded
* as explained above.  Thus, x[1][2] is the number of copies of
* the phenotype Aa/BB (1/0).
*
* n is the sample size and should equal the sum of x.

* g is a 4X4 matrix of gamete frequencies.
* g[0][3] is the frequency of the genotype that combines gamete 0 (AB) with gamete 3 (ab).

* p is a 3X3 matrix of genotype frequencies, recoded as described for the input matrix x.

*/
void esem_step(douv &hh, douv&h, vector<douv>&x)
{
	vector<douv> g(4, douv(4, 0)), p(3, douv(3, 0));


	double n;
	for (size_t i = 0; i<4; i++)
	{
		g[i][i] = h[i] * h[i];
		for (size_t j = 0; j<i; j++)
		{
			g[i][j] = 2 * h[i] * h[j];
		}
	}
	p[0][0] = g[0][0];
	p[0][1] = g[1][0];
	p[0][2] = g[1][1];

	p[1][0] = g[2][0];
	p[1][1] = g[3][0] + g[2][1];
	p[1][2] = g[3][1];

	p[2][0] = g[2][2];
	p[2][1] = g[3][2];
	p[2][2] = g[3][3];

	hh[0] = 2 * x[0][0] + x[0][1] + x[1][0] + x[1][1] * g[3][0] / p[1][1];
	hh[1] = x[0][1] + 2 * x[0][2] + x[1][1] * g[2][1] / p[1][1] + x[1][2];
	hh[2] = x[1][0] + x[1][1] * g[2][1] / p[1][1] + 2 * x[2][0] + x[2][1];
	hh[3] = x[1][1] * g[3][0] / p[1][1] + x[1][2] + x[2][1] + 2 * x[2][2];

	/* convert gamete counts counts to relative frequencies*/
	n = hh[0] + hh[1] + hh[2] + hh[3];
	hh[0] /= n;
	hh[1] /= n;
	hh[2] /= n;
	hh[3] /= n;
}
// y is a vector of the first marker's genotype values, coded as 0, 1, and 2
//
// z is a vector of the first marker's genotype values

void count_genotypes(vector<douv> &x, intv& y, intv&z)
{
	for (size_t i = 0; i<y.size(); i++)
	{
		if (y[i] != 5 && z[i] != 5) //
		{
			x[y[i]][z[i]] += 1;
		}
	}
}
//calculate square of correlation coefficient
double esem_r(intv&y, intv &z, douv& h)
{
	double pA, pB, qA, qB, D, denom, r;

	pA = h[0] + h[1];
	pB = h[0] + h[2];
	qA = 1.0 - pA;
	qB = 1.0 - pB;
	D = h[0] * h[3] - h[1] * h[2];
	denom = pA*qA*pB*qB;
	if (denom != 0)
		r = D*D / denom;
	else
		r = 0.0;
	return r;
}
//
double esem(intv&y, intv&z, int &itern)
{
	vector<douv> x(3, douv(3, 0));
	double dh;
	douv hh(4, 0);
	douv h(4, 0.25); // initialize the haplotype frequencies as 0.25 for each


	count_genotypes(x, y, z);

	int i;
	for (i = 0; i < MAX_ITR; i++)
	{
		esem_step(hh, h, x); // return h estimated by EM
		dh = 0;
		for (size_t j = 0; j < 4; j++)
		{
			dh += fabs(h[j] - hh[j]);
			h[j] = hh[j];
		}
		if (dh <= TOL) break;
	}

	itern = i;
	return esem_r(y, z, h);
}

//-----------------------------------------------------------------------------------------//

void read(const string& file_str, vector<intv>&markers, strv&ID)
{
	ifstream fs(file_str.c_str());

	if (!fs)
	{
		std::cout << "Cannot open!" << endl;
		return;
	}

	string str;
	strv splits;
	intv individual;
	int n = 0;
	while (!fs.eof())
	{
		getline(fs, str);

		if (str.size() > 0)
		{
			splitstr(splits, str, ",");
			ID.push_back(splits[0] );

			for (size_t i = 1; i < splits.size(); i++) //starting point for data loading
			{
				individual.push_back(atoi(splits[i].c_str()));
			}

			markers.push_back(individual);
			n++;

			splits.clear();
			individual.clear();
		}
	}

}



//---------------------------------------------------------//
// Functions to check adjacency and get weight
// no use right now
bool comp(wedge& e1, int target)
{
	return e1.b < target;
}
double get_weight(wedgevm& adj, int a, int b)
{
	if (adj.find(a) == adj.end())
		return -1;

	wedgev::iterator i;
	i = lower_bound(adj[a].begin(), adj[a].end(), b, comp);

	//cout << i->b << "\t" << i->w << endl;

	if (i == adj[a].end() || i->b != b)
		return -1;
	else
		return i->w;
}
//---------------------------------------------------------//


//-----------------------------------------------------------------------------//
// Build and output the UNDERLINED association adjacency lists 
void asso_adj(vector<intv>& markers, const string& file_str)
{
	//When the matrix is too large, write to disk
	int mxitern = -1;
	double total_itern = 0, caln = 0;
	ofstream fout(file_str.c_str());

	int i, j, totaln = (int)markers.size();
	for (i = 0; i < totaln; i++)
	{
		intv iterns(totaln, -1);
		douv asso_array(totaln, -1);
		fout << i << ": ";

		//------------------------------------------------------//
#pragma omp parallel for
		for (j = i + 1; j < totaln; j++)
		{
			int itern;
			asso_array[j] = esem(markers[i], markers[j], itern); // call LD calculation function
			iterns[j] = itern;
		}
		//------------------------------------------------------//

		for (j = i + 1; j < totaln; j++)
		{
			if (asso_array[j] > LEAST_THRE)
				fout << j << "," << asso_array[j] << " ";

			if (iterns[j] > mxitern)
				mxitern = iterns[j];

			total_itern += iterns[j];
			caln++;
		}
		fout << endl;
	}

	std::cout << endl << "Max interation: " << mxitern << endl;
	std::cout << "Average interation: " << total_itern << "/" << caln
		<< " = " << total_itern / caln << endl << endl;

	fout.close();
}
int asso_adj(vector<intv>& markers, wedgevm& adj)
{
	// Load to meory if fitting 
	// return total count of nodes
	wedge e;  int mxitern = -1;
	double total_itern = 0, caln = 0;

	int i, j, totaln = (int)markers.size();
	for (i = 0; i < totaln; i++)
	{
		intv iterns(totaln, -1);
		douv asso_array(totaln, -1);

		//------------------------------------------------------//
#pragma omp parallel for
		for (j = i + 1; j < totaln; j++)
		{
			int itern;
			asso_array[j] = esem(markers[i], markers[j], itern);
			iterns[j] = itern;
		}
		//------------------------------------------------------//

		for (j = i + 1; j < totaln; j++)
		{
			if (asso_array[j] > LEAST_THRE)
			{
				e.w = asso_array[j];

				if (adj.find(i) == adj.end())
				{
					adj[i] = wedgev();
				}
				e.b = j;
				adj[i].push_back(e);

				if (adj.find(j) == adj.end())
				{
					adj[j] = wedgev();
				}
				e.b = i;
				adj[j].push_back(e);
			}

			if (iterns[j] > mxitern)
				mxitern = iterns[j];

			total_itern += iterns[j];
			caln++;
		}
	}

	std::cout << endl << "Max interation: " << mxitern << endl;
	std::cout << "Average interation: " << total_itern << "/" << caln
		<< " = " << total_itern / caln << endl << endl;

	return totaln;
}
void print_adj_full(wedgevm& adj, int totaln, const string& file_str)
{
	ofstream fout(file_str.c_str());

	int ison = 0;
	wedgev::iterator j;

	for (int i = 0; i < totaln; i++)
	{
		fout << i << ": ";

		if (adj.find(i) != adj.end())
		{
			for (j = adj[i].begin(); j != adj[i].end(); j++)
			{
				fout << j->b << "," << j->w << " ";
			}
		}
		else
		{
			ison++;
		}
		fout << endl;
	}

	fout << "Isolated node count: " << ison << endl;
	fout.close();
}
void print_adj_upper(wedgevm& adj, int totaln, const string& file_str)
{
	ofstream fout(file_str.c_str());

	int ison = 0;
	wedgev::iterator j;

	for (int i = 0; i < totaln; i++)
	{
		fout << i << ": ";

		if (adj.find(i) != adj.end())
		{
			for (j = adj[i].begin(); j != adj[i].end(); j++)
			{
				if (j->b > i)
				{
					fout << j->b << "," << j->w << " ";
				}
			}
		}
		else
		{
			ison++;
		}
		fout << endl;
	}

	fout << "Isolated node count: " << ison << endl;
	fout.close();
}
void build_adj_mem(const string& fname, wedgevm& adj, domsol& ds)
{
	strv ID;
	vector<intv> markers;

	string ipath = "input/infile/" + fname + ".csv";
	string opath = "output/" + INPUT + "/" + fname + "_adj.txt";

	std::cout << "Reading " + fname + "..." << endl;
	read(ipath, markers, ID);
	std::cout << "Reading is done!" << endl << endl
		<< "Calculating association matrix (put in memory)" << endl;

	double t = omp_get_wtime();

	ds.totaln = asso_adj(markers, adj);

	//print_adj_full(adj, n, opath);
	print_adj_upper(adj, ds.totaln, opath);

	t = omp_get_wtime() - t;

	std::cout << "Calculation is done!" << endl
		<< "The running time is " << t << " seconds" << endl << endl;
}
void build_adj_disk(const string& fname)
{
	strv ID;
	vector<intv> markers;

	string ipath = "input/infile/" + fname + ".csv";
	string opath = "output/" + INPUT + "/" + fname + "_adj.txt";

	std::cout << "Reading " + fname + "..." << endl;
	read(ipath, markers, ID);
	std::cout << "Reading is done!" << endl << endl
		<< "Calculating association matrix (write to disk)" << endl;

	double t = omp_get_wtime();
	asso_adj(markers, opath);
	t = omp_get_wtime() - t;

	std::cout << "Calculation is done!" << endl;
	std::cout << "The running time is " << t << " seconds" << endl << endl;
}
//----------------------------------------------------------------------------//


//------------------------------------------------------------------------------------------------------//
void get_graph(wedgevm &adj, intv& prior_sol, double thre, intsm &g, intm& nid2xid, intv& xid2nid)
{// Get prior dom and build graph for non-isolated nodes only
	std::cout << "Getting the sub-graph with threshold " << thre << endl;
	wedgevm::iterator i;
	wedgev::iterator j;
	intsm::iterator k;

	for (i = adj.begin(); i != adj.end(); i++)
	{
		if (prior_sol[i->first]>0)	// >0: get prior dom obtained
		{
			for (j = i->second.begin(); j != i->second.end(); j++)
			{
				if (prior_sol[j->b]>0 && j->w > thre) //lazy without eps for thre !!!!!
				{	// this 'if' block allows not creating isolated node
					if (g.find(i->first) == g.end())
					{
						g[i->first] = ints();
					}
					g[i->first].insert(j->b);

					if (g.find(j->b) == g.end())
					{
						g[j->b] = ints();
					}
					g[j->b].insert(i->first);
				}
			}
		}
	}

	int subg_n = 0;
	for (k = g.begin(); k != g.end(); k++)
	{
		nid2xid[k->first] = subg_n;
		xid2nid.push_back(k->first);
		subg_n++;
	}

	std::cout << "Done with getting the sub-graph" << endl << endl;
}
void print_g_upper(intsm& g, int totaln, const string& file_str)
{
	ofstream fout(file_str.c_str());

	ints::iterator j;
	for (int i = 0; i < totaln; i++)
	{
		fout << i << ": ";

		if (g.find(i) != g.end())
		{
			for (j = g[i].begin(); j != g[i].end(); j++)
			{
				if (*j > i)
				{
					fout << *j << " ";
				}
			}
		}
		fout << endl;
	}

	fout.close();
}
void print_g_full(intsm& g, int totaln, const string& file_str)
{
	ofstream fout(file_str.c_str());

	ints::iterator j;
	for (int i = 0; i < totaln; i++)
	{
		fout << i << ": ";

		if (g.find(i) != g.end())
		{
			for (j = g[i].begin(); j != g[i].end(); j++)
			{
				fout << *j << " ";
			}
		}
		fout << endl;
	}

	fout.close();
}
//------------------------------------------------------------------------------------------------------//

//-------------------------------------------------------------------------------------------------------------------------------//
// solve using k-domination model
GRBVar* grb_vars(int n, GRBModel &model, double obj_coeff, double lb = 0, double ub = GRB_INFINITY, char type = GRB_CONTINUOUS)
{
	// delete after use
	GRBVar* vars = new GRBVar[n];

	for (int i = 0; i<n; i++)
		vars[i] = model.addVar(lb, ub, obj_coeff, type);

	return vars;
}
void dom(intsm& g, intm& nid2xid, intv& xid2nid, int k, domsol &ds)
{
	std::cout << endl << "IP solving with k = " << k << "..." << endl;

	int i, nid, xid;
	int subg_n = (int)g.size();
	int totaln = ds.totaln;

	ints::iterator j;

	GRBEnv env = GRBEnv();
	env.set(GRB_DoubleParam_TimeLimit, TIMELIMIT);
	//env.set(GRB_IntParam_Threads, 4);

	GRBModel model = GRBModel(env);
	GRBVar *x = grb_vars(subg_n, model, 1, 0, 1, GRB_BINARY);
	model.update();

	for (i = 0; i < subg_n; i++)
	{
		nid = xid2nid[i];
		if (g[nid].size() < k)
		{
			model.addConstr(x[i] == 1);
		}
		else
		{
			GRBLinExpr lhs = k*x[i];
			for (j = g[nid].begin(); j != g[nid].end(); j++)
			{
				xid = nid2xid[*j];
				lhs += x[xid];
			}
			model.addConstr(lhs >= k);
		}
	}

	model.optimize();

	ds.pron = subg_n;

	intv prior_allsol = ds.allsol;
	ds.allsol.assign(totaln, 0);

	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
		ds.isopt = 1;
	else
		ds.isopt = 0;

	ds.soln = model.get(GRB_DoubleAttr_ObjVal);

	for (i = 0; i < subg_n; i++)
	{
		if (x[i].get(GRB_DoubleAttr_X) > 0.5)
		{
			nid = xid2nid[i];
			ds.allsol[nid] = 1;
		}
	}

	ds.ison = 0;
	for (i = 0; i < totaln; i++)
	{
		if (prior_allsol[i]>0 && g.find(i) == g.end())
		{
			ds.allsol[i] = 2;
			ds.ison++;
		}
	}

	std::cout << endl << "IP is done" << endl << endl << endl;
}
void print_domsol(const string& fname, domsol& ds)
{
	int i;
	string opath = "output/" + INPUT + "/" + fname + "_sol.txt";

	ofstream fout(opath.c_str());

	fout << "Original size: " << ds.totaln << endl
		<< "Problem size: " << ds.pron << endl
		<< "Optimization status: " << ds.isopt << endl
		<< "Solution size: " << ds.soln << endl << endl
		<< "Solution node IDs: " << endl;

	for (i = 0; i<ds.totaln; i++)
	{
		if (ds.allsol[i] == 1)
			fout << i << ", ";
	}
	fout << endl << endl << "Isolated node IDs (" << ds.ison << "): " << endl;

	for (i = 0; i<ds.totaln; i++)
	{
		if (ds.allsol[i] == 2)
			fout << i << ", ";
	}
	fout << endl << endl;

	fout.close();
}
//-------------------------------------------------------------------------------------------------------------------------------//

void solve_one(const string& fname)
{
	wedgevm adj; domsol ds;
	build_adj_mem(fname, adj, ds);

	int i, j;

	int thren = sizeof(THRE) / sizeof(*THRE);
	int kn = sizeof(K) / sizeof(*K);
	for (i = 0; i<kn; i++)
	{
		ds.allsol.assign(ds.totaln, 3);	// make program iteratively solve

		int k = atoi(K[i].c_str());

		for (j = 0; j<thren; j++)
		{
			double thre = atof(THRE[j].c_str());

			intsm g; intm nid2xid; intv xid2nid;
			get_graph(adj, ds.allsol, thre, g, nid2xid, xid2nid);

			dom(g, nid2xid, xid2nid, k, ds);

			string final_fname = fname + "_" + K[i] + "_" + THRE[j];

			print_domsol(final_fname, ds);
		}
	}
}


int main()
{
	strv fnames;
	get_input(fnames);

	strv::iterator i;
	for (i = fnames.begin(); i != fnames.end(); i++)
	{
		if (INMEM)
		{
			solve_one(*i);
		}
		else
		{
			build_adj_disk(*i);
		}
	}
	system("pause");

	return 0;
}
