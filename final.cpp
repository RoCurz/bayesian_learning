#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <bits/stdc++.h>

// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory
using namespace std;

// Our graph consists of a list of nodes where each node is represented as follows:
class Graph_Node
{

private:
    string Node_Name;       // Variable name
    vector<int> Children;   // Children of a particular node - these are index of nodes in graph.
    vector<string> Parents; // Parents of a particular node- note these are names of parents
    int nvalues;            // Number of categories a variable represented by this node can take
    vector<string> values;  // Categories of possible values
    vector<double> CPT;     // conditional probability table as a 1-d array . Look for BIF format to understand its meaning

public:
    // Constructor- a node is initialised with its name and its categories
    Graph_Node(string name, int n, vector<string> val)
    {
        Node_Name = name;

        nvalues = n;
        values = val;
    }
    string get_name()
    {
        return Node_Name;
    }
    vector<int> get_children()
    {
        return Children;
    }
    vector<string> get_Parents()
    {
        return Parents;
    }
    vector<double> get_CPT()
    {
        return CPT;
    }
    int get_nvalues()
    {
        return nvalues;
    }
    vector<string> get_values()
    {
        return values;
    }
    void set_CPT(vector<double> new_CPT)
    {
        CPT.clear();
        CPT = new_CPT;
    }
    void set_Parents(vector<string> Parent_Nodes)
    {
        Parents.clear();
        Parents = Parent_Nodes;
    }
    // add another node in a graph as a child of this node
    int add_child(int new_child_index)
    {
        for (int i = 0; i < Children.size(); i++)
        {
            if (Children[i] == new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }

    // -------------------------------------------------------- //
    int get_value_ind(string value)
    {
        for (int i = 0; i < values.size(); i++)
        {
            if (values[i] == value)
                return i;
        }
        return -1;
    }
};

// The whole network represted as a list of nodes
class network
{

    list<Graph_Node> Pres_Graph;

public:
    int addNode(Graph_Node node)
    {
        Pres_Graph.push_back(node);
        return 0;
    }

    int netSize()
    {
        return Pres_Graph.size();
    }
    // get the index of node with a given name
    int get_index(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        int count = 0;
        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
        {
            if (listIt->get_name().compare(val_name) == 0)
                return count;
            count++;
        }
        return -1;
    }
    // get the node at nth index
    list<Graph_Node>::iterator get_nth_node(int n)
    {
        list<Graph_Node>::iterator listIt;
        int count = 0;
        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
        {
            if (count == n)
                return listIt;
            count++;
        }
        return listIt;
    }
    // get the iterator of a node with a given name
    list<Graph_Node>::iterator search_node(string val_name)
    {
        list<Graph_Node>::iterator listIt;
        for (listIt = Pres_Graph.begin(); listIt != Pres_Graph.end(); listIt++)
        {
            if (listIt->get_name().compare(val_name) == 0)
                return listIt;
        }

        cout << "node not found\n";
        return listIt;
    }
};

network read_network(string filename)
{
    network Alarm;
    string line;
    int find = 0;
    ifstream myfile(filename);
    string temp;
    string name;
    vector<string> values;

    if (myfile.is_open())
    {
        while (!myfile.eof())
        {
            stringstream ss;
            getline(myfile, line);

            ss.str(line);
            ss >> temp;

            if (temp.compare("variable") == 0)
            {

                ss >> name;
                getline(myfile, line);

                stringstream ss2;
                ss2.str(line);
                for (int i = 0; i < 4; i++)
                {

                    ss2 >> temp;
                }
                values.clear();
                while (temp.compare("};") != 0)
                {
                    values.push_back(temp);

                    ss2 >> temp;
                }
                Graph_Node new_node(name, values.size(), values);
                int pos = Alarm.addNode(new_node);
            }
            else if (temp.compare("probability") == 0)
            {

                ss >> temp;
                ss >> temp;

                list<Graph_Node>::iterator listIt;
                list<Graph_Node>::iterator listIt1;
                listIt = Alarm.search_node(temp);
                int index = Alarm.get_index(temp);
                ss >> temp;
                values.clear();
                while (temp.compare(")") != 0)
                {
                    listIt1 = Alarm.search_node(temp);
                    listIt1->add_child(index);
                    values.push_back(temp);

                    ss >> temp;
                }
                listIt->set_Parents(values);
                getline(myfile, line);
                stringstream ss2;

                ss2.str(line);
                ss2 >> temp;

                ss2 >> temp;

                vector<double> curr_CPT;
                string::size_type sz;
                while (temp.compare(";") != 0)
                {

                    curr_CPT.push_back(atof(temp.c_str()));

                    ss2 >> temp;
                }

                listIt->set_CPT(curr_CPT);
            }
            else
            {
            }
        }

        if (find == 1)
            myfile.close();
    }

    return Alarm;
}
vector<vector<int>> graph_data;
vector<int> missing;
vector<vector<int>> datapoints;
vector<double> probabilities;
// ofstream outf("output.txt");
// read data from records.dat and store them in above variables
void read_data(string datafile, network &Alarm)
{
    ifstream myfile(datafile);
    string val, line;
    while (!myfile.eof())
    {
        stringstream ss;
        getline(myfile, line);
        ss.str(line);
        int ind = 0;
        bool flag = false;
        vector<int> temp;
        while (ss >> val)
        {
            if (val == "\"?\"")
            {
                missing.push_back(ind);
                flag = true;
            }
            temp.push_back(Alarm.get_nth_node(ind)->get_value_ind(val));
            ind++;
        }
        if (!flag)
        {
            missing.push_back(-1);
        }
        graph_data.push_back(temp);
        if (!flag)
        {
            datapoints.push_back(temp);
        }
        else
        {
            int missing_ind = missing[missing.size() - 1];
            int N = Alarm.get_nth_node(missing_ind)->get_nvalues();
            for (int i = 0; i < N; i++)
            {
                // temp[missing_ind] = i;
                datapoints.push_back(temp);
                datapoints[datapoints.size() - 1][missing_ind] = i;
            }
        }
    }
    myfile.close();
}

void initialise_cpt(network &Alarm)
{
    for (int i = 0; i < Alarm.netSize(); i++)
    {
        int siz = Alarm.get_nth_node(i)->get_CPT().size();
        vector<double> new_CPT(siz, 0.0);
        for (int j = 0; j < siz; j++)
        {
            // number of values for i-th node
            int N = Alarm.get_nth_node(i)->get_nvalues();
            new_CPT[j] = 1.0 / N;
        }
        Alarm.get_nth_node(i)->set_CPT(new_CPT);
    }
}

int find_cpt_ind(vector<int> &val, vector<int> &sizes)
{
    if (val.size() == 0)
    {
        return 0;
    }
    int idx = 0, b = 1;
    int M = sizes.size();
    for (int i = M - 1; i >= 0; i--)
    {
        idx = idx + b * val[i];
        b = b * sizes[i];
    }
    return idx;
}

// calculate the expectation of missing values
void expectation_maximization(network &Alarm)
{
    //expectation 
    probabilities.clear();
    for (int i = 0; i < graph_data.size(); i++)
    {
        int missing_ind = missing[i];
        if (missing_ind == -1)
        {
            probabilities[i] = 1;
        }
        else
        {
            int N = Alarm.get_nth_node(missing_ind)->get_nvalues();
            double numerator, denominator = 0.0;
            vector<double> all_numerators;
            for (int t = 0; t < N; t++)
            {
                numerator = 1.0;
                vector<int> temp(graph_data[i].begin(), graph_data[i].end());
                temp[missing_ind] = t;
                vector<int> val, sizes;
                vector<int> children = Alarm.get_nth_node(missing_ind)->get_children();
                for (int j = 0; j < children.size(); j++)
                {
                    val.clear(), sizes.clear();
                    int child = children[j];
                    int curr_size = Alarm.get_nth_node(child)->get_nvalues();
                    val.push_back(graph_data[i][child]);
                    sizes.push_back(curr_size);
                    vector<string> parents = Alarm.get_nth_node(child)->get_Parents();
                    for (int k = 0; k < parents.size(); k++)
                    {
                        string parent = parents[k];
                        int parent_ind = Alarm.get_index(parent);
                        val.push_back(temp[parent_ind]);
                        sizes.push_back(Alarm.get_nth_node(parent_ind)->get_nvalues());
                    }
                    numerator = numerator * Alarm.get_nth_node(child)->get_CPT()[find_cpt_ind(val, sizes)];
                    // cout<<Alarm.get_nth_node(child)->get_CPT().size()<<" "<<find_cpt_ind(val, sizes)<<"\n";
                }
                denominator += numerator;
                val.clear(), sizes.clear();
                val.push_back(t);
                sizes.push_back(N);
                vector<string> parents = Alarm.get_nth_node(missing_ind)->get_Parents();
                for (int k = 0; k < parents.size(); k++)
                {
                    string parent = parents[k];
                    // cout<<parent<<"\n";
                    int parent_ind = Alarm.get_index(parent);
                    val.push_back(temp[parent_ind]);
                    sizes.push_back(Alarm.get_nth_node(parent_ind)->get_nvalues());
                }
                numerator = numerator * Alarm.get_nth_node(missing_ind)->get_CPT()[find_cpt_ind(val, sizes)];
                all_numerators.push_back(numerator);
                // all_numerators.push_back(numerator);
            }
            for (int j = 0; j < all_numerators.size(); j++)
            {
                probabilities.push_back(all_numerators[j] / denominator);
                // if(isnan(denominator)){
                //     cout<<denominator<<"\n";
                // }
            }
        }
    }

    // maximization
    // double max_diff = 0.0;
    for (int i = 0; i < Alarm.netSize(); i++)
    {
        vector<int> val, val_ind_map, sizes;
        val_ind_map.push_back(i);
        sizes.push_back(Alarm.get_nth_node(i)->get_nvalues());
        vector<string> parents = Alarm.get_nth_node(i)->get_Parents();
        for (int j = 0; j < parents.size(); j++)
        {
            string parent = parents[j];
            int parent_ind = Alarm.get_index(parent);
            val_ind_map.push_back(parent_ind);
            sizes.push_back(Alarm.get_nth_node(parent_ind)->get_nvalues());
        }
        vector<double> curr_cpt = Alarm.get_nth_node(i)->get_CPT();
        int MOD = curr_cpt.size() / (Alarm.get_nth_node(i)->get_nvalues());
        vector<double> denominators(MOD, 0.0), numerators(curr_cpt.size(), 0.0);
        // start editing the cpt
        for (int j = 0; j < datapoints.size(); j++)
        {
            val.clear();
            int index;
            for (int k = 0; k < val_ind_map.size(); k++)
            {
                val.push_back(datapoints[j][val_ind_map[k]]);
            }
            index = find_cpt_ind(val, sizes);
            // int curr_ind = find_cpt_ind(datapoints[j], val_ind_map);
            denominators[index % MOD] += probabilities[j];
            numerators[index] += probabilities[j];
        }
        double temp;
        for (int j = 0; j < curr_cpt.size(); j++)
        {
            // if(denominators[j%MOD]+0.001*Alarm.get_nth_node(i)->get_nvalues() < 0.001){
            //     cout<<denominators[j%MOD]+0.001*Alarm.get_nth_node(i)->get_nvalues()<<"\n";
            // }
            temp = (numerators[j] + 0.001) / (denominators[j % MOD] + 0.001 * Alarm.get_nth_node(i)->get_nvalues());
            // max_diff = max(max_diff, abs(temp - curr_cpt[j]));
            curr_cpt[j] = temp;
            // cout<<temp<<"\n";
            if (curr_cpt[j] < 0.0001)
            {
                curr_cpt[j] = 0.0001; // Addition for smoothing and avoiding a zero
            }

            // need to update previous cpt with curr updated cpt
            Alarm.get_nth_node(i)->set_CPT(curr_cpt);
        }
    }
    // if (max_diff < 0.00001)
    //     return true;
    // return false;
}

// write cpt to a file
void write_network(network &Alarm, string filename)
{
    ifstream ref_file(filename);
    ofstream outfile;
    outfile.open("solved_alarm.bif");
    outfile << setprecision(4) << fixed;
    string line;
    string temp;
    string name;
    int counter = 0;
    if (ref_file.is_open())
    {
        while (!ref_file.eof())
        {
            stringstream ss;
            getline(ref_file, line);
            ss.str(line);
            ss >> temp;
            if (temp.compare("probability") == 0)
            {
                outfile << line << "\n";
                getline(ref_file, line);
                outfile << "\ttable ";
                // cpt of index counter
                vector<double> CPT = Alarm.get_nth_node(counter)->get_CPT();
                for (int i = 0; i < CPT.size(); i++)
                    outfile << CPT[i] << " ";
                outfile << " ;\n";
                counter++;
            }
            else
            {
                if (line.size() != 0)
                    outfile << line << "\n";
                else
                    outfile << line;
            }
        }
        ref_file.close();
        outfile.close();
    }
}

int main(int argc, char **argv)
{
    // start time in chrono
    auto start = chrono::high_resolution_clock::now();

    network Alarm;
    Alarm = read_network(argv[1]);
    initialise_cpt(Alarm);
    read_data(argv[2], Alarm);
    while (true)
    {
        expectation_maximization(Alarm);
        auto curr_time = chrono::high_resolution_clock::now();
        if ((chrono::duration_cast<chrono::milliseconds>(curr_time-start).count())>=115*(1000))
            break;
    }
    write_network(Alarm, argv[1]);

    cout << "Perfect! Hurrah! \n";
    // end time
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "Time taken: " << static_cast<double>(duration.count())/1000.0 << " seconds\n";
}
