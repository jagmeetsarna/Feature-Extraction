//#include <GLEW/include/GL/glew.h>
#include <GL/glut.h>
#include <GLFW/glfw3.h>
#include <ANN/ANN.h>
#include <filesystem>
#include <iostream>
#include <io.h>
#include <filesystem>
#include <conio.h>
#include <map>

#include <fstream>
#include <filesystem>
#include <tuple>
#include <stdio.h>
#include <string>
#include "..\header\Grid.h"
#include "..\header\FilterItem.h"
#include "..\header\OFFReader.h"
#include "..\header\Renderer.h"

using namespace std;
namespace fs = std::filesystem;
void loadFilter();

string fileName;

vector<pair<string, vector<float>>> feature_vectors;
vector<Grid> result_grids;

int drawing_style = 0;
const int FILTER_SIZE = 250;

int N = 1000000;                            // Number of computations per feature

ANNkd_tree* tree;                           // kd search tree for ANN


float D1min = FLT_MAX, D2min = FLT_MAX, D3min = FLT_MAX, D4min = FLT_MAX, A3min = FLT_MAX, SAmin = FLT_MAX, COmin = FLT_MAX, BBVmin = FLT_MAX, DIAmin = FLT_MAX, ECCmin = FLT_MAX;
float D1max = FLT_MIN, D2max = FLT_MIN, D3max = FLT_MIN, D4max = FLT_MIN, A3max = FLT_MIN, SAmax = FLT_MIN, COmax = FLT_MIN, BBVmax = FLT_MIN, DIAmax = FLT_MIN, ECCmax = FLT_MIN;
float SAavg = 0, COavg = 0, BBVavg = 0, DIAavg = 0, ECCavg = 0, SAsd = 0, COsd = 0, BBVsd = 0, DIAsd = 0, ECCsd = 0;
int bins = 14;

float SAval1, SAval2, COval1, COval2, BBVval1, BBVval2, DIAval1, DIAval2, ECCval1, ECCval2;

float w_f = 0.05, w_h = 0.19;

// Dictionary for the number of shapes in the whole data base
map<string, int> db_count;
int num_entries = 0;


Grid* grid = 0;
Renderer renderer;

FilterItem* fis;

int index = 0;

void mkdir(const char* dir)
{
    size_t size = strlen(dir) + 1;
    wchar_t* ndb = new wchar_t[size];
    size_t out;
    mbstowcs_s(&out, ndb, size, dir, size - 1);
    _wmkdir(ndb);
    delete[] ndb;
}

void draw()												                //Render the 3D mesh (GLUT callback function)
{
    renderer.draw(*grid);
}

bool sortbysec(const pair<string, float>& a, const pair<string, float>& b)
{
    return (a.second < b.second);
}

void scanFolder(string location)
{
    string folder;
    int i = 0;
    for (const auto& entry : fs::directory_iterator(location))
    {
        folder = entry.path().string();
        for (const auto& fl : fs::directory_iterator(folder + "/"))
        {
            string file = fl.path().string();
            string suffix = ".off";
            if (!(file.size() >= suffix.size() && 0 == file.compare(file.size() - suffix.size(), suffix.size(), suffix)))
                continue;
            FilterItem fi = scanFile(file);
            int a = folder.find("/");
            if (a <= folder.size())
                fi.cls = folder.substr(folder.find("/") + 1);
            else
                fi.cls = folder.substr(folder.find("\\") + 1);
            fis[i] = fi;
            i++;
        }
    }
    cout << "Scanning complete!" << endl;
}


void keyboard(unsigned char c, int, int)					//Callback for keyboard events:
{
    switch (c)
    {
    case ' ':											    // space:   Toggle through the various drawing styles of the mesh renderer
    {
        break;
    }
    }
}

float normalize(float value, float min, float max)
{
    return (value - min) / (max - min);
}

float standardize(float val, float avg, float sd)
{
    return (val - avg) / (sd * 2);
}

void loadFilter()
{
    fstream filtin;
    string line;
    filtin.open("FilterOutput_after.csv", ios::in);
    if (filtin)
    {
        getline(filtin, line);
        getline(filtin, line);
        while (!filtin.eof())
        {
            getline(filtin, line);
            vector<string> vec = split(line, ',');
            if (size(vec) == 0)
                continue;
            FilterItem fi;
            fi.cls = vec[1];
            fi.numFaces = stoi(vec[2]);
            fi.numVertices = stoi(vec[3]);
            fi.typeOfFace = vec[4];
            fi.minX = stof(vec[5]);
            fi.maxX = stof(vec[6]);
            fi.minY = stof(vec[7]);
            fi.maxY = stof(vec[8]);
            fi.minZ = stof(vec[9]);
            fi.maxZ = stof(vec[10]);
            fi.bX = stof(vec[11]);
            fi.bY = stof(vec[12]);
            fi.bZ = stof(vec[13]);
            fi.path = vec[14];
            fis[stoi(vec[0])] = fi;
        }
        filtin.close();
    }
    else
    {
        std::cout << "No previous filter output found" << endl;
    }
}

void loadDB(string file) 
{

    fstream fin;
    fin.open(file, ios::in);

    vector<string> row;
    vector<float> features;
    ANNpointArray data_points;
    string line, word, temp, name, shape_class;
    float val;

    ANNpoint* pa; // an array of points
    ANNdist* ANNdistArray; // an array of squared distances

    getline(fin, line);

    row.clear();
    getline(fin, line);
    stringstream s(line);
    while (getline(s, word, ';')) {
        row.push_back(word);
    }

    SAval1 = stof(row[0]);
    SAval2 = stof(row[1]);
    COval1 = stof(row[2]);
    COval2 = stof(row[3]);
    BBVval1 = stof(row[4]);
    BBVval2 = stof(row[5]);
    DIAval1 = stof(row[6]);
    DIAval2 = stof(row[7]);
    ECCval1 = stof(row[8]);
    ECCval2 = stof(row[9]);

    row.clear();
    s.clear();
    getline(fin, line);
    stringstream s2(line);
    while (getline(s2, word, ';')) {
        row.push_back(word);
    }

    D1min = stof(row[0]);
    D1max = stof(row[1]);
    D2min = stof(row[2]);
    D2max = stof(row[3]);
    D3min = stof(row[4]);
    D3max = stof(row[5]);
    D4min = stof(row[6]);
    D4max = stof(row[7]);
    A3min = stof(row[8]);
    A3max = stof(row[9]);

    getline(fin, line);

    while (getline(fin, line)) 
    {
        row.clear();
        features.clear();
        stringstream s(line);

        while (getline(s, word, ';')) {
            row.push_back(word);
        }
        name = row[0];
        int index = name.find_last_of("\\/");
        string shape_class = name.substr(index + 1);

        if (db_count.count(shape_class) == 0) {
            db_count.insert({ shape_class, 0 });
        }
        db_count[shape_class] += 1;

        for (int i = 1; i < row.size(); i++) {                  // CHANGE i TO 2 WHEN ACCOUNTING FOR SHAPE
            // cout << row[i] << endl;
            features.push_back(stof(row[i]));
        }

        pair<string, vector<float>> p = make_pair(name, features);
        feature_vectors.push_back(p);
        num_entries += 1;
    }

    int maxPoints = feature_vectors.size();
    int dims = row.size() - 1;
    data_points = annAllocPts(maxPoints, dims);

    for (int i = 0; i < maxPoints; i++) {
        ANNpoint p = annAllocPt(dims);
        for (int j = 0; j < dims; j++) {
            p[j] = feature_vectors[i].second[j];
        }
        data_points[i] = p;
    }

    tree = new ANNkd_tree(data_points, maxPoints, dims);
}

float eucleanDist(vector<float> s1, vector<float> s2)
{
    float distance = 0;

    for (int i = 0; i < s1.size(); i++)
    {
        distance += pow((s1[i] - s2[i]), 2.0);
    }
   
    return distance;
}

float crossBinDist(vector<float> s1, vector<float> s2) {
    float sum = 0;
    float distance = 0;
    float w1 = 0.15;
    float w2 = 0.7;

    for (int i = 0; i < s1.size(); i++) {
        for (int j = 0; j < s2.size(); j++) {

            float val = pow((s1[i] - s2[j]), 2.0);
            if (i == j) {
                distance += pow((s1[i] - s2[j]), 2.0) * w2;
            }
            else if (abs(i - j) == 1) {
                distance += pow((s1[i] - s2[j]), 2.0) * w1;
            }
        }
    }

    /*for (int i = 0; i < s1.size(); i++)
    {
        distance += pow((s1[i] - s2[i]), 2);
    }*/
    return distance;
}

void featureExtractNormalized(string input)
{
    vector<Grid*> grids;
    vector<string> names;
    vector<string> shapes;

    for (const auto& entry : fs::directory_iterator(input))
    {
        string folder = entry.path().string();
        for (const auto& fl : fs::directory_iterator(folder + "/"))
        {
            string file = fl.path().string();
            string suffix = ".off";
            if (!(file.size() >= suffix.size() && 0 == file.compare(file.size() - suffix.size(), suffix.size(), suffix)))
                continue;
            if (file.find("/") <= size(file))
            {
                names.push_back(file.substr(file.find_last_of("/") + 1));
            }
            else
            {
                names.push_back(file.substr(file.find_last_of("\\") + 1));
            }

            grid = std::get<0>(openFile(file));
            
            grid->calculateD1();
            if (grid->D1min < D1min)
                D1min = grid->D1min;
            if (grid->D1max > D1max)
                D1max = grid->D1max;

            grid->calculateD2(1000000);
            if (grid->D2min < D2min)
                D2min = grid->D2min;
            if (grid->D2max > D2max)
                D2max = grid->D2max;

            grid->calculateD3(1000000);
            if (grid->D3min < D3min)
                D3min = grid->D3min;
            if (grid->D3max > D3max)
                D3max = grid->D3max;

            grid->calculateD4(1000000);
            if (grid->D4min < D4min)
                D4min = grid->D4min;
            if (grid->D4max > D4max)
                D4max = grid->D4max;

            grid->calculateA3(1000000);
            if (grid->A3min < A3min)
                A3min = grid->A3min;
            if (grid->A3max > A3max)
                A3max = grid->A3max;

            float SA = grid->calculateSurfaceArea();
            if (SA < SAmin)
                SAmin = SA;
            if (SA > SAmax)
                SAmax = SA;

            float CO = grid->calculateSphericity();
            if (CO < COmin)
                COmin = CO;
            if (CO > COmax)
                COmax = CO;

            float BBV = grid->calculateBoundingBoxVol();
            if (BBV < BBVmin)
                BBVmin = BBV;
            if (BBV > BBVmax)
                BBVmax = BBV;

            float DIAM = grid->calculateDiameter();
            if (DIAM < DIAmin)
                DIAmin = DIAM;
            if (DIAM > DIAmax)
                DIAmax = DIAM;

            float ECC = grid->calculateEccentricity();
            if (ECC < ECCmin)
                ECCmin = ECC;
            if (ECC > ECCmax)
                ECCmax = ECC;

            grids.push_back(grid);
        }
    }

    fstream filtout;
    filtout.open("output.csv", ios::out);
    filtout << "sep=;" << endl;
    filtout << SAmin << ";" << SAmax << ";" << COmin << ";" << COmax << ";" << BBVmin << ";" << BBVmax << ";" << DIAmin << ";" << DIAmax << ";" << ECCmin << ";" << ECCmax << endl;
    filtout << D1min << ";" << D1max << ";" << D2min << ";" << D2max << ";" << D3min << ";" << D3max << ";" << D4min << ";" << D4max << ";" << A3min << ";" << A3max << endl;
    filtout << "name;Surface Area;Compactness;Bounding Box Volume;Diameter;Eccentricity;D1-1;D1-2;D1-3;D1-4;D1-5;D1-6;D1-7;D1-8;D1-9;D1-10;D1-11;D1-12;D1-13;D1-14;D2-1;D2-2;D2-3;D2-4;D2-5;D2-6;D2-7;D2-8;D2-9;D2-10;D2-11;D2-12;D2-13;D2-14;D3-1;D3-2;D3-3;D3-4;D3-5;D3-6;D3-7;D3-8;D3-9;D3-10;D3-11;D3-12;D3-13;D3-14;D4-1;D4-2;D4-3;D4-4;D4-5;D4-6;D4-7;D4-8;D4-9;D4-10;D4-11;D4-12;D4-13;D4-14;A3-1;A3-2;A3-3;A3-4;A3-5;A3-6;A3-7;A3-8;A3-9;A3-10;A3-11;A3-12;A3-13;A3-14" << endl;

    for (int i = 0; i < grids.size(); i++)
    {
        vector<float> D1hist = grids[i]->getD1hist(D1min, D1max, bins);
        vector<float> D2hist = grids[i]->getD2hist(D2min, D2max, bins);
        vector<float> D3hist = grids[i]->getD3hist(D3min, D3max, bins);
        vector<float> D4hist = grids[i]->getD4hist(D4min, D4max, bins);
        vector<float> A3hist = grids[i]->getA3hist(A3min, A3max, bins);

        filtout << names[i] << ";";

        filtout << grids[i]->calculateSurfaceArea() << ";";
        filtout << grids[i]->calculateSphericity() << ";";
        filtout << grids[i]->calculateBoundingBoxVol() << ";";
        filtout << grids[i]->calculateDiameter() << ";";
        filtout << grids[i]->calculateEccentricity() << ";";

        for (int j = 0; j < bins; j++)
        {
            filtout << D1hist[j] << ";";
        }
        for (int j = 0; j < bins; j++)
        {
            filtout << D2hist[j] << ";";
        }
        for (int j = 0; j < bins; j++)
        {
            filtout << D3hist[j] << ";";
        }
        for (int j = 0; j < bins; j++)
        {
            filtout << D4hist[j] << ";";
        }
        for (int j = 0; j < bins; j++)
        {
            filtout << A3hist[j] << ";";
        }

        filtout << endl;
    }
    filtout.close();
}

void featureExtractStandardized(string input)
{
    vector<Grid*> grids;
    vector<string> names;
    vector<string> shapes;
    vector<float> SAs, COs, BBVs, DIAs, ECCs;

    for (const auto& entry : fs::directory_iterator(input))
    {
        string folder = entry.path().string();
        for (const auto& fl : fs::directory_iterator(folder + "/"))
        {
            //string shape = folder.substr(3);
            //shapes.push_back(shape);
            string file = fl.path().string();
            string suffix = ".off";
            if (!(file.size() >= suffix.size() && 0 == file.compare(file.size() - suffix.size(), suffix.size(), suffix)))
                continue;
            if (file.find("/") <= size(file))
            {
                names.push_back(file.substr(file.find_last_of("/") + 1));
            }
            else
            {
                names.push_back(file.substr(file.find_last_of("\\") + 1));
            }
            string shape = folder.substr(3);
            shapes.push_back(shape);
            grid = std::get<0>(openFile(file));

            grid->calculateD1();
            if (grid->D1min < D1min)
                D1min = grid->D1min;
            if (grid->D1max > D1max)
                D1max = grid->D1max;

            grid->calculateD2(1000000);
            if (grid->D2min < D2min)
                D2min = grid->D2min;
            if (grid->D2max > D2max)
                D2max = grid->D2max;

            grid->calculateD3(1000000);
            if (grid->D3min < D3min)
                D3min = grid->D3min;
            if (grid->D3max > D3max)
                D3max = grid->D3max;

            grid->calculateD4(1000000);
            if (grid->D4min < D4min)
                D4min = grid->D4min;
            if (grid->D4max > D4max)
                D4max = grid->D4max;

            grid->calculateA3(1000000);
            if (grid->A3min < A3min)
                A3min = grid->A3min;
            if (grid->A3max > A3max)
                A3max = grid->A3max;

            SAs.push_back(grid->calculateSurfaceArea());
            COs.push_back(grid->calculateSphericity());
            BBVs.push_back(grid->calculateBoundingBoxVol());
            DIAs.push_back(grid->calculateDiameter());
            ECCs.push_back(grid->calculateEccentricity());

            grids.push_back(grid);
        }
    }

    for (int i = 0; i < grids.size(); i++)
    {
        SAavg += SAs[i];
        COavg += COs[i];
        BBVavg += BBVs[i];
        DIAavg += DIAs[i];
        ECCavg += ECCs[i];
    }

    SAavg /= grids.size();
    COavg /= grids.size();
    BBVavg /= grids.size();
    DIAavg /= grids.size();
    ECCavg /= grids.size();

    for (int i = 0; i < grids.size(); i++)
    {
        SAsd += pow(SAs[i] - SAavg, 2);
        COsd += pow(COs[i] - COavg, 2);
        BBVsd += pow(BBVs[i] - BBVavg, 2);
        DIAsd += pow(DIAs[i] - DIAavg, 2);
        ECCsd += pow(ECCs[i] - ECCavg, 2);
    }

    SAsd = sqrt(SAsd / grids.size());
    COsd = sqrt(COsd / grids.size());
    BBVsd = sqrt(BBVsd / grids.size());
    DIAsd = sqrt(DIAsd / grids.size());
    ECCsd = sqrt(ECCsd / grids.size());

    fstream filtout;
    filtout.open("outputStand.csv", ios::out);
    filtout << "sep=;" << endl;
    filtout << SAavg << ";" << SAsd << ";" << COavg << ";" << COsd << ";" << BBVavg << ";" << BBVsd << ";" << DIAavg << ";" << DIAsd << ";" << ECCavg << ";" << ECCsd << endl;
    filtout << D1min << ";" << D1max << ";" << D2min << ";" << D2max << ";" << D3min << ";" << D3max << ";" << D4min << ";" << D4max << ";" << A3min << ";" << A3max << endl;
    filtout << "name;Surface Area;Compactness;Bounding Box Volume;Diameter;Eccentricity;D1-1;D1-2;D1-3;D1-4;D1-5;D1-6;D1-7;D1-8;D1-9;D1-10;D1-11;D1-12;D1-13;D1-14;D2-1;D2-2;D2-3;D2-4;D2-5;D2-6;D2-7;D2-8;D2-9;D2-10;D2-11;D2-12;D2-13;D2-14;D3-1;D3-2;D3-3;D3-4;D3-5;D3-6;D3-7;D3-8;D3-9;D3-10;D3-11;D3-12;D3-13;D3-14;D4-1;D4-2;D4-3;D4-4;D4-5;D4-6;D4-7;D4-8;D4-9;D4-10;D4-11;D4-12;D4-13;D4-14;A3-1;A3-2;A3-3;A3-4;A3-5;A3-6;A3-7;A3-8;A3-9;A3-10;A3-11;A3-12;A3-13;A3-14" << endl;

    for (int i = 0; i < grids.size(); i++)
    {
        vector<float> D1hist = grids[i]->getD1hist(D1min, D1max, bins);
        vector<float> D2hist = grids[i]->getD2hist(D2min, D2max, bins);
        vector<float> D3hist = grids[i]->getD3hist(D3min, D3max, bins);
        vector<float> D4hist = grids[i]->getD4hist(D4min, D4max, bins);
        vector<float> A3hist = grids[i]->getA3hist(A3min, A3max, bins);

        filtout << shapes[i] + "/" + names[i] << ";";

        filtout << standardize(SAs[i], SAavg, SAsd) << ";";
        filtout << standardize(COs[i], COavg, COsd) << ";";
        filtout << standardize(BBVs[i], BBVavg, BBVsd) << ";";  
        filtout << standardize(DIAs[i], DIAavg, DIAsd) << ";";
        filtout << standardize(ECCs[i], ECCavg, ECCsd) << ";";

        for (int j = 0; j < bins; j++)
        {
            filtout << D1hist[j] << ";";
        }
        for (int j = 0; j < bins; j++)
        {
            filtout << D2hist[j] << ";";
        }
        for (int j = 0; j < bins; j++)
        {
            filtout << D3hist[j] << ";";
        }
        for (int j = 0; j < bins; j++)
        {
            filtout << D4hist[j] << ";";
        }
        for (int j = 0; j < bins; j++)
        {
            filtout << A3hist[j] << ";";
        }

        filtout << endl;
    }
    filtout.close();
}

float cosineSimilarity(vector<float> s1, vector<float> s2)
{
    float dot = 0.0, denom_a = 0.0, denom_b = 0.0;
    
    for (int i = 0; i < s1.size(); ++i) 
    {
        dot += s1[i] * s2[i];
        denom_a += s1[i] * s1[i];
        denom_b += s2[i] * s2[i];
    }
    
    return dot / (sqrt(denom_a) * sqrt(denom_b));
}

vector<pair<string, float>> startNewQuery(string file_name, int K, bool ann_flag) {

    int numF = feature_vectors.size();

    vector<float> query_vector_f;
    vector<float> query_vector_h1, query_vector_h2, query_vector_h3, query_vector_h4, query_vector_h5;
    vector<float> ann_vector;

    tuple tup = openFile(file_name);
    Grid* query_grid = get<0>(tup);
    ANNpoint query_point = annAllocPt(feature_vectors[0].second.size());

    /*float s = query_grid->calculateSurfaceArea();
    float r = query_grid->calculateSphericity();
    float b = query_grid->calculateBoundingBoxVol();
    float d = query_grid->calculateDiameter();
    float e = query_grid->calculateEccentricity();*/

    float s = standardize(query_grid->calculateSurfaceArea(), SAval1, SAval2);
    float r = standardize(query_grid->calculateSphericity(), COval1, COval2);
    float b = standardize(query_grid->calculateBoundingBoxVol(), BBVval1, BBVval2);
    float d = standardize(query_grid->calculateDiameter(), DIAval1, DIAval2);
    float e = standardize(query_grid->calculateEccentricity(), ECCval1, ECCval2);

    query_grid->calculateD1();
    query_grid->calculateD2(N);
    query_grid->calculateD3(N);
    query_grid->calculateD4(N);
    query_grid->calculateA3(N);

    vector<float> D1hist = query_grid->getD1hist(D1min, D1max, bins);
    vector<float> D2hist = query_grid->getD2hist(D2min, D2max, bins);
    vector<float> D3hist = query_grid->getD3hist(D3min, D3max, bins);
    vector<float> D4hist = query_grid->getD4hist(D4min, D4max, bins);
    vector<float> A3hist = query_grid->getA3hist(A3min, A3max, bins);

    query_vector_f.push_back(s);
    ann_vector.push_back(s);
    query_vector_f.push_back(r);
    ann_vector.push_back(r);
    query_vector_f.push_back(b);
    ann_vector.push_back(b);
    query_vector_f.push_back(d);
    ann_vector.push_back(d);
    query_vector_f.push_back(e);
    ann_vector.push_back(e);

    for (int i = 0; i < bins; i++) {
        query_vector_h1.push_back(D1hist[i]);
        ann_vector.push_back(D1hist[i]);
    }
    for (int i = 0; i < bins; i++) {
        query_vector_h2.push_back(D2hist[i]);
        ann_vector.push_back(D2hist[i]);
    }
    for (int i = 0; i < bins; i++) {
        query_vector_h3.push_back(D3hist[i]);
        ann_vector.push_back(D3hist[i]);
    }
    for (int i = 0; i < bins; i++) {
        query_vector_h4.push_back(D4hist[i]);
        ann_vector.push_back(D4hist[i]);
    }
    for (int i = 0; i < bins; i++) {
        query_vector_h5.push_back(A3hist[i]);
        ann_vector.push_back(A3hist[i]);
    }

    vector<pair<string, float>> result;
    float distance = 0;

    for (int i = 0; i < numF; i++) {
        float distance = 0;

        vector<float> vec1 = get<1>(feature_vectors[i]);

        vector<float> vecf(vec1.begin(), vec1.begin() + 5);
        vector<float> vec_h1(vec1.begin() + 5 , vec1.begin() + 5 + bins);
        vector<float> vec_h2(vec1.begin() + 5 + bins , vec1.begin() + 5 + bins*2);
        vector<float> vec_h3(vec1.begin() + 5 + bins*2, vec1.begin() + 5 + bins*3);
        vector<float> vec_h4(vec1.begin() + 5 + bins*3, vec1.begin() + 5 + bins*4);
        vector<float> vec_h5(vec1.begin() + 5 + bins*4, vec1.end());


        //distance += cosineSimilarity(vecf, query_vector_f);

        distance = eucleanDist(vecf, query_vector_f) * w_f;


        distance += crossBinDist(vec_h1, query_vector_h1) * w_h;
        distance += crossBinDist(vec_h2, query_vector_h2) * w_h;
        distance += crossBinDist(vec_h3, query_vector_h3) * w_h;
        distance += crossBinDist(vec_h4, query_vector_h4) * w_h;
        distance += crossBinDist(vec_h5, query_vector_h5) * w_h;

        distance = sqrt(distance);                                   // To correctly compute Euclidean


        string name = get<0>(feature_vectors[i]);

        pair<string, float> p2 = make_pair(name, distance);
        result.push_back(p2);
    }

    if (ann_flag) {
        //ANN QUERY
        for (int i = 0; i < ann_vector.size(); i++) {
            query_point[i] = ann_vector[i];
        }
        ANNidx* nnIdx = new ANNidx[K];
        ANNdist* dists = new ANNdist[K];

        tree->annkSearch(query_point, K, nnIdx, dists, 0);

        cout << "#############" << endl;
        cout << "CLOSEST SHAPES USING ANN: " << endl;
        for (int i = 0; i < K; i++) {
            int index = nnIdx[i];
            cout << result[index].first << endl;
            cout << "distance: ";
            cout << dists[i] << endl;
            cout << endl;
        }
        cout << "#############" << endl;
    }
    sort(result.begin(), result.end(), sortbysec);
    vector<pair<string, float>> output(result.begin(), result.begin() + K);
    grid = get<0>(tup);
    //delete query_grid;
    return output;
}

void performEvaluation(int K) {

    loadDB("outputStand.csv");

    vector<vector<pair<string, float>>> results;
    vector<pair<string, float>> shape_nums;
    vector<float> total_accuracies;
    vector<float> total_TPRs;
    vector<float> total_PPVs;
    vector<string> class_names;


    string folder;
    float total_acc = 0.0;
    float total_TPR = 0.0;
    float total_PPV = 0.0;
    int total_queries = 0;
    int total_TP = 0;
    int total_TN = 0;
    int total_FP = 0;
    int total_FN = 0;
    for (const auto& entry : fs::directory_iterator("Evaluation_DB"))
    {
        folder = entry.path().string();
        int index = folder.find_last_of("\\/");
        string query_shape = folder.substr(index + 1);
        class_names.push_back(query_shape);
        //cout << folder << endl;
        float current_acc = 0.0;
        float current_PPV = 0.0;
        float current_TPR = 0.0;
        int current_queries = 0;
        total_queries += 1;
        for (const auto& fl : fs::directory_iterator(folder + "/"))
        {
            int current_TP = 0;
            int current_TN = 0;
            int current_FP = 0;
            int current_FN = 0;
            string file = fl.path().string();
            string name = file;
            string suffix = ".off";
            if (!(file.size() >= suffix.size() && 0 == file.compare(file.size() - suffix.size(), suffix.size(), suffix)))
                continue;
            cout << file << endl;
            vector<pair<string, float>> result = startNewQuery(file, K, false);
            current_queries += 1;

            for (int i = 0; i < 10; i++) {
                int index = result[i].first.find_last_of("\\/");
                string result_shape = result[i].first.substr(index + 1);
                cout << query_shape + ", " + result_shape << endl;;
                if (query_shape == result_shape) {
                    current_TP += 1;
                    total_TP += 1;
                }
                else {
                    current_FP += 1;
                    total_FP += 1;
                }
            }
            current_FN = db_count[query_shape] - current_TP;
            current_TN = num_entries - db_count[query_shape] - current_FN;
            current_acc += float(current_TP + current_TN) / float(num_entries);
            current_TPR += float(current_TP) / float(db_count[query_shape]);
            current_PPV += float(current_TP) / float(K);
        }
        current_acc /= current_queries;
        current_TPR /= current_queries;
        current_PPV /= current_queries;
        total_accuracies.push_back(current_acc);
        total_TPRs.push_back(current_TPR);
        total_PPVs.push_back(current_PPV);
        total_acc += current_acc;
        total_TPR += current_TPR;
        total_PPV += current_PPV;
    }

    total_acc /= total_queries;
    total_TPR /= total_queries;
    total_PPV /= total_queries;
    cout << total_acc << endl;
    cout << total_TPR << endl;
    cout << total_PPV << endl;

    fstream evout;
    evout.open("evaluationOutput.csv", ios::out);
    evout << "sep=;" << endl;
    evout << "shape class;average accuracy;average TPR;average PPV" << endl;
    for (int i = 0; i < total_accuracies.size(); i++) {
        evout << class_names[i] << ";";
        evout << total_accuracies[i] << ";";
        evout << total_TPRs[i] << ";";
        evout << total_PPVs[i] << endl;
    }

    evout << endl;
    evout << "overall accuracy;overall TPR;overall PPV" << endl;
    evout << total_acc << ";" << total_TPR << ";" << total_PPV << endl;


}

void displayQueryResult(int argc, char* argv[], vector<pair<string, float>> result) {

    glutInit(&argc, argv);								                //Initialize the GLUT toolkit
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    //Ask GLUT to create next windows with a RGB framebuffer and a Z-buffer too
    glutInitWindowSize(500, 500);							            //Tell GLUT how large are the windows we want to create next
    glutCreateWindow(fileName.c_str());	                                //Create our window

    //glutMouseFunc(mouseclick);							                //Bind the mouse click and mouse drag (click-and-move) events to callbacks. This allows us
    //glutMotionFunc(mousemotion);
    glutKeyboardFunc(keyboard);
    glutDisplayFunc(draw);
    //glutReshapeFunc(viewing);
    glutMainLoop();

}

int main(int argc, char* argv[])
{

    cout << "Press 'q' to start a new query." << endl;

    cout << "Press 'n' to load a new data base without normalization." << endl;

    cout << "Press 's' to load a new data base with standardization." << endl;

    cout << "Press 'e' to perform an evaluation of the system." << endl;

    char input = _getch();

    if (input == 'q') {

        string db_file;
        loadDB("outputStand.csv");
        //loadDB("outputq.csv");
        cout << "Please specify the query file::" << endl;
        string file_name;
        cin >> file_name;
        vector<pair<string, float>> result = startNewQuery(file_name, 10, true);

        for (int i = 0; i < result.size(); i++) {
            cout << result[i].first << endl;
        }

        cout << "#############" << endl;
        cout << "CLOSEST SHAPES USING CUSTOM METRIC: " << endl;
        for (int i = 0; i < 10; i++) {
            cout << result[i].first << endl;
            cout << "distance: ";
            cout << result[i].second << endl;
            cout << endl;
        }
        cout << "#############" << endl;

        displayQueryResult(argc, argv, result);
    }
    else if (input == 'n') {

        cout << "Please specify the database folder" << endl;
        string finput;
        cin >> finput;
        featureExtractNormalized(finput);
    }
    else if (input == 's') {

        cout << "Please specify the database folder" << endl;
        string finput;
        cin >> finput;
        featureExtractStandardized(finput);
    }
    else if (input == 'e') {
        performEvaluation(10);
    }

    
    return 0;
}
