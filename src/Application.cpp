//#include <GL/glew.h>
//#include <GL/glut.h>
#include <GLFW/glfw3.h>
#include <filesystem>
#include <iostream>
#include <io.h>
#include <filesystem>
#include <conio.h>

#include <fstream>
#include <filesystem>
#include <tuple>
#include <stdio.h>
#include <string>
#include "..\header\Grid.h"
#include "..\header\FilterItem.h"
#include "..\header\OFFReader.h"

using namespace std;
namespace fs = std::filesystem;
void loadFilter();

string fileName;

vector<tuple<string, string, vector<float>>> feature_vectors;

int drawing_style = 0;
const int FILTER_SIZE = 250;

int N = 1000000;                            // Number of computations per feature


Grid* grid = 0;
//Renderer renderer;

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

bool sortbysec(const pair<pair<string, string>, float>& a, const pair<pair<string, string>, float>& b)
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
        cout << folder << endl;
        for (const auto& fl : fs::directory_iterator(folder + "/"))
        {
            string file = fl.path().string();
            string suffix = ".off";
            if (!(file.size() >= suffix.size() && 0 == file.compare(file.size() - suffix.size(), suffix.size(), suffix)))
                continue;
            cout << file << endl;
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
    /*case 'o':
    {
        cout << "Please enter a name for the file" << endl;
        string fileName;
        cin >> fileName;
        outputFilter(fileName);
        break;
    }*/
    case 's':
    {
        cout << "Please enter a folder to scan" << endl;
        string fileName;
        cin >> fileName;
        scanFolder(fileName);
        break;
    }
    case 'l':
    {
        cout << "Loading from output file" << endl;
        loadFilter();
        break;
    }
    }
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

void loadDB() {

    fstream fin;

    string file = "output.csv";
    fin.open(file, ios::in);

    vector<string> row;
    vector<float> features;
    tuple<string, string, vector<float>> tup;
    string line, word, temp, name, shape_class;
    float val;
    getline(fin, line);
    getline(fin, line);

    while (getline(fin, line)) {
        row.clear();
        features.clear();

        stringstream s(line);

        while (getline(s, word, ';')) {
            row.push_back(word);
        }
        name = row[0];
        shape_class = row[1];

        for (int i = 1; i < row.size(); i++) {                  // CHANGE i TO 2 WHEN ACCOUNTING FOR SHAPE
            // cout << row[i] << endl;
            features.push_back(stof(row[i]));
        }

        tuple tup = make_tuple(name, shape_class, features);
        feature_vectors.push_back(tup);
    }
}

float eucleanDist(vector<float> s1, vector<float> s2)
{
    float Sum = 0;
    float distance = 0;

    for (int i = 0; i < s1.size(); i++)
    {
        Sum += pow((s1[i] - s2[i]), 2.0);
        distance += sqrt(Sum);
    }
    return distance;
}

float matchFeatures(vector<float> vec1, vector<float> vec2) {

    float distance = 0;


    return distance;

}

void startNewQuery() {

    cout << "Please specify the query file::" << endl;
    string file_name;
    cin >> file_name;

    vector<float> query_vector;

    tuple tup = openFile(file_name);
    Grid* query_grid = get<0>(tup);

    float s = query_grid->calculateSurfaceArea();
    float v = query_grid->calculateVolume();
    float r = query_grid->calculateSphericity();
    float b = query_grid->calculateBoundingBoxVol();
    float d = query_grid->calculateDiameter();
    float e = query_grid->calculateEccentricity();

    query_grid->calculateD1();
    query_grid->calculateD2(N);
    query_grid->calculateD3(N);
    query_grid->calculateD4(N);
    query_grid->calculateA3(N);

    /*query_vector.push_back(s);
    query_vector.push_back(v);
    query_vector.push_back(query_grid->compactness);
    query_vector.push_back(r);
    query_vector.push_back(b);
    query_vector.push_back(d);
    query_vector.push_back(e);*/

    for (int i = 0; i < size(query_grid->D1hist); i++) {
        query_vector.push_back(query_grid->D1hist[i]);
    }
    for (int i = 0; i < size(query_grid->D2hist); i++) {
        query_vector.push_back(query_grid->D2hist[i]);
    }
    for (int i = 0; i < size(query_grid->D3hist); i++) {
        query_vector.push_back(query_grid->D3hist[i]);
    }
    for (int i = 0; i < size(query_grid->D4hist); i++) {
        query_vector.push_back(query_grid->D4hist[i]);
    }
    for (int i = 0; i < size(query_grid->A3hist); i++) {
        query_vector.push_back(query_grid->A3hist[i]);
    }

    vector<pair<pair<string, string>, float>> result;
    float distance;
    int numF = feature_vectors.size();

    for (int i = 0; i < numF; i++) {

        cout << i << endl;

        vector<float> vec1 = get<2>(feature_vectors[i]);
        distance = eucleanDist(vec1, query_vector);

        string name = get<0>(feature_vectors[i]);
        string shape = get<1>(feature_vectors[i]);

        cout << name + shape << endl;
        cout << distance << endl;

        pair<string, string> p1 = make_pair(name, shape);
        pair<pair<string, string>, float> p2 = make_pair(p1, distance);
        result.push_back(p2);
    }


    // COMPARE query_vector WITH feature_vectors HERE

    sort(result.begin(), result.end(), sortbysec);

    cout << "#############" << endl;
    cout << "CLOSEST 5 SHAPES: " << endl;
    for (int i = 0; i < 5; i++) {
        cout << result[i].first.first + ",  " + result[i].first.second << endl;
        cout << "distance: ";
        cout << result[i].second << endl;
        cout << endl;
    }
    cout << "#############" << endl;


}

int main(int argc, char* argv[])
{

    cout << "Loading default data base..." << endl;
    loadDB();

    cout << "Press q to start a new query." << endl;

    cout << "Press l to load a new data base." << endl;

    char input = _getch();

    if (input == 'q') {
        startNewQuery();
    }
    else if (input == 'l') {

        cout << "Please specify the new database folder" << endl;
        cin >> input;
    }

    
    return 0;
}
