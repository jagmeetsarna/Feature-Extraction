//#include <GL/glew.h>
//#include <GL/glut.h>
#include <GLFW/glfw3.h>
#include <filesystem>
#include <iostream>
#include <io.h>
#include <filesystem>

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

int drawing_style = 0;
const int FILTER_SIZE = 250;


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

int main(int argc, char* argv[])
{
    string input;

    for (int i = 0; i < 3; i++) {

        string input;
        cout << "Please specify the file you want to view:" << endl;
        cin >> input;

        std::tuple<Grid*, FilterItem> tup = openFile(input);
        grid = std::get<0>(tup);


       /* cout << "Surface Area: ";
        cout << grid->calculateSurfaceArea() << endl;

        cout << "Volume: ";
        cout << grid->calculateVolume() << endl;

        //cout << "Compactness: ";
        float compactness;
        float sphericity;
        sphericity = grid->calculateSphericity();
        //cout << grid.compactness << endl;

        cout << "Sphericity: ";
        cout << sphericity << endl;

        cout << "Bounding Box Volume: ";
        cout << grid->calculateBoundingBoxVol() << endl;

        cout << "Diameter: ";
        cout << grid->calculateDiameter() << endl;


        cout << "Eccentricity: ";
        cout << grid->calculateEccentricity() << endl;*/

        cout << "D1: " << endl;
        grid->calculateD1();
        vector<float> hist1 = grid->getD1hist();
        for (int i = 0; i < 10; i++) {
            cout << hist1[i] << endl;
        }

        cout << "D2: " << endl;
        grid->calculateD2(1000000);
        vector<float> hist2 = grid->getD2hist();
        for (int i = 0; i < 10; i++) {
            cout << hist2[i] << endl;
        }

        cout << "D3: " << endl;
        grid->calculateD3(1000000);
        vector<float> hist3 = grid->getD3hist();
        for (int i = 0; i < 10; i++) {
            cout << hist3[i] << endl;
        }

    }
    
    return 0;
}
