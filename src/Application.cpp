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

    cout << "Please specify the folder to extract the features from" << endl;
    cin >> input;

    fstream filtout;
    filtout.open("output.csv", ios::out);
    filtout << "sep=;" << endl;
    filtout << "D1-1;D1-2;D1-3;D1-4;D1-5;D1-6;D1-7;D1-8;D1-9;D1-10;D2-1;D2-2;D2-3;D2-4;D2-5;D2-6;D2-7;D2-8;D2-9;D2-10;D3-1;D3-2;D3-3;D3-4;D3-5;D3-6;D3-7;D3-8;D3-9;D3-10;D4-1;D4-2;D4-3;D4-4;D4-5;D4-6;D4-7;D4-8;D4-9;D4-10;A3-1;A3-2;A3-3;A3-4;A3-5;A3-6;A3-7;A3-8;A3-9;A3-10";

    for (const auto& entry : fs::directory_iterator(input))
    {
        string folder = entry.path().string();
        cout << folder << endl;
        for (const auto& fl : fs::directory_iterator(folder + "/"))
        {
            string file = fl.path().string();
            string suffix = ".off";
            if (!(file.size() >= suffix.size() && 0 == file.compare(file.size() - suffix.size(), suffix.size(), suffix)))
                continue;
            cout << file << endl;

            

            if (file.find("/") <= size(file))
            {
                filtout << file.substr(file.find_last_of("/") + 1);
            }
            else
            {
                filtout << file.substr(file.find_last_of("\\") + 1);
            }

            grid = std::get<0>(openFile(file));
            
            grid->calculateD1();
            vector<float> hist1 = grid->getD1hist();
            for (int i = 0; i < 10; i++)
                filtout << hist1[i] << ";";

            grid->calculateD2(1000000);
            vector<float> hist2 = grid->getD2hist();
            for (int i = 0; i < 10; i++)
                filtout << hist2[i] << ";";

            grid->calculateD3(1000000);
            vector<float> hist3 = grid->getD3hist();
            for (int i = 0; i < 10; i++)
                filtout << hist3[i] << ";";

            grid->calculateD4(1000000);
            vector<float> hist4 = grid->getD4hist();
            for (int i = 0; i < 10; i++)
                filtout << hist4[i] << ";";

            grid->calculateA3(1000000);
            vector<float> hist5 = grid->getA3hist();
            for (int i = 0; i < 10; i++)
                filtout << hist5[i] << ";";

            filtout << endl;
            delete(grid);
        }
    }

    filtout.close();

    for (int i = 0; i < 3; i++) {

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

        cout << "D4: " << endl;
        grid->calculateD4(1000000);
        vector<float> hist4 = grid->getD4hist();
        for (int i = 0; i < 10; i++) {
            cout << hist4[i] << endl;
        }
        
        cout << "A3: " << endl;
        grid->calculateA3(1000000);
        vector<float> hist5 = grid->getA3hist();
        for (int i = 0; i < 10; i++) {
         cout << hist5[i] << endl;
     }
        

    }
    
    return 0;
}
