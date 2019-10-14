#include <QCoreApplication>
#include <iostream>
#include <string>

using namespace std;

#include "dob.h"
#include "RobotArm/robotarm.h"

//void load_data(string file_name, vector<double> *data, string delimiter)
//{
//    data->clear();
//    FILE *fp_in;
//    const int buffer = 1000000;
//    char *ptr, basic[buffer];
//    fp_in = fopen(file_name.c_str(), "r");
//    while (fgets(basic, buffer, fp_in) != nullptr)
//    {
//        ptr = strtok(basic, delimiter.c_str());
//        while (ptr != nullptr) {
//            data->push_back(atof(ptr));
//            ptr = strtok(nullptr, delimiter.c_str());
//        }
//    }
//    fclose(fp_in);
//}

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    vector<double> q_data;
    load_data("/home/hajun/Desktop/DOB/recurdyn_result_collision.txt", &q_data, "\t");

    uint row = 2001, col = 8;

//    DOB *dob = new DOB(1, 1);

    double qi, qi_dot, Ta;
    int collision = 0;
    RobotArm robot(2, 1);

    for(uint indx = 0; indx < row; indx++){
        qi = q_data[indx*col + 2];
        qi_dot = q_data[indx*col + 6];
        Ta = 2;

    //    robot.run_dynamics();
        robot.run_dynamics_with_DOB(&qi, &qi_dot, &Ta);

//        dob->run(&qi, &qi_dot, &Ta, &collision);
    }

//    delete dob;

//    return a.exec();
    return 0;
}
