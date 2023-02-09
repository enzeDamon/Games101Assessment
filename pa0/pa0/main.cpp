#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>

int main(){
    //  pi
    float pi = std::acos(-1);
    // Point (2, 1)
    Eigen::Vector3f point(2, 1, 1);
    // rotation matrix
    Eigen::Matrix3f transform;
    // rotate by 45 degree
    transform << std::cos(45.0/180.0 * pi), -std::sin(45.0/180.0 * pi), 0, std::sin(45.0/180.0 * pi), std::cos(45.0/180.0 * pi), 0, 0, 0, 1;
    std::cout << "rotation matrix is like" << std:: endl;
    std::cout << transform << std::endl;
    Eigen::Vector3f temp;
    temp = transform * point;
    // rotated result is stored in temp
    std::cout << "rotated point is like:" << std::endl;
    std::cout << temp << std::endl ;
    // get in the shifting matrix
    transform << 1, 0, 1, 0, 1, 2, 0, 0, 1;
    std ::cout << "the final result is:" << std ::endl;
    std::cout << transform * temp << std::endl;

    // Example of matrix
    // std::cout << "Example of matrix \n";
    // matrix definition
    // Eigen::Matrix3f i,j;
    // i << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0; 
    // j << 2.0, 3.0, 1.0, 4.0, 6.0, 5.0, 9.0, 7.0, 8.0;
    // matrix output
    // std::cout << "Example of output \n";
    // std::cout << i << std::endl;
    // matrix add i + j
    // matrix scalar multiply i * 2.0
    // matrix multiply i * j
    // matrix multiply vector i * v

    return 0;
}