//#include "Poly.h"
//#include <iostream>
//
//
//Polynom::Polynom(std::vector<std::vector<template>> coef)
//{
//    for(int i = 0; i < coef.size(); i++){
//        this->coef.push_back(coef[i]); // [1 2 5]
//    }
//}
//
//
//void Polynom::PrintPolynom(){
//    std::cout<<"[";
//    for(int i = 0; i < this->coef.size(); i++)
//    {
//        std::cout<<this->coef[i]<<(i == this->coef.size() -1 ? "" : " ");
//    }
//    std::cout<<"]"<<std::endl;
//}
//
//
//Polynom Polynom::operator+(const Polynom& p)
//{
//    int maxSize = std::max(this->coef.size(), p.coef.size());
//    std::vector<double> newPolyCoef(maxSize);
//    for(int i = 0; i < maxSize; i++)
//    {
//        double p1 = this->coef.size() > i ? this->coef[i] : 0;
//        double p2 = p.coef.size() > i ? p.coef[i] : 0;
//        newPolyCoef[i] = p1 + p2;
//    }
//    return Polynom(newPolyCoef);
//}
//
//
//Polynom Polynom::operator-(const Polynom& p){
//    int maxSize = std::max(this->coef.size(), p.coef.size());
//    std::vector<double> newPolyCoef(maxSize);
//    for(int i = 0; i < maxSize; i++)
//    {
//        double p1 = this->coef.size() > i ? this->coef[i] : 0;
//        double p2 = p.coef.size() > i ? p.coef[i] : 0;
//        newPolyCoef[i] = p1 - p2;
//    }
//    return Polynom(newPolyCoef);
//}
//
//
//Polynom Polynom::operator*(const Polynom& p){
//    int Size = this->coef.size() + p.coef.size() - 1;
//    std::vector<double> newPolyCoef(Size, 0);
//    for(int i = 0; i < Size; i++)
//    {
//        for(int j = 0; j <= i; j++)
//        {
//            newPolyCoef[i + j] += this->coef[i] * p.coef[j];
//        }
//    }
//    return Polynom(newPolyCoef);
//}
//
//
//Polynom Polynom::diff()
//{
//    std::vector<double> newPolyCoef(this->coef.size()-1, 0);
//    for(int i = 0; i < this->coef.size() - 1; i++){
//        newPolyCoef[i] = this->coef[i+1] * (i+1);
//    }
//    return Polynom(newPolyCoef);
//}
//
//
//Polynom Polynom::intg()
//{
//    std::vector<double> newPolyCoef(this->coef.size()+1, 0);
//    for(int i = 1; i < this->coef.size()+1; i++){
//        newPolyCoef[i] = this->coef[i-1] / i;
//    }
//    return Polynom(newPolyCoef);
//
//}
//
//
//
//
//
//
//
