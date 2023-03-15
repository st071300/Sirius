#pragma once
#include <iostream>
#include <vector>
#include <stack>
#include <string>
#include <complex>

using namespace std;

template <class T>
class Polynom
{
public:
    std::vector<T> coef;
    Polynom() = default;
    Polynom(const std::vector<T>& coef)
    {
        {
            this->coef=coef;
        }
    }

    Polynom(int n)
    {
        coef.resize(n);
        for (int i = 0; i < coef.size(); ++i) {
            coef[i]=0;
        }
    };

    void PrintPolynom(std::ostream& os){

        os<<"[";
        for(int i = 0; i < this->coef.size(); i++){
            os<<this->coef[i]<<(i == this->coef.size() -1 ? "" : ", ");
        }
        os<<"]"<<std::endl;
    }

    Polynom operator+(const Polynom<T>& p)
    {
        int maxSize = std::max(this->coef.size(), p.coef.size());
        Polynom<T> newPoly(maxSize);

        if(this->coef.size() >= p.coef.size())
        {
            for (int i = 0; i < this->coef.size(); ++i) {
                newPoly.coef[i] = this->coef[i];
            }
            for (int i = 0; i < p.coef.size(); ++i) {
                newPoly.coef[i] += p.coef[i];
            }
        }
        else
        {
            for (int i = 0; i < p.coef.size(); ++i) {
                newPoly.coef[i] = p.coef[i];
            }
            for (int i = 0; i < this->coef.size(); ++i) {
                newPoly.coef[i] += this->coef[i];
            }
        }

        return Polynom(newPoly);
    }

    Polynom operator-(const Polynom<T>& p)
    {
        int maxSize = std::max(this->coef.size(), p.coef.size());
        Polynom<T> newPoly(maxSize);

        if(this->coef.size() >= p.coef.size())
        {
            for (int i = 0; i < this->coef.size(); ++i) {
                newPoly.coef[i] = this->coef[i];
            }
            for (int i = 0; i < p.coef.size(); ++i) {
                newPoly.coef[i] -= p.coef[i];
            }
        }
        else
        {
            for (int i = 0; i < p.coef.size(); ++i) {
                newPoly.coef[i] = p.coef[i];
            }
            for (int i = 0; i < this->coef.size(); ++i) {
                newPoly.coef[i] -= this->coef[i];
            }
        }

        return Polynom(newPoly);
    }

    Polynom operator*(const Polynom<T>& p)
    {
        int Size = this->coef.size() + p.coef.size() - 1;
        Polynom<T> newPoly(Size);
        std::vector<double> newPolyCoef(Size, 0);
        for(int i = 0; i < coef.size(); i++)
        {
            for(int j = 0; j <= p.coef.size(); j++)
            {
                newPolyCoef[i + j] += this->coef[i] * p.coef[j];
            }
        }
        return Polynom(newPolyCoef);
    }

};




