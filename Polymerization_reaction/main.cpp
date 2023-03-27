#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <vector>
#include <random>
#include <iostream>

#define NSTEPS 100000
#define eps  5
#define sigma 1
#define m 1.0
#define L  10
#define r_c  4.0
#define p  0.001
#define lambda  10
#define KB  0.00831
#define tau 0.001
#define T  300

struct Atom
{
    float x, y, z;
    float vx, vy, vz;
    float fx, fy, fz;
    char name;// A B C
    bool real;
};

float transferPBC(float x)
{
    if (x < 0)
    {
        return x + L;
    }
    else if (x > L)
    {
        return x - L;
    }
    return x;
}

void saveFrame(const char* filename, const char* modifier, std::vector<Atom> atoms)
{
    FILE* out = fopen(filename, modifier);
    fprintf(out, "%ld\nA+B->A+C\n", atoms.size());
    for (int i = 0; i < atoms.size(); i++)
    {
        fprintf(out, "%c\t%f\t%f\t%f\n",
                atoms[i].name,
                atoms[i].x*10.0,
                atoms[i].y*10.0,
                atoms[i].z*10.0);
    }
    fclose(out);
}

float distance(Atom atom1, Atom atom2) {

    float dist = (atom1.x - atom2.x) * (atom1.x - atom2.x) +
                 (atom1.y - atom2.y) * (atom1.y - atom2.y) +
                 (atom1.z - atom2.z) * (atom1.z - atom2.z);
    return sqrtf(dist);
}

int main(int argc, char* argv[])
{
//    float currentTemperature = T;
    int N_b=100, N_c = 100, N_a = 10;
    int N = N_b + N_a + N_c;
    int N_first = N_b + N_a ;
    std::vector<Atom> atoms(N);
    std::random_device randomDevice;
    std::mt19937 randomGenerator(randomDevice());
    std::uniform_real_distribution<> distributionX(0, L);
    std::uniform_real_distribution<> distributionKSI(0, 1);
    std::normal_distribution<> distributionV(0.0, sqrtf(KB*T));
    std::normal_distribution<> distributionR(0.0, 1.0);

    Atom null_atom;

    null_atom.name = 'B';
    null_atom.real = false;

    null_atom.x = 0.0;
    null_atom.y = 0.0;
    null_atom.z = 0.0;

    null_atom.vx = 0.0;
    null_atom.vy = 0.0;
    null_atom.vz = 0.0;

    null_atom.fx = 0.0;
    null_atom.fy = 0.0;
    null_atom.fz = 0.0;


    for (int i = 0; i < N; i++)
    {


        if (i <= 10)
        {
            atoms[i].name = 'A';
            atoms[i].real = true;
        }
        if (i > 10 and i < 110)
        {
            atoms[i].name = 'B';
            atoms[i].real = true;
        }
        if (i >= 110 and i < 211)
        {
            atoms[i].name = 'C';
            atoms[i].real = false;
        }

        if(i<110){
            atoms[i].x = distributionX(randomGenerator);
            atoms[i].y = distributionX(randomGenerator);
            atoms[i].z = distributionX(randomGenerator);

            atoms[i].vx = distributionV(randomGenerator);
            atoms[i].vy = distributionV(randomGenerator);
            atoms[i].vz = distributionV(randomGenerator);

            atoms[i].fx = 0.0;
            atoms[i].fy = 0.0;
            atoms[i].fz = 0.0;
        }
        else
        {
            atoms[i].x = 0.0;
            atoms[i].y = 0.0;
            atoms[i].z = 0.0;

            atoms[i].vx = 0.0;
            atoms[i].vy = 0.0;
            atoms[i].vz = 0.0;

            atoms[i].fx = 0.0;
            atoms[i].fy = 0.0;
            atoms[i].fz = 0.0;
        }

    }

    saveFrame("ABC.xyz", "w", atoms);
//    double temperature = 0.0;
//    int nTemperature = 0;

    for (int n = 0; n < NSTEPS; n++)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if(atoms[i].real and atoms[j].real) {
                    if ((atoms[i].name == 'B' and atoms[j].name == 'A') or
                        (atoms[i].name == 'A' and atoms[j].name == 'B')) {

                        float dx = atoms[i].x - atoms[j].x;
                        float dy = atoms[i].y - atoms[j].y;
                        float dz = atoms[i].z - atoms[j].z;

                        dx -= rint(dx / L) * L;
                        dy -= rint(dy / L) * L;
                        dz -= rint(dz / L) * L;

                        float dr2 = dx * dx + dy * dy + dz * dz;
                        float dr_n = sqrtf(dr2);

                        if (dr_n < r_c) {
                            float ksi = distributionKSI(randomGenerator);
                            if (ksi < p) {
                                if (atoms[i].name == 'B') {

                                    atoms[i + N_b] = atoms[i];

                                    atoms[i + N_b].name = 'C';

                                    atoms[i] = null_atom;

                                    std::cout << "i" << atoms[i].x << " " << atoms[i].y << " " <<atoms[i].z << std::endl;
                                    std::cout << "i+N_b" << atoms[i+N_b].x << " " << atoms[i+N_b].y << " " <<atoms[i+N_b].z << std::endl;



                                }
                                if (atoms[j].name == 'B') {

                                    atoms[i + N_b] = atoms[i];

                                    atoms[i + N_b].name = 'C';

                                    atoms[i] = null_atom;

                                }

                            }
                        }

                    }
                }
            }
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i; j++)
            {
                if(atoms[i].real and atoms[j].real)
                {

                    float C1, C2;

                    if (atoms[i].name == 'C' and atoms[j].name == 'C') {
                        C1 = 1;
                        C2 = -2;
                    } else {
                        C1 = 0;
                        C2 = 1;
                    }

                    float dx = atoms[i].x - atoms[j].x;
                    float dy = atoms[i].y - atoms[j].y;
                    float dz = atoms[i].z - atoms[j].z;

                    dx -= rint(dx / L) * L;
                    dy -= rint(dy / L) * L;
                    dz -= rint(dz / L) * L;

                    float dr2 = dx * dx + dy * dy + dz * dz;
                    float dr_n = sqrtf(dr2);

                    float df = eps * (4 * C1 * pow(sigma, 4) * pow(dr_n, -6) + 2 * C2 * pow(sigma, 2) * pow(dr_n, -4));

                    atoms[i].fx += df * dx;
                    atoms[i].fy += df * dy;
                    atoms[i].fz += df * dz;

                    atoms[j].fx -= df * dx;
                    atoms[j].fy -= df * dy;
                    atoms[j].fz -= df * dz;

                }
            }
        }




        for (int i = 0; i < N; i++){
            if(atoms[i].real){


            atoms[i].vx = (1/((m/tau)+(lambda*m/2)))*(atoms[i].fx + atoms[i].vx * ((m/tau)-(lambda*m/2)) + sqrtf(2*KB*T*lambda*m/tau) * distributionR(randomGenerator));
            atoms[i].vy = (1/((m/tau)+(lambda*m/2)))*(atoms[i].fy + atoms[i].vy * ((m/tau)-(lambda*m/2)) + sqrtf(2*KB*T*lambda*m/tau) * distributionR(randomGenerator));
            atoms[i].vz = (1/((m/tau)+(lambda*m/2)))*(atoms[i].fz + atoms[i].vz * ((m/tau)-(lambda*m/2)) + sqrtf(2*KB*T*lambda*m/tau) * distributionR(randomGenerator));

            atoms[i].x = atoms[i].x + tau * atoms[i].vx;
            atoms[i].y = atoms[i].y + tau * atoms[i].vy;
            atoms[i].z = atoms[i].z + tau * atoms[i].vz;

            atoms[i].x = transferPBC(atoms[i].x);
            atoms[i].y = transferPBC(atoms[i].y);
            atoms[i].z = transferPBC(atoms[i].z);

            atoms[i].fx = 0.0;
            atoms[i].fy = 0.0;
            atoms[i].fz = 0.0;

//            float v2 = atoms[i].vx*atoms[i].vx +
//                       atoms[i].vy*atoms[i].vy +
//                       atoms[i].vz*atoms[i].vz;
//            temperature += m * v2;
//            nTemperature ++;

        }}

        if (n % 100 == 0)
        {
//            currentTemperature = (temperature/(3.0*KB))/nTemperature;
//            printf("%d\t%f\n", n, currentTemperature);
//            temperature = 0.0;
//            nTemperature = 0;
            saveFrame("ABC.xyz", "a", atoms);
        }
    }
}