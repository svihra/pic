#include <iostream>
#include <stdio.h>
#include <deque>
#include <math.h>

using namespace std;

const static uint m_sizeX = 100;
const static uint m_sizeY = 100;
const static float m_step = 1e-3;
const static float m_eps_0 = 8.854e-12;
const static float error = 1e-10;

struct t_particle {
    float pos[2];
    float vel[3];
    uint node[2];
};



// ################################################################################
// FIELDS CODE
// ################################################################################

//
// Initialize specified 2D field to 0
void initField(float field[m_sizeX][m_sizeY])
{
    for (uint x = 0; x < m_sizeX; x++)
    {
        for (uint y = 0; y < m_sizeY; y++)
        {
            field[x][y] = 0;
        }
    }

    cout << "Field initialized" << endl;
}

//
// Save specified 2D field to a file
void saveField(FILE* file, float field[m_sizeX][m_sizeY])
{
    for (uint x = 0; x < m_sizeX; x++)
    {
        for (uint y = 0; y < m_sizeY; y++)
        {
            fprintf(file, "%u %u %f \n",x,y,field[x][y]);
        }
    }
}

//
// Calculate phi using phi and rho
void calculateFields(float phi[m_sizeX][m_sizeY], float rho[m_sizeX][m_sizeY])
{
    bool loopFlag;
    float phi_diff;
    float omega = 1.56;
    uint iter = 0;

    while (loopFlag)
    {
        loopFlag = false;
        for (uint x = 0; x < m_sizeX; x++)
        {
            for (uint y = 0; y < m_sizeY; y++)
            {
                if (x == 0 || y == 0 || x == m_sizeX-1 || y == m_sizeY-1)
                {
                    phi[x][y] = 0;
                }
                else
                {
                    phi_diff = phi[x][y];
                    phi[x][y] = (1-omega)*phi[x][y] + omega*(phi[x+1][y] + phi[x-1][y] + phi[x][y+1] + phi[x][y-1] + (rho[x][y]*m_step*m_step/m_eps_0))/4.0;
                    if (phi_diff - phi[x][y] > error || phi_diff - phi[x][y] < -error)
                    {
                        loopFlag = true;
                    }
                }
            }
        }
        if (++iter%1000 == 0)
            cout << "iteration: " << iter << endl;
    }
}

// ################################################################################
// PARTICLES CODE
// ################################################################################

//
// Calculate charge to grid
void chargeToGrid(deque<t_particle > *particles, float rho_p, float rho[m_sizeX][m_sizeY])
{
    initField(rho);
    uint n_particles = particles->size();
    uint *x, *y;
    for (uint iter = 0; iter < n_particles; iter++)
    {
        x = &particles->at(iter).node[0];
        y = &particles->at(iter).node[1];

        rho[(*x)][(*y)] += ( rho_p
                        * ((m_step * (1 + (*x))) - particles->at(iter).pos[0]) // Delta x - dx
                        * ((m_step * (1 + (*y))) - particles->at(iter).pos[1]) // Delta y - dy
                        ) / (m_step * m_step);

        cout << "particle: " << iter << ", at position: " << *x << " " << *y << endl;
    }
}

// ################################################################################
// MAIN CODE
// ################################################################################
int main()
{
    FILE* file_phi = fopen("./phi.txt","w");
    FILE* file_rho = fopen("./rho.txt","w");
    float rho[m_sizeX][m_sizeY];
    float phi[m_sizeX][m_sizeY];

    initField(rho);
    initField(phi);

    deque<t_particle > electrons;
    t_particle electron;
    float electron_charge = 1.602e-10;

    for (uint y = 0; y < m_sizeX; y++)
    {
        uint x = 50;
//        for (uint y = 0; y < m_sizeY; y++)
//        {
//            if (pow((x - 50),2) + pow((y - 50),2) < 100)
            {
                electron.node[0] = x;
                electron.node[1] = y;
                electron.pos[0] = x*m_step;
                electron.pos[1] = y*m_step;
                electrons.push_back(electron);
            }
//        }
    }

    chargeToGrid(&electrons,electron_charge,rho);
    calculateFields(phi,rho);

    saveField(file_phi,phi);
    saveField(file_rho,rho);

    return 0;
}

