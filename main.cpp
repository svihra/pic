#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <ctime>
#include <cstdlib>
#include <math.h>

using namespace std;

// general stats
const static double m_error = 1e-27;

const static double m_widthX = 1.5e-3;                                        // X width in m
const static double m_widthY = 1.5e-3;                                        // Y width in m
const static unsigned long m_density = 1e12;                                // electrons per meter^2
const static unsigned long m_nParticles = m_density * m_widthX * m_widthY;  // computed number of particles
const static uint m_nMacroParticles = 100000;                               // used number of particles
const static uint m_multStats = m_nParticles/m_nMacroParticles;             // multiplication factor for stats (real/used particles)
const static uint m_nIterations = 1;
const static uint m_nIterationsEmpty = 1;

// used constants
const static double m_eps_0    = 8.854e-12;
const static double m_boltzman = 1.38064852e-23;    // m2 kg s-2 K-1
const static double m_charge   = 1.602e-19;         // C
const static double m_magB[3]  = {5,-5,0};          // B
const static double m_magBSize = sqrt(pow(m_magB[0],2)+pow(m_magB[1],2)+pow(m_magB[2],2));

const static double m_elMass  = 9.1094e-31;          // kg
const static double m_elMC    = -1.758820024e11;     // C/kg
const static double m_elTemp  = 1e7;                 // K 1e7
const static double m_elVel   = sqrt((m_boltzman*m_elTemp)/(m_multStats*m_elMass));//1 // m s-1
const static double m_elOmega = m_magBSize * m_elMC  / (2.0 * M_PI);
const static double m_elLarm  = m_elVel / m_elOmega;

const static double m_ionMass  = 1.6726e-27;         // kg
const static double m_ionMC    = 9.57883322e7;       // C/kg
const static double m_ionTemp  = 1e7;                // K
const static double m_ionVel   = sqrt((m_boltzman*m_ionTemp)/(m_multStats*m_ionMass));//7e-2 // m s-1
const static double m_ionOmega = m_magBSize * m_ionMC  / (2.0 * M_PI);
const static double m_ionLarm  = m_ionVel / m_ionOmega;

// used scales
const static int m_sizeX = 300;
const static int m_sizeY = 300;
const static double dx = m_widthX/m_sizeX;
const static double dy = m_widthY/m_sizeY;

const static double dt = (dx/m_elVel) < (-0.2/m_elOmega) ? (dx/m_elVel) : (-0.2/m_elOmega);   // dt is either one spatial step or fifth of cyclotron frequency

// used geometry
bool m_geometry[m_sizeX][m_sizeY];

// defined structures
struct t_particle {
    double pos[2];
    double vel[3];
    int node[2];
};

struct t_magField {
    double tan[3];
    double sin[3];
};

// ################################################################################
// HELPFUL MATH CODE
// ################################################################################

//
// Power of sum of squares of 3D vector (3 elements in array)
double VectorPow(double a[3], double power)
{
    return pow((pow(a[0],2)+pow(a[1],2)+pow(a[2],2)),power);
}

//
// Vector cross product
void VectorCross(double a[3], double b[3], double result[3])
{
    result[0] = (a[1]*b[2]) - (a[2]*b[1]);
    result[1] = (a[2]*b[0]) - (a[0]*b[2]);
    result[2] = (a[0]*b[1]) - (a[1]*b[0]);
}

// ################################################################################
// FIELDS CODE
// ################################################################################

//
// Initialize specified 2D field to 0
void InitField(double field[m_sizeX][m_sizeY])
{
    for (int x = 0; x < m_sizeX; x++)
    {
        for (int y = 0; y < m_sizeY; y++)
        {
            field[x][y] = 0;
        }
    }

    cout << "Field initialized." << endl;
}

//
// Initialize geometry
void InitGeometry()
{
    for (int x = 0; x < m_sizeX; x++)
    {
        for (int y = 0; y < m_sizeY; y++)
        {
            if (x == 0 || y == 0 || x == m_sizeX-1 || y == m_sizeY-1)
                m_geometry[x][y] = true;
            else if (y < m_sizeY / 4 && (x < m_sizeX / 3 || x > m_sizeX - (m_sizeX / 3) ))
                m_geometry[x][y] = true;
            else
                m_geometry[x][y] = false;
        }
    }

    cout << "Geometry initialized." << endl;
}

//
// Initialize magnetic field
void InitMagField(t_magField* magFieldConst)
{
    double magBMagnitudeSquared = VectorPow((double*) m_magB, 1);

    // tan for mag field
    if (magBMagnitudeSquared != 0)
    {
        magFieldConst->tan[0] = m_magB[0] * tan(m_elMC * magBMagnitudeSquared * dt / 2.0) / magBMagnitudeSquared;
        magFieldConst->tan[1] = m_magB[1] * tan(m_elMC * magBMagnitudeSquared * dt / 2.0) / magBMagnitudeSquared;
        magFieldConst->tan[2] = 0;
    }
    else
    {
        magFieldConst->tan[0] = 0;
        magFieldConst->tan[1] = 0;
        magFieldConst->tan[2] = 0;
    }

    double magTanElMagnitudeSquared = VectorPow(magFieldConst->tan,1);

    // sin for mag field
    magFieldConst->sin[0] = 2 * magFieldConst->tan[0] / (1 + magTanElMagnitudeSquared);
    magFieldConst->sin[1] = 2 * magFieldConst->tan[1] / (1 + magTanElMagnitudeSquared);
    magFieldConst->sin[2] = 0;
}

//
// Save specified 2D field to a file
void SaveField(ofstream *file, double field[m_sizeX][m_sizeY])
{
    for (int x = 0; x < m_sizeX; x++)
    {
        for (int y = 0; y < m_sizeY; y++)
        {
            *file << x << " " << y << " " << field[x][y] << endl;
        }
    }

    cout << "Field saved to a file." << endl;
}

int SaveField(uint iter, double field[m_sizeX][m_sizeY], string name)
{
    ofstream file;
    stringstream fileName;
    fileName << "./output/fields/" << name << "_" << iter << ".txt";
    file.open(fileName.str().c_str());
    if (!file)
        return -1;
    for (int x = 0; x < m_sizeX; x++)
    {
        for (int y = 0; y < m_sizeY; y++)
        {
            file << x << " " << y << " " << field[x][y] << endl;
        }
    }
    file.close();

    cout << "Field " << name << " " << iter << " saved to a file." << endl;
    return 0;
}

int SaveField(uint iter, double field[m_sizeX][m_sizeY][3], string name)
{
    ofstream file;
    stringstream fileName;
    fileName << "./output/fields/" << name << "_" << iter << ".txt";
    file.open(fileName.str().c_str());
    if (!file)
        return -1;
    for (int x = 0; x < m_sizeX; x++)
    {
        for (int y = 0; y < m_sizeY; y++)
        {
            file << x << " " << y << " " << field[x][y][0] << " " << field[x][y][1] << " " << field[x][y][2] << endl;
        }
    }
    file.close();

    cout << "Field " << name << " " << iter << " saved to a file." << endl;
    return 0;
}

//
// Calculate electric potential using phi and rho
void CalculatePotential(double phi[m_sizeX][m_sizeY], double rho[m_sizeX][m_sizeY])
{
    cout << "Calculating potential:" << endl;
    bool loopFlag = true;
    double phi_diff, err;
    double omega = 1.56;
    uint iter = 0;

    while (loopFlag)
    {
        err = 0;
        loopFlag = false;
        for (int x = 0; x < m_sizeX; x++)
        {
            for (int y = 0; y < m_sizeY; y++)
            {
                if (m_geometry[x][y])
                {
                    phi[x][y] = 0;
                }
                else
                {
                    phi_diff = phi[x][y];
                    phi[x][y] = (1-omega)*phi[x][y]
                              + omega*(phi[x+1][y] + phi[x-1][y] + phi[x][y+1] + phi[x][y-1] + (rho[x][y]*pow(dx,2)/m_eps_0))/4.0;//pow(dx,2)/m_eps_0))/4.0;
                    err += pow(phi_diff-phi[x][y],2);
                }
            }
        }
        if (err > m_error || iter < 100)
            loopFlag = true;
        if (++iter%1000 == 0)
            cout << " - iteration: " << iter << endl;
    }
    cout << "Finished calculating potential, total iterations: " << iter << endl;
}

void CalculateField(double elField[m_sizeX][m_sizeY][3], double phi[m_sizeX][m_sizeY])
{
    cout << "Calculating field:" << endl;

    for (int x = 0; x < m_sizeX; x++)
    {
        for (int y = 0; y < m_sizeY; y++)
        {
            if (x != m_sizeX-1 && y != m_sizeY-1)
            {
                elField[x][y][0] = (phi[x][y] - phi[x+1][y])/dx;
                elField[x][y][1] = (phi[x][y] - phi[x][y+1])/dy;
            }
            else
            {
                if (x != m_sizeX-1)
                    elField[x][y][0] = (phi[x][y] - phi[x+1][y])/dx;
                else
                    elField[x][y][0] = (phi[x-1][y] - phi[x][y])/dx;

                if (y != m_sizeY-1)
                    elField[x][y][1] = (phi[x][y] - phi[x][y+1])/dy;
                else
                    elField[x][y][1] = (phi[x][y-1] - phi[x][y])/dy;
            }
        }
    }

    cout << "Finished calculating field. "<< elField[m_sizeX/2][m_sizeY/2][0] << endl;
}

// ################################################################################
// PARTICLES CODE
// ################################################################################

//
// Calculate charge to grid
void ChargeToGrid(deque<t_particle > *particles, double rho_p, double rho[m_sizeX][m_sizeY])
{
    uint n_particles = particles->size();
    int *x, *y;
    for (uint iter = 0; iter < n_particles; iter++)
    {
        x = &particles->at(iter).node[0];
        y = &particles->at(iter).node[1];

        rho[(*x)][(*y)] += ( rho_p
                        * ((dx * (1 + (*x))) - particles->at(iter).pos[0]) // Delta x - dx
                        * ((dy * (1 + (*y))) - particles->at(iter).pos[1]) // Delta y - dy
                        ) / (dx * dy);

//        cout << "particle: " << iter << ", at position: " << *x << " " << *y << endl;
    }

    cout << "Charge calculated to the grid." << endl;
}

//
// Save particles to a file
void SaveParticles(ofstream *file, deque<t_particle > *particles, int charge)
{
    uint n_particles = particles->size();
    for (uint iter = 0; iter < n_particles; iter++)
    {
        *file << particles->at(iter).pos[0] << " " << particles->at(iter).pos[1] << " " << charge << " " << particles->at(iter).vel[0] << " " << particles->at(iter).vel[1] << " " << particles->at(iter).vel[2] << endl;
    }

    cout << "Particles position and velocities saved to a file." << endl;
}

//
// Save particles to a file according to iteration number
int SaveParticles(uint iter, deque<t_particle > *particles, int charge, string name)
{
    ofstream file;
    stringstream fileName;
    fileName << "./output/particles/" << name << "_" << iter << ".txt";
    file.open(fileName.str().c_str());
    if (!file)
        return -1;

    double field[m_sizeX][m_sizeY];
    InitField(field);

    uint n_particles = particles->size();
    for (uint iter = 0; iter < n_particles; iter++)
    {
        field[particles->at(iter).node[0]][particles->at(iter).node[1]] += charge;
        file << particles->at(iter).pos[0] << " " << particles->at(iter).pos[1] << " " << charge << " " << particles->at(iter).vel[0] << " " << particles->at(iter).vel[1] << " " << particles->at(iter).vel[2] << endl;
    }
    file.close();

    SaveField(iter,field,name);

    cout << "Particles position and velocities run " << iter << " saved to a file." << endl;
    return 0;
}

//
// Generate particles
void GenerateParticles(deque<t_particle > *particles, double vel)
{
    t_particle particle;

    double U1, U2, U3, U4;
    for (uint iter = 0; iter < m_nMacroParticles; iter++)
    {
        // position
        particle.pos[0] = m_widthX*((double) rand())/(RAND_MAX);
        particle.pos[1] = m_widthY*((double) rand())/(RAND_MAX);
        particle.node[0] = (int) (particle.pos[0]/dx);
        particle.node[1] = (int) (particle.pos[1]/dy);

        // velocity
        U1 = ((double) rand())/RAND_MAX;
        U2 = ((double) rand())/RAND_MAX;
        U3 = ((double) rand())/RAND_MAX;
        U4 = ((double) rand())/RAND_MAX;
        particle.vel[0] = ((vel/sqrt(3))*sqrt(-2 * log(U1)) * cos(2 * M_PI * U2));
        particle.vel[1] = ((vel/sqrt(3))*sqrt(-2 * log(U1)) * sin(2 * M_PI * U2));
        particle.vel[2] = ((vel/sqrt(3))*sqrt(-2 * log(U3)) * cos(2 * M_PI * U4));

        particles->push_back(particle);
    }
}

//
// Collisions of particles
void EvaluateCollisions(deque<t_particle > *particles, deque<t_particle > *particlesPassed)
{
    bool passed = false;
    uint nPassed = 0;

    int *x, *y;
    uint n_particles = particles->size();
    for (uint iter = 0; iter < n_particles; iter++)
    {
        x = &particles->at(iter).node[0];
        y = &particles->at(iter).node[1];

        passed = false;
        particles->at(iter).pos[0] += particles->at(iter).vel[0] * dt;
        particles->at(iter).pos[1] += particles->at(iter).vel[1] * dt;

        while (particles->at(iter).pos[0] > m_widthX)
        {
            particles->at(iter).pos[0] -= m_widthX;
            passed = true;
        }

        while (particles->at(iter).pos[0] < 0)
        {
            particles->at(iter).pos[0] += m_widthX;
            passed = true;
        }

        while (particles->at(iter).pos[1] > m_widthY)
        {
            particles->at(iter).pos[1] -= m_widthY;
            passed = false;
        }

        while (particles->at(iter).pos[1] < 0)
        {
            particles->at(iter).pos[1] += m_widthY;
            passed = true;
        }

        *x = (int) (particles->at(iter).pos[0] / dx);
        *y = (int) (particles->at(iter).pos[1] / dy);

        if (*x > m_sizeX-1 || *x < 0 || *y > m_sizeY-1 || *y < 0 )
        {
            particles->erase(particles->begin() + iter);
            continue;
        }

        // BUG
        // particles can enter to geometry, even though the charge will be 0
        // and they will cease to exist after one step, however, if the velocity
        // is large enough, they may emerge from corner of geometry
        if (passed)
        {
            nPassed++;
            t_particle particle = particles->at(iter);
            particlesPassed->push_front(particle);
        }
    }

    cout << "Collisions evaluated, " << nPassed << " particles left the volume." << endl;
}

//
// Evaluate forces and update particles position/velocities
// magnetic field done by birsdale method, used global variable magTanEl and magSinEl
void EvaluateForces(deque<t_particle > *particles, double particleMC, double elField[m_sizeX][m_sizeY][3], t_magField *magConst)
{
    long n_particles = particles->size();
    double vMin[3], vPrime[3], vMax[3];

    int *x, *y;
    for (long iter = n_particles-1; iter >= 0; iter--)
    {
        x = &particles->at(iter).node[0];
        y = &particles->at(iter).node[1];

        particles->at(iter).pos[0] += particles->at(iter).vel[0] * dt;
        particles->at(iter).pos[1] += particles->at(iter).vel[1] * dt;

        *x = (int) (particles->at(iter).pos[0] / dx);
        *y = (int) (particles->at(iter).pos[1] / dy);

        //
        // erase particle if out of boundaries or in defined geometry
        if (*x > m_sizeX-1 || *x < 0 || *y > m_sizeY-1 || *y < 0 )
        {
            particles->erase(particles->begin() + iter);
            continue;
        }
        else if (*x != m_sizeX-1 && *x != 0 && *y != m_sizeY-1 && *y != 0 && m_geometry[*x][*y]) //
        {
            particles->erase(particles->begin() + iter);
            continue;
        }

        for (int i = 0; i < 3; i++)
        {
            vMin[i] = particles->at(iter).vel[i] + (0.5*dt*particleMC*elField[*x][*y][i]);
            vMax[i] = 0;
        }

        VectorCross(vMin,magConst->tan,vPrime);

        for (int i = 0; i < 3; i++)
        {
            vPrime[i] += vMin[i];
        }

        VectorCross(vPrime,magConst->sin,vMax);

        for (int i = 0; i < 3; i++)
        {
            vMax[i] += vMin[i];
            particles->at(iter).vel[i] = vMax[i] + (0.5*dt*particleMC*elField[*x][*y][i]);
        }
    }
    cout << "Forces evaluated. " << 0.5*dt*particleMC*elField[m_sizeX/2][m_sizeX/2][0] << endl;
}

// ################################################################################
// INFO CODE
// ################################################################################

bool PrintInfo()
{
    cout << "##############################################" << endl;
    cout << "######## Starting the PIC simulation #########" << endl;
    cout << "##############################################" << endl;
    cout << endl;
    cout << "##############################################" << endl;
    cout << "################ Constants ###################" << endl;
    cout << "##############################################" << endl;
    cout << "       eps0  = " << m_eps_0    << endl;
    cout << "       kb    = " << m_boltzman << endl;
    cout << "       c     =" << m_charge   << endl;
    cout << "       Bx    = " << m_magB[0] << endl;
    cout << "       By    = " << m_magB[1] << endl;
    cout << "       Bz    = " << m_magB[2] << endl;
    cout << "       Bsize = " << m_magBSize << endl;
    cout << endl;
    cout << "##############################################" << endl;
    cout << "################ Electrons ###################" << endl;
    cout << "##############################################" << endl;
    cout << "       m     = " << m_elMass   << endl;
    cout << "       c/m   = " << m_elMC     << endl;
    cout << "       T     = " << m_elTemp   << endl;
    cout << "       v     = " << m_elVel    << endl;
    cout << "       omega = " << m_elOmega  << endl;
    cout << "       rl    = " << m_elLarm   << endl;
    cout << endl;
    cout << "##############################################" << endl;
    cout << "################# Ions #######################" << endl;
    cout << "##############################################" << endl;
    cout << "        m     = " << m_ionMass  << endl;
    cout << "        c/m   = " << m_ionMC    << endl;
    cout << "        T     = " << m_ionTemp  << endl;
    cout << "        v     = " << m_ionVel   << endl;
    cout << "        omega = " << m_ionOmega << endl;
    cout << "        rl    = " << m_ionLarm  << endl;
    cout << endl;
    cout << "##############################################" << endl;
    cout << "################ Parameters ##################" << endl;
    cout << "##############################################" << endl;
    cout << "        error = " << m_error            << endl;
    cout << "         x    = " << m_widthX           << endl;
    cout << "        Nx    = " << m_sizeX            << endl;
    cout << "        dx    = " << dx                 << endl;
    cout << "         y    = " << m_widthY           << endl;
    cout << "        Ny    = " << m_sizeY            << endl;
    cout << "        dy    = " << dy                 << endl;
    cout << "        dt    = " << dt                 << endl;
    cout << "        n0    = " << m_density          << endl;
    cout << "        N     = " << m_nParticles       << endl;
    cout << "        NReal = " << m_nMacroParticles  << endl;
    cout << "        Ndiff = " << m_multStats        << endl;
    cout << "        STEP  = " << m_nIterations      << endl;
    cout << "        step  = " << m_nIterationsEmpty << endl;
    cout << "##############################################" << endl;
    cout << endl;
    cout << "Iterations needed: " << endl;
    cout << "    electron to traverse: " << m_widthX  / (m_elVel * dt) << endl;
    cout << "    electron to rotate:   " << 1 / (m_elOmega * dt) << endl;
    cout << "    ion to traverse: " << m_widthX  / (m_ionVel * dt) << endl;
    cout << "    ion to rotate:   " << 1 / (m_ionOmega * dt) << endl;
    cout << endl;

    char proceed;
    while (proceed != 'y' && proceed != 'n')
    {
        cout << "Are the input parameters correct? (y/n): ";
        cin >> proceed;
    }
    if (proceed == 'y')
        return true;
    else
        return false;
}

// ################################################################################
// MAIN CODE
// ################################################################################

int main()
{
    //
    // Print constants and info
    if (!PrintInfo())
        return -1;

    //
    // Randomize input
    srand(time(NULL));

    //
    // Define all needed fields
    double rho[m_sizeX][m_sizeY];
    double phi[m_sizeX][m_sizeY];
    double elField[m_sizeX][m_sizeY][3];
    t_magField magElConst, magIonConst;

    //
    // Initialize fields
    InitField(rho);
    InitField(phi);

    //
    // Initialize magnetic field constants for forces evaluation
    InitMagField(&magElConst);
    InitMagField(&magIonConst);

    //
    // Initialize geometry
    InitGeometry();

    //
    // Generate particle sources and monitored particles
    deque<t_particle > electronSource;
    GenerateParticles(&electronSource, m_elVel);
    deque<t_particle > electrons;

    deque<t_particle > ionSource;
    GenerateParticles(&ionSource, m_ionVel);
    deque<t_particle > ions;

    double electronCharge = (-1) * m_charge * m_multStats;
    double ionCharge = m_charge * m_multStats;

    //
    // Main part of the code
    for (uint ionMove = 0; ionMove < m_nIterations; ionMove++)
    {
        cout << "-----electrons---------------------------------" << endl;
        cout << "nx: " << electronSource.begin()->node[0] << ", ny: " << electronSource.begin()->node[1] << endl;
        cout << " x: " << electronSource.begin()->pos[0] << ", y: "   << electronSource.begin()->pos[1] << endl;
        cout << "vx: " << electronSource.begin()->vel[0] << ", vy: "  << electronSource.begin()->vel[1] << endl;
        cout << "number of particles: " << electrons.size() << endl;
        cout << "-----ions--------------------------------------" << endl;
        cout << "nx: " << ionSource.begin()->node[0] << ", ny: "<< ionSource.begin()->node[1] << endl;
        cout << " x: " << ionSource.begin()->pos[0] << ", y: "  << ionSource.begin()->pos[1] << endl;
        cout << "vx: " << ionSource.begin()->vel[0] << ", vy: " << ionSource.begin()->vel[1] << endl;
        cout << "number of particles: " << ions.size() << endl;

        for (uint elMove = 0; elMove < m_nIterationsEmpty; elMove++)
        {
            EvaluateCollisions(&ionSource,&ions);
            EvaluateCollisions(&electronSource,&electrons);
            InitField(rho);
            ChargeToGrid(&electrons,electronCharge,rho);
            ChargeToGrid(&ions,ionCharge,rho);
            InitField(phi);
            CalculatePotential(phi,rho);
            CalculateField(elField,phi);
            EvaluateForces(&ions, m_ionMC, elField, &magIonConst);
            EvaluateForces(&electrons, m_elMC, elField, &magElConst);
        }

        SaveField(ionMove,phi,"phi");
        SaveField(ionMove,rho,"rho");
        SaveField(ionMove,elField,"E");
        SaveParticles(ionMove,&electrons,-1,"el");
        SaveParticles(ionMove,&ions,1,"ion");
        cout << "##############################################" << endl;
        cout << "##### " << (100*(ionMove+1))/m_nIterations <<"% completed ###########################" << endl;
        cout << "##############################################" << endl;
    }
    // -------------------------

    //
    // Erase everything not needed
    electrons.clear();
    // -------------------------

    return 0;
}

