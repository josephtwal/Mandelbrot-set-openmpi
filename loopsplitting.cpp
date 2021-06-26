#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <time.h>
#include <mpi.h>

using namespace std;
const int imgX = 3000;      //  horizontal image resolution
const int imgY = 3000;      //  vertical image resolution
const int iter_n = 3000;    // max. iteration number
const double yMin = -0.135; // Mandelbrot scene y - range
const double yMax = -0.115;
const double xMin = -0.79; // Mandelbrot scene x -range
const double xMax = -0.77;
int img_array[imgX][imgY] = {0}; // our MAndelbrot set values array

// convergation function - base of the Mandelbrot set value generation
// it will get two parameters (x and y coordinates) and will give an iteration count in return
int converges(double cx, double cy)
{
    int n = 0;
    double zx = 0;
    double new_zx = 0;
    double zy = 0;
    // we iterate until max. iteration count iter_n, or until z^2 (complex!) runs over 4 - this means, our series will run to infinity, so it's not part of the set
    while ((n < iter_n) && (zx * zx + zy * zy) < 4)
    {
        // z * z => new_zx = (zx*zx - zy*zy)  new_zy = (zx*zy + zx*zy)   // we work with complex numbers
        // z*z + c = zx^2 - zy^2 +cx   +  i(zx*zy*2 + cy)
        new_zx = zx * zx - zy * zy + cx;
        zy = 2 * zx * zy + cy;
        zx = new_zx;
        n++;
    }
    return n;
}

int main(int argc, char **argv)
{
    //variables for MPI communication:
    int id, nproc;
    MPI_Status status;
    int id_from;

    double resX = 0;  // Resolution of our iteration steps
    double resY = 0;  // this will be calculated by (yMax-yMin) / imgY later..
    double cx = xMin; // calculation will start at this point and we will change this dynamically
    double cy = yMin;
    double s_time, e_time;
    ofstream myfile;
    char *fileName1; // we will show some timing data

    // Initialize MPI:
    MPI_Init(&argc, &argv);
    // Get my rank:
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    // Get the total number of processors:
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    MPI_Barrier(MPI_COMM_WORLD); //for precize timing

    if (id == 0)
    { //Master

        if (argc < 2)
        {
            cout << "Usage:" << endl;
            cout << argv[0] << " out_file.ppm" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            exit(1);
        }

        
        string filename1(argv[1]);
        // filename1 = "mandel.ppm";
        fileName1 = (char *)filename1.c_str();

        s_time = MPI_Wtime(); // we get a time value at this point

        //recieve
        for (int j = 1; j < nproc; ++j)
        {
            MPI_Recv(&id_from, 1, MPI_INT, MPI_ANY_SOURCE, 1,
                     MPI_COMM_WORLD, &status);
            //it is important to recive from ANY_SOURCE,
            //as process execution may differ greatly

            for (int i = id_from - 1; i < imgX; i += (nproc - 1))
                MPI_Recv(&img_array[i][0], imgY, MPI_INT, id_from,
                         2, MPI_COMM_WORLD, &status);
        }

        // we get another time at this point, so we can calculate the elapsed time for the calculation
        e_time = MPI_Wtime();

        cout << "Time elapsed during calculation: " << e_time - s_time << " secs." << endl;
        ;

        // file IO
        myfile.open(fileName1);
        myfile << "P3\r\n";
        myfile << imgX;
        myfile << " ";
        myfile << imgY;
        myfile << "\r\n";
        myfile << "255\r\n";

        // we have to colour our dataset. Actually, the members of the Mandelbrot set are used to be the same colour (black?) and have from point of visualisations view no interest.
        // the outer points are represented by assigning colours to their iteration steps and this generates vivid forms and colors
        for (int i = 0; i < imgX; i++)
        {
            for (int j = 0; j < imgY; j++)
            {

                if ((img_array[i][j] < 256)) // we go from black to red in this range
                {
                    myfile << img_array[i][j] << " 0 0";// (int)(84*pow(img_array[i][j],0.2)) << " 0 0"; //myfile << img_array[i][j] << " 0 0";
                }
                else if (img_array[i][j] < 512) // we go from red to yellow in this range
                {
                    myfile << "255 " << img_array[i][j] - 256 << " 0";
                }
                else if (img_array[i][j] < 768) // we go from yellow to white in this range
                {
                    myfile << "255 255 " << img_array[i][j] - 512;
                }
                 // we could refine our palette for more resolution, more iteration-step images
 		      else if( img_array[i][j] < 1024)
 		      {
 		      myfile << 1024-img_array[i][j] << " 255 255"; 
 		      }
 		      else if( img_array[i][j] < 1280)
 		      {
 		      myfile << "0 "  << 1280-img_array[i][j] << " 255"; 
 		      }
 		      else if( img_array[i][j] < 1536)
 		      {
 		      myfile << "0 0 "  << 1536-img_array[i][j]; 
 		      }
 		  
                else // everything else is black
                {
                    myfile << "0 0 0 ";
                }

                myfile << " ";
            }
            myfile << "\r\n";
        }

        myfile.close(); // we close our file

        e_time = MPI_Wtime(); // and give another elapsed time info (IO included)

        cout << "Time elapsed total: " << e_time - s_time << " secs \r\n";
    }
    else
    { //Slave
        //because the slaves numbered 1..nproc:
        //id -> id-1
        //nproc -> nproc-1

        //prepare the step resolution
        resX = (xMax - xMin) / imgX;
        resY = (yMax - yMin) / imgY;

        //because of the Master-Slave execution process 0 does no work!

        cx = cx + resX * (id - 1); //jump to the right cx
        //beware, round-off error!
        //this value is slightly different from
        //adding resX startval times to cx
        //this will cause unspottable difference to the picture

        for (int i = id - 1; i < imgX; i += (nproc - 1))
        {
            cy = yMin; // at every new step in X direction,
            //we start at the first Y value
            for (int j = 0; j < imgY; j++)
            {
                img_array[i][j] = converges(cx, cy);
                cy = cy + resY;
            }

            cx = cx + resX * (nproc - 1); //the steps are by nproc!
                                          //beware, round-off error!
        }

        MPI_Send(&id, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        for (int i = id - 1; i < imgX; i += (nproc - 1))
        {
            MPI_Send(&img_array[i][0], imgY, MPI_INT, 0, 2, MPI_COMM_WORLD);
        }
    }

    // Terminate MPI:
    MPI_Finalize();

    return 0;
}