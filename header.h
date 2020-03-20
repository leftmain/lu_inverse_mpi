#ifndef HEADER_H
#define HEADER_H

//#define DEBUG
//#define RES_TEST
//#define CHECK

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <sys/time.h>
#include <limits>
#include <algorithm>

double
get_earth_time ();

class MPI_Auto_Timer
{
  private:
    double *time = nullptr;

  public:
    MPI_Auto_Timer (double *t)
      {
        time = t;
        *time = 0. - MPI_Wtime ();
      }
    ~MPI_Auto_Timer ()
      {
        *time += MPI_Wtime ();
      }
};

template <typename T>
class Deleter
{
  private:
    std::vector<T> allocated_memory;
    
  public:
    ~Deleter ()
      {
        for (T i : allocated_memory)
          {
            if (i)
              delete [] i;
          }
      }
    void add (T new_mem)
      {
        allocated_memory.push_back (new_mem);
      }
};

class Closer
{
  private:
    std::vector<FILE *> fp;

  public:
    ~Closer ()
      {
        for (FILE *i : fp)
          {
            if (i)
              fclose (i);
          }
      }
    void add (FILE *new_fp)
      {
        fp.push_back (new_fp);
      }
};

enum Error
{
  ALL_RIGHT = 0,
  CANNOT_OPEN = 1,
  CANNOT_READ = 2,
  MEMORY_ERROR = 3,
  METHOD_ERROR = 4,
};

#endif
