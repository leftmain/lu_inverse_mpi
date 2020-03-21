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
#include <sys/resource.h>
#include <limits>
#include <algorithm>

double
get_earth_time ();

double
get_process_time ();


struct Time
{
  double global_time = 0.;
  double process_time = 0.;
  double decomposition_time = 0.;
  double invert_time = 0.;
  double mult_time = 0.;

  double residual_global_time = 0.;
  double residual_process_time = 0.;
};


class ProcessTimer
{
  private:
    double *time = nullptr;
    bool stopped = false;

  public:
    ProcessTimer (double &t)
      {
        time = &t;
        *time = 0. - get_process_time ();
      }
    ~ProcessTimer ()
      {
        if (!stopped)
          *time += get_process_time ();
      }
    void start ()
      {
        *time = 0. - get_process_time ();
      }
    void stop ()
      {
        *time += get_process_time ();
        stopped = true;
      }
};


class GlobalTimer
{
  private:
    double *time = nullptr;
    bool stopped = false;

  public:
    GlobalTimer (double &t)
      {
        time = &t;
        *time = 0. - MPI_Wtime ();
      }
    ~GlobalTimer ()
      {
        if (!stopped)
          {
            *time += MPI_Wtime ();
          }
      }
    void start ()
      {
        *time = 0. - MPI_Wtime ();
      }
    void stop ()
      {
        *time += MPI_Wtime ();
        stopped = true;
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
