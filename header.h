#ifndef HEADER_H
#define HEADER_H

//#define DEBUG
#define TEST

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>

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
};

#endif
