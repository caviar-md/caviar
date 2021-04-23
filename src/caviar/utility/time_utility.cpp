
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// This file is part of the CAVIAR package.
//
// The CAVIAR package is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the CAVIAR distribution.
//
//========================================================================

//  Windows
#ifdef _WIN32
#include <Windows.h>

namespace caviar {
  
double get_wall_time(){
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)){
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)){
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}

double get_cpu_time(){
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
            ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}

}

//  Posix/Linux
#else

#include <time.h>
#include <sys/time.h>

namespace caviar {

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

}
#endif

/*
void example_for_openmp (){
    {
      //  Start Timers
      double wall0 = get_wall_time();
      double cpu0  = get_cpu_time();

      //  Perform some computation.
      double sum = 0;
  //#pragma omp parallel for reduction(+ : sum)
      for (long long i = 1; i < 100000; i++){
          sum += log((double)i);
      }

      //  Stop timers
      double wall1 = get_wall_time();
      double cpu1  = get_cpu_time();

      cout << "Wall Time = " << wall1 - wall0 << endl;
      cout << "CPU Time  = " << cpu1  - cpu0  << endl;

      //  Prevent Code Elimination
      cout << endl;
      cout << "Sum = " << sum << endl;
    } 
    //----------------------------------------------
    {
      //  Start Timers
      double wall0 = get_wall_time();
      double cpu0  = get_cpu_time();

      //  Perform some computation.
      double sum = 0;
    //#pragma omp parallel for reduction(+ : sum)
      for (long long i = 1; i < 100000; i++){
          sum += log((double)i);
      }

      //  Stop timers
      double wall1 = get_wall_time();
      double cpu1  = get_cpu_time();

      cout << "Wall Time = " << wall1 - wall0 << endl;
      cout << "CPU Time  = " << cpu1  - cpu0  << endl;

      //  Prevent Code Elimination
      cout << endl;
      cout << "Sum = " << sum << endl;
    }

}
*/
