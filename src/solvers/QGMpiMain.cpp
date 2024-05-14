#include <mpi.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>

#include "MinotaurConfig.h"
#include "Environment.h"
#include "Option.h"
#include "Problem.h"
#include "Types.h"
#include "QGMpi.h"

using namespace Minotaur;

std::string share_problem_file(std::string);

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int mpirank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

  EnvPtr env = (EnvPtr) new Environment();
  int err = 0;
  QGMpi qg(env);
  ProblemPtr p = 0;
  std::string fname, dname;
 
  qg.doSetup();
  env->startTimer(err);
  if (err) {
    goto CLEANUP;
  }

  // Parse command line for options set by the user.
  env->readOptions(argc, argv);
  
  if (0!=qg.showInfo()) {
    goto CLEANUP;
  }

  dname = env->getOptions()->findString("debug_sol")->getValue();
  fname = env->getOptions()->findString("problem_file")->getValue();
  if (""==fname) {
    qg.showHelp();
    goto CLEANUP;
  }

  // modify so that only rank 0 reads the problem and sends it to all others.
  // read the file and broadcast its contents to all others

  fname = share_problem_file(fname);

  p = qg.readProblem(fname, dname, "mqg", err);

  if (mpirank != 0)
  {
    std::filesystem::remove(fname);
  }

  if (err) {
    goto CLEANUP;
  }

  err = qg.solve(p);
  if (err) {
    goto CLEANUP;
  }

CLEANUP:
  MPI_Finalize();
  if (p) {
    delete p;
  }
  delete env;

  return 0;
}

std::string share_problem_file(std::string fname)
{
  int mpirank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  std::string ret;

  if (mpirank == 0)
  {
    std::ifstream file_in(fname);
    std::ostringstream buffer;
    buffer << file_in.rdbuf();
    std::string fstr = buffer.str();

    int fsize = fstr.size();
    MPI_Bcast(&fsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&fstr[0], fsize, MPI_CHAR, 0, MPI_COMM_WORLD);

    ret = fname;
  }
  else
  {
    int fsize;
    MPI_Bcast(&fsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::string fstr(fsize, '*');
    MPI_Bcast(&fstr[0], fsize, MPI_CHAR, 0, MPI_COMM_WORLD);

    int i = fname.size();
    while(i - 1 >= 0 && fname[i - 1] != '/')
      i--;
    
    while(fname[i] != '.')
    {
      ret.push_back(fname[i]);
      i++;
    }

    ret += "_copy" + std::to_string(mpirank);

    while (i < fname.size())
    {
      ret.push_back(fname[i]);
      i++;
    }

    std::ofstream out(ret);
    out << fstr;
    out.close();
  }

  return ret;
}
