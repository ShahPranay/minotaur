#include <mpi.h>
#include <stdio.h>

#include "MinotaurConfig.h"
#include "Environment.h"
#include "Option.h"
#include "Problem.h"
#include "Types.h"
#include "QGMpi.h"

using namespace Minotaur;

int main(int argc, char** argv)
{
  MPI_Init(NULL, NULL);

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
  p = qg.readProblem(fname, dname, "mqg", err);
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

