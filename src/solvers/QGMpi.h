#ifndef QGMPI_H
#define QGMPI_H

#include "QG.h"

namespace Minotaur {

  class QGMpi : public QG {
    public:
      using QG::QG;

      virtual int solve(ProblemPtr p);
      virtual int writeBnbStatus_(BranchAndBound *bab);
  };
}

#endif
