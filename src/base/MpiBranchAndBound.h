#ifndef MINOTAURMPIBNB_H
#define MINOTAURMPIBNB_H

#include "Types.h"
#include "BranchAndBound.h"

namespace Minotaur {
  class MpiBranchAndBound : public BranchAndBound {
    public:
      MpiBranchAndBound(EnvPtr env, ProblemPtr p, int &cur_rank, int &comm_world_size);

      /// Destroy.
      virtual ~MpiBranchAndBound();

      void solve() override;

    private:
      void collectData_();
      bool shouldBalanceLoad_();
      NodePtr LoadBalance_(NodePtr current_node);

      unsigned mpirank_, comm_world_size_;
      bool all_finished_;
  };
}
#endif
