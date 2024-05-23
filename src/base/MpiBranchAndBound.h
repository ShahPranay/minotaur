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
      int getStatusFlag();

      void solve() override;
      void showStatus_(bool current_uncounted, bool last_line) override;
      bool shouldStop_() override;

    private:
      void collectData_();
      bool shouldBalanceLoad_();
      void sendToAll_(double value);
      NodePtr LoadBalance_(NodePtr current_node);

      unsigned mpirank_, comm_world_size_;
      bool all_finished_;

      Timer *lb_timer_;
  };
}
#endif
