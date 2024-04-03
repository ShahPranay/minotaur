#include <cmath>

#include "mpi.h"
#include "MpiBranchAndBound.h"
#include "Serializer.h"

using namespace Minotaur;

MpiBranchAndBound::MpiBranchAndBound(EnvPtr env, ProblemPtr p, int &cur_rank, int &comm_world_size) :
  BranchAndBound(env, p), mpirank_(cur_rank), comm_world_size_(comm_world_size), all_finished_(false)
{

}

MpiBranchAndBound::~MpiBranchAndBound()
{

}

void MpiBranchAndBound::collectData_()
{

}

bool MpiBranchAndBound::shouldBalanceLoad_()
{
  if(tm_->getActiveNodes() == 0)
    return true;

  if(mpirank_ == 0)
  {
    if(tm_->getActiveNodes() >= 2)
      return true;
    else return false;
  }
  else return false;
}

// returns new current_node
NodePtr MpiBranchAndBound::LoadBalance_()
{
  int nodesCnt = 0, myContrib = tm_->getActiveNodes();

  MPI_Allreduce(&myContrib, &nodesCnt, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(nodesCnt == 0)
  {
    all_finished_ = true;
    return nullptr;
  }

  if(mpirank_ == 0)
  {
    NodePtr curnode = tm_->getCandidate();
    while(curnode)
    {
      Serializer nodesr;
      nodesr.writeNode(curnode);
      std::string msg = nodesr.get_string();
      std::cout << "Message size = " << msg.size() << std::endl;
      // need to maintain buffers and requests for Isend.
      MPI_Send((void *) &(msg[0]), msg.size(), MPI_CHAR, 1, 0, MPI_COMM_WORLD);

      tm_->removeActiveNode(curnode);
      curnode = tm_->getCandidate();
    }

    return nullptr;
  }
  else
  {
    std::cout << "Num Nodes Recv = " << nodesCnt << "\n";

    for(int i = 0; i < nodesCnt; i++)
    {
      std::string tmpbuf(200, '*');
      MPI_Status status;
      MPI_Recv((void *) &(tmpbuf[0]), 200, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);

      DeSerializer nodedes(tmpbuf);
      NodePtr newnode = nodedes.readNode(nodeRlxr_->getRelaxation());
      tm_->insertRecvCandidate(newnode);
    }
    return tm_->getCandidate();
  }
}

void MpiBranchAndBound::solve()
{
  bool should_dive = false, dived_prev = false;
  bool should_prune = false;
  NodePtr current_node = NodePtr();
  NodePtr new_node = NodePtr();
  Branches branches = 0;
  WarmStartPtr ws;
  RelaxationPtr rel = RelaxationPtr();
  // initialize timer
  timer_->start();
  logger_->msgStream(LogInfo) << me_ << "starting branch-and-bound"
    << std::endl;

  // get problem size and statistics to detect problem type.
  problem_->calculateSize();
#if SPEW
  problem_->writeSize(logger_->msgStream(LogExtraInfo));
#endif

  // initialize statistics
  if (stats_) {
    delete stats_;
  }
  stats_ = new BabStats();

  // initialize solution pool
  // TODO: use user options to set the pool size. For now it is 1.
  solPool_ = (SolutionPoolPtr) new SolutionPool(env_, problem_, 1);

  // call heuristics before the root, if needed 
  for (HeurVector::iterator it=preHeurs_.begin(); it!=preHeurs_.end(); ++it) {
    (*it)->solve(current_node, rel, solPool_);
  }
  tm_->setUb(solPool_->getBestSolutionValue());

  if (mpirank_ == 0){
    // do the root
    current_node = processRoot_(&should_prune, &dived_prev);

    // stop if done
    if (!current_node) {
      tm_->updateLb();
      if (tm_->getUb() <= -INFINITY) {
        status_ = SolvedUnbounded;
      } else  if (tm_->getUb() < INFINITY) {
        status_ = SolvedOptimal; 
      } else {
        status_ = SolvedInfeasible; 
      }
#if SPEW
      logger_->msgStream(LogDebug) << me_ << "stopping after root node "
        << std::endl;
#endif
      // TODO:signal all to stop
    } else if (shouldStop_()) {
      tm_->updateLb();
      // TODO:signal all to stop
    } else {
#if SPEW
      logger_->msgStream(LogDebug) << std::setprecision(8)
        << me_ << "lb = " << tm_->updateLb() << std::endl
        << me_ << "ub = " << tm_->getUb() << std::endl;
#endif
    }
  }
  else{
    if (options_->createRoot == true) {
      nodeRlxr_->createRootRelaxation(current_node, should_prune)->setProblem(problem_);
    }
  }

  // solve root outside the loop. save the useful information.
  while(!all_finished_) {
    collectData_();

    if (shouldBalanceLoad_())
      current_node = LoadBalance_();

    if(!current_node)
      continue;

#if SPEW
    logger_->msgStream(LogDebug1) << me_ << "processing node "
      << current_node->getId() << std::endl
      << me_ << "depth = " << current_node->getDepth() << std::endl
      << me_ << "did we dive = " << dived_prev << std::endl;
#endif

    should_dive = false;
    rel = nodeRlxr_->createNodeRelaxation(current_node, dived_prev, 
        should_prune);
    nodePrcssr_->process(current_node, rel, solPool_);

    ++stats_->nodesProc;
#if SPEW
    logger_->msgStream(LogDebug1) << me_ << "node lower bound = " << 
      current_node->getLb() << std::endl;
#endif

    if (nodePrcssr_->foundNewSolution()) {
      tm_->setUb(solPool_->getBestSolutionValue());
    }

    should_prune = shouldPrune_(current_node);
    if (should_prune) {
#if SPEW
      logger_->msgStream(LogDebug1) << me_ << "node pruned" << 
        std::endl;
#endif
      nodeRlxr_->reset(current_node, false);
      tm_->pruneNode(current_node);
      if (!dived_prev) {
        tm_->removeActiveNode(current_node);
      }
      new_node = tm_->getCandidate();
      dived_prev = false;
    } else {
#if SPEW
      logger_->msgStream(LogDebug1) << me_ << "branching" << 
        std::endl;
#endif
      branches = nodePrcssr_->getBranches();

      ws = nodePrcssr_->getWarmStart();
      if (!dived_prev) {
        tm_->removeActiveNode(current_node);
      }
      should_dive = tm_->shouldDive();

      new_node = tm_->branch(branches, current_node, ws);
      assert((should_dive && new_node) || (!should_dive && !new_node));
      if (should_dive) {
        dived_prev = true;
      } else {
        nodeRlxr_->reset(current_node, false);
        new_node = tm_->getCandidate(); // Can be NULL. The branches that were
                                        // created could have large lb and tm 
                                        // might have eliminated them.
        dived_prev = false;
      }
    }
    current_node = new_node;

    showStatus_(should_dive, false);

    // TODO:stop if done
    /* if (shouldStop_()) { */
    /*   tm_->updateLb(); */
    /*   // broadcast to other compute nodes */
    /*   break; */
    /* } else { */
#if SPEW
    logger_->msgStream(LogDebug) << std::setprecision(8)
      << me_ << "lb = " << tm_->getLb() << std::endl 
      << me_ << "ub = " << tm_->getUb() << std::endl;
#endif
    /* } */
  }

  // TODO:modify to execute only when shouldStop_() did not quit solver
  if (true){
    tm_->updateLb();
    if (tm_->getUb() <= -INFINITY) {
      status_ = SolvedUnbounded;
    } else if (tm_->getUb() < INFINITY) {
      status_ = SolvedOptimal; // TODO: get the right status
    } else {
      status_ = SolvedInfeasible; // TODO: get the right status
    }
#if SPEW
    logger_->msgStream(LogDebug) << me_ << "all nodes have "
      << "been processed" << std::endl;
#endif
  }

  showStatus_(false, true);
  //logger_->msgStream(LogError) << " " << std::endl;
  logger_->msgStream(LogError) << "----------------------------------------------------------------------------------------------" << std::endl;
  logger_->msgStream(LogError) << " " << std::endl;

  logger_->msgStream(LogInfo) << me_ << "stopping branch-and-bound"
    << std::endl
    << me_ << "nodes processed = " << stats_->nodesProc << std::endl
    << me_ << "nodes created   = " << tm_->getSize() << std::endl;
  stats_->timeUsed = timer_->query();
  timer_->stop();
}
