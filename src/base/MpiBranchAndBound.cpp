#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

#include "mpi.h"
#include "MpiBranchAndBound.h"
#include "Serializer.h"

using namespace Minotaur;
using std::cout, std::endl;

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
    if(tm_->getActiveNodes() >= comm_world_size_)
      return true;
    else return false;
  }
  else return false;
}

// returns new current_node
NodePtr MpiBranchAndBound::LoadBalance_()
{
  constexpr unsigned MIN_NODES_PER_RANK = 1;
  constexpr double MAX_LB = 1e18;
  unsigned num_send_nodes = MIN_NODES_PER_RANK * comm_world_size_, num_tot_nodes = num_send_nodes * comm_world_size_;

  std::vector<double> myLbs, worldLbs(num_tot_nodes, 9);
  std::vector<NodePtr> poppedNodes;

  NodePtr curnode = tm_->getCandidate();
  for(unsigned i = 0; curnode && i < num_send_nodes; i++)
  {
    myLbs.push_back(curnode->getLb());

    poppedNodes.push_back(curnode);
    tm_->removeActiveNode(curnode);
    curnode = tm_->getCandidate();
  }

  // change default value
  myLbs.resize(num_send_nodes, MAX_LB);
  /* if (mpirank_ <= 1) */
  /* { */
  /*   for ( double tmp : myLbs ) */
  /*   { */
  /*     cout << tmp << " "; */
  /*   } */
  /*   cout << endl; */
  /* } */

  MPI_Allgather(&myLbs[0], num_send_nodes, MPI_DOUBLE, &worldLbs[0], num_send_nodes, MPI_DOUBLE, MPI_COMM_WORLD);

  /* if (mpirank_ == 0) */
  /* { */
  /*   for ( double tmp : worldLbs ) */
  /*   { */
  /*     cout << tmp << " "; */
  /*   } */
  /*   cout << endl; */
  /* } */

  struct NodeInfo
  {
    unsigned owner_rank, local_index;
    double Lb;

    public:
    NodeInfo(unsigned owr, unsigned li, double lb) : owner_rank(owr), local_index(li), Lb(lb) {  }
  };

  auto cmpNodeInfo = [](const NodeInfo &a, const NodeInfo &b) -> bool 
  {
    return a.Lb < b.Lb;
  };

  std::vector<NodeInfo> worldNodeInfos;

  for (unsigned blockrank = 0; blockrank < comm_world_size_; blockrank++)
  {
    for (unsigned local_index = 0; local_index < num_send_nodes; local_index++)
    {
      worldNodeInfos.push_back(NodeInfo(blockrank, local_index, worldLbs[blockrank * num_send_nodes + local_index]));
    }
  }

  std::sort(worldNodeInfos.begin(), worldNodeInfos.end(), cmpNodeInfo);

  if (mpirank_ == 1)
  {
    for (auto curinfo : worldNodeInfos)
    {
      cout << curinfo.owner_rank << ", " << curinfo.Lb << endl;
    }
  }

  if ( worldNodeInfos[0].Lb == MAX_LB ){
    all_finished_ = true;
    return nullptr;
  }

  for (unsigned i = 0; i < num_tot_nodes; i++)
  {
    unsigned receiver_rank = i % comm_world_size_;
    NodeInfo &curinfo = worldNodeInfos[i];

    if ( curinfo.Lb == MAX_LB )
      break;

    if ( mpirank_ == receiver_rank )
    {
      if ( curinfo.owner_rank == mpirank_ )
      {
        tm_->insertRecvCandidate(poppedNodes[curinfo.local_index]);
      }
      else 
      {
        std::string tmpbuf(200, '*');
        MPI_Status status;
        MPI_Recv((void *) &(tmpbuf[0]), 200, MPI_CHAR, curinfo.owner_rank, 0, MPI_COMM_WORLD, &status);

        DeSerializer nodedes(tmpbuf);
        NodePtr newnode = nodedes.readNode(nodeRlxr_->getRelaxation());

        cout << mpirank_ << ": Recv node ID: " << newnode->getId() << " with Lb: " << newnode->getLb() << endl;
        tm_->insertRecvCandidate(newnode);
      }
    }
    else if ( worldNodeInfos[i].owner_rank == mpirank_ )
    {
      Serializer nodesr;
      nodesr.writeNode(poppedNodes[curinfo.local_index]);
      std::string msg = nodesr.get_string();
      /* std::cout << "Message size = " << msg.size() << std::endl; */
      // need to maintain buffers and requests for Isend.
      MPI_Send((void *) &(msg[0]), msg.size(), MPI_CHAR, receiver_rank, 0, MPI_COMM_WORLD);
    }
  }

  return tm_->getCandidate();
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

  unsigned itr_cnt = 0;

  // solve root outside the loop. save the useful information.
  while(!all_finished_ && itr_cnt++ < 5) {
    collectData_();

    /* if (shouldBalanceLoad_()) */
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
