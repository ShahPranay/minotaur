#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>

#include "mpi.h"
#include "MpiBranchAndBound.h"
#include "Serializer.h"

using namespace Minotaur;
using std::cout, std::endl;

MpiBranchAndBound::MpiBranchAndBound(EnvPtr env, ProblemPtr p, int &cur_rank, int &comm_world_size) :
  BranchAndBound(env, p), mpirank_(cur_rank), comm_world_size_(comm_world_size), all_finished_(false), lb_frequency_(comm_world_size)
{
  /*lb_timer_ = env->getNewTimer();*/
}

MpiBranchAndBound::~MpiBranchAndBound()
{

}

int MpiBranchAndBound::getStatusFlag()
{
  if (status_ == TimeLimitReached)
    return 1;
  else if (status_ == IterationLimitReached)
    return 2;
  else if (status_ == SolLimitReached)
    return 4;
  else return 0;
}

void MpiBranchAndBound::collectData_()
{
  int ismsg = 0;
  double value;
  while (true) {
    MPI_Status status;
    MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &ismsg, &status);
    if (!ismsg)
      return;
    MPI_Recv(&value, 1, MPI_DOUBLE, status.MPI_SOURCE, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    /* cout << "Recieved sol value = " << value << endl; */
    if (value < tm_->getUb())
      tm_->setUb(value);
  }
}

bool MpiBranchAndBound::shouldBalanceLoad_()
{
  /*if (lb_timer_->query() < lb_frequency_)*/
  /*  return false;*/
  return true;
}

void MpiBranchAndBound::updateLbFrequency_(std::vector<double> &worldLbs)
{
  lb_frequency_ = 50 * comm_world_size_;
  /*int cnt = 0;*/
  /*double sum = 0, sqsum = 0;*/
}

static SolveStatus flagToStatus(int flag)
{
  if (flag & 1)
    return TimeLimitReached;
  else if (flag & 2)
    return IterationLimitReached;
  else if (flag & 4)
    return SolLimitReached;
  else return NotStarted;
}

// returns new current_node
NodePtr MpiBranchAndBound::LoadBalance_(NodePtr current_node)
{
  constexpr unsigned MIN_NODES_PER_RANK = 50;
  constexpr double MAX_LB = INFINITY;
  unsigned num_send_nodes = MIN_NODES_PER_RANK * comm_world_size_, num_tot_nodes = num_send_nodes * comm_world_size_;

  int statusflag = getStatusFlag();
  MPI_Allreduce(MPI_IN_PLACE, &statusflag, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
  if (statusflag != 0)
  {
    status_ = flagToStatus(statusflag);
    all_finished_ = true;
    return nullptr;
  }

  std::vector<double> myLbs, worldLbs(num_tot_nodes, 9);
  std::vector<NodePtr> poppedNodes;
  NodePtr curnode = current_node;
  for(unsigned i = 0; curnode && i < num_send_nodes; i++)
  {
    myLbs.push_back(curnode->getLb());

    poppedNodes.push_back(curnode);
    tm_->removeActiveNode(curnode);
    curnode = tm_->getCandidate();
  }

  myLbs.resize(num_send_nodes, MAX_LB);

  MPI_Allgather(&myLbs[0], num_send_nodes, MPI_DOUBLE, &worldLbs[0], num_send_nodes, MPI_DOUBLE, MPI_COMM_WORLD);

  updateLbFrequency_(worldLbs);

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
        /* cout << "Rank " << mpirank_ << "; ReInserting NodeID: " << poppedNodes[curinfo.local_index]->getId() << endl; */
        tm_->insertPoppedCandidate(poppedNodes[curinfo.local_index]);
      }
      else 
      {
        MPI_Status status;
        MPI_Probe(curinfo.owner_rank, 0, MPI_COMM_WORLD, &status);
        
        int count;
        MPI_Get_count(&status, MPI_CHAR, &count);

        /*cout << count << endl;*/

        std::string tmpbuf(count, '*');
        MPI_Recv((void *) &(tmpbuf[0]), count, MPI_CHAR, curinfo.owner_rank, 0, MPI_COMM_WORLD, &status);

        DeSerializer nodedes(tmpbuf);
        NodePtr newnode = nodedes.readNode(nodeRlxr_->getRelaxation());

        /* cout << "Rank " << mpirank_ << "; Recv node ID: " << newnode->getId() << " with Lb: " << newnode->getLb() << endl; */
        tm_->insertRecvCandidate(newnode);
      }
    }
    else if ( worldNodeInfos[i].owner_rank == mpirank_ )
    {
      Serializer nodesr;
      nodesr.writeNode(poppedNodes[curinfo.local_index]);
      tm_->pruneNode(poppedNodes[curinfo.local_index]); // remove from current tree
      std::string msg = nodesr.get_string();
      /* std::cout << "Message size = " << msg.size() << std::endl; */
      // need to maintain buffers and requests for Isend.
      MPI_Send((void *) &(msg[0]), msg.size(), MPI_CHAR, receiver_rank, 0, MPI_COMM_WORLD);
      tm_->sentNode();
    }
  }

  /*lb_timer_->stop();*/
  /*lb_timer_->start();*/
  /* cout << endl; */

  return tm_->getCandidate();
}

void MpiBranchAndBound::sendToAll_(double value)
{
  /* cout << "Sendin sol value = " << value << endl; */
  for (unsigned r = 0; r < comm_world_size_; r++)
  {
    if (r == mpirank_)
      continue;
    MPI_Request req;
    MPI_Isend(&value, 1, MPI_DOUBLE, r, 1, MPI_COMM_WORLD, &req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
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

  if (mpirank_ == 0) {
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
  } else {
    if (options_->createRoot == true) {
      nodeRlxr_->createRootRelaxation(current_node, should_prune)->setProblem(problem_);
    }
  }

  /*lb_timer_->start();*/
  int times_balanced = 0, itr_since_last_balance = 0;

  // solve root outside the loop. save the useful information.
  while (!all_finished_) {
    collectData_();

    if (shouldStop_() || !current_node || itr_since_last_balance == lb_frequency_) {
      times_balanced += 1;
      current_node = LoadBalance_(current_node);
      itr_since_last_balance = 0;
    }

    if (!current_node)
      continue;

    itr_since_last_balance++;

    /* cout << "Rank " << mpirank_ << "; Processing; No. ActiveNodes = " << tm_->getActiveNodes() << "; Current NodeID = " << current_node->getId() << "; Lb = " << current_node->getLb() << endl; */

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

    if (nodePrcssr_->foundNewSolution()) 
    {
      double solval = solPool_->getBestSolutionValue();
      if (solval < tm_->getUb())
        tm_->setUb(solval);
      sendToAll_(solval);
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

    /* if(!current_node) */
      /* cout << "Rank " << mpirank_ << ": Afterprocessing; ActiveNodes = " << tm_->getActiveNodes() << "Node pruned!" << endl; */

    if (mpirank_ == 0)
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

  double tmpub = tm_->getUb();
  MPI_Allreduce(MPI_IN_PLACE, &tmpub, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  tm_->setUb(tmpub);

  // TODO:modify to execute only when shouldStop_() did not quit solver
  if (!shouldStop_()){
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

  // get all data to rank 1
  if (mpirank_ == 0)
  {
    showStatus_(false, true);
    logger_->msgStream(LogError) << "----------------------------------------------------------------------------------------------" << std::endl;
    logger_->msgStream(LogError) << " " << std::endl;

    std::vector<unsigned> allstats(2 * comm_world_size_), curstats(2);
    curstats[0] = stats_->nodesProc;
    curstats[1] = tm_->getSize();
    MPI_Gather(&curstats[0], 2, MPI_UNSIGNED, &allstats[0], 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    std::ostringstream out;
    out << me_ << "stopping branch-and-bound" << std::endl;
    out << "-------------------------------------------------" << std::endl;
    out << std::setw(15) << "MPI Rank"
      << std::setw(15) << "Nodes Proc"
      << std::setw(15) << "Nodes Created" << std::endl;
    out << "-------------------------------------------------" << std::endl;

    for (unsigned i = 0; i < comm_world_size_; i++)
    {
      out << std::setw(15) << i
        << std::setw(15) << allstats[2 * i]
        << std::setw(15) << allstats[2 * i + 1] << std::endl;
    }
    out << "-------------------------------------------------" << std::endl;
    out << "Number of times balanced = " << times_balanced << std::endl;
    logger_->msgStream(LogInfo) << out.str();
  }
  else
  {
    std::vector<unsigned> curstats(2);
    curstats[0] = stats_->nodesProc;
    curstats[1] = tm_->getSize();
    MPI_Gather(&curstats[0], 2, MPI_UNSIGNED, &curstats[0], 2, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  }

  stats_->timeUsed = timer_->query();
  timer_->stop();
  /*lb_timer_->stop();*/

  MPI_Barrier(MPI_COMM_WORLD);
}

void MpiBranchAndBound::showStatus_(bool current_uncounted, bool last_line)
{
  static bool header = false;
  static bool firstRow = true; // Add a flag for the first row

  UInt off = 0;
  if (current_uncounted) {
    off = 1;
  }

  std::streamsize defaultPrecision = logger_->msgStream(LogError).precision(); // Store default precision

  if (!header) {
    logger_->msgStream(LogError) << " " << std::endl;	  
    logger_->msgStream(LogError) << "----------------------------------------------------------------------------------------------" << std::endl;
    logger_->msgStream(LogError) << std::setw(7) << "Cpu(s)"
      << std::setw(10) << "Wall(s)"
      << std::setw(16) << "LB"
      << std::setw(13) << "UB"
      << std::setw(12) << "Gap%"
      << std::setw(15) << "   Nodes-Proc"
      << std::setw(14) << "   Nodes-Rem"
      << std::setw(7) << "#Sol"
      << std::endl;
    logger_->msgStream(LogError) << "----------------------------------------------------------------------------------------------" << std::endl;
    header = true;
  }

  if (firstRow) {
    // Print the initial row with all values set to zero
    logger_->msgStream(LogError) << std::setw(6) << "0"
      << std::setw(10) << "0"
      << std::setw(17) << "-inf"
      << std::setw(13) << "inf"
      << std::setw(12) << "inf"
      << std::setw(15) << "0"
      << std::setw(14) << "0"
      << std::setw(7) << "0"
      << std::endl;
    firstRow = false;
  }

  if (timer_->query() - stats_->updateTime > options_->logInterval || last_line == true) {
    logger_->msgStream(LogError).precision(defaultPrecision); // Reset precision to default
    logger_->msgStream(LogError) << std::setw(6) << std::fixed << std::setprecision(0) << timer_->query()
      << std::setw(10) << std::fixed << std::setprecision(0) << timer_->wQuery() // Print wall time
      << std::setw(17) << std::setprecision(4) << std::scientific << tm_->updateLb()
      << std::setw(13) << std::setprecision(4) << std::scientific << tm_->getUb()
      << std::setw(12) << std::setprecision(2) << std::scientific << tm_->getPerGap()
      << std::setw(15) << tm_->getSize() - tm_->getActiveNodes() - off
      << std::setw(14) << tm_->getActiveNodes() + off
      << std::setw(7) << std::setprecision(5) << solPool_->getNumSolsFound()
      << std::endl;
    stats_->updateTime = timer_->query();
  }
}

bool MpiBranchAndBound::shouldStop_()
{
  bool stop_bnb = false;

  /*if ( tm_->getPerGap() <= options_->perGapLimit) {*/
  /*  stop_bnb = true;*/
  /*  status_ = SolvedGapLimit;*/
  /*}*/

  if (timer_->query() > options_->timeLimit) {
    stop_bnb = true;
    status_ = TimeLimitReached;
  } else if (stats_->nodesProc >= options_->nodeLimit) {
    stop_bnb = true;
    status_ = IterationLimitReached;
  } else if (solPool_->getNumSolsFound()>=options_->solLimit) { 
    stop_bnb = true;
    status_ = SolLimitReached;
  }

  return stop_bnb;
}
