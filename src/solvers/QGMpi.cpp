#include <mpi.h>
#include "QGMpi.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <iostream>
#include <sstream>

#include "AMPLHessian.h"
#include "AMPLInterface.h"
#include "AMPLJacobian.h"
#include "EngineFactory.h"
#include "Environment.h"
#include "Handler.h"
#include "IntVarHandler.h"
#include "LPEngine.h"
#include "LexicoBrancher.h"
#include "LinearHandler.h"
#include "Logger.h"
#include "MaxVioBrancher.h"
#include "MinotaurConfig.h"
#include "NLPEngine.h"
#include "NlPresHandler.h"
#include "NodeIncRelaxer.h"
#include "Objective.h"
#include "Option.h"
#include "PCBProcessor.h"
#include "Presolver.h"
#include "Problem.h"
#include "QGHandler.h"
#include "QPEngine.h"
#include "RCHandler.h"
#include "Reader.h"
#include "ReliabilityBrancher.h"
#include "SamplingHeur.h"
#include "StrongBrancher.h"
#include "Timer.h"
#include "TransSep.h"
#include "Types.h"
#include "MpiBranchAndBound.h"

using namespace Minotaur;

int QGMpi::solve(ProblemPtr p)
{
  OptionDBPtr options = env_->getOptions();
  Timer* timer = env_->getNewTimer();
  EnginePtr nlp_e = 0;
  LPEnginePtr lp_e = 0; // lp engine
  VarVector* orig_v = 0;
  MpiBranchAndBound* bab = 0;
  PresolverPtr pres = 0;
  BrancherPtr br = BrancherPtr(); // NULL
  PCBProcessorPtr nproc;
  NodeIncRelaxerPtr nr;

  // handlers
  HandlerVector handlers;
  IntVarHandlerPtr v_hand;
  LinearHandlerPtr l_hand;
  QGHandlerPtr qg_hand;
  RCHandlerPtr rc_hand;

  oinst_ = p;
  oinst_->calculateSize();
  int err = 0;

  int mpirank, comm_world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  std::ostringstream out;

  timer->start();
  if(options->findBool("display_problem")->getValue() == true) {
    oinst_->write(env_->getLogger()->msgStream(LogNone), 12);
  }

  if(options->findBool("display_size")->getValue() == true) {
    oinst_->writeSize(env_->getLogger()->msgStream(LogNone));
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "Starting constraint classification\n";
    oinst_->classifyCon(false);
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "Finished constraint classification\n";
  }

  // setup the jacobian and hessian
  if(false == options->findBool("use_native_cgraph")->getValue()) {
    JacobianPtr jac = new MINOTAUR_AMPL::AMPLJacobian(iface_);
    oinst_->setJacobian(jac);

    // create the hessian
    HessianOfLagPtr hess = new MINOTAUR_AMPL::AMPLHessian(iface_);
    oinst_->setHessian(hess);

    // set initial point
    oinst_->setInitialPoint(iface_->getInitialPoint(),
                            oinst_->getNumVars() - iface_->getNumDefs());
  }

  if(oinst_->getObjective() &&
     oinst_->getObjective()->getObjectiveType() == Maximize) {
    objSense_ = -1.0;
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "objective sense: maximize (will be converted to Minimize)"
        << std::endl;
  } else {
    objSense_ = 1.0;
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "objective sense: minimize" << std::endl;
  }

  // get presolver.
  orig_v = new VarVector(oinst_->varsBegin(), oinst_->varsEnd());
  pres = presolve_(handlers);
  for(HandlerVector::iterator it = handlers.begin(); it != handlers.end();
      ++it) {
    delete(*it);
  }
  handlers.clear();
  status_ = pres->getStatus();
  if(Finished != status_ && NotStarted != status_) {
    env_->getLogger()->msgStream(LogInfo)
        << me_ << "status of presolve: " << getSolveStatusString(status_)
        << std::endl;
    writeSol_(env_, orig_v, pres, SolutionPtr(), status_, iface_);
    goto CLEANUP;
  }

  // transform to exploit separability
  sepDetection();

  // create engines for solving LPs and NLPs
  err = getEngines_(&nlp_e, &lp_e);
  if(err) {
    goto CLEANUP;
  }

  if(true == options->findBool("use_native_cgraph")->getValue()) {
    oinst_->setNativeDer();
  }
  if(options->findBool("rc_fix")->getValue() && 0) {
    rc_hand = (RCHandlerPtr) new RCHandler(env_);
    rc_hand->setModFlags(false, true);
    handlers.push_back(rc_hand);
    assert(rc_hand);
  }
  // Initialize the handlers for branch-and-cut
  l_hand = (LinearHandlerPtr) new LinearHandler(env_, oinst_);
  l_hand->setModFlags(false, true);
  handlers.push_back(l_hand);
  assert(l_hand);

  v_hand = (IntVarHandlerPtr) new IntVarHandler(env_, oinst_);
  v_hand->setModFlags(false, true);
  handlers.push_back(v_hand);
  assert(v_hand);

  qg_hand = (QGHandlerPtr) new QGHandler(env_, oinst_, nlp_e);
  qg_hand->setModFlags(false, true);

  handlers.push_back(qg_hand);
  assert(qg_hand);

  // report name
  env_->getLogger()->msgStream(LogExtraInfo)
      << me_ << "handlers used:" << std::endl;
  for(HandlerIterator h = handlers.begin(); h != handlers.end(); ++h) {
    env_->getLogger()->msgStream(LogExtraInfo)
        << me_ << (*h)->getName() << std::endl;
  }

  // Only store bound-changes of relaxation (not problem)
  nr = (NodeIncRelaxerPtr) new NodeIncRelaxer(env_, handlers);
  nr->setModFlag(false);
  nr->setEngine(lp_e);
  nproc = (PCBProcessorPtr) new PCBProcessor(env_, lp_e, handlers);
  if(env_->getOptions()->findString("brancher")->getValue() == "rel" ||
     env_->getOptions()->findString("brancher")->getValue() == "weak" ||
     env_->getOptions()->findString("brancher")->getValue() == "relstronger") {
    ReliabilityBrancherPtr rel_br =
        (ReliabilityBrancherPtr) new ReliabilityBrancher(env_, handlers);
    rel_br->setEngine(lp_e);
    nproc->setBrancher(rel_br);
    br = rel_br;
  } else if(env_->getOptions()->findString("brancher")->getValue() ==
            "maxvio") {
    MaxVioBrancherPtr mbr =
        (MaxVioBrancherPtr) new MaxVioBrancher(env_, handlers);
    nproc->setBrancher(mbr);
    br = mbr;
  } else if(env_->getOptions()->findString("brancher")->getValue() == "lex") {
    LexicoBrancherPtr lbr =
        (LexicoBrancherPtr) new LexicoBrancher(env_, handlers);
    br = lbr;
  } else if(env_->getOptions()
                ->findString("brancher")
                ->getValue()
                .find("strong") != std::string::npos) {
    StrongBrancherPtr str_br =
        (StrongBrancherPtr) new StrongBrancher(env_, handlers);
    str_br->setEngine(lp_e);
    if(env_->getOptions()->findString("brancher")->getValue() == "stronger") {
      str_br->doStronger();
      str_br->setProblem(oinst_);
    }
    br = str_br;
  }
  nproc->setBrancher(br);
  env_->getLogger()->msgStream(LogExtraInfo)
      << me_ << "brancher used = " << br->getName() << std::endl;

  bab = new MpiBranchAndBound(env_, oinst_, mpirank, comm_world_size);
  bab->setNodeRelaxer(nr);
  bab->setNodeProcessor(nproc);
  bab->shouldCreateRoot(true);

  if(env_->getOptions()->findBool("samplingheur")->getValue() == true) {
    SamplingHeurPtr s_heur = (SamplingHeurPtr) new SamplingHeur(env_, oinst_);
    bab->addPreRootHeur(s_heur);
  }

  // start solving
  bab->solve();

  status_ = bab->getStatus();
  lb_ = bab->getLb();
  ub_ = bab->getUb();
  sol_ = bab->getSolution();

  out << "----------------------------------------------------------------------------------------------\n";
  out << "Rank " << mpirank << ":\n\n";

  bab->writeStats(out);
  nlp_e->writeStats(out);
  lp_e->writeStats(out);

  for(HandlerVector::iterator it = handlers.begin(); it != handlers.end();
      ++it) {
    (*it)->writeStats(out);
  }

  /* out << "----------------------------------------------------------------------------------------------\n"; */

  env_->getLogger()->msgStream(LogExtraInfo) << out.str();

  err = writeSol_(env_, orig_v, pres, sol_, status_, iface_);
  if(err) {
    goto CLEANUP;
  }
  if (mpirank == 0)
    err = writeBnbStatus_(bab);

CLEANUP:
  for(HandlerVector::iterator it = handlers.begin(); it != handlers.end();
      ++it) {
    delete(*it);
  }
  if(lp_e) {
    delete lp_e;
  }
  if(nlp_e) {
    delete nlp_e;
  }
  if(iface_ && ownIface_) {
    delete iface_;
    iface_ = 0;
  }
  if(pres) {
    delete pres;
  }
  if(bab) {
    if(bab->getNodeRelaxer()) {
      delete bab->getNodeRelaxer();
    }
    if(bab->getNodeProcessor()) {
      delete bab->getNodeProcessor();
    }
    delete bab;
  }
  if(orig_v) {
    delete orig_v;
  }
  if(timer) {
    delete timer;
  }
  oinst_ = 0;
  return err;
}


