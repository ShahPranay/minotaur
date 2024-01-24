//
//     Minotaur -- It's only 1/2 bull
//
//     (C)opyright 2008 - 2024 The Minotaur Team.
//

/**
 * \file IpoptEngine.cpp
 * \brief Define IpoptEngine class for solving NLPs using IPOPT.
 * \author Ashutosh Mahajan, Argonne National Laboratory
 */

#include <cmath>
#include <iostream>

// undefine some package names
// #include "MinotaurDeconfig.h" // no need for now.
#include "coin-or/IpTNLP.hpp"
#include "coin-or/IpIpoptApplication.hpp"
#include "coin-or/IpIpoptCalculatedQuantities.hpp"
#include "coin-or/IpSolveStatistics.hpp"
#undef F77_FUNC_
#undef F77_FUNC

#include "MinotaurConfig.h"
#include "Constraint.h"
#include "Environment.h"
#include "HessianOfLag.h"
#include "IpoptEngine.h"
#include "IpoptEngineTnlp.h"
#include "Jacobian.h"
#include "Logger.h"
#include "Objective.h"
#include "Option.h"
#include "Problem.h"
#include "ProblemSize.h"
#include "Timer.h"
#include "Variable.h"

//#define SPEW 1

namespace Minotaur {

const std::string IpoptEngine::me_ = "IpoptEngine: ";
      
// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //
IpoptSolution::IpoptSolution()
: Solution(),
  dualXLow_(0),
  dualXUp_(0)
{
}


IpoptSolution::IpoptSolution(const double *x, double objval, 
                             ProblemPtr problem)
: Solution(objval, x, problem),
  dualXLow_(0),
  dualXUp_(0)
{

}


IpoptSolution::IpoptSolution(ConstIpoptSolPtr sol)
{
  n_ = sol->n_;
  m_ = sol->m_;

  if (sol->x_) {
    x_ = new double[n_];
    std::copy(sol->x_, sol->x_+n_, x_);
  } else {
    x_ = 0;
  }

  if (sol->dualCons_) {
    dualCons_ = new double[m_];
    std::copy(sol->dualCons_, sol->dualCons_+m_, dualCons_);
  } else {
    dualCons_ = 0;
  }

  if (sol->dualX_) {
    dualX_ = new double[n_];
    std::copy(sol->dualX_, sol->dualX_+n_, dualX_);
  } else {
    dualX_ = 0;
  }

  if (sol->dualXLow_) {
    dualXLow_ = new double[n_];
    std::copy(sol->dualXLow_, sol->dualXLow_+n_, dualXLow_);
  } else {
    dualXLow_ = 0;
  }

  if (sol->dualXUp_) {
    dualXUp_ = new double[n_];
    std::copy(sol->dualXUp_, sol->dualXUp_+n_, dualXUp_);
  } else {
    dualXUp_ = 0;
  }

  consViol_ = INFINITY;
  objValue_ = INFINITY;
}


IpoptSolution::~IpoptSolution()
{
  if (dualXLow_) {
    delete [] dualXLow_;
  }
  if (dualXUp_) {
    delete [] dualXUp_;
  }
  // dualX_ is freed in base class dtor.
}


void IpoptSolution::setDualOfVars(const double *lower, const double *upper)
{
  if (lower && upper) {
    double *l, *u, *d;
    if (!dualXLow_) {
      dualXLow_ = new double[n_];
      dualXUp_ = new double[n_];
    }
    std::copy(lower, lower+n_, dualXLow_);
    std::copy(upper, upper+n_, dualXUp_);
    if (!dualX_) {
      dualX_ = new double[n_];
    }
    d = dualX_; 
    l = dualXLow_;
    u = dualXUp_;
    for (UInt i=0; i<n_; ++i, ++d, ++l, ++u) {
      *d = *l + *u;
    }
  }
}


void IpoptSolution::write(std::ostream &out) const
{
  const double *d;
  out << "Number of variables = "   << n_ << std::endl
      << "Number of constraints = " << m_ << std::endl
      << "Primal values:"           << std::endl;

  d=x_;
  if (problem_) {
    for (VariableConstIterator it=problem_->varsBegin(); 
        it!=problem_->varsEnd(); ++it, ++d) {
      out << (*it)->getName() << "   " << *d << std::endl;
    }
  } else {
    for (UInt i=0; i<n_; ++i) {
      out << x_[i] << std::endl;
    }
  }

  out << "Dual values of constraint multipliers:" << std::endl;
  d=dualCons_;
  for (UInt i=0; i<m_; ++i, ++d) {
    out << *d << std::endl;
  }

  out << "Dual values of variable lower bound multipliers:" << std::endl;
  d = dualXLow_;
  for (UInt i=0; i<n_; ++i, ++d) {
    out << *d << std::endl;
  }
  out << "Dual values of variable upper bound multipliers:" << std::endl;
  d = dualXUp_;
  for (UInt i=0; i<n_; ++i, ++d) {
    out << *d << std::endl;
  }
}

// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

IpoptWarmStart::IpoptWarmStart()
: sol_(0)
{
}


/// Copy constructor. Creates a full copy, not just copies pointers.
IpoptWarmStart::IpoptWarmStart(ConstIpoptWarmStartPtr ws)
{
  if (ws && ws->sol_) {
    sol_ = (IpoptSolPtr) new IpoptSolution(ws->sol_);
  } else {
    sol_ = 0;
  }
}


IpoptWarmStart::~IpoptWarmStart()
{
  if (sol_) {
    delete sol_;
  }
}


IpoptSolPtr IpoptWarmStart::getPoint()
{
  return sol_;
}


bool IpoptWarmStart::hasInfo()
{
  // return true if warm start carries a starting solution
  if (sol_ && sol_->getPrimal()) {
    return true;
  } else {
    return false;
  }
}


void IpoptWarmStart::setPoint(IpoptSolPtr sol)
{
  sol_ = sol;
}


void IpoptWarmStart::write(std::ostream &out) const
{
  out << "Ipopt warm start information:"  << std::endl;
  sol_->write(out);
}


// ----------------------------------------------------------------------- //
// ----------------------------------------------------------------------- //

IpoptEngine::IpoptEngine(EnvPtr env)
: bndChanged_(false),
  consChanged_(false),
  env_(env),
  etol_(1e-7),
  justLoaded_(false),
  myapp_(0),
  mynlp_(0),
  sol_(0),
  strBr_(false),
  timer_(0),
  ws_(0), 
  wsSb_(0) 
{
#if defined(USE_IPOPT)  
  problem_ = 0;
  logger_ = env_->getLogger();
  myapp_ = new Ipopt::IpoptApplication();
  setOptionsForRepeatedSolve();

  status_ = EngineError;
  if (env->getOptions()->findBool("use_warmstart")->getValue()==true) {
    prepareWs_ = true;
    useWs_ = true;
  } else {
    prepareWs_ = false;
    useWs_ = false;
  }
  timer_ = env->getNewTimer();

  stats_ = new IpoptStats();
  stats_->calls    = 0;
  stats_->opt      = 0;
  stats_->reopt    = 0;
  stats_->strCalls = 0;
  stats_->time     = 0;
  stats_->ptime    = 0;
  stats_->strTime  = 0;
  stats_->iters    = 0;
  stats_->strIters = 0;
  
#else 
  assert(!"ipopt engine can only be called when compiled with ipopt!")
#endif
}


IpoptEngine::~IpoptEngine()
{
  if (wsSb_) {
    delete wsSb_;
  }
  if (ws_) {
    delete ws_; // deletes sol_ as well
    sol_ = 0;
  } else if (sol_) {
    delete sol_;
  }
  if (timer_) {
    delete timer_;
  }
  if (myapp_) {
    delete myapp_; // deletes mynlp_ as well
  }
  if (stats_) {
    if (0==stats_->opt && 0==stats_->reopt && mynlp_) {
      delete mynlp_;
    }
    delete stats_;
  }
  if (problem_) {
    problem_->unsetEngine();
    problem_ = 0;
  }
  return;
}


void IpoptEngine::addConstraint(ConstraintPtr)
{
  consChanged_ = true;
}


void IpoptEngine::changeBound(ConstraintPtr, BoundType, double)
{
  bndChanged_ = true;
}


void IpoptEngine::changeBound(VariablePtr, BoundType, double)
{
  bndChanged_ = true;
}


void IpoptEngine::changeBound(VariablePtr, double, double)
{
  bndChanged_ = true;
}


void IpoptEngine::changeConstraint(ConstraintPtr, LinearFunctionPtr, 
                                   double , double)
{
  consChanged_ = true;
}


void IpoptEngine::changeConstraint(ConstraintPtr, NonlinearFunctionPtr)
{
  consChanged_ = true;
}


void IpoptEngine::changeObj(FunctionPtr, double)
{
  consChanged_ = true;
}


void IpoptEngine::clear() 
{
  if (ws_) {
    delete ws_; ws_ = 0; // deletes sol_ as well
    sol_ = 0;
  } else if (sol_) {
    delete sol_; sol_ = 0;
  }
  if (wsSb_) {
    delete wsSb_; wsSb_ = 0;
  }
  if (problem_) {
    problem_->unsetEngine();
    problem_ = 0;
  }
}


void IpoptEngine::disableStrBrSetup()
{
  prepareWs_ = useWs_;
  strBr_     = false;
  if (wsSb_) {
    delete wsSb_; wsSb_ = 0;
  }
}


EnginePtr IpoptEngine::emptyCopy()
{
  return (IpoptEnginePtr) new IpoptEngine(env_);
}


void IpoptEngine::enableStrBrSetup()
{
  if (ws_) {
    wsSb_ = new IpoptWarmStart(ws_);
  }
  prepareWs_ = false;
  strBr_     = true;
}


std::string IpoptEngine::getName() const
{
  return "IPOPT";
}


ConstSolutionPtr IpoptEngine::getSolution()
{
  if (sol_) {
    return sol_;
  } 
  return SolutionPtr(); // NULL
}


ConstWarmStartPtr IpoptEngine::getWarmStart()
{
  return ws_;
}


WarmStartPtr IpoptEngine::getWarmStartCopy()
{
  IpoptWarmStartPtr ws;
  if (ws_) {
    ws = (IpoptWarmStartPtr) new IpoptWarmStart(ws_); // copy
  } else {
    ws = 0; // NULL
  }
  return ws;
}


double IpoptEngine::getSolutionValue() 
{
  return mynlp_->getSolutionValue();
}


EngineStatus IpoptEngine::getStatus()
{
  return status_;
}


void IpoptEngine::load(ProblemPtr problem)
{
  if (problem_) {
    problem_->unsetEngine();
  }
  problem_ = problem;

  if (sol_) {
    delete sol_;
  }
  sol_ = new IpoptSolution(0, INFINITY, problem);
  mynlp_ = new Ipopt::IpoptFunInterface(env_, problem, sol_);
  //Ipopt::ApplicationReturnStatus status;
  //status = myapp_->Initialize();
  myapp_->Initialize();

  // check if warm start needs to be saved
  if (prepareWs_) {
    ws_ = new IpoptWarmStart();
    ws_->setPoint(sol_);
  }

  // set the status of the engine to unknown.
  status_ = EngineUnknownStatus;
  bndChanged_ = true;
  consChanged_ = true;
  problem->calculateSize();
  setOptionsForProb_();
  problem->setEngine(this);
  justLoaded_ = true;
}


void IpoptEngine::loadFromWarmStart(const WarmStartPtr ws)
{
  if (ws) {
    // dynamic cast can't seem to be avoided.
    ConstIpoptWarmStartPtr ws2 = dynamic_cast <const IpoptWarmStart*> (ws);

    // now create a full copy.
    if (ws_) {
      delete ws_; // sol_ gets freed
      ws_ = 0; sol_ = 0;
    }
    ws_ = new IpoptWarmStart(ws2);
    sol_ = ws_->getPoint();
    if (!useWs_) {
      logger_->msgStream(LogInfo) << me_ << "setWarmStart() method is called "
        << "but warm-start is not enabled." << std::endl;
    }
  } else {
    ws_ = 0;
  }
}


void IpoptEngine::negateObj()
{
  consChanged_ = true;
}


bool IpoptEngine::presolve_()
{
  bool should_stop = false;
  VariablePtr v;
  double diff;
  bool all_fixed = true;
  int e=0;

  status_ = EngineUnknownStatus;
  // visit each variable and check bounds. Stop if bad bounds are found or if
  // all variables are fixed.
  for (VariableConstIterator it=problem_->varsBegin(); it!=problem_->varsEnd(); 
      ++it) {
    v = *it;
    diff = v->getUb() - v->getLb();
    if (diff < -etol_) {
      status_ = ProvenLocalInfeasible;
      return true;
    }
    if (fabs(diff)>etol_) {
      all_fixed = false;
    }
  }

  if (all_fixed == true) {
    double obj_value;
    double *x;
    double vio = -INFINITY;
    double act;

    x = new double[problem_->getNumVars()];
    should_stop = true;
    for (VariableConstIterator it=problem_->varsBegin();
        it!=problem_->varsEnd(); ++it) {
      v = *it;
      x[v->getIndex()] = v->getUb();
    }
    obj_value = problem_->getObjValue(x, &e);

    // check if constraints are violated.
    for (ConstraintConstIterator it=problem_->consBegin(); 
        it!=problem_->consEnd(); ++it) {
      act = (*it)->getActivity(x, &e);
      vio = std::max(act-(*it)->getUb(), (*it)->getLb()-act);
      vio = std::max(vio, 0.);
      if (vio>etol_ || e!=0) {
        break;
      }
    }

   // sol_ must be set even if the problem is infeasible because QG handler
   // may need this point for cut generation etc.
    sol_->setPrimal(x); 
    if (vio>etol_) {
      status_ = ProvenLocalInfeasible;
    } else {
      sol_->setObjValue(obj_value);
      mynlp_->setSolution(sol_);
      status_ = ProvenLocalOptimal;
    }
    delete [] x;
  }
  return should_stop;
}


void IpoptEngine::removeCons(std::vector<ConstraintPtr> &)
{
  consChanged_ = true;
}


void IpoptEngine::resetIterationLimit()
{
  myapp_->Options()->SetIntegerValue("max_iter", maxIterLimit_);
}


void IpoptEngine::setIterationLimit(int limit)
{
  if (limit<1) {
    limit = maxIterLimit_;
  }
  myapp_->Options()->SetIntegerValue("max_iter", limit);
}


void IpoptEngine::setOptionsForProb_()
{
  if (problem_->isQP() || problem_->isLinear()) {
    myapp_->Options()->SetStringValue("hessian_constant","yes", true, true);
    myapp_->Options()->SetStringValue("jac_c_constant","yes", true, true);
    myapp_->Options()->SetStringValue("jac_d_constant","yes", true, true);
  } else if (problem_->getSize()->cons == problem_->getSize()->linCons) {
    myapp_->Options()->SetStringValue("jac_c_constant","yes", true, true);
    myapp_->Options()->SetStringValue("jac_d_constant","yes", true, true);
  } else {
    myapp_->Options()->SetStringValue("hessian_constant","no", true, true);
    myapp_->Options()->SetStringValue("jac_c_constant","no", true, true);
    myapp_->Options()->SetStringValue("jac_d_constant","no", true, true);
  }
}


void IpoptEngine::setOptionsForRepeatedSolve()
{
  if (myapp_) {
    myapp_->Options()->SetIntegerValue("print_level", 0);
    //myapp_->Options()->SetNumericValue("tol", 1e-7);
    //myapp_->Options()->SetIntegerValue("max_iter", 30);
    //myapp_->Options()->SetStringValue("mu_strategy", "adaptive");
    //myapp_->Options()->SetStringValue("output_file", "ipopt.out");
    // myapp_->Options()->SetStringValue("hessian_approximation", "limited-memory");
    //myapp_->Initialize("");
    //myapp_->Options()->SetNumericValue("dual_inf_tol", 1e-7);    
    //myapp_->Options()->SetNumericValue("constr_viol_tol", 1e-7);
    //myapp_->Options()->SetNumericValue("compl_inf_tol", 1e-12);
    //myapp_->Options()->SetNumericValue("gamma_phi", 1e-8, true, true);
    //myapp_->Options()->SetNumericValue("gamma_theta", 1e-4, true, true);
    // myapp_->Options()->SetNumericValue("required_infeasibility_reduction",
    // 0.3, true, true);
    myapp_->Options()->SetStringValue("expect_infeasible_problem","yes", true, 
                                      true);
    //myapp_->Options()->SetNumericValue("expect_infeasible_problem_ctol",
    //1e-1, true, true);
    //myapp_->Options()->SetNumericValue("expect_infeasible_problem_ytol",1e-1, 
    //true, true);
    myapp_->Options()->SetStringValue("mu_strategy", "adaptive", true, true);
    myapp_->Options()->SetStringValue("mu_oracle","probing", true, true);
    //myapp_->Options()->SetStringValue("fixed_variable_treatment", 
    //"make_parameter");
    //myapp_->Options()->SetNumericValue("bound_relax_factor", 1e-5);
  }
}


void IpoptEngine::setOptionsForSingleSolve()
{
  if (myapp_) {
    // myapp_->Options()->SetStringValue("hessian_approximation", "limited-memory");
    myapp_->Options()->SetIntegerValue("print_level", 0);
    myapp_->Options()->SetStringValue("expect_infeasible_problem","yes", true, 
                                      true);
    myapp_->Options()->SetStringValue("mu_strategy", "adaptive", true, true);
    myapp_->Options()->SetStringValue("mu_oracle","probing", true, 
                                      true);
    //myapp_->Options()->SetNumericValue("required_infeasibility_reduction",
    //0.3, true, true);
  }
}


EngineStatus IpoptEngine::solve()
{
  Ipopt::ApplicationReturnStatus status = Ipopt::Internal_Error; 
  bool should_stop;

  stats_->calls += 1;
  if (!(bndChanged_ || consChanged_ || true)) {
    return status_;
  }
  // Check if warm start is enabled. If so, load the information and
  // resolve. Otherwise solve from scratch.
  if (useWs_ &&
      ((false==strBr_ && ws_ && ws_->hasInfo()) ||
       (strBr_ && wsSb_ && wsSb_->hasInfo()))) {
    myapp_->Options()->SetStringValue("warm_start_init_point", "yes");
    if (true == strBr_) {
      mynlp_->copySolution(wsSb_->getPoint());
    } else {
      mynlp_->setSolution(ws_->getPoint());
    }
  } else {
    myapp_->Options()->SetStringValue("warm_start_init_point", "no");
  }

  if (consChanged_) {
    problem_->prepareForSolve(); // reset hessian and jacobian structures.
    setOptionsForProb_();
  }
  // Ask Ipopt to solve the problem. In order to do that, one should convert
  // mynlp_ from the derived class IpoptFunInterface (that we define) to
  // Ipopt::TNLP
  status_ = EngineUnknownStatus;
  //Ipopt::SmartPtr<Ipopt::TNLP> base = Ipopt::SmartPtr<Ipopt::TNLP>
  //(&(*mynlp_));
  timer_->start();
  should_stop = presolve_();
  stats_->ptime += timer_->query();
  if (should_stop) {
    if (status_== ProvenLocalOptimal) {
      status = Ipopt::Solve_Succeeded;
    } else if (status_ == ProvenLocalInfeasible) {
      status = Ipopt::Infeasible_Problem_Detected;
    }
  } else if (((ws_ && ws_->hasInfo()) || (strBr_ && wsSb_ && wsSb_->hasInfo()))
             && false==justLoaded_) {
    status = myapp_->ReOptimizeTNLP(mynlp_);
    justLoaded_ = false;
    stats_->reopt += 1;
  } else {
    status = myapp_->OptimizeTNLP(mynlp_);
    justLoaded_ = false;
    stats_->opt += 1;
  }

#if SPEW
  logger_->msgStream(LogDebug) << me_ << "time taken = " << timer_->query() 
    << std::endl;
  logger_->msgStream(LogDebug) << me_ << "Ipopt's status = " << status 
    << std::endl;
#endif

  // See IpReturnCodes_inc.h for the list
  switch (status) {
   case Ipopt::Solve_Succeeded : 
     status_ = ProvenLocalOptimal;
     sol_ = mynlp_->getSolution();
     break;
   case Ipopt::Solved_To_Acceptable_Level :
     status_ = ProvenLocalOptimal;
     sol_ = mynlp_->getSolution();
     break;
   case Ipopt::Infeasible_Problem_Detected :
     status_ = ProvenLocalInfeasible;
     break;
   case Ipopt::Maximum_Iterations_Exceeded :
     status_ = EngineIterationLimit;
     break;
   case Ipopt::Maximum_CpuTime_Exceeded :
     status_ = EngineIterationLimit;
     break;
   case Ipopt::Restoration_Failed :  // don't know what else to do.
     logger_->msgStream(LogInfo) << me_ << "restoration failed, "
       << "assuming local infeasible." << std::endl;
     status_ = ProvenLocalInfeasible;
     sol_ = mynlp_->getSolution();
     break;
   case Ipopt::Diverging_Iterates :
     status_ = ProvenUnbounded;
     break;
   case Ipopt::Search_Direction_Becomes_Too_Small :
     assert(!"Ipopt: search direction becomes too small.");
     break;
   case Ipopt::User_Requested_Stop:
     assert(!"Ipopt: user requested stop.");
     break;
   case Ipopt::Feasible_Point_Found:
     assert(!"Ipopt: feasible point found.");
     break;
   case Ipopt::Error_In_Step_Computation:
   case Ipopt::Not_Enough_Degrees_Of_Freedom:
   case Ipopt::Invalid_Problem_Definition:
   case Ipopt::Invalid_Option:
   case Ipopt::Invalid_Number_Detected:
   case Ipopt::Unrecoverable_Exception:
   case Ipopt::NonIpopt_Exception_Thrown:
   case Ipopt::Insufficient_Memory:
   case Ipopt::Internal_Error:
   default:
     logger_->msgStream(LogNone) << me_ << "error reported." << std::endl;
     status_ = EngineError;
  }
  if (prepareWs_) {
    // save warm start information
    ws_->setPoint(sol_);
  }

  Ipopt::SmartPtr<Ipopt::SolveStatistics> stats = myapp_->Statistics();
  // sometimes, (e.g. when all variables are fixed) ipopt does not solve
  // anything returns an empty Statistics object.
  UInt iters = (Ipopt::IsValid(stats)) ? stats->IterationCount() : 0;
#if SPEW
  logger_->msgStream(LogDebug) << me_ << "solve number = " << stats_->calls
    << std::endl;
  logger_->msgStream(LogDebug) << me_ << "number of iterations = " << iters 
    << std::endl;
  logger_->msgStream(LogDebug) << me_ << "status = " << getStatusString() 
    << std::endl;
  logger_->msgStream(LogDebug) << me_ << "obj = ";
  if (sol_) {
    logger_->msgStream(LogDebug) << mynlp_->getSolutionValue() << std::endl;
  } else {
    logger_->msgStream(LogDebug) << 1e40 << std::endl;
  }
#endif 
  if (true == strBr_) {
    stats_->strCalls += 1;
    stats_->strTime  += timer_->query();
    stats_->strIters += iters;
  } 
  stats_->time  += timer_->query();
  stats_->iters += iters;
  timer_->stop();

  bndChanged_ = false;
  consChanged_ = false;
  return status_;
}
	  

void IpoptEngine::fillStats(std::vector<double> &nlpStats)
{
  if (nlpStats.size()) {
    nlpStats[0] += stats_->calls;
    nlpStats[1] += stats_->opt;
    nlpStats[2] += stats_->reopt;
    nlpStats[3] += stats_->strCalls;
    nlpStats[4] += stats_->time;
    nlpStats[5] += stats_->ptime;
    nlpStats[6] += stats_->strTime;
    nlpStats[7] += stats_->iters;
    nlpStats[8] += stats_->strIters;
  }
}


void IpoptEngine::writeStats(std::ostream &out) const 
{
  if (stats_) {
    std::string me = "Ipopt: ";
    out << me << "total calls            = " << stats_->calls << std::endl
      << me << "calls to Optimize      = " << stats_->opt << std::endl
      << me << "calls to ReOptimize    = " << stats_->reopt << std::endl
      << me << "strong branching calls = " << stats_->strCalls << std::endl
      << me << "total time in solving  = " << stats_->time  << std::endl
      << me << "total time in presolve = " << stats_->ptime  << std::endl
      << me << "time in str branching  = " << stats_->strTime << std::endl
      << me << "total iterations       = " << stats_->iters << std::endl
      << me << "strong br iterations   = " << stats_->strIters << std::endl;
  }
}

} // namespace Minotaur

// 
// ############# IpOptFunInterface ###################
//
  
/* Constructor. */
namespace Ipopt{

IpoptFunInterface::IpoptFunInterface(Minotaur::EnvPtr env,
                                     Minotaur::ProblemPtr problem, 
                                     Minotaur::IpoptSolPtr sol)
: bOff_(1e-9),
  bTol_(1e-6),
  problem_(problem),
  sol_(sol)
{
  evalWithinBnds_ = env->getOptions()->findBool("eval_within_bnds")->getValue();
  logger_ = env->getLogger();
}


IpoptFunInterface::~IpoptFunInterface()
{
  if (sol_) {
    sol_ = 0;
  }
  if (problem_) {
    problem_ = 0;
  }
}

void IpoptFunInterface::copySolution(Minotaur::IpoptSolPtr sol)
{
  sol_->setPrimal(sol->getPrimal());
}

#ifdef NDEBUG
bool IpoptFunInterface::get_bounds_info(Index , Number* x_l, Number* x_u, 
                                        Index , Number* g_l, Number* g_u)
#else
bool IpoptFunInterface::get_bounds_info(Index n, Number* x_l, Number* x_u, 
                                        Index m, Number* g_l, Number* g_u)
#endif
{
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  double u,l;
  assert(n == (int) problem_->getNumVars());
  assert(m == (int) problem_->getNumCons());

  // variable bounds
  Minotaur::VariablePtr vPtr;
  Minotaur::VariableConstIterator vIter;

  for (vIter=problem_->varsBegin(); vIter!=problem_->varsEnd(); 
      ++vIter,++x_l,++x_u) {
    vPtr = *vIter;
    l = vPtr->getLb();
    u = vPtr->getUb();
    if (fabs(u-l)<bTol_) {
      u=l+bOff_;
    }
    *x_l = l;
    *x_u = u;
  }

  // constraint bounds
  Minotaur::ConstraintPtr cPtr;
  Minotaur::ConstraintConstIterator cIter;
  for (cIter=problem_->consBegin(); cIter!=problem_->consEnd(); 
      ++cIter,++g_l,++g_u) {
    cPtr = *cIter;
    *g_l = cPtr->getLb();
    *g_u = cPtr->getUb();
  }

  return true;
}


bool IpoptFunInterface::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                     Index& nnz_h_lag, 
                                     IndexStyleEnum& index_style)
{
  n = problem_->getNumVars();
  m = problem_->getNumCons();

  // nonzeros in jacobian. 
  nnz_jac_g = problem_->getNumJacNnzs();

  // ... and in hessian. 
  nnz_h_lag = problem_->getNumHessNnzs(); 

  // We use the standard c index style for row/col entries
  index_style = TNLP::C_STYLE;

  return true; 
}


#ifdef NDEBUG
bool IpoptFunInterface::get_starting_point(Index n, bool , Number* x,
                                           bool init_z, Number* z_L, 
                                           Number* z_U, Index m, 
                                           bool init_lambda, Number* lambda)
#else
bool IpoptFunInterface::get_starting_point(Index n, bool init_x, Number* x,
                                           bool init_z, Number* z_L, 
                                           Number* z_U, Index m, 
                                           bool init_lambda, Number* lambda)
#endif
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the others if
  // you wish.

  assert(init_x == true);

  if (init_z == false || init_lambda == false) {
    double* xp = x;
#if SPEW
    logger_->msgStream(Minotaur::LogDebug2) << "IpoptFunInterface: "
                                  << "setting initial point using the values "
                                  << " available in the Problem object."
                                  << std::endl;
#endif
    for (Minotaur::VariableConstIterator vit=problem_->varsBegin();
         vit!=problem_->varsEnd(); ++vit, ++xp) {
      *xp = (*vit)->getInitVal();
    }
  } else {
    const double *initial_point;
#if SPEW
    logger_->msgStream(Minotaur::LogDebug2) << "IpoptFunInterface: "
                                  << "setting initial point using the values "
                                  << " available from previous solve."
                                  << std::endl;
#endif
    // we should have some starting information in sol_.
    assert(sol_);
    initial_point = sol_->getPrimal();
    if (initial_point) {
      std::copy(initial_point, initial_point+n, x);
    } else {
      std::fill(x, x+n, 0);
    }

    initial_point = sol_->getDualOfCons();
    if (initial_point) {
    std::copy(initial_point, initial_point+m, lambda);
    } else {
      std::fill(lambda, lambda+m, 0);
    }

    initial_point = sol_->getLowerDualOfVars();
    if (initial_point) {
      std::copy(initial_point, initial_point+n, z_L);
    } else {
      std::fill(z_L, z_L+n, 0);
    }

    initial_point = sol_->getUpperDualOfVars();
    if (initial_point) {
    std::copy(initial_point, initial_point+n, z_U);
    } else {
      std::fill(z_U, z_U+n, 0);
    }
  }

  return true;
}


bool IpoptFunInterface::eval_f(Index n, const Number* x, bool, Number& obj_value)
{
  int error = 0;
  double *ex = 0; // it will be assigned to xx if evalWithinBnds_ is true,
                  // otherwise to x.
  double *xx = 0;

  if (evalWithinBnds_) {
    xx = pullXToBnds_(x, n);
    ex = xx;
  } else {
    ex = (double *) x;
  }
  // return the value of the objective function
  obj_value = problem_->getObjValue(ex, &error);
  return (0==error);
}


bool IpoptFunInterface::eval_g(Index n, const Number* x, bool, Index, Number* g)
{
  // return the value (activity) of the constraints: g(x)
  Minotaur::ConstraintConstIterator cIter;
  Minotaur::ConstraintPtr cPtr;
  Minotaur::UInt i=0;
  int error = 0, e = 0;
  double *ex = 0; // it will be assigned to xx if evalWithinBnds_ is true,
                  // otherwise to x.
  double *xx = 0;

  if (evalWithinBnds_) {
    xx = pullXToBnds_(x, n);
    ex = xx;
  } else {
    ex = (double *) x;
  }
  for (cIter=problem_->consBegin(); cIter!=problem_->consEnd(); ++cIter) {
    cPtr = *cIter;
    e = 0;
    g[i] = cPtr->getActivity(ex, &e);
    if (e!=0) {
      error=1;
#if SPEW
      logger_->msgStream(Minotaur::LogDebug)
        << "IpoptFunInterface: error in evaluating constraint\n";
      cPtr->write(logger_->msgStream(Minotaur::LogDebug));
#endif
    }
    ++i;
  }

  if (xx) {
    delete [] xx ;
  }

  return (0==error);
}


bool IpoptFunInterface::eval_grad_f(Index n, const Number* x, bool, 
                                    Number* grad_f)
{
  // return the gradient of the objective function grad_{x} f(x)

  int error = 0;
  Minotaur::ObjectivePtr o;
  double *ex = NULL; // it will be assigned to xx if evalWithinBnds_ is true,
                     // otherwise to x.
  double *xx = NULL;

  std::fill(grad_f, grad_f+n, 0);

  if (evalWithinBnds_) {
    xx = pullXToBnds_(x, n);
    ex = xx;
  } else {
    ex = (double *) x;
  }

  o = problem_->getObjective();
  if (o) {
    o->evalGradient(ex, (double *)grad_f, &error);
  }

  if (xx) {
    delete [] xx ;
  }

  //for (int i=0; i<n; ++i) {
  //  std::cout << "grad obj [" << i << "] = " << grad_f[i] << std::endl;
  //}
  return (0==error);
}


bool IpoptFunInterface::eval_h(Index n, const Number* x, bool, Number obj_factor, 
                               Index, const Number* lambda, bool, Index,
                               Index* iRow, Index* jCol, Number* values)
{
  int error = 0;
  double *ex = NULL; // it will be assigned to xx if evalWithinBnds_ is true,
                     // otherwise to x.
  double *xx = NULL;

  if (x==0 && lambda==0 && values==0) {
    problem_->getHessian()->fillRowColIndices((Minotaur::UInt *)iRow,
        (Minotaur::UInt *)jCol);
  } else if (x!=0 && lambda!=0 && values!=0) {
    if (evalWithinBnds_) {
      xx = pullXToBnds_(x, n);
      ex = xx;
    } else {
      ex = (double *) x;
    }
    problem_->getHessian()->fillRowColValues(ex, 
        (double) obj_factor, (double *) lambda, 
        (double *)values, &error);
    //std::cout << "error = " << error << std::endl;
    //for (int i=0; i<problem_->getNumVars(); ++i) {
    //  std::cout << std::setprecision(8) << "x["<<i<<"] = "<<x[i] << std::endl;
    //}
    //for (int i=0; i<problem_->getHessian()->getNumNz(); ++i) {
    //  std::cout << std::setprecision(8) << "h["<<i<<"] = "<<values[i] << std::endl;
    //}
    if (xx) {
      delete [] xx ;
    }
  } else {
    assert (!"one of x, lambda and values is NULL!");
  }
  return (0==error); 
}


bool IpoptFunInterface::eval_jac_g(Index n, const Number* x, bool, Index, 
                                   Index, Index* iRow, Index *jCol, 
                                   Number* values)
{
  int error = 0;
  double *ex = NULL; // it will be assigned to xx if evalWithinBnds_ is true,
                     // otherwise to x.
  double *xx = NULL;


  if (values == 0) {
    // return the structure of the jacobian of the constraints
    problem_->getJacobian()->fillRowColIndices((Minotaur::UInt *) iRow,
        (Minotaur::UInt *) jCol);
  }
  else {
    if (evalWithinBnds_) {
      xx = pullXToBnds_(x, n);
      ex = xx;
    } else {
      ex = (double *) x;
    }
    // return the values of the jacobian of the constraints
    problem_->getJacobian()->fillRowColValues((double *) ex, 
        (double *) values, &error);
    if (xx) {
      delete [] xx ;
    }
  }
  if (error!=0) {
    logger_->msgStream(Minotaur::LogError) << "IpoptFunInterface: " 
      << "error in evaluating jacobian. Setting eval_within_bnds to True."
      << std::endl;
    evalWithinBnds_ = true;
  }

  return (0==error);
}


void IpoptFunInterface::finalize_solution(SolverReturn, Index, const Number* x, 
                                          const Number* z_L, const Number* z_U,
                                          Index, const Number*,
                                          const Number* lambda,
                                          Number obj_value, const IpoptData*, 
                                          IpoptCalculatedQuantities* /* q */)
{
  sol_->setPrimal(x);
  sol_->setObjValue(obj_value);
  sol_->setDualOfCons(lambda);
  sol_->setDualOfVars(z_L, z_U);
  //std::cout << "vio = " << q->curr_constraint_violation() << std::endl;
}


double IpoptFunInterface::getSolutionValue() const
{
  return sol_->getObjValue();
}


double* IpoptFunInterface::pullXToBnds_(const Number* x, Index n)
{
  double *xx = new double[n];
  for (int ii=0; ii<n; ++ii) {
    if (x[ii]<problem_->getVariable(ii)->getLb()) {
      xx[ii] = problem_->getVariable(ii)->getLb();
    } else if (x[ii]>problem_->getVariable(ii)->getUb()) {
      xx[ii] = problem_->getVariable(ii)->getUb();
    } else {
      xx[ii] = x[ii];
    }
  }
  return xx;
}

} // namespace Ipopt

// Local Variables: 
// mode: c++ 
// eval: (c-set-style "k&r") 
// eval: (c-set-offset 'innamespace 0) 
// eval: (setq c-basic-offset 2) 
// eval: (setq fill-column 78) 
// eval: (auto-fill-mode 1) 
// eval: (setq column-number-mode 1) 
// eval: (setq indent-tabs-mode nil) 
// End:
