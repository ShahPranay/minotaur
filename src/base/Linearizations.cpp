// 
//     MINOTAUR -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2017 The MINOTAUR Team.
// 

/** 
 * \file Linearizations.cpp
 * \Briefly define a class for adding linearizations in linearization
 * based methods. Added for problems with nonlinear constraints
 * \Author Meenarli Sharma, Indian Institute of Technology Bombay
 */


#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>

#include "MinotaurConfig.h"

#include "CNode.h"
#include "Constraint.h"
#include "Engine.h"
#include "Environment.h"
#include "Function.h"
#include "Logger.h"
#include "Node.h"
#include "NonlinearFunction.h"
#include "QuadraticFunction.h"
#include "Objective.h"
#include "Operations.h"
#include "Option.h"
#include "ProblemSize.h"
#include "Linearizations.h"
#include "Relaxation.h"
#include "Solution.h"
#include "SolutionPool.h"
#include "VarBoundMod.h"
#include "Variable.h"

using namespace Minotaur;

typedef std::vector<ConstraintPtr>::const_iterator CCIter;
const std::string Linearizations::me_ = "Linearizations: ";

Linearizations::Linearizations(EnvPtr env, RelaxationPtr rel,
                               ProblemPtr minlp, std::vector<ConstraintPtr> nlCons)
: env_(env),
  rel_(rel),
  minlp_(minlp),
  //lpe_(EnginePtr()),
  nlpe_(EnginePtr()),
  solC_(NULL),
  nlpx_(NULL),
  varPtrs_(0)
{
  nlCons_ = nlCons;
  logger_ = env->getLogger();
  //MS: set the option for root_LinSchemes
  rs1_ = env_->getOptions()->findDouble("root_linScheme1")->getValue();
  rs2Per_ = env_->getOptions()->findDouble("root_linScheme2_per")->getValue();
  rs2NbhSize_ = env_->getOptions()->findDouble("root_linScheme2_nbhSize")->getValue();
  rs3_ = env_->getOptions()->findInt("root_linScheme3")->getValue();
  rgs1_ = env_->getOptions()->findBool("root_genLinScheme1")->getValue();
  //rgs2_ = env_->getOptions()->findBool("root_genLinScheme2")->getValue();
  rsg2Per_ = env_->getOptions()->findDouble("root_linGenScheme2_per")->getValue();
  intTol_ = env_->getOptions()->findDouble("int_tol")->getValue();
  solAbsTol_ = env_->getOptions()->findDouble("feasAbs_tol")->getValue();
  solRelTol_ = env_->getOptions()->findDouble("feasRel_tol")->getValue();
  objATol_ = env_->getOptions()->findDouble("solAbs_tol")->getValue();
  objRTol_ = env_->getOptions()->findDouble("solRel_tol")->getValue();

  stats_ = new LinStats();
  stats_->rs1Cuts = 0;
  stats_->rs2Cuts = 0;
  stats_->rs3Cuts = 0;
  stats_->rgs1Cuts = 0;
  stats_->rgs2Cuts = 0;
 }


Linearizations::~Linearizations()
{ 
  if (stats_) {
    delete stats_;
  }
  if (solC_) {
    delete [] solC_;
    solC_ = 0;
  }
  env_ = 0;
  rel_ = 0;
  minlp_ = 0;
  nlCons_.clear();
}


bool Linearizations::addCutAtRoot_(double *x, ConstraintPtr con, UInt &newConId)
{
  int error = 0; 
  FunctionPtr f;
  double c, act, cUb;
  ConstraintPtr newcon;
  std::stringstream sstm;
  LinearFunctionPtr lf = LinearFunctionPtr();

  act = con->getActivity(x, &error);
  if (error == 0) {
    f = con->getFunction();
    linearAt_(f, act, x, &c, &lf, &error);
    if (error == 0) {
      cUb = con->getUb();
      ++(stats_->rs1Cuts);
      sstm << "_OAcut_" << stats_->rs1Cuts << "_AtRoot";
      f = (FunctionPtr) new Function(lf);
      newcon = rel_->newConstraint(f, -INFINITY, cUb-c, sstm.str());
      newConId = newcon->getIndex();
      sstm.str("");
      return true;
    }
  }	else {
    logger_->msgStream(LogError) << me_ << "Constraint" <<  con->getName() <<
      " is not defined at this point." << std::endl;
  }
  return false;
}

//MS: delete unused functions
bool Linearizations::cutAtLineSearchPt_(const double *xIn, const double *xOut,
                                        double* xNew, ConstraintPtr con)
{
  double nlpact;
  bool lsPtFound;
  lsPtFound = lineSearchPt_(xIn, xOut, xNew, con, nlpact);
  if (lsPtFound) {
    int error = 0;
    std::stringstream sstm;
    //ConstraintPtr newcon;
    LinearFunctionPtr lf = 0;
    FunctionPtr f = con->getFunction();
    double c, cUb = con->getUb();
    if (error == 0) {
      linearAt_(f, nlpact, xNew, &c, &lf, &error);
      ++(stats_->rs3Cuts);
      sstm << "_OAcut_" << stats_->rs3Cuts;
      f = (FunctionPtr) new Function(lf);
      rel_->newConstraint(f, -INFINITY, cUb-c, sstm.str());
      //newcon = rel_->newConstraint(f, -INFINITY, cUb-c, sstm.str());
    }	else {
      return false;
    }
   } else {
     return false;   
   }
  return true;
}


bool Linearizations::linPart_(double *b1, UInt lVarIdx, ConstraintPtr con,
                           double lVarCoeff, double act)
{
  int error = 0;
  double nlTerm = 0; 
  QuadraticFunctionPtr qf = con->getQuadraticFunction();
  NonlinearFunctionPtr nlf = con->getNonlinearFunction();

  if (nlf) {
    nlTerm = nlf->eval(b1, &error); 
  } 

  if (error == 0) {
    if (qf) {
      nlTerm = nlTerm + qf->eval(b1); 
    }
    b1[lVarIdx] = (con->getUb() - nlTerm - act)/lVarCoeff;    
    return true;
  }
  return false;
}


bool Linearizations::addNewCut_(double *b1, ConstraintPtr con,
                           UInt &newConId)
{
  bool found;  
  found = addCutAtRoot_(b1, con, newConId);
  if (found) {
    return true;
  } 
  return false;
}


void Linearizations::findCenter()
{
  // Center is found if the feasible region is compact and has
  // non-empty inetrior. Otherwise, not.
  double lb, ub;
  bool NLPFeas = false;
  VariablePtr vPtr, v;
  ConstraintPtr con;
  FunctionType fType;
  std::vector<ConstraintPtr > cp;
  ProblemPtr inst_C = minlp_->clone();
  //inst_C->write(std::cout);
  //UInt numVars = minlp_->getNumVars();
  double *sol1;

  FunctionPtr fnewc;
  LinearFunctionPtr lfc = (LinearFunctionPtr) new LinearFunction();
  vPtr = inst_C->newVariable(-INFINITY, 0, Continuous, "eta", VarHand);
  vPtr->setFunType_(Nonlinear);
  inst_C->removeObjective();
  lfc->addTerm(vPtr, 1.0);
  fnewc = (FunctionPtr) new Function(lfc);
  inst_C->newObjective(fnewc, 0.0, Minimize);

  for (ConstraintConstIterator it=inst_C->consBegin(); it!=inst_C->consEnd();
     ++it) {
    con = *it;
    lb = con->getLb();
    ub = con->getUb();
    fType = con->getFunctionType();
    if (fType == Constant || fType == Linear) {
      continue;
    } else {
      if (con->getLinearFunction()) {
        lfc = con->getLinearFunction()->clone();
        lfc->addTerm(vPtr, -1.0);
      } else {
        lfc = (LinearFunctionPtr) new LinearFunction();
        lfc->addTerm(vPtr, -1.0);
      }
    }
    inst_C->changeConstraint(con, lfc, lb, ub);
  }
  
  //inst_C->write(std::cout);
  inst_C->prepareForSolve();
  nlpe_->load(inst_C);
  solveNLP_();

  if (solC_) {
    NLPFeas = true;
    sol1 = new double[minlp_->getNumVars()];
    std::copy(solC_, solC_+minlp_->getNumVars(), sol1);
  } 
//else {
    //return;
  //}
  
  
 // Solving more restricted proiblem to find center 
  for (ConstraintConstIterator it=inst_C->consBegin(); it!=inst_C->consEnd();
     ++it) {
    con = *it;
    lb = con->getLb();
    ub = con->getUb();
    fType = con->getFunctionType();
    if (fType == Constant) {
      inst_C->markDelete(con);
      continue;
    } else if (fType == Linear)  {
      if (lb != -INFINITY && ub != INFINITY) {
        if (lb == ub) {
          continue;       
        }
        cp.push_back(con);
        inst_C->markDelete(con);
        continue;
      } else if (lb != -INFINITY) {
        ub = INFINITY;
        lfc = con->getLinearFunction()->clone();
        lfc->addTerm(vPtr, 1.0);
      } else if (ub != INFINITY ) {
        lb = -INFINITY;
        lfc = con->getLinearFunction()->clone();
        lfc->addTerm(vPtr, -1.0);
      } else {
        inst_C->markDelete(con);
        continue;
      }
    } else {
      continue;
    }
    inst_C->changeConstraint(con, lfc, lb, ub);
  }  

  for (UInt i = 0; i < cp.size(); ++i) {
    con = cp[i];
    lb = con->getLb(), ub = con->getUb();
    lfc = con->getLinearFunction()->clone();
    lfc->addTerm(vPtr, 1.0);
    fnewc = (FunctionPtr) new Function(lfc);
    inst_C->newConstraint(fnewc, lb, INFINITY);

    lfc = con->getLinearFunction()->clone();
    lfc->addTerm(vPtr, -1.0);
    fnewc = (FunctionPtr) new Function(lfc);
    inst_C->newConstraint(fnewc, -INFINITY, ub);
  }
  cp.clear();
  inst_C->delMarkedCons();

  for (VariableConstIterator vit=inst_C->varsBegin(); vit!=inst_C->varsEnd()-1;
       ++vit) {
    v = *vit;
    lb = v->getLb(), ub = v->getUb();
    if (lb == ub) {
      continue;
    }

    if (lb != -INFINITY) {
      lfc = (LinearFunctionPtr) new LinearFunction();
      lfc->addTerm(vPtr, 1.0);
      lfc->addTerm(v, 1.0);
      fnewc = (FunctionPtr) new Function(lfc);
      inst_C->newConstraint(fnewc, lb, INFINITY);
    }
    
    if (ub != INFINITY) {
      lfc = (LinearFunctionPtr) new LinearFunction();
      lfc->addTerm(vPtr, -1.0);
      lfc->addTerm(v, 1.0);
      fnewc = (FunctionPtr) new Function(lfc);
      inst_C->newConstraint(fnewc, -INFINITY, ub);
    }
  }
  inst_C->prepareForSolve();
  
  nlpe_->clear();
  nlpe_->load(inst_C);
  solveNLP_();
  
  //inst_C->write(std::cout);
  if (solC_ == 0) {
    if (NLPFeas) {
      solC_ = new double[minlp_->getNumVars()];
      std::copy(sol1, sol1 + minlp_->getNumVars(), solC_);
    }
  }
  
  if (NLPFeas) {
    delete [] sol1;
    sol1 = 0;
  }
  
  delete nlpe_;
  nlpe_ = 0;
  delete inst_C;
  inst_C = 0;
  return;
}

void Linearizations::solveNLP_()
{ 
  EngineStatus nlpStatus = nlpe_->solve();

  switch(nlpStatus) {
  case (ProvenOptimal):
  case (ProvenLocalOptimal):
    if (solC_) {
      delete [] solC_;
      solC_ = 0;
    }
    //std::cout << "Center " << std::setprecision(6) << nlpe_->getSolution()->getObjValue();
    //exit(1);
    if (fabs(nlpe_->getSolution()->getObjValue()) > solAbsTol_) {
      const double *temp = nlpe_->getSolution()->getPrimal();
      solC_ = new double[minlp_->getNumVars()];
      std::copy(temp, temp+minlp_->getNumVars(), solC_);
    } 
    break;
  case (EngineIterationLimit):
  case (ProvenInfeasible):
  case (ProvenLocalInfeasible): 
  case (ProvenObjectiveCutOff):
    break;
  case (FailedFeas):
  case (EngineError):
  case (FailedInfeas):
  case (ProvenUnbounded):
  case (ProvenFailedCQFeas):
  case (EngineUnknownStatus):
  case (ProvenFailedCQInfeas):
  default:
    logger_->msgStream(LogError) << me_ << "NLP engine status = "
      << nlpe_->getStatusString() << std::endl;
    break;
  }

  return;
}


bool Linearizations::findIntersectPt_(std::vector<UInt > newConsId, VariablePtr vl,
                                VariablePtr vnl, double * iP)
{
  ConstraintPtr con;
  LinearFunctionPtr lf;
  con = rel_->getConstraint(newConsId[0]);
  lf = con->getLinearFunction();
  //std::cout << "linear and nonlinear  pointer " << vl << " " << vnl << std::endl;
  double a = lf->getWeight(vl);
  double b = lf->getWeight(vnl);
  double e = con->getUb();
  
  con = rel_->getConstraint(newConsId[1]);
  lf = con->getLinearFunction();
  double c = lf->getWeight(vl);
  double d = lf->getWeight(vnl);
  double f = con->getUb();

  /* we solve the linear system
   * ax+by=e
   * cx+dy=f
   * where, x is iP[1] and y is iP[0]
   */
  double determinant = a*d - b*c; 
  if(determinant != 0) {
    iP[1] = (e*d - b*f)/determinant;
    iP[0] = (a*f - e*c)/determinant;
  } else {
    std::cout << "Cramer equations system: determinant is zero\n"
              "there are either no solutions or many solutions exist.\n";
    return false; 
  }
  return true;
}


void Linearizations::insertNewPt_(UInt j, UInt i, std::vector<double > & xc,
                             std::vector<double> & yc, ConstraintPtr newcon,
                             VariablePtr vl, VariablePtr vnl, bool & shouldCont)
{
  double f = newcon->getUb();
  LinearFunctionPtr lf = newcon->getLinearFunction();
  
  double d = lf->getWeight(vl);
  double c = lf->getWeight(vnl);
  double x1 = xc[j], y1 = yc[j], x2 = xc[i], y2 = yc[i], x, y;

  // point of intersection of newcon with the lin from j and j-1
  double a = y1-y2;
  double b = x2-x1;
  double e = y1*(x2-x1) - x1*(y2-y1);
  double determinant = a*d - b*c; 
  if(determinant != 0) {
    x = (e*d - b*f)/determinant;
    y = (a*f - e*c)/determinant;
    xc.insert(xc.begin()+j,x);
    yc.insert(yc.begin()+j,y);
  } else {
    std::cout << "Cramer equations system: determinant is zero\n"
              "there are either no solutions or many solutions exist.\n"; 
    shouldCont = false;
  }
}


std::vector<UInt > Linearizations::isFeas_(double *x, 
                                           std::vector<UInt > varConsPos,
                                           bool &foundActive, bool &foundVio)
{
  UInt j;
  int error = 0;
  ConstraintPtr c;
  double act, cUb;
  std::vector<UInt > consPos;

  for (UInt i = 0; i < varConsPos.size(); ++i) {
    j = varConsPos[i];
    c = nlCons_[j];
    act = c->getActivity(x, &error);
    if (error == 0) {
      cUb = c->getUb();
      if ((act > cUb + solAbsTol_) &&
          (cUb == 0 || act > cUb + fabs(cUb)*solRelTol_)) {
        if (!foundVio) {
          foundVio = true;
          if (consPos.size() != 0) {
            consPos.clear();          
          }         
        }
        consPos.push_back(j);
      } else if ((fabs(act-cUb) <= solAbsTol_) ||
            (cUb != 0 && fabs(act- cUb) <= fabs(cUb)*solRelTol_)) {
        if (!foundVio) {
          if (!foundActive) {
            foundActive = true;         
            if (consPos.size() != 0) {
              consPos.clear();          
            }         
          }
          consPos.push_back(j);
        }
      } 
    }
  }
  return consPos;
}


void Linearizations::linearAt_(FunctionPtr f, double fval, const double *x,
                          double *c, LinearFunctionPtr *lf, int *error)
{
  int n = rel_->getNumVars();
  double *a = new double[n];
  VariableConstIterator vbeg = rel_->varsBegin(), vend = rel_->varsEnd();
  const double linCoeffTol =
    env_->getOptions()->findDouble("conCoeff_tol")->getValue();

  std::fill(a, a+n, 0.);
  f->evalGradient(x, a, error);
  
  if (*error==0) {
    *lf = (LinearFunctionPtr) new LinearFunction(a, vbeg, vend, linCoeffTol);
    *c  = fval - InnerProduct(x, a, minlp_->getNumVars());
  } else {
    logger_->msgStream(LogError) << me_ <<"gradient not defined at this point."
      << std::endl;
#if SPEW
    logger_->msgStream(LogDebug) << me_ <<"gradient not defined at this point."
      << std::endl;
#endif
  }
  delete [] a;
  return;
}


//MS: rootLinScheme3_ - individual cons - working
bool Linearizations::lineSearchPt_(const double *xIn, const double *xOut,
                                   double* x, ConstraintPtr con, double &nlpact)
{
  int error = 0;
  bool ptFound = true, isLS = true;
  double cUb = con->getUb();
  FunctionPtr f = con->getFunction();
  
  nlpact = f->eval(xOut, &error);
  ////MS: later try extending the point towards outside of the feasible
  //// region and adding cut there
  if (error == 0) {
    if (nlpact <= cUb+solAbsTol_ || (cUb != 0 && 
                                     nlpact <= cUb+fabs(cUb)*solRelTol_)) {
      isLS = false;
    }
  } else {
    return false;  
  }

  UInt numVars =  minlp_->getNumVars();
  double *xl, *xu;
  if (isLS) {
    xl = new double[numVars];
    xu = new double[numVars];
    std::copy(xIn,xIn+numVars,xl);
    std::copy(xOut,xOut+numVars,xu);
  } else {
    ptFound = false;
  }

  while (isLS) {
    //MS: all vars or only in the nonlinear ones
    for (UInt i = 0 ; i < minlp_->getNumVars(); ++i) {
      x[i] = 0.5*(xl[i] + xu[i]);
    }
    nlpact = f->eval(x, &error);
    if (error == 0) {
      if (nlpact > cUb + solAbsTol_ && (cUb==0 ||
                                      nlpact > cUb + fabs(cUb)*solRelTol_)) {
        std::copy(x,x+numVars,xu);
      } else {
        if (fabs(cUb-nlpact) <= solAbsTol_ || 
            (cUb!=0 && fabs(cUb-nlpact) <= fabs(cUb)*solRelTol_)) {
          break;
        } else {
          std::copy(x,x+numVars,xl);
        }
      }
    } else {
      //MS: think later what can be done here
      ptFound = false;
      break;
    }
  }

  if (isLS) {
    delete [] xl;
    delete [] xu;
  }
  return ptFound;
}



void Linearizations::rootLinearizations(const double * nlpx)
{
  ConstraintPtr con;
  bool isFound = false;
  UInt nVarIdx, lVarIdx;
  double lVarCoeff = 0, nVarCoeff = 0;
    
  nlpx_ = nlpx;
 
  if (rs1_ || rs2Per_) { 
    for (CCIter it = nlCons_.begin(); it != nlCons_.end(); ++it) {
      con = *it;
      lVarIdx = 0; lVarCoeff = 0; nVarCoeff = 0;
      // constraints with only one var in the nonlinear part
      isFound = uniVarNlFunc_(con, lVarCoeff, lVarIdx, nVarIdx, nVarCoeff); 
      //MS: see if this if-else has to be changed
      if (isFound == false) {
        continue;    
      } else {
        if (rs1_ > 0) {
          rootLinScheme1_(con, lVarCoeff, lVarIdx, nVarIdx, nVarCoeff);
        }
        if (rs2Per_ > 0) { // there is a default neighborhood
          //nlpx_ = nlpx;
          rootLinScheme2_(con, lVarCoeff, lVarIdx, nVarIdx);  
        }
      }
    }
    // To be written here for objective
    //if (oNl_) {
    
    //}
  }
  /// General scheme at root
  // Option for general scheme
  //isFound = false;
  if (rgs1_ || rsg2Per_) {
    //nlpe_->clear();
    //findCenter_();
    //findCenter_(isFound);
    if (solC_) {
      // populate varPtrs_ with index of variables in nonlinear constraints
      varsInNonlinCons_();
      if (varPtrs_.size() == 0) {
        return;      
      }
      // General scheme using center and  positive spanning vectors
      if (rgs1_) {
        rootLinGenScheme1_();
      }

      // General scheme using nonlinear solution and positive spanning vectors
      if (rsg2Per_) {
        rootLinGenScheme2_();
      }
    }
  }

  return;
}


void Linearizations::rootLinGenScheme2_()
{
  int error, firstnnz= -1;
  VariablePtr v;
  FunctionPtr f;
  ConstraintPtr con;
  //bool isFound = false;
  std::vector<UInt > varIdx;
  UInt fixIdx, vIdx;
  std::vector<double* > nlconsGrad;
  std::vector<UInt > varConsPos;
  std::vector<UInt > varPos;
  std::vector<int > alphaSign;
  double vbnd, rhs = 0.0, fixCoeff, coeff;
  int n = minlp_->getNumVars(), nr = rel_->getNumVars();
  
  double *xOut = new double[n]; 
  std::copy(nlpx_, nlpx_+n, xOut);
    
  // rhs of the plane passing through nlpx_ and has solC-nlpx_ as normal;
  // fixIdx is the first variable in varPtrs_ with nonzero coeff on this plane
  for (UInt i = 0; i < varPtrs_.size(); ++i) {
    v = varPtrs_[i];
    vIdx = v->getIndex();
    rhs = rhs + nlpx_[vIdx]*(solC_[vIdx] - nlpx_[vIdx]);
    if ((firstnnz == -1) && (solC_[vIdx] - nlpx_[vIdx]) != 0) {
      firstnnz = i;
      fixIdx = vIdx; //MS: is needed?
      //isFound = true;              
    }
  }

  //if (!isFound || rhs == 0) 
  if (firstnnz == -1 || rhs == 0) {
    return;  
  }
   
 //// Gradient of constraints at nonlinear solution nlpx_
  for (UInt i = 0; i < nlCons_.size(); ++i) {
    error = 0;
    con = nlCons_[i];
    double *grad = new double[nr];
    std::fill(grad, grad+nr, 0);
    
    f = con->getFunction();
    f->evalGradient(nlpx_, grad, &error);
    
    if (error == 0) {
      nlconsGrad.push_back(grad);
    } else {
      nlconsGrad.push_back(0);
    }
  }

  // coefficient of the fix var in any direction
  fixCoeff = rhs/(solC_[fixIdx] - nlpx_[fixIdx]);
     
  // linearly independent directions for each var
  for (UInt i = 0; i < varPtrs_.size(); ++i) {
    v = varPtrs_[i];
    vIdx = v->getIndex();
    if (int(i) < firstnnz) {
      // Coeff of var is zero in the hyperplane expression
      // Nonlinear constraints containing var v
      varPos.push_back(i);
      varCons_(varPos, varConsPos);
      varPos.clear();
      
      if (varConsPos.size() == 0) {
        continue;      
      }
      
      vbnd = v->getUb();
      varIdx.push_back(vIdx);
      alphaSign.push_back(1);
      exploreDir_(xOut, varIdx, varConsPos, nlconsGrad, alphaSign, vbnd);
      
      //// reverse the direction
      vbnd = v->getLb();
      alphaSign[0] = -1;
      exploreDir_(xOut, varIdx, varConsPos, nlconsGrad, alphaSign, vbnd);
    } else if (int(i) > firstnnz) {
      // in this case coeff of v could be nonzero in the hyperplane expression 
      // list of constraints containing var v
      varPos.push_back(i);
      varPos.push_back(firstnnz);
      varCons_(varPos, varConsPos);
      varPos.clear();
      
      if (varConsPos.size() == 0) {
        continue;      
      }
      
      if ((solC_[vIdx] - nlpx_[vIdx]) == 0) {
        // in this case coeff of v is zero in hyperplane expression 
        vbnd = v->getUb();
        varIdx.push_back(vIdx);
        alphaSign.push_back(1);
        exploreDir_(xOut, varIdx, varConsPos, nlconsGrad, alphaSign, vbnd);
        //reverse the direction
        vbnd = v->getLb();
        alphaSign[0] = -1;
        exploreDir_(xOut, varIdx, varConsPos, nlconsGrad, alphaSign, vbnd);
      } else {
        coeff = rhs/(solC_[vIdx] - nlpx_[vIdx]);
        boundingVar_(vbnd, i, varPtrs_[firstnnz], coeff,-fixCoeff, alphaSign, varIdx);
        exploreDir_(xOut, varIdx, varConsPos, nlconsGrad, alphaSign, vbnd);
        varIdx.clear();
        alphaSign.clear();
       
        //// reverse the direction
        boundingVar_(vbnd, i, varPtrs_[firstnnz], -coeff, fixCoeff, alphaSign, varIdx);
        exploreDir_(xOut, varIdx, varConsPos, nlconsGrad, alphaSign, vbnd);
      }
    } else {
      continue;    
    }
    varIdx.clear();
    alphaSign.clear();
    varConsPos.clear();
  }

  for (UInt i = 0; i < nlCons_.size(); ++i) {
    if (nlconsGrad[i]) {
      delete [] nlconsGrad[i];
      nlconsGrad[i] = 0;
    }
  }
  delete [] xOut;
  return;
}

void Linearizations::exploreDir_(double *xOut, std::vector<UInt > varIdx,
                                 std::vector<UInt > varConsPos, 
                                std::vector<double *> nlconsGrad,
                                std::vector<int> alpha,
                                double &vbnd)
{
  bool isFound = false;
  UInt idx = varIdx[0], idx1;
  int nr = rel_->getNumVars();
  std::vector<double* > lastGrad;
  std::vector<double > unitVec;
  lastGrad.resize(nlCons_.size(), 0); //MS: Can be made efficient.
  double delta = 0.5, val = fabs(vbnd - nlpx_[idx]), nbhSize = 20;    //MS: parametrize nbhSize;
 
  if (val < 1) {
    delta = fabs(vbnd - nlpx_[idx]);
  } else if (val > nbhSize) {
    vbnd = xOut[idx] + alpha[0]*nbhSize;
  }

  if (delta != 0) {
    val = 0;
    for (UInt i = 0; i < varIdx.size(); ++i) {
      idx1 = varIdx[i];
      val = val + pow(nlpx_[idx1], 2);   
    }
    val = sqrt(val);

    for (UInt i = 0; i < varIdx.size(); ++i) {
      idx1 = varIdx[i];
      unitVec.push_back(fabs(nlpx_[idx1])/val);
      xOut[idx1] = xOut[idx1] + alpha[i]*delta*unitVec[i];
    }
    
    for (UInt i = 0; i < varConsPos.size(); ++i) {
      idx1 = varConsPos[i];
      if (nlconsGrad[idx1]) {
        double * grad = new double[nr];
        std::copy(nlconsGrad[idx1], nlconsGrad[idx1]+nr, grad);
        lastGrad[idx1] = grad;
      }
    }

    while (alpha[0]*xOut[idx] <= alpha[0]*vbnd) {
      isFound = findLinPoint_(xOut, varConsPos, lastGrad);
      // double the stepsize delta if linearization at this point is not useful
      if (!isFound) {
        delta = 2*delta;
      }
      for (UInt i = 0; i < varIdx.size(); ++i) {
        idx1 = varIdx[i];
        xOut[idx1] = xOut[idx1] + alpha[i]*delta*unitVec[i];
      }
    }

    for (UInt i = 0; i < varConsPos.size(); ++i) {
      idx1 = varConsPos[i];   
      if (lastGrad[idx1]) {
        delete [] lastGrad[idx1];
        lastGrad[idx1] = 0;
      }
    }

    for (UInt i = 0; i < varIdx.size(); ++i) {
      idx1 = varIdx[i];
      xOut[idx1] = nlpx_[idx1];
    }
 
  }
  return;
} 


bool Linearizations::findLinPoint_(double *xOut, 
                                   std::vector<UInt > varConsPos,
                                   std::vector<double* > & lastGrad) 
{
  double lambda = 0.5;
  bool foundActive = false, foundVio = false, cutAdded = false;
  std::vector<UInt> vConsPos = isFeas_(xOut, varConsPos, foundActive, foundVio);
        
  if (foundVio) {
    // point outside feasible region. Do line search
    cutAdded = boundaryPt_(xOut, vConsPos, lastGrad);
  } else if (foundActive) {
    //point on the feasible region. Directly add linearizations
    cutAdded = genLin_(xOut, vConsPos, lastGrad);
  } else {
    double bnd;
    VariablePtr v;
    bool isCont = true;
    UInt n =  minlp_->getNumVars(), vIdx;
    
    double* x = new double[n];
    std::copy(solC_, solC_+n, x);
    
    // point inside feasible region. Move along xOut - x^C from x^C
    while (isCont) {
     
      for (UInt i = 0 ; i < n; ++i) {
        x[i] = x[i] + lambda*(xOut[i] - solC_[i]);
      }
      
      vConsPos = isFeas_(x, varConsPos, foundActive, foundVio);
      if (foundVio) {
        // point outside feasible region. Do line search
        cutAdded = boundaryPt_(x, vConsPos, lastGrad);
        break;
      } else if (foundActive) {
        //point on the feasible region. Directly add linearizations
        cutAdded = genLin_(x, vConsPos, lastGrad);
        break;
      }
      vConsPos.clear();

      for (UInt i = 0; i < varPtrs_.size(); ++i) {
        v = varPtrs_[i];
        vIdx = v->getIndex();
        bnd = v->getLb();
        if ((x[vIdx] < bnd - solAbsTol_) || 
            (bnd!=0 && (x[vIdx] < bnd - bnd*solRelTol_))) {
          isCont = false;
          break;        
        } 
        bnd = v->getUb();
        if ((x[vIdx] > bnd + solAbsTol_) || 
            (bnd!=0 && (x[vIdx] > bnd + bnd*solRelTol_))) {
          isCont = false;
          break;        
        }
      }
    }
    delete [] x;
  }
  return cutAdded;
}


void Linearizations::varCons_(std::vector<UInt > varPos,
                              std::vector<UInt > & varConsPos)
{
  VariablePtr v;
  ConstraintPtr con;
  QuadraticFunctionPtr qf;
  NonlinearFunctionPtr nlf;

  for (UInt j = 0; j < nlCons_.size(); ++j) {
    con = nlCons_[j];
    qf = con->getQuadraticFunction();
    nlf = con->getNonlinearFunction();
    for (UInt i = 0; i < varPos.size(); ++i) {
      v = varPtrs_[i];
      if (nlf) {
        if (nlf->hasVar(v)) {
          varConsPos.push_back(j);
          continue;
        }
      }
      if (qf) {
        if (qf->hasVar(v)) {
          varConsPos.push_back(j);
        }
      }
    }
  }
  return;
}

//void Linearizations::rootLinGenScheme2_()
//{
  //VariablePtr v;
  //bool isFound = false;
  //std::vector<UInt > varIdx;
  //UInt fixIdx, vIdx, pos = 0, firstnnz;
  //std::vector<double > alphaSign;
  //int n = minlp_->getNumVars();
  //double varbound, rhs = 0.0, coeff, fixCoeff;   
  //double *xOut = new double[n]; 
  //double *lastDir = new double[n]; 
  
  //std::fill(lastDir, lastDir+n, 0);
  //std::copy(nlpx_, nlpx_+n, xOut);
    
  ////// variable to be fixed in finding search direction
  //for (UInt i = 0; i < varPtrs_.size(); ++i) {
    //v = varPtrs_[i];
    //vIdx = v->getIndex();
    //rhs = rhs + nlpx_[vIdx]*(solC_[vIdx] - nlpx_[vIdx]);
    //if ((isFound == false) && (solC_[vIdx] - nlpx_[vIdx]) != 0) {
      //firstnnz = i;
      //fixIdx = vIdx;
      //isFound = true;              
    //}
  //}

  //if (!isFound) {
    //return;  
  //}
  
  ////// in case rhs is 0, shift the plane to find search directions 
  //if (rhs == 0) {
    //rhs = 1;  
  //}

  //// coefficient of the fix var in any direction
  //fixCoeff = rhs/(solC_[fixIdx] - nlpx_[fixIdx]);
     
  //// coordinate direction for each variable 
  //for (UInt i = 0; i < varPtrs_.size(); ++i) {
    //v = varPtrs_[i];
    //vIdx = v->getIndex();
    //if (i < firstnnz) {
      //// in this case coeff of v is zero in hyperplane expression 
      //varbound = v->getUb();
      //alphaSign.push_back(1);
      //varIdx.push_back(vIdx);
      //search_(varbound, vIdx, nlpx_[vIdx], varIdx, xOut, alphaSign, 0, 0);
      ////// reverse the direction
      //alphaSign[0] = -1;
      //varbound = v->getLb();
      //xOut[vIdx] = nlpx_[vIdx];
      //search_(varbound, vIdx, nlpx_[vIdx], varIdx, xOut, alphaSign, 0, 0);
    //} else if (i > firstnnz) {
      //// in this case coeff of v could be nonzero in the hyperplane expression 
      //coeff = rhs/(solC_[vIdx] - nlpx_[vIdx]);
      //lastDir[vIdx] = -coeff;
      //lastDir[fixIdx] = lastDir[fixIdx] + fixCoeff;
      //boundingVar_(varbound,i,varPtrs_[firstnnz],coeff,-fixCoeff,alphaSign,varIdx);
      //search_(varbound,varIdx[0],nlpx_[varIdx[0]],varIdx,xOut,alphaSign,0,0);
      ////// reverse direction
      //xOut[vIdx] = nlpx_[vIdx];
      //xOut[fixIdx] = nlpx_[fixIdx];
      //varIdx.clear();
      //alphaSign.clear();
      //boundingVar_(varbound, i,varPtrs_[firstnnz],-coeff, fixCoeff, alphaSign, varIdx);
      //search_(varbound,varIdx[0],nlpx_[varIdx[0]],varIdx,xOut,alphaSign,0,0);
      //xOut[fixIdx] = nlpx_[fixIdx];
    //} else {
      //continue;    
    //}
    //varIdx.clear();
    //alphaSign.clear();
    //xOut[vIdx] = nlpx_[vIdx];
  //}

  ////// last direction in positive spanning set 
  //varIdx.resize(0, 0);
  //alphaSign.resize(n, 0);
  //boundingVar_(varbound, pos, lastDir, alphaSign);
  //vIdx = varPtrs_[pos]->getIndex();
  //search_(varbound, vIdx, nlpx_[vIdx] , varIdx, xOut, alphaSign, pos, 1);

  ////// reverse the direction
  //for (UInt i = 0; i < varPtrs_.size(); ++i) {
    //lastDir[i] = -lastDir[i];
  //}
  //std::copy(nlpx_, nlpx_+n, xOut);
  //boundingVar_(varbound, pos, lastDir, alphaSign);
  //vIdx = varPtrs_[pos]->getIndex();
  //search_(varbound, vIdx, nlpx_[vIdx], varIdx, xOut, alphaSign, pos, 1);
  //delete [] lastDir;

  //delete [] xOut;
  //return;
//}


void Linearizations::search_(double varbound, UInt vIdx, double val,  
                             std::vector<UInt > varIdx, double *xOut,
                             std::vector<double > &alphaSign, UInt pos,
                             bool isAllOne)
{
  double alpha;
  //std::cout << "pos alphaSign" << pos << " " << alphaSign.size() << std::endl;
  setStepSize_(varbound, alpha, vIdx, val, alphaSign[pos]);
  //std::cout <<"Stepsize step: varUb, alpha, val " << varbound << " " << alpha << " " << val << "\n";
  
  if (alpha == 0) {
    return;  
  }

  if (isAllOne) {
    UInt idx;
    for (UInt i = 0; i < varPtrs_.size(); ++i) {
      idx = varPtrs_[i]->getIndex();
      alphaSign[i] = alpha*alphaSign[i];  
      xOut[idx] = xOut[idx] + alphaSign[i];
    } 
  } else {
    for (UInt i = 0; i < alphaSign.size(); ++i) {
      alphaSign[i] = alpha*alphaSign[i];  
      xOut[varIdx[i]] = xOut[varIdx[i]] + alphaSign[i];
    }
  }
  foundLinPt_(vIdx, varIdx, pos, alphaSign, varbound, xOut, isAllOne);
  return;
}

// determine which out of varPtrs_ is bounding 
void Linearizations::boundingVar_(double &varbound, 
                                  UInt &pos, 
                                  double *lastDir,
                                  std::vector<double > &alphaSign)
{
  VariablePtr v;
  UInt idx;
  double diff, minDiff = INFINITY, bound;
    
  for (UInt i = 0; i < varPtrs_.size(); ++i) {
    v = varPtrs_[i];
    idx = v->getIndex();
    if (lastDir[i] < 0) {
      bound = v->getLb();
      diff = nlpx_[idx] - bound;
      alphaSign[i] = -1;
    } else if (lastDir[i] > 0) {
      bound = v->getUb();
      diff = bound - nlpx_[idx];
      alphaSign[i] = 1;
    } else {
      alphaSign[i] = 0;
      continue;    
    }
    if (diff < minDiff) {
      pos = i;
      varbound = bound;    
    }
  }
  return;
}

void Linearizations::boundingVar_(double &varbound,
                                  UInt i, VariablePtr fixVar, double coeff,
                                 double fixCoeff, 
                                 std::vector<int > &alphaSign, 
                                 std::vector<UInt > &varIdx)
{
  VariablePtr v = varPtrs_[i];
  double diffCurrent, diffFix, val1, val2;
  UInt vIdx = v->getIndex(), fixIdx = fixVar->getIndex();

  if (coeff < 0) {
    val1 = v->getLb();
    diffCurrent = nlpx_[vIdx] - val1;
    if (fixCoeff < 0) {
      val2 = fixVar->getLb();
      alphaSign.push_back(-1); alphaSign.push_back(-1); // one for each
      diffFix = nlpx_[fixIdx] - val2;
      if (diffCurrent < diffFix) {
        varbound = val1;
        varIdx.push_back(vIdx);
        varIdx.push_back(fixIdx);
      } else {
        varbound = val2;
        varIdx.push_back(fixIdx);
        varIdx.push_back(vIdx);
      }
    } else {
      val2 = fixVar->getUb();
      diffFix = val2 - nlpx_[fixIdx];
      if (diffCurrent < diffFix) {
        varbound = val1;
        varIdx.push_back(vIdx); alphaSign.push_back(-1);
        varIdx.push_back(fixIdx); alphaSign.push_back(1);
      } else {
        varbound = val2;
        varIdx.push_back(fixIdx); alphaSign.push_back(1);
        varIdx.push_back(vIdx); alphaSign.push_back(-1);
      }
    }
  } else {
    val1 = v->getUb();
    diffCurrent = val1 - nlpx_[vIdx];
    if (fixCoeff < 0) {
      val2 = fixVar->getLb();
      diffFix = nlpx_[fixIdx] - val2;
      if (diffCurrent < diffFix) {
        varbound = val1;
        varIdx.push_back(vIdx); alphaSign.push_back(1);
        varIdx.push_back(fixIdx); alphaSign.push_back(-1);
      } else {
        varbound = val2;
        varIdx.push_back(fixIdx); alphaSign.push_back(-1);
        varIdx.push_back(vIdx); alphaSign.push_back(1);
      }
    } else {
      val2 = fixVar->getUb();
      diffFix = val2 - nlpx_[fixIdx];
      alphaSign.push_back(1); alphaSign.push_back(1); // one for each
      if (diffCurrent < diffFix) {
        varbound = val1;
        varIdx.push_back(vIdx);
        varIdx.push_back(fixIdx);
      } else {
        varbound = val2;
        varIdx.push_back(fixIdx);
        varIdx.push_back(vIdx);
      }
    } 
  }
  return;
}
// determine which out of vIdx and fixIdx is bounding 
//void Linearizations::boundingVar_(double &varbound,
                                  //UInt vpos, VariablePtr fixVar, double coeff,
                                 //double fixCoeff, 
                                 //std::vector<double > &alphaSign, 
                                 //std::vector<UInt > &varIdx)
//{
  //VariablePtr v = varPtrs_[vpos];
  //double diffCurrent, diffFix, val1, val2;
  //UInt vIdx = v->getIndex(), fixIdx = fixVar->getIndex();

  //if (coeff < 0) {
    //val1 = v->getLb();
    //diffCurrent = nlpx_[vIdx] - val1;
    //if (fixCoeff < 0) {
      //val2 = fixVar->getLb();
      //alphaSign.push_back(-1); alphaSign.push_back(-1); // one for each
      //diffFix = nlpx_[fixIdx] - val2;
      //if (diffCurrent < diffFix) {
        //varbound = val1;
        //varIdx.push_back(vIdx);
        //varIdx.push_back(fixIdx);
      //} else {
        //varbound = val2;
        //varIdx.push_back(fixIdx);
        //varIdx.push_back(vIdx);
      //}
    //} else {
      //val2 = fixVar->getUb();
      //diffFix = val2 - nlpx_[fixIdx];
      //if (diffCurrent < diffFix) {
        //varbound = val1;
        //varIdx.push_back(vIdx); alphaSign.push_back(-1);
        //varIdx.push_back(fixIdx); alphaSign.push_back(1);
      //} else {
        //varbound = val2;
        //varIdx.push_back(fixIdx); alphaSign.push_back(1);
        //varIdx.push_back(vIdx); alphaSign.push_back(-1);
      //}
    //}
  //} else {
    //val1 = v->getUb();
    //diffCurrent = val1 - nlpx_[vIdx];
    //if (fixCoeff < 0) {
      //val2 = fixVar->getLb();
      //diffFix = nlpx_[fixIdx] - val2;
      //if (diffCurrent < diffFix) {
        //varbound = val1;
        //varIdx.push_back(vIdx); alphaSign.push_back(1);
        //varIdx.push_back(fixIdx); alphaSign.push_back(-1);
      //} else {
        //varbound = val2;
        //varIdx.push_back(fixIdx); alphaSign.push_back(-1);
        //varIdx.push_back(vIdx); alphaSign.push_back(1);
      //}
    //} else {
      //val2 = fixVar->getUb();
      //diffFix = val2 - nlpx_[fixIdx];
      //alphaSign.push_back(1); alphaSign.push_back(1); // one for each
      //if (diffCurrent < diffFix) {
        //varbound = val1;
        //varIdx.push_back(vIdx);
        //varIdx.push_back(fixIdx);
      //} else {
        //varbound = val2;
        //varIdx.push_back(fixIdx);
        //varIdx.push_back(vIdx);
      //}
    //} 
  //}
  //return;
//}


void Linearizations::rootLinGenScheme1_()
{
  //VariablePtr v;
  //UInt lPos, uPos, vIdx;
  //std::vector<UInt > varIdx;  
  //int n = minlp_->getNumVars();
  //double *xOut = new double[n];
  //std::vector<double > alphaSign;
  //double varbound, vLb = INFINITY, vUb = INFINITY;   

  //std::copy(solC_, solC_ + n, xOut);

  //// coordinate direction along each variable in varPtrs_.
  //for (UInt i = 0; i < varPtrs_.size(); ++i) {
    //v = varPtrs_[i];
    //vIdx = v->getIndex();
    //// determine indices of variables that are bounding in the last
    //// direction (and its opposite) of search and bound values 
    //varbound = v->getUb() - solC_[vIdx];
    //if (varbound < vUb) {
      //uPos = i;
      //vUb = varbound;
    //}
    //varbound = solC_[vIdx] - v->getLb();
    //if (varbound < vLb) {
      //lPos = i;
      //vLb = varbound;
    //}

    //// coordinate direction for each variable 
    //varbound = v->getUb();
    //alphaSign.push_back(1);
    //varIdx.push_back(vIdx);
    //search_(varbound, vIdx, solC_[vIdx], varIdx, xOut, alphaSign, 0, 0);
    
    /////reverse search direction if previous direction was unsuccessful 
    //alphaSign[0] = -1;
    //varbound = v->getLb();
    //xOut[vIdx] = solC_[vIdx];
    //search_(varbound, vIdx, solC_[vIdx], varIdx, xOut, alphaSign, 0, 0);
    //varIdx.clear();
    //alphaSign.clear();
    //xOut[vIdx] = solC_[vIdx];
  //}
  
   //// last direction in positive spanning set 
  ////if (vLb == INFINITY) {
    ////lPos = 0; 
    ////vIdx = varPtrs_[0]->getIndex(); 
  ////} else {
    ////vLb = varPtrs_[lPos]->getLb();
    ////vIdx = varPtrs_[lPos]->getIndex(); 
  ////}
  ////for (UInt i = 0; i < varPtrs_.size(); ++i) {
    ////alphaSign.push_back(-1);
  ////}
  ////search_(vLb, vIdx, solC_[vIdx], varIdx, xOut, alphaSign, lPos, 1);
 
  ////if (vUb == INFINITY) {
    ////uPos = 0;  
    ////vIdx = varPtrs_[0]->getIndex(); 
  ////} else {
    ////vUb = varPtrs_[uPos]->getUb();
    ////vIdx = varPtrs_[uPos]->getIndex(); 
  ////}
  ////std::copy(solC_, solC_+minlp_->getNumVars(), xOut);
  ////std::fill(alphaSign.begin(), alphaSign.end(),1);
  ////search_(vUb, vIdx, solC_[vIdx], varIdx, xOut, alphaSign, uPos, 1);
  
  //delete [] xOut;
  //return;
}


void Linearizations::setStepSize_(double &varbound, double &alpha,
                                   UInt vIdx, double val, double boundSign)
{ 
  if (varbound != boundSign*INFINITY) {
    alpha = fabs(varbound-val);
  } else {
    if (fabs(nlpx_[vIdx] - solC_[vIdx]) != 0) {
      alpha = fabs(nlpx_[vIdx] - solC_[vIdx]);
    } else {
      alpha = fabs(val) + 4;
    }
    varbound = val + boundSign*(fabs(10*val) + 10); // parameter here
  }
  alpha = 0.25*alpha;
  return;
}
          

void Linearizations::varsInNonlinCons_()
{
  VariablePtr v;
  FunctionType type;
  for (VariableConstIterator vit = minlp_->varsBegin(); 
       vit != minlp_->varsEnd(); ++vit) {
    v = *vit;
    type = v->getFunType();
    if (!(type == Linear || type == Constant)) {
      //std::cout << v->getName() << std::endl;
      varPtrs_.push_back(v);
    }
  }
  //exit(1);
  return;
}


void Linearizations::foundLinPt_(UInt vIdx, std::vector<UInt> varIdx, 
                                 UInt pos,
                                 std::vector<double> alphaSign, double varbound,
                                 double *xOut, bool isAllOne)
{
  int aSign = 1, error = 0;
  double act, cUb;
  ConstraintPtr con;
  std::vector<ConstraintPtr > vioCons;
  
  /* find constraints violated at xOut. If no constraint is violated then
   * move further along the given direction, If all linear constraints are
   * violated then stop and return. 
   */
  if (alphaSign[pos] < 0) {
    aSign = -1;
  }
  while (true) {
    for (CCIter it = nlCons_.begin(); it != nlCons_.end(); ++it) {
      con = *it;
      cUb = con->getUb();
      act = con->getActivity(xOut, &error);
      if (error == 0) {
        if ((act > cUb + solAbsTol_) &&
            (cUb == 0 || act > cUb + fabs(cUb)*solRelTol_)) { // violated cons
          vioCons.push_back(con);
        }
      }   
    }
    if (vioCons.size() == 0) {
      newPoint_(isAllOne, varIdx, xOut, alphaSign);  
      //if (alphaSign[pos]*(xOut[vIdx]-varbound) > 0) 
      if (aSign*(xOut[vIdx]-varbound) > 0) {
        return;
      }
    } else {
      break;    
    }
  }
  
  /* find point on boundary along the direction and add cuts*/
  UInt n =  minlp_->getNumVars();
  double* xnew = new double[n];
  findBoundaryPt_(xOut, solC_, xnew, vioCons);

  delete [] xnew;
  return;
}

void Linearizations::newPoint_(bool isAllOne,
                               std::vector<UInt> varIdx, double *xOut,
                               std::vector<double> alphaSign)
{
  UInt idx;
  if (isAllOne) {
    for (UInt i = 0; i < varPtrs_.size(); ++i) {
      idx = varPtrs_[i]->getIndex();
      xOut[idx] = xOut[idx] + alphaSign[i];
    }
  } else {
    for (UInt i = 0; i < varIdx.size(); ++i) {
      idx  = varIdx[i];
      xOut[idx] = xOut[idx] + alphaSign[i];
    }
  }
  return;
}

bool Linearizations::boundaryPt_(const double *x,
                                 std::vector<UInt > &vioConsPos,
                                 std::vector<double* > &lastGrad)
{
  ConstraintPtr con;
  int error = 0, repPt = 0;
  UInt numVars =  minlp_->getNumVars();
  bool firstVio, firstActive, cutAdded = false;
  double act, cUb, repSol, repSolOld = 0, lambda1 = 0.5, lambda2 = 0.5;
 
  double* xnew = new double[numVars];
  double* xl = new double[numVars];
  double* xu = new double[numVars];

  std::copy(x, x+numVars, xu);
  std::copy(solC_, solC_+numVars, xl);

  while (true) {
    for (UInt i = 0 ; i < numVars; ++i) {
      xnew[i] = lambda1*xl[i] + lambda2*xu[i];
    }
    firstVio = false, firstActive = false;
    repSol = 0;

    for (UInt k = 0; k < vioConsPos.size(); ) {
      con = nlCons_[vioConsPos[k]];
      cUb = con->getUb();
      act = con->getActivity(xnew, &error);
      repSol = repSol + act;
      if (error != 0) {
        delete [] xnew;
        delete [] xl;
        delete [] xu;
        return false;
      }
      if ((act > cUb + solAbsTol_) &&
          (cUb == 0 || act > cUb + fabs(cUb)*solRelTol_)) { // violated
        if (!firstVio) {
          firstVio = true;
          if (k != 0) {
            vioConsPos.erase(vioConsPos.begin(), vioConsPos.begin() + k);
            k = 0;
          }
        }
        ++k;
      } else if ((fabs(act-cUb) <= solAbsTol_) ||
            (cUb != 0 && fabs(act- cUb) <= fabs(cUb)*solRelTol_)) { // active
        if (firstVio) {
          vioConsPos.erase(vioConsPos.begin() + k);
        } else {
          if (!firstActive) {
            firstActive = true;         
            if (k != 0) {
              vioConsPos.erase(vioConsPos.begin(),vioConsPos.begin() + k);
              k = 0;
            }
          }
          //activeConsAct.push_back(act);
          ++k;
        }
      } else {
        if (firstVio || firstActive) {
           vioConsPos.erase(vioConsPos.begin() + k);
        } else {
           ++k;
        }   
      }
    }
    //MS: Implement this login in Scheme 3 as well
    if (repSol == repSolOld) {
      ++repPt;    
    } else {
      repPt = 0;
      repSolOld = repSol;    
    }

    if (repPt == 5) {
      firstVio = false;
      firstActive = true;    
    }

    if (!firstVio) {
      if (!firstActive) {
        std::copy(xnew,xnew+numVars,xl);
      } else {
        // add linearization to active nonlinear constraints only if the
        // linearizations are far apart
        cutAdded = genLin_(xnew, vioConsPos, lastGrad);
        break;
      }
    } else {
      std::copy(xnew,xnew+numVars,xu);
    } 
  }
  
  delete [] xnew;
  delete [] xl;
  delete [] xu;
  return cutAdded;
}

bool Linearizations::genLin_(double *x, std::vector<UInt > vioConsPos,
                                     std::vector<double *> &lastGrad)
{
  UInt cIdx;
  FunctionPtr f;
  ConstraintPtr con;
  std::stringstream sstm;
  LinearFunctionPtr lf = 0;
  bool isCont, isFound = false;
  int error, n = rel_->getNumVars();
  double angle, PI = 3.14159265, d, c, m1, m2, cUb, act;
  VariableConstIterator vbeg = rel_->varsBegin(), vend = rel_->varsEnd();
  const double linCoeffTol =
    env_->getOptions()->findDouble("conCoeff_tol")->getValue();

  for (UInt j = 0; j < vioConsPos.size(); ++j) {
    cIdx = vioConsPos[j];
    angle = 0;
    error = 0;
    isCont = false;
    double *a = new double[n];
    std::fill(a, a+n, 0.);
    
    con = nlCons_[cIdx];
    f = con->getFunction();
    f->evalGradient(x, a, &error);
    
    if (error == 0) {
      if (lastGrad[cIdx] == NULL) {
        isCont = true;
      } else {
        // compute angle
        d = InnerProduct(a, lastGrad[cIdx], n);
        m1 = sqrt(InnerProduct(a, a, n));
        m2 = sqrt(InnerProduct(lastGrad[cIdx], lastGrad[cIdx], n));
        angle  = acos(d/(m1*m2))*180/PI;
      }
    } else {
      delete [] a;
      a = 0;
      continue;          
    }
    if (fabs(angle) >= rsg2Per_ || isCont) {
      cUb = con->getUb();
      act = con->getActivity(x, &error);
      if (error == 0) {
        lf = (LinearFunctionPtr) new LinearFunction(a, vbeg, vend, linCoeffTol);
        c  = act - InnerProduct(x, a, minlp_->getNumVars());
        ++(stats_->rgs2Cuts);
        sstm << "_OACutRoot_" << stats_->rgs2Cuts;
        f = (FunctionPtr) new Function(lf);
        rel_->newConstraint(f, -INFINITY, cUb-c, sstm.str());
        //newcon = rel_->newConstraint(f, -INFINITY, cUb-c, sstm.str());
        sstm.str("");
        if (lastGrad[cIdx]) {
          delete [] lastGrad[cIdx];
          lastGrad[cIdx] = 0;
        }
        lastGrad[cIdx] = a;
        isFound = true;
      } else {
        delete [] a;      
        a = 0;
      }
    } else {
      delete [] a;
      a = 0;
    }
  }
  return isFound;
}


bool Linearizations::findBoundaryPt_(const double *xOut, const double *xIn,
                                     double *x, 
                                     std::vector<ConstraintPtr> &vioCons)
{
  int error = 0;
  double act, cUb; 
  ConstraintPtr con;
  bool firstVio, firstActive;
  std::vector<double > activeConsAct;
  UInt numVars =  minlp_->getNumVars();
 
  double* xl = new double[numVars];
  double* xu = new double[numVars];

  std::copy(xOut, xOut+numVars, xu);
  std::copy(xIn, xIn+numVars, xl);
 
  while (true) { 
    for (UInt i = 0 ; i < numVars; ++i) {
      x[i] = 0.5*(xl[i] + xu[i]);
    }
    firstVio = false, firstActive = false;

    for (UInt k = 0; k < vioCons.size(); ) {
      con = vioCons[k];
      cUb = con->getUb();
      act = con->getActivity(x, &error);
      if (error != 0) {
        delete [] xl;
        delete [] xu;
        return false;
      }
      if ((act > cUb + solAbsTol_) &&
          (cUb == 0 || act > cUb + fabs(cUb)*solRelTol_)) { // violated
        if (!firstVio) {
          firstVio = true;
          if (k != 0) {
            vioCons.erase(vioCons.begin(), vioCons.begin() + k);
            k = 0;
          }
        }
        //MS: seems some error
        ++k;
      } else if ((fabs(act-cUb) <= solAbsTol_) ||
            (cUb != 0 && fabs(act- cUb) <= fabs(cUb)*solRelTol_)) { // active
        if (firstVio) {
          vioCons.erase(vioCons.begin() + k);
        } else {
          if (!firstActive) {
            firstActive = true;         
            if (k != 0) {
              vioCons.erase(vioCons.begin(),vioCons.begin() + k);
              k = 0;
            }
          }
          activeConsAct.push_back(act);
          ++k;
        }
      } else {
        if (firstVio || firstActive) {
           vioCons.erase(vioCons.begin() + k);
        } else {
           ++k;
        }   
      }
    }

    if (!firstVio) {
      if (!firstActive) {
        std::copy(x,x+numVars,xl);
      } else {
        // add linearization to active nonlinear constraints 
        double c;
        FunctionPtr f;
        std::stringstream sstm;
        LinearFunctionPtr lf = 0;
        for (UInt j = 0; j < vioCons.size(); ++j) {
          con = vioCons[j];
          f = con->getFunction();
          linearAt_(f, activeConsAct[j], x, &c, &lf, &error);
          if (error == 0) {
            cUb = con->getUb();
            if (rgs1_) {
              ++(stats_->rgs1Cuts);
              sstm << "_OACutRoot_" << stats_->rgs1Cuts;
            } else if (rsg2Per_) {
              ++(stats_->rgs2Cuts);
              sstm << "_OACutRoot_" << stats_->rgs2Cuts;
            } else if (rs3_) {
              ++(stats_->rs3Cuts);
              sstm << "_OACutRoot_" << stats_->rs3Cuts; 
            } else {
              sstm << "_OACutRoot_"; 
              // Later: print message here.            
            }
            f = (FunctionPtr) new Function(lf);
            rel_->newConstraint(f, -INFINITY, cUb-c, sstm.str());
            //newcon = rel_->newConstraint(f, -INFINITY, cUb-c, sstm.str());
            sstm.str("");
          }
        }
        delete [] xl;
        delete [] xu;
        return true;
      }
    } else {
      std::copy(x,x+numVars,xu);
    } 
  }
  
  delete [] xl;
  delete [] xu;
  return false;    
}
 

void Linearizations::rootLinScheme1_(ConstraintPtr con, double lVarCoeff,
                            UInt lVarIdx, UInt nVarIdx, double nVarCoeff)
{
  double iP[2]; // intersection point
  UInt newConId;
  bool shouldCont;
  ConstraintPtr newcon;
  std::vector<UInt > newConsId;
  VariablePtr vnl = NULL, vl = NULL;
  std::vector<double > linVioVal, xc, yc; // xc and yc  nonlinear and lin var
  int i, error = 0, n = rel_->getNumVars();
  double act, cUb, y1, y2, vLb, vUb, maxVio, stopCond, consUb; 
  double *b1 = new double[n];

  std::fill(b1, b1+n, 0.);
  vl = rel_->getVariable(lVarIdx);
  vnl = rel_->getVariable(nVarIdx);
  
  vLb = vnl->getLb();
  vUb = vnl->getUb();

  if (vLb == -INFINITY) {
    if (vUb == INFINITY) {
      vLb = -50;
      vUb = 50;
    } else {
      vLb = vUb - 100;
    }
  } else {
    if (vUb == INFINITY) {
      vUb = vLb + 100;
    } 
  }
    
  b1[nVarIdx] = vLb;
  act = nVarCoeff*vLb;
  shouldCont = linPart_(b1, lVarIdx, con, lVarCoeff, act);  
  if (shouldCont) {
    shouldCont = addNewCut_(b1, con, newConId);
    if (shouldCont) {
      y1 = b1[lVarIdx];
      newConsId.push_back(newConId); 
    } else {
      delete [] b1;
      return;    
    }
  } else {
    delete [] b1;
    return;    
  }

  // upper bound of var in nonlinear cons
  b1[nVarIdx] = vUb;
  act = nVarCoeff*vUb;
  shouldCont = linPart_(b1, lVarIdx, con, lVarCoeff, act);  
  if (shouldCont) {
    shouldCont = addNewCut_(b1, con, newConId);
    if (shouldCont) {
      y2 = b1[lVarIdx];
      newConsId.push_back(newConId); 
    } else {
      delete [] b1;
      return;    
    }
  } else {
    delete [] b1;
    return;    
  }

  shouldCont = findIntersectPt_(newConsId, vl, vnl, iP);
  if (shouldCont == false) {
    delete [] b1;
    return;    
  }

  // populate points and their violation in cyclic order
  xc.push_back(vLb);
  yc.push_back(y1);
  linVioVal.push_back(0);
      
  b1[nVarIdx] = iP[0], b1[lVarIdx] = iP[1];
  act = con->getActivity(b1, &error);
  if (error != 0) {
    delete [] b1;
    return;    
  }  
  consUb = con->getUb();
  act = std::max(act-consUb, 0.0); 
  xc.push_back(iP[0]), yc.push_back(iP[1]), linVioVal.push_back(act);
  
  xc.push_back(vUb), yc.push_back(y2), linVioVal.push_back(0);

  // starting from intersection point
  i = 1;
  maxVio = linVioVal[i];
  if (fabs(consUb) > solAbsTol_) { 
    stopCond = consUb*rs1_/100;    
  } else {
    stopCond = maxVio*rs1_/100;    
  }

  if ((stopCond < solAbsTol_) || 
      (consUb!=0 && stopCond < fabs(consUb)*solRelTol_ )) { 
    delete [] b1;
    return;
  }

  while (maxVio >= stopCond) { 
    //add a new cut at the point indexed i
    b1[nVarIdx] = xc[i];
    shouldCont = addNewCut_(b1, con, newConId);
    if (shouldCont) {
      newcon = rel_->getConstraint(newConId);
      cUb = newcon->getUb();
    } else {
      break;
    }
    
    // Move right and determine first point that satisfy the newcon
    for (UInt j = i+1; j < xc.size(); ) {
      b1[nVarIdx] = xc[j], b1[lVarIdx] = yc[j];
      act = newcon->getActivity(b1, &error);
      if (error == 0) {
        if ((act < cUb + solAbsTol_) ||
            (cUb =! 0 && act < cUb + fabs(cUb)*solRelTol_)) {
            //insert new point just before index j
          insertNewPt_(j, j-1, xc, yc, newcon, vl, vnl, shouldCont); 
          b1[nVarIdx] = xc[j], b1[lVarIdx] = yc[j];
          act = con->getActivity(b1, &error);
          if (error != 0) {
            shouldCont = false;    
          } else {
            act = std::max(act-consUb, 0.0); 
            linVioVal.insert(linVioVal.begin()+j,act);
          }
          break;
        } else {
          // delete point if violates newcon
          xc.erase(xc.begin() + j);
          yc.erase(yc.begin() + j);        
          linVioVal.erase(linVioVal.begin() + j);        
        }  
      }
    }
    if (shouldCont == false) {
      break;
    }
    int j = i-1;
    while (j >= 0) {
      b1[nVarIdx] = xc[j], b1[lVarIdx] = yc[j];
      act = newcon->getActivity(b1, &error);
      if (error == 0) {
        if ((act < cUb + solAbsTol_) ||
            (cUb =! 0 && act < cUb + fabs(cUb)*solRelTol_)) {
          insertNewPt_(j+1, j, xc, yc, newcon, vl, vnl, shouldCont); 
          b1[nVarIdx] = xc[j+1], b1[lVarIdx] = yc[j+1];
          act = con->getActivity(b1, &error);
          if (error != 0) {
            shouldCont = false;    
          } else {
            act = std::max(act-consUb, 0.0); 
            linVioVal.insert(linVioVal.begin()+j+1,act);
            xc.erase(xc.begin()+j+2);
            yc.erase(yc.begin()+j+2);        
            linVioVal.erase(linVioVal.begin()+j+2);        
          }
          break;
        }  else {
          // delete point if violates newcon
          xc.erase(xc.begin() + j);
          yc.erase(yc.begin() + j);        
          linVioVal.erase(linVioVal.begin() + j);        
          --j;
        }  
      }
    }
    if (shouldCont == false) {
      break;
    }
    maxVio = *(std::max_element(linVioVal.begin(), linVioVal.end()));
    if ((maxVio < solAbsTol_) || 
        (consUb!=0 && maxVio < fabs(consUb)*solRelTol_ )) { 
      break;
    }
    i = std::max_element(linVioVal.begin(), linVioVal.end())-linVioVal.begin();     
  }
  delete [] b1;
  return;
}


void Linearizations::rootLinScheme2_(ConstraintPtr con,
                                     double lVarCoeff,
                                     UInt lVarIdx, UInt nVarIdx)
{
  int error = 0;
  FunctionPtr f;
  VariablePtr vnl;
  UInt n = minlp_->getNumVars();
  double lastSlope, delta, nlpSlope, nbhSize;
  
  vnl = rel_->getVariable(nVarIdx);
  
  double* npt = new double[n];
  std::fill(npt, npt+n, 0.);
  
  double *grad = new double[n];
  std::fill(grad, grad+n, 0.);
  
  f = con->getFunction();
  f->evalGradient(nlpx_, grad, &error);

  if (error != 0) {
    return;
  }
  
  nlpSlope = -1*(grad[nVarIdx]/lVarCoeff);
  lastSlope = nlpSlope;       // nlpSlope is going to be used later on as well
  
  nbhSize = std::max(vnl->getLb(), nlpx_[nVarIdx] - rs2NbhSize_);    
  //nbhSize = vnl->getLb();    
  if (nlpx_[nVarIdx] - nbhSize >= 1) {
    delta = 0.5;  
  } else {
    delta = nlpx_[nVarIdx] - nbhSize;  
  }

  npt[nVarIdx] = nlpx_[nVarIdx] - delta;
     
  if (delta != 0) {
    while (npt[nVarIdx] >= nbhSize) {
      grad[nVarIdx] = 0; grad[lVarIdx] = 0;
      rScheme2Cut_(con, delta, lVarCoeff, lastSlope, nVarIdx, npt, grad);
      npt[nVarIdx] =  npt[nVarIdx] - delta;
    }
  }
  
  nbhSize = std::min(vnl->getUb(), nlpx_[nVarIdx] + rs2NbhSize_);
  //nbhSize = vnl->getUb();    
  if (nbhSize - nlpx_[nVarIdx] >= 1) {
    delta = 0.5;  //MS: can be a parameter  
  } else {
    delta = nbhSize - nlpx_[nVarIdx];  
  }

  lastSlope = nlpSlope;
  npt[nVarIdx] = nlpx_[nVarIdx] + delta;

  if (delta != 0) {
    while (npt[nVarIdx] <= nbhSize) {
      grad[nVarIdx] = 0; grad[lVarIdx] = 0;
      rScheme2Cut_(con, delta, lVarCoeff, lastSlope, nVarIdx, npt, grad);
      npt[nVarIdx] =  npt[nVarIdx] + delta;
    }
  }
  delete [] grad;
  delete [] npt;
  return;
}


void Linearizations::rScheme2Cut_(ConstraintPtr con, double &delta,
                                double lVarCoeff, double &lastSlope,
                                UInt nVarIdx, double * npt, double * grad)
{
  int error = 0;
  FunctionPtr f = con->getFunction();
  double newSlope, angle, tanTheta, PI = 3.14159265;
  
  f->evalGradient(npt, grad, &error);
  if (error != 0) {
    return;
  } 
  
  newSlope = -1*(grad[nVarIdx]/lVarCoeff);
  tanTheta = (newSlope-lastSlope)/(1+newSlope*lastSlope);
  angle = atan (tanTheta) * 180 / PI;

 // MS: old stuff beased on gradient comparison 
  //if ((lastSlope == 0 && newSlope == 0) ||
      //(lastSlope != 0 && fabs((newSlope-lastSlope)/lastSlope)*100 <  rs2Per_)) {
    //delta = 2*delta;
    //return;
  //}

  // Add new linearization if angle between the two lines is at least rs2Per_
  if (fabs(angle) <  rs2Per_) {
    delta = 2*delta;
    return;
  }

  lastSlope = newSlope;
  //ConstraintPtr newcon;
  std::stringstream sstm;
  LinearFunctionPtr lf = 0;
  VariableConstIterator vbeg = rel_->varsBegin(), vend = rel_->varsEnd();
  const double linCoeffTol =
    env_->getOptions()->findDouble("conCoeff_tol")->getValue();
  double c, cUb = con->getUb(), act = con->getActivity(npt, &error);
  
  lf = (LinearFunctionPtr) new LinearFunction(grad, vbeg, vend, linCoeffTol);
  c  = act - InnerProduct(npt, grad, minlp_->getNumVars());

  ++stats_->rs2Cuts;
  sstm << "_OAcut_" << stats_->rs2Cuts << "_AtRoot";
  f = (FunctionPtr) new Function(lf);
  rel_->newConstraint(f, -INFINITY, cUb-c, sstm.str());
  //newcon = rel_->newConstraint(f, -INFINITY, cUb-c, sstm.str());
  //newcon->write(std::cout);
  return;
}


bool Linearizations::shouldStop_(EngineStatus eStatus)
{
  bool shouldStop = false;
  switch (eStatus) {
   case (FailedInfeas):
     logger_->msgStream(LogInfo) << me_ << "failed to converge "
     << "(infeasible) in root" << std::endl;
     shouldStop = true;
     break;
   case (ProvenFailedCQInfeas):
     logger_->msgStream(LogInfo) << me_ << "constraint qualification "
                                        << "violated in root "
                                        << std::endl;
   case (ProvenInfeasible):
   case (ProvenLocalInfeasible):
     shouldStop = true;
     break;
   case (ProvenObjectiveCutOff):
     shouldStop = true;
     break;
   case (ProvenUnbounded):
     shouldStop = false;
     logger_->msgStream(LogDebug2) << me_ << "problem relaxation is "
                                   << "unbounded!" << std::endl;
     assert(!"Relaxation unbounded."); 
     break;
   case (FailedFeas):
     logger_->msgStream(LogInfo) << me_ << "Failed to converge " 
                                 << "(feasible) in root " << std::endl;
     break;
   case (ProvenFailedCQFeas):
     logger_->msgStream(LogInfo) << me_ << "constraint qualification "
                                 << "violated in root" << std::endl;
     break;
   case (EngineIterationLimit):
     logger_->msgStream(LogInfo) << me_ << "engine hit iteration limit, "
                                 << "continuing in root" << std::endl;
     // continue with this node by following ProvenLocalOptimal case.
   case (ProvenLocalOptimal):
   case (ProvenOptimal):
     break;
   case (EngineError):
     shouldStop = true;
     break;
   default:
     break;
  }
  return shouldStop;
}

//MS: add esh all from LP solution - working
void Linearizations::rootLinScheme3(EnginePtr lpe, VariablePtr objVar,
                                    SeparationStatus *status)
{
 //// ESH to all nonlinear constraints 
  if (solC_ == 0) {
    return;  
  }
  //int error = 0;
  //UInt numOldCuts;
  //double act, cUb;
  //ConstraintPtr con;
  //const double *lpx;
  //double *x = new double[minlp_->getNumVars()]; 
  
  //for (UInt i = 1; i <= rs3_; ++i) {
    //numOldCuts = stats_->rs3Cuts;
    //lpx = lpe->getSolution()->getPrimal();
    ////lpe->getSolution()->writePrimal(std::cout);
    ////std::cout << lpe->getSolution()->getObjValue() << std::endl;
    //for (CCIter it = nlCons_.begin(); it!=nlCons_.end(); ++it) {
      //con = *it;
      //cUb = con->getUb();
      //act = con->getActivity(lpx, &error);
      //if (error == 0) {
        //if ((act > cUb + solAbsTol_) &&
            //(cUb == 0 || act > cUb + fabs(cUb)*solRelTol_)) {
          //cutAtLineSearchPt_(solC_, lpx, x, con);
        //}
      //} else {
        //logger_->msgStream(LogError) << me_ << "Constraint" <<  con->getName()
          //<< " is not defined at this point." << std::endl;
      //}
    //}
    //if (numOldCuts < stats_->rs3Cuts) {
      //lpe->solve();
      //if (shouldStop_(lpe->getStatus())) {
        //break;    
      //}
    //} else {
      //break;
    //}
  //}

  //if (stats_->rs3Cuts > 0) {
    //*status = SepaResolve;
  //}
  //delete [] x;
  //return;


  //// ESH only at the boundary point
  int error = 0;
  FunctionPtr f;
  UInt numOldCuts;
  ConstraintPtr con;
  double nlpact, cUb;
  const double *lpx;
  std::vector<ConstraintPtr > vioCons;
  double *xnew = new double[minlp_->getNumVars()];
  
  for (UInt i = 1; i <= rs3_; ++i) {
    lpx = lpe->getSolution()->getPrimal();
    for (CCIter it=nlCons_.begin(); it!=nlCons_.end(); ++it) {
      con = *it;
      f = con->getFunction();
      nlpact = f->eval(lpx, &error);
      if (error == 0) {
        cUb = con->getUb();
        if ((nlpact > cUb + solAbsTol_) &&
            (cUb == 0 || nlpact > cUb + fabs(cUb)*solRelTol_)) {
          vioCons.push_back(con);
        }
      }
    }
    if (vioCons.size() == 0) {
      break;    
    }
    numOldCuts = stats_->rs3Cuts;
    findBoundaryPt_(lpx, solC_, xnew, vioCons);
    if (numOldCuts < stats_->rs3Cuts) {
      lpe->solve();
      if (shouldStop_(lpe->getStatus())) {
        break;    
      }
    } else {
      break;
    }
    vioCons.clear();
  }

  if (stats_->rs3Cuts > 0) {
    *status = SepaResolve;
  }
  
  delete [] xnew;
  return;
}


bool Linearizations::uniVarNlFunc_(ConstraintPtr con, double &lVarCoeff,
                            UInt & lVarIdx, UInt & nVarIdx, double &nVarCoeff)
{
  double coeff;
  bool foundVar = false, foundNVar = false;
  UInt nlTerms = 0, qTerms = 0, idx;
  LinearFunctionPtr lf = con->getLinearFunction();
  QuadraticFunctionPtr qf = con->getQuadraticFunction();
  NonlinearFunctionPtr nlf = con->getNonlinearFunction();
  const double linCoeffTol =
    env_->getOptions()->findDouble("conCoeff_tol")->getValue();
  if (nlf) {
    nlTerms = nlf->numVars();
    if (nlTerms != 1) {
      return false;    
    }
    nVarIdx= (*(nlf->varsBegin()))->getIndex(); //index of var in nonlinear term
  }

  if (qf) {
    qTerms = qf->getNumVars();
    if (qTerms != 0) {
      if (qTerms > 1) {
        return false;    
      }
      if (nlTerms > 1) {
        if (nVarIdx != ((qf->varsBegin())->first)->getIndex()) {
          return false;      
        }
      } else {
        nVarIdx = ((qf->varsBegin())->first)->getIndex(); 
      }
    }
  }
  
  if (lf) {
    for(VariableGroupConstIterator vit = lf->termsBegin();
        vit != lf->termsEnd(); ++vit) {
      coeff = vit->second;
      idx = (vit->first)->getIndex();
      if (idx == nVarIdx) {
        foundNVar = true;
        nVarCoeff = coeff;
        continue;      
      }
      if (fabs(coeff) > linCoeffTol && foundVar == false) {
        lVarIdx = idx;
        foundVar = true;
        lVarCoeff = coeff;
      }
      if (foundVar && foundNVar) {
        break;      
      }
    }
  }
  
  if (foundVar) {
    return true;
  }
  return false;
}


double Linearizations::maxVio(const double *x, int &index)
{
  ConstraintPtr c;
  int error=0, i = 0;
  double act, cUb, vio = 0.0, max = -INFINITY;

  for (CCIter it=nlCons_.begin(); it!=nlCons_.end(); ++it, ++i) {
    c = *it;
    act = c->getActivity(x, &error);
    if (error == 0) {
      cUb = c->getUb();
      if (act > cUb+solAbsTol_ && (cUb == 0 ||
                                   act > cUb + fabs(cUb)*solRelTol_)) {
        //if (fabs(cUb) > solAbsTol_ && fabs(cUb) > solRelTol_) {
        if (fabs(cUb) > solAbsTol_) {
          vio = 100*(act - cUb)/fabs(cUb);      
        } else {
          vio = act - cUb;
        }
        if (vio > max) {
          max = vio;          
          index = i; 
        }
      }      
    }
  }
  return max;
}


void Linearizations::writeStats(std::ostream &out) const
{
  out
    << me_ << "number of cuts in root scheme 1      = "
    << stats_->rs1Cuts << std::endl
    << me_ << "number of cuts in root scheme 2      = "
    << stats_->rs2Cuts << std::endl
    << me_ << "number of cuts in root scheme 3      = "
    << stats_->rs3Cuts << std::endl
    << me_ << "number of cuts in root gen. scheme 1 = "
    << stats_->rgs1Cuts << std::endl
    << me_ << "number of cuts in root gen. scheme 2 = "
    << stats_->rgs2Cuts << std::endl;

  return;
}


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
