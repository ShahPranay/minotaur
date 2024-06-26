// 
//     Minotaur -- It's only 1/2 bull
// 
//     (C)opyright 2009 - 2024 The Minotaur Team.
// 

#include <cmath>
#include <string>
#include <iostream>

#include <Types.h>
#include "MinotaurConfig.h"
#include "Environment.h"
#include "LinearFunction.h"
#include "Branch.h"
#include "Node.h"
#include "Serializer.h"
#include "Variable.h"
#include "ProblemSize.h"
#include "SerializeUT.h"

CPPUNIT_TEST_SUITE_REGISTRATION(SerializeUT);
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(SerializeUT, "SerializeUT");

using namespace Minotaur;

void SerializeUT::setUp()
{
  env_ = new Environment();
  instance_ = new Problem(env_);

  std::vector<VariablePtr> vars;

  vars.push_back(instance_->newVariable(0.0, INFINITY, Integer));
  vars.push_back(instance_->newVariable(0.0, 3.0, Integer));

  BranchPtr br = new Branch();
  ModificationPtr mod1 = new VarBoundMod(vars[1], BoundType::Upper, 2);
  br->addRMod(mod1);

  par_ = new Node();
  mynode_ = new Node(par_, br); 
  ModificationPtr mod2 = new VarBoundMod(vars[0], BoundType::Lower, -1);

  mynode_->addRMod(mod1);
  mynode_->setId(69);
  
}

void SerializeUT::test_serialize_reconstruct()
{
  Serializer cereal;

  cereal.writeNode(mynode_);
  std::string cpy = cereal.get_string();

  DeSerializer des(cpy);

  NodePtr newnode = des.readNode(instance_);

  CPPUNIT_ASSERT((*newnode) == (*mynode_));

}

void SerializeUT::tearDown()
{
  delete instance_;
  delete env_;
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
