#include "Node.h"
#include "Branch.h"
#include "VarBoundMod.h"
#include "Problem.h"
#include "Serializer.h"
#include <iostream>

using namespace Minotaur;

/**
 * Serializing Methods
 */

// change so that only Arithmetic types are allowed
template <typename T>
void Serializer::writeArith(const T& var)
{
  _stream.write((char *) &var, sizeof(T));
}

std::string Serializer::get_string()
{
  return _stream.str();
}

void Serializer::writeNode(NodePtr node)
{
  writeBranch(node->getBranch());

  writeArith(node->getId());
  writeArith(node->getLb());

  
  writeMods(node->modsrBegin(), node->modsrEnd());
}

void Serializer::writeBranch(BranchPtr branch)
{
  writeMods(branch->rModsBegin(), branch->rModsEnd());
  writeArith(branch->getActivity());
}

void Serializer::writeVbm(VarBoundModPtr vbm)
{
  writeArith((vbm->getVar())->getId());
  writeArith((short) vbm->getLU());
  writeArith(vbm->getNewVal());
  writeArith(vbm->getOldVal());
}

// assumes all the modifications are VarBoundMod. Modify to incorporate other modifications too.
void Serializer::writeMods(ModificationConstIterator begin, ModificationConstIterator end)
{
  writeArith((size_t) (end - begin));
  for(auto ptr = begin; ptr != end; ptr++)
  {
    writeVbm((VarBoundModPtr) *ptr);
  }
}


/**
 * DeSerializing Methods
 */

// change template type too
// template <typename T>
// T read() {}
template <typename T>
T DeSerializer::readArith()
{
  T var;
  _stream.read((char *) &var, sizeof(T));
  return var;
}

NodePtr DeSerializer::readNode(ProblemPtr prob, NodePtr root)
{
  BranchPtr cur_branch = readBranch(prob);

  NodePtr node = new Node(root, cur_branch);


  node->setId(readArith<UInt>());
  node->setLb(readArith<double>());


  std::vector<ModificationPtr> rmods = readMods(prob);


  for(auto i = 0; i < rmods.size(); i++)
    node->addRMod(rmods[i]);

  return node;
}

BranchPtr DeSerializer::readBranch(ProblemPtr prob)
{
  BranchPtr br = new Branch();

  std::vector<ModificationPtr> rmods = readMods(prob);
  for(auto i = 0; i < rmods.size(); i++)
    br->addRMod(rmods[i]);

  br->setActivity(readArith<double>());

  return br;
}

std::vector<ModificationPtr> DeSerializer::readMods(ProblemPtr prob)
{
  std::vector<ModificationPtr> vec;


  vec.resize(readArith<size_t>());

  for(auto i = 0; i < vec.size(); i++)
  {
    VarBoundModPtr vbm = readVarBoundMod(prob);
    vec[i] = (ModificationPtr) vbm;
  }

  return vec;
}

VarBoundModPtr DeSerializer::readVarBoundMod(ProblemPtr prob)
{
  
  UInt varid = readArith<UInt>();
  short bndtype = readArith<short>();

  VariablePtr var = prob->getVariable(varid);
  BoundType bt = static_cast<BoundType>(bndtype);

  VarBoundModPtr vbm = new VarBoundMod(var, bt, readArith<double>());

  vbm->setOldVal(readArith<double>());

  return vbm;
}
