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
  // need to remove/ update this
  writeArith(node->getId());
  writeArith(node->getLb());

  std::vector<NodePtr> ancestors;
  NodePtr curNode = node;
  while(curNode)
  {
    ancestors.push_back(curNode);
    curNode = curNode->getParent();
  }

  std::vector<ModificationPtr> allmods;
  
  for(int i = ancestors.size() - 1; i >= 0; i--)
  {
    BranchPtr curbranch = ancestors[i]->getBranch();
    if(curbranch)
    {
      for(auto itr = curbranch->rModsBegin(); itr != curbranch->rModsEnd(); itr++)
        allmods.push_back(*itr);
    }
     
    for(auto itr = ancestors[i]->modsrBegin(); itr != ancestors[i]->modsrEnd(); itr++)
      allmods.push_back(*itr);
  }

  writeMods(allmods.begin(), allmods.end());
}

void Serializer::writeBranch(BranchPtr branch)
{
  writeMods(branch->rModsBegin(), branch->rModsEnd());
  writeArith(branch->getActivity());
}

static std::pair<int, short> getKey(VarBoundModPtr vbm)
{
  return { 
    (vbm->getVar())->getIndex(),
      (short) vbm->getLU()
  };
}

static std::pair<double, double> getVals(VarBoundModPtr vbm)
{
  return {
    vbm->getOldVal(),
      vbm->getNewVal()
  };
}

void Serializer::writeVbm(VarBoundModPtr vbm)
{
  writeArith((vbm->getVar())->getIndex());
  writeArith((short) vbm->getLU());
  writeArith(vbm->getNewVal());
  writeArith(vbm->getOldVal());
}

// assumes all the modifications are VarBoundMod. Modify to incorporate other modifications too.
void Serializer::writeMods(ModificationConstIterator begin, ModificationConstIterator end)
{
  std::map<std::pair<int, short>, std::pair<double, double>> vbm_map; 
  for(auto ptr = begin; ptr != end; ptr++)
  {
    auto key = getKey((VarBoundModPtr) *ptr);
    auto vals = getVals((VarBoundModPtr) *ptr);
    if (vbm_map.find(key) == vbm_map.end())
      vbm_map[key] = vals;
    else
    {
      vals.first = vbm_map[key].first;
      vbm_map[key] = vals;
    }
  }
  writeArith(vbm_map.size());
  for (auto &p: vbm_map)
  {
    writeArith(p.first.first);
    writeArith(p.first.second);
    writeArith(p.second.first);
    writeArith(p.second.second);
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

NodePtr DeSerializer::readNode(ProblemPtr prob)
{
  NodePtr node;
  
  node = new Node();

  node->setId(readArith<UInt>());
  node->setLb(readArith<double>());

  std::vector<ModificationPtr> rmods = readMods(prob);

  for(size_t i = 0; i < rmods.size(); i++)
    node->addRMod(rmods[i]);

  return node;
}

BranchPtr DeSerializer::readBranch(ProblemPtr prob)
{
  BranchPtr br = new Branch();

  std::vector<ModificationPtr> rmods = readMods(prob);
  for(size_t i = 0; i < rmods.size(); i++)
    br->addRMod(rmods[i]);

  br->setActivity(readArith<double>());

  return br;
}

std::vector<ModificationPtr> DeSerializer::readMods(ProblemPtr prob)
{
  std::vector<ModificationPtr> vec;


  vec.resize(readArith<size_t>());

  for(size_t i = 0; i < vec.size(); i++)
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

  double oldval = readArith<double>(), newval = readArith<double>();

  VarBoundModPtr vbm = new VarBoundMod(var, bt, newval);
  vbm->setOldVal(oldval);

  return vbm;
}
