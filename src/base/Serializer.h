
/**
 * \file Serializer.h
 * \brief Declare classes for serializing and deserialzing objects to strings / char arrays and vice versa.
 * \author Pranay Shah, Indian Institute of Technology, Delhi
 */

#ifndef MINOTAURSERIALIZER_H
#define MINOTAURSERIALIZER_H

#include <vector>
#include <sstream>
#include <string>

#include "Types.h" 

namespace Minotaur {
  // c++ 20 feature
  /* template <typename T> */
  /*   concept Arithmetic = std::is_arithmetic_v<T>; */

  class Serializer {
    public:
      Serializer() {  };

      ~Serializer() {  };

      // change so that only Arithmetic types are allowed
      template<typename T>
        void writeArith(const T& var);

      void writeNode(NodePtr node);
      void writeBranch(BranchPtr branch);
      void writeVbm(VarBoundModPtr vbm);
      void writeMods(ModificationConstIterator begin, ModificationConstIterator end);

      std::string get_string();

    private:
      std::ostringstream _stream;
  };

  class DeSerializer {
    public:
      DeSerializer(const std::string& str): _stream(str) {  }; 

      ~DeSerializer() {  };

      template<typename T>
        T readArith();

      NodePtr readNode (ProblemPtr prob, NodePtr root);
      BranchPtr readBranch (ProblemPtr prob);
      std::vector<ModificationPtr> readMods (ProblemPtr prob);
      VarBoundModPtr readVarBoundMod (ProblemPtr prob);
    private:
      std::istringstream _stream;
  };
}
#endif
