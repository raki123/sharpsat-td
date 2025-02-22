/*
 * stack.h
 *
 *  Created on: Aug 23, 2012
 *      Author: Marc Thurley
 */

#ifndef STACK_H_
#define STACK_H_

#include <gmpxx.h>
#include <vector>
#include <unordered_map>

template<class T_num>
class StackLevel {
public:
  /// active Component, once initialized, it should not change
  const unsigned super_component_ = 0;
  // branch
  bool active_branch_ = false;

  // offset in the literal stack where to store set lits
  const unsigned literal_stack_ofs_ = 0;

  //  Solutioncount
  T_num branch_model_count_[2] = {T_num::Zero(), T_num::Zero()};
  bool branch_found_unsat_[2] = {false,false};

  /// remaining Components

  // the start offset in the component stack for
  // the remaining components in this decision level
  // all remaining components can hence be found in
  // [remaining_components_ofs_, "nextLevel".remaining_components_begin_)
  const unsigned remaining_components_ofs_ = 0;

  // boundary of the stack marking which components still need to be processed
  // all components to be processed can be found in
  // [remaining_components_ofs_, unprocessed_components_end_)
  // also, all processed, can be found
  // in [unprocessed_components_end_, component_stack.size())
  unsigned unprocessed_components_end_ = 0;

  unordered_map<LiteralID, T_num> dec_weights;

  bool hasUnprocessedComponents() {
    assert(unprocessed_components_end_ >= remaining_components_ofs_);
    return unprocessed_components_end_ > remaining_components_ofs_;
  }
  void nextUnprocessedComponent() {
    assert(unprocessed_components_end_ > remaining_components_ofs_);
    unprocessed_components_end_--;
  }

  void resetRemainingComps() {
    unprocessed_components_end_ = remaining_components_ofs_;
  }
  unsigned super_component() {
    return super_component_;
  }
  unsigned remaining_components_ofs() {
    return remaining_components_ofs_;
  }
  void set_unprocessed_components_end(unsigned end) {
    unprocessed_components_end_ = end;
    assert(remaining_components_ofs_ <= unprocessed_components_end_);
  }

  StackLevel(unsigned super_comp, unsigned lit_stack_ofs,
      unsigned comp_stack_ofs) :
      super_component_(super_comp),
      literal_stack_ofs_(lit_stack_ofs),
      remaining_components_ofs_(comp_stack_ofs),
      unprocessed_components_end_(comp_stack_ofs) {
    assert(super_comp < comp_stack_ofs);
  }

  unsigned currentRemainingComponent() {
    assert(remaining_components_ofs_ <= unprocessed_components_end_ - 1);
    return unprocessed_components_end_ - 1;
  }
  bool isSecondBranch() {
    return active_branch_;
  }

  void changeBranch() {
    active_branch_ = true;
  }

  bool anotherCompProcessible() {
    return (!branch_found_unsat()) && hasUnprocessedComponents();
  }

  unsigned literal_stack_ofs() {
    return literal_stack_ofs_;
  }

  void includeSolution(const T_num& solutions);

  bool branch_found_unsat() const;
  void mark_branch_unsat() {
    branch_found_unsat_[active_branch_] = true;
  }

//  void set_both_branches_unsat(){
//	  branch_found_unsat_[0] =
//			  branch_found_unsat_[1] = true;
//	  branch_model_count_[0] = branch_model_count_[1] = 0;
//	  active_branch_ = 1;
//  }
  const T_num getTotalModelCount() const;
};

template<typename T_num>
inline bool StackLevel<T_num>::branch_found_unsat() const {
  return branch_found_unsat_[active_branch_];
}


template <typename T_num>
inline void StackLevel<T_num>::includeSolution(const T_num& solutions) {
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_].IsAlgZero());
    return;
  }
  if (solutions.IsAlgZero()) {
    branch_found_unsat_[active_branch_] = true;
  }
  if (branch_model_count_[active_branch_].IsAlgZero()) {
    T_num factor = T_num::One();
    for(auto it : dec_weights) {
      factor *= it.second;
    }
    branch_model_count_[active_branch_] = solutions * factor;
  }
  else {
    branch_model_count_[active_branch_] *= solutions;
  }
}


template <>
inline void StackLevel<dDNNFNode>::includeSolution(const dDNNFNode& solutions) {
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_].IsAlgZero());
    return;
  }
  if (solutions.IsAlgZero()) {
    branch_found_unsat_[active_branch_] = true;
  }
  if (branch_model_count_[active_branch_].IsAlgZero()) {
    bool has = solutions.id != 0;
    int nr_relevant = 0;
    for(auto it : dec_weights) {
      if(it.second.id != 1) {
        nr_relevant++;
        has &= (it.second.id != 0);
      }
    }
    if(has) {
      if(nr_relevant > 0 && solutions.id != 1) {
        dDNNFNode::buffer.push_back('A');
        dDNNFNode::buffer.push_back(' ');
        dDNNFNode::WriteID(nr_relevant + 1);
        dDNNFNode::buffer.push_back(' ');
        dDNNFNode::WriteID(solutions.id);
        for(auto it : dec_weights) {
          if(it.second.id != 1) {
            dDNNFNode::buffer.push_back(' ');
            dDNNFNode::WriteID(it.second.id);
          }
        }
        dDNNFNode::buffer.push_back('\n');
        dDNNFNode::edges += 1 + nr_relevant;
        branch_model_count_[active_branch_].id = dDNNFNode::nodes++;
      } else if(nr_relevant > 1) {
        dDNNFNode::buffer.push_back('A');
        dDNNFNode::buffer.push_back(' ');
        dDNNFNode::WriteID(nr_relevant);
        for(auto it : dec_weights) {
          if(it.second.id != 1) {
            dDNNFNode::buffer.push_back(' ');
            dDNNFNode::WriteID(it.second.id);
          }
        }
        dDNNFNode::buffer.push_back('\n');
        dDNNFNode::edges += nr_relevant;
        branch_model_count_[active_branch_].id = dDNNFNode::nodes++;
      } else if(nr_relevant == 1) {
        for(auto it : dec_weights) {
          if(it.second.id != 1) {
            branch_model_count_[active_branch_].id = it.second.id;
          }
        }
      } else if(solutions.id != 1) {
        branch_model_count_[active_branch_].id = solutions.id;
      } else {
        branch_model_count_[active_branch_].id = 1;
      }
    }
  }
  else {
    branch_model_count_[active_branch_] *= solutions;
  }
}

template <>
inline void StackLevel<instantdDNNFNode>::includeSolution(const instantdDNNFNode& solutions) {
  if (branch_found_unsat_[active_branch_]) {
    assert(branch_model_count_[active_branch_].IsAlgZero());
    return;
  }
  if (solutions.IsAlgZero()) {
    branch_found_unsat_[active_branch_] = true;
  }
  if (branch_model_count_[active_branch_].IsAlgZero()) {
    bool has = solutions.id != 0;
    int nr_relevant = 0;
    for(auto it : dec_weights) {
      if(it.second.id != 1) {
        nr_relevant++;
        has &= (it.second.id != 0);
      }
    }
    if(has) {
      if(nr_relevant > 0 && solutions.id != 1) {
        *instantdDNNFNode::out << "A " << nr_relevant + 1 << " " << solutions.id;
        for(auto it : dec_weights) {
          if(it.second.id != 1) {
            *instantdDNNFNode::out << " " << it.second.id;
          }
        }
        *instantdDNNFNode::out << "\n";
        instantdDNNFNode::edges += 1 + nr_relevant;
        branch_model_count_[active_branch_].id = instantdDNNFNode::nodes++;
      } else if(nr_relevant > 1) {
        *instantdDNNFNode::out << "A " << nr_relevant;
        for(auto it : dec_weights) {
          if(it.second.id != 1) {
            *instantdDNNFNode::out << " " << it.second.id;
          }
        }
        *instantdDNNFNode::out << "\n";
        instantdDNNFNode::edges += nr_relevant;
        branch_model_count_[active_branch_].id = instantdDNNFNode::nodes++;
      } else if(nr_relevant == 1) {
        for(auto it : dec_weights) {
          if(it.second.id != 1) {
            branch_model_count_[active_branch_].id = it.second.id;
          }
        }
      } else if(solutions.id != 1) {
        branch_model_count_[active_branch_].id = solutions.id;
      } else {
        branch_model_count_[active_branch_].id = 1;
      }
    }
  }
  else {
    branch_model_count_[active_branch_] *= solutions;
  }
}

template <typename T_num>
inline const T_num StackLevel<T_num>::getTotalModelCount() const {
  return branch_model_count_[0] + branch_model_count_[1];
}

template<class T_num>
class DecisionStack: public vector<StackLevel<T_num>> {
public:

  StackLevel<T_num> &top() {
    assert(vector<StackLevel<T_num>>::size() > 0);
    return vector<StackLevel<T_num>>::back();
  }
  int get_decision_level() const {
    assert(vector<StackLevel<T_num>>::size() > 0);
    return vector<StackLevel<T_num>>::size() - 1;
  } // 0 means pre-1st-decision

};



#endif /* STACK_H_ */
