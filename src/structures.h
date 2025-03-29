/*
 * structures.h
 *
 *  Created on: Jun 25, 2012
 *      Author: Marc Thurley
 */

#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <vector>
#include <set>
#include <iostream>
#include <memory>
#include <cmath>
#include <limits>
#include <cassert>
#include "primitive_types.h"
#include <gmpxx.h>
#include "mpfr/mpreal.h"
using namespace std;

struct SDouble
{
public:
  SDouble()
  {
    n = 0;
    has = false;
  }
  SDouble(const SDouble &other)
  {
    n = other.n;
    has = other.has;
  }
  SDouble &operator=(const SDouble &other)
  {
    n = other.n;
    has = other.has;
    return *this;
  }
  bool operator==(const SDouble &other) const
  {
    return n == other.n && has == other.has;
  }
  bool IsAlgZero() const
  {
    return !has;
  }
  SDouble operator*(SDouble other) const
  {
    SDouble ret = other;
    ret.n *= n;
    ret.has &= has;
    return ret;
  }
  SDouble operator+(SDouble other) const
  {
    SDouble ret = other;
    ret.n += n;
    ret.has |= has;
    return ret;
  }
  SDouble &operator*=(const SDouble &other)
  {
    n *= other.n;
    has &= other.has;
    return *this;
  }
  size_t InternalSize() const
  {
    return 0;
  }
  double Get() const
  {
    return n;
  }
  static SDouble Zero()
  {
    SDouble ret;
    return ret;
  }
  static SDouble One()
  {
    SDouble ret;
    ret.n = 1;
    ret.has = true;
    return ret;
  }
  static SDouble FromString(string s)
  {
    SDouble ret;
    ret.n = stod(s);
    ret.has = ret.n != 0;
    return ret;
  }
  static bool Idempotent(const SDouble d)
  {
    return !d.has;
  }

private:
  double n = 0;
  bool has = false;
};

template <typename Number>
struct Complex
{
public:
  Complex()
  {
    real = 0;
    imag = 0;
    has = false;
  }
  Complex(const Complex &other)
  {
    real = other.real;
    imag = other.imag;
    has = other.has;
  }
  Complex &operator=(const Complex &other)
  {
    real = other.real;
    imag = other.imag;
    has = other.has;
    return *this;
  }
  bool operator==(const Complex &other) const
  {
    return real == other.real && imag == other.imag && has == other.has;
  }
  bool IsAlgZero() const
  {
    return !has;
  }
  Complex operator*(Complex other) const
  {
    Complex ret = other;
    ret.real = real * other.real - imag * other.imag;
    ret.imag = real * other.imag + other.real * imag;
    ret.has &= has;
    return ret;
  }
  Complex operator+(Complex other) const
  {
    Complex ret = other;
    ret.real += real;
    ret.imag += imag;
    ret.has |= has;
    return ret;
  }
  Complex &operator*=(const Complex &other)
  {
    auto old_real = real;
    real = real * other.real - imag * other.imag;
    imag = old_real * other.imag + other.real * imag;
    has &= other.has;
    return *this;
  }
  size_t InternalSize() const
  {
    return 0;
  }
  Number Get(size_t i) const
  {
    return i == 0 ? real : imag;
  }
  static Complex Zero()
  {
    Complex ret;
    return ret;
  }
  static Complex One()
  {
    Complex ret;
    ret.real = 1;
    ret.imag = 0;
    ret.has = true;
    return ret;
  }
  static Complex FromString(string s)
  {
    Complex ret;
    size_t start = 0;
    size_t end = s.find(';', 1);
    ret.real = stod(s.substr(start, end - start));

    start = end + 1;
    end = s.find(';', end + 1);
    ret.imag = stod(s.substr(start, end - start));

    ret.has = ret.real != 0 || ret.imag != 0;
    return ret;
  }
  static bool Idempotent(const Complex d)
  {
    return !d.has;
  }

private:
  Number real = 0;
  Number imag = 0;
  bool has = false;
};

struct instantdDNNFNode
{
public:
  static uint64_t nodes;
  static uint64_t edges;
  static ostream *out;
  instantdDNNFNode()
  {
    id = 0;
  }
  instantdDNNFNode(const instantdDNNFNode &other)
  {
    id = other.id;
  }
  instantdDNNFNode &operator=(const instantdDNNFNode &other)
  {
    id = other.id;
    return *this;
  }

  bool operator==(const instantdDNNFNode &other) const
  {
    return id == other.id;
  }

  bool IsAlgZero() const
  {
    return id == 0;
  }

  instantdDNNFNode operator*(instantdDNNFNode other) const
  {
    if (other.IsAlgZero() || id == 1)
    {
      return other;
    }
    if (IsAlgZero() || other.id == 1)
    {
      return *this;
    }
    *out << "A 2 " << other.id << " " << id << "\n";
    instantdDNNFNode ret;
    edges += 2;
    ret.id = nodes++;
    return ret;
  }
  instantdDNNFNode operator+(instantdDNNFNode other) const
  {
    if (IsAlgZero() || other.id == 1)
    {
      return other;
    }
    if (other.IsAlgZero() || id == 1)
    {
      return *this;
    }
    *out << "O 0 2 " << other.id << " " << id << "\n";
    instantdDNNFNode ret;
    edges += 2;
    ret.id = nodes++;
    return ret;
  }
  instantdDNNFNode &operator*=(const instantdDNNFNode &other)
  {
    if (other.IsAlgZero() || id == 1)
    {
      id = other.id;
      return *this;
    }
    if (IsAlgZero() || other.id == 1)
    {
      return *this;
    }
    *out << "A 2 " << other.id << " " << id << "\n";
    edges += 2;
    id = nodes++;
    return *this;
  }

  size_t InternalSize() const
  {
    return 0;
  }
  static instantdDNNFNode Zero()
  {
    instantdDNNFNode ret;
    ret.id = 0;
    return ret;
  }
  static instantdDNNFNode One()
  {
    instantdDNNFNode ret;
    ret.id = 1;
    return ret;
  }
  static instantdDNNFNode FromString(string s)
  {
    *out << "L " << s << "\n";
    instantdDNNFNode ret;
    ret.id = nodes++;
    return ret;
  }
  static bool Idempotent(const instantdDNNFNode d)
  {
    return d.id < 2;
  }
  uint64_t id;
};

uint64_t instantdDNNFNode::nodes = 2;
uint64_t instantdDNNFNode::edges = 0;
ostream *instantdDNNFNode::out;

struct dDNNFNode
{
public:
  static uint64_t nodes;
  static uint64_t edges;
  static vector<uint8_t> buffer;
  dDNNFNode()
  {
    id = 0;
  }
  dDNNFNode(const dDNNFNode &other)
  {
    id = other.id;
  }
  dDNNFNode &operator=(const dDNNFNode &other)
  {
    id = other.id;
    return *this;
  }

  bool operator==(const dDNNFNode &other) const
  {
    return id == other.id;
  }

  bool IsAlgZero() const
  {
    return id == 0;
  }

  dDNNFNode operator*(dDNNFNode other) const
  {
    if (other.IsAlgZero() || id == 1)
    {
      return other;
    }
    if (IsAlgZero() || other.id == 1)
    {
      return *this;
    }
    buffer.push_back('A');
    buffer.push_back(' ');
    buffer.push_back('2');
    buffer.push_back(' ');
    WriteID(other.id);
    buffer.push_back(' ');
    WriteID(id);
    buffer.push_back('\n');
    dDNNFNode ret;
    edges += 2;
    ret.id = nodes++;
    return ret;
  }
  dDNNFNode operator+(dDNNFNode other) const
  {
    if (IsAlgZero() || other.id == 1)
    {
      return other;
    }
    if (other.IsAlgZero() || id == 1)
    {
      return *this;
    }

    buffer.push_back('O');
    buffer.push_back(' ');
    buffer.push_back('0');
    buffer.push_back(' ');
    buffer.push_back('2');
    buffer.push_back(' ');
    WriteID(other.id);
    buffer.push_back(' ');
    WriteID(id);
    buffer.push_back('\n');
    dDNNFNode ret;
    edges += 2;
    ret.id = nodes++;
    return ret;
  }
  dDNNFNode &operator*=(const dDNNFNode &other)
  {
    if (other.IsAlgZero() || id == 1)
    {
      id = other.id;
      return *this;
    }
    if (IsAlgZero() || other.id == 1)
    {
      return *this;
    }
    buffer.push_back('A');
    buffer.push_back(' ');
    buffer.push_back('2');
    buffer.push_back(' ');
    WriteID(other.id);
    buffer.push_back(' ');
    WriteID(id);
    buffer.push_back('\n');
    edges += 2;
    id = nodes++;
    return *this;
  }

  size_t InternalSize() const
  {
    return 0;
  }
  static dDNNFNode Zero()
  {
    dDNNFNode ret;
    ret.id = 0;
    return ret;
  }
  static dDNNFNode One()
  {
    dDNNFNode ret;
    ret.id = 1;
    return ret;
  }
  static dDNNFNode FromString(string s)
  {
    int64_t lit = stoll(s.c_str());
    if (lit > 0)
    {
      buffer.push_back('L');
      buffer.push_back(' ');
    }
    else
    {
      buffer.push_back('L');
      buffer.push_back(' ');
      buffer.push_back('-');
    }
    WriteID(abs(lit));
    buffer.push_back('\n');
    dDNNFNode ret;
    ret.id = nodes++;
    return ret;
  }
  static bool Idempotent(const dDNNFNode d)
  {
    return d.id < 2;
  }
  uint64_t id;

  // buffer bookkeeping
  static void WriteID(uint64_t to_write);
};

uint64_t dDNNFNode::nodes = 2;
uint64_t dDNNFNode::edges = 0;
vector<uint8_t> dDNNFNode::buffer;

void dDNNFNode::WriteID(uint64_t to_write)
{
  int8_t idx = 18;
  uint64_t div = 1000000000000000000u;
  while (to_write / div == 0 && idx > 0)
  {
    idx--;
    div /= 10;
  }
  while (idx >= 0)
  {
    uint8_t nibble = to_write / div;
    to_write -= nibble * div;
    buffer.push_back(nibble + 48);
    idx--;
    div /= 10;
  }
  assert(to_write == 0);
}

struct Mmpr
{
public:
  static size_t N;
  Mmpr()
  {
    n.resize(N);
    has = false;
  }
  Mmpr(const Mmpr &other)
  {
    n.resize(N);
    for (size_t i = 0; i < N; i++)
    {
      n[i] = other.n[i];
    }
    has = other.has;
  }
  Mmpr &operator=(const Mmpr &other)
  {
    for (size_t i = 0; i < N; i++)
    {
      n[i] = other.n[i];
    }
    has = other.has;
    return *this;
  }
  bool operator==(const Mmpr &other) const
  {
    for (size_t i = 0; i < N; i++)
    {
      if (n[i] != other.n[i])
      {
        return false;
      }
    }
    return has == other.has;
  }
  bool IsAlgZero() const
  {
    return !has;
  }
  Mmpr operator*(Mmpr other) const
  {
    Mmpr ret(other);
    for (size_t i = 0; i < N; i++)
    {
      ret.n[i] *= n[i];
    }
    ret.has = has && other.has;
    return ret;
  }
  Mmpr operator+(Mmpr other) const
  {
    Mmpr ret(other);
    for (size_t i = 0; i < N; i++)
    {
      ret.n[i] += n[i];
    }
    ret.has = has || other.has;
    return ret;
  }
  Mmpr &operator*=(const Mmpr &other)
  {
    for (size_t i = 0; i < N; i++)
    {
      n[i] *= other.n[i];
    }
    has &= other.has;
    return *this;
  }
  size_t InternalSize() const
  {
    return 0;
  }
  mpfr::mpreal Get(size_t i) const
  {
    assert(i < N);
    return n[i];
  }
  static Mmpr Zero()
  {
    Mmpr ret;
    for (size_t i = 0; i < N; i++)
    {
      ret.n[i] = 0;
    }
    ret.has = false;
    return ret;
  }
  static Mmpr One()
  {
    Mmpr ret;
    for (size_t i = 0; i < N; i++)
    {
      ret.n[i] = 1;
    }
    ret.has = true;
    return ret;
  }
  static Mmpr FromString(string s)
  {
    Mmpr ret;
    size_t start = 0;
    size_t end = 0;
    for (size_t i = 0; i < N - 1; i++)
    {
      end = s.find(';', end + 1);
      ret.n[i] = stod(s.substr(start, end - start));
      start = end + 1;
    }
    ret.n[N - 1] = stod(s.substr(start));
    ret.has = false;
    for (size_t i = 0; i < N; i++)
    {
      ret.has |= (ret.n[i] != 0);
    }
    return ret;
  }
  static bool Idempotent(const Mmpr d)
  {
    return d.IsAlgZero();
  }

private:
  std::vector<mpfr::mpreal> n;
  bool has;
};
size_t Mmpr::N = 0;

struct MDouble
{
public:
  static size_t N;
  MDouble()
  {
    n.resize(N);
    has = false;
  }
  MDouble(const MDouble &other)
  {
    n.resize(N);
    for (size_t i = 0; i < N; i++)
    {
      n[i] = other.n[i];
    }
    has = other.has;
  }
  MDouble &operator=(const MDouble &other)
  {
    for (size_t i = 0; i < N; i++)
    {
      n[i] = other.n[i];
    }
    has = other.has;
    return *this;
  }
  bool operator==(const MDouble &other) const
  {
    for (size_t i = 0; i < N; i++)
    {
      if (n[i] != other.n[i])
      {
        return false;
      }
    }
    return has == other.has;
  }
  bool IsAlgZero() const
  {
    return !has;
  }
  MDouble operator*(MDouble other) const
  {
    MDouble ret(other);
    for (size_t i = 0; i < N; i++)
    {
      ret.n[i] *= n[i];
    }
    ret.has = has && other.has;
    return ret;
  }
  MDouble operator+(MDouble other) const
  {
    MDouble ret(other);
    for (size_t i = 0; i < N; i++)
    {
      ret.n[i] += n[i];
    }
    ret.has = has || other.has;
    return ret;
  }
  MDouble &operator*=(const MDouble &other)
  {
    for (size_t i = 0; i < N; i++)
    {
      n[i] *= other.n[i];
    }
    has &= other.has;
    return *this;
  }
  size_t InternalSize() const
  {
    return 0;
  }
  double Get(size_t i) const
  {
    assert(i < N);
    return n[i];
  }
  static MDouble Zero()
  {
    MDouble ret;
    for (size_t i = 0; i < N; i++)
    {
      ret.n[i] = 0;
    }
    ret.has = false;
    return ret;
  }
  static MDouble One()
  {
    MDouble ret;
    for (size_t i = 0; i < N; i++)
    {
      ret.n[i] = 1;
    }
    ret.has = true;
    return ret;
  }
  static MDouble FromString(string s)
  {
    MDouble ret;
    size_t start = 0;
    size_t end = 0;
    for (size_t i = 0; i < N - 1; i++)
    {
      end = s.find(';', end + 1);
      ret.n[i] = stod(s.substr(start, end - start));
      start = end + 1;
    }
    ret.n[N - 1] = stod(s.substr(start));
    ret.has = false;
    for (size_t i = 0; i < N; i++)
    {
      ret.has |= (ret.n[i] != 0);
    }
    return ret;
  }
  static bool Idempotent(const MDouble d)
  {
    return d.IsAlgZero();
  }

private:
  std::vector<double> n;
  bool has;
};
size_t MDouble::N = 0;

struct Smpr
{
public:
  Smpr()
  {
    n = 0;
    has = false;
  }
  Smpr(const Smpr &other)
  {
    n = other.n;
    has = other.has;
  }
  Smpr &operator=(const Smpr &other)
  {
    n = other.n;
    has = other.has;
    return *this;
  }
  bool operator==(const Smpr &other) const
  {
    return n == other.n && has == other.has;
  }
  bool IsAlgZero() const
  {
    return !has;
  }
  Smpr operator*(Smpr other) const
  {
    Smpr ret = other;
    ret.n *= n;
    ret.has &= has;
    return ret;
  }
  Smpr operator+(Smpr other) const
  {
    Smpr ret = other;
    ret.n += n;
    ret.has |= has;
    return ret;
  }
  Smpr &operator*=(const Smpr &other)
  {
    n *= other.n;
    has &= other.has;
    return *this;
  }
  size_t InternalSize() const
  {
    return 0;
  }
  mpfr::mpreal Get() const
  {
    return n;
  }
  static Smpr Zero()
  {
    Smpr ret;
    return ret;
  }
  static Smpr One()
  {
    Smpr ret;
    ret.n = 1;
    ret.has = true;
    return ret;
  }
  static Smpr FromString(string s)
  {
    Smpr ret;
    ret.n = stod(s);
    ret.has = ret.n != 0;
    return ret;
  }
  static bool Idempotent(const Smpr d)
  {
    return d.IsAlgZero();
  }

private:
  mpfr::mpreal n = 0;
  bool has = false;
};

struct Smpz
{
public:
  Smpz()
  {
    n = 0;
    has = false;
  }
  Smpz(const Smpz &other)
  {
    n = other.n;
    has = other.has;
  }
  Smpz &operator=(const Smpz &other)
  {
    n = other.n;
    has = other.has;
    return *this;
  }
  bool operator==(const Smpz &other) const
  {
    return n == other.n && has == other.has;
  }
  bool IsAlgZero() const
  {
    return !has;
  }
  Smpz operator*(Smpz other) const
  {
    Smpz ret = other;
    ret.n *= n;
    ret.has &= has;
    return ret;
  }
  Smpz operator+(Smpz other) const
  {
    Smpz ret = other;
    ret.n += n;
    ret.has |= has;
    return ret;
  }
  Smpz &operator*=(const Smpz &other)
  {
    n *= other.n;
    has &= other.has;
    return *this;
  }
  size_t InternalSize() const
  {
    return n.get_mpz_t()->_mp_alloc * sizeof(mp_limb_t);
  }
  mpz_class Get() const
  {
    return n;
  }
  static Smpz Zero()
  {
    Smpz ret;
    return ret;
  }
  static Smpz One()
  {
    Smpz ret;
    ret.n = 1;
    ret.has = true;
    return ret;
  }
  static Smpz FromString(string s)
  {
    Smpz ret;
    ret.n = stoi(s);
    ret.has = ret.n != 0;
    return ret;
  }
  static bool Idempotent(const Smpz d)
  {
    return d.IsAlgZero();
  }

private:
  mpz_class n = 0;
  bool has = false;
};

#define INVALID_DL -1

typedef unsigned char TriValue;
#define F_TRI 0
#define T_TRI 1
#define X_TRI 2

class LiteralID
{
public:
  LiteralID()
  {
    value_ = 0;
  }
  LiteralID(int lit)
  {
    value_ = (abs(lit) << 1) + (unsigned)(lit > 0);
  }

  LiteralID(VariableIndex var, bool sign)
  {
    value_ = (var << 1) + (unsigned)sign;
  }

  VariableIndex var() const
  {
    return (value_ >> 1);
  }

  int toInt() const
  {
    return ((int)value_ >> 1) * ((sign()) ? 1 : -1);
  }

  void inc() { ++value_; }

  void copyRaw(unsigned int v)
  {
    value_ = v;
  }

  bool sign() const
  {
    return (bool)(value_ & 0x01);
  }

  bool operator!=(const LiteralID &rL2) const
  {
    return value_ != rL2.value_;
  }

  bool operator==(const LiteralID &rL2) const
  {
    return value_ == rL2.value_;
  }

  const LiteralID neg() const
  {
    return LiteralID(var(), !sign());
  }

  void print() const
  {
    cout << (sign() ? " " : "-") << var() << " ";
  }

  unsigned raw() const { return value_; }

private:
  unsigned value_;

  template <class _T>
  friend class LiteralIndexedVector;
  friend struct std::hash<LiteralID>;
};

namespace std
{
  template <>
  struct hash<LiteralID>
  {
    size_t operator()(const LiteralID &p) const
    {
      return p.value_;
    }
  };
}

static const LiteralID NOT_A_LIT(0, false);
#define SENTINEL_LIT NOT_A_LIT

class Literal
{
public:
  vector<LiteralID> binary_links_ = vector<LiteralID>(1, SENTINEL_LIT);
  vector<ClauseOfs> watch_list_ = vector<ClauseOfs>(1, SENTINEL_CL);
  float activity_score_ = 0.0f;

  void increaseActivity(unsigned u = 1)
  {
    activity_score_ += u;
  }

  void removeWatchLinkTo(ClauseOfs clause_ofs)
  {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
      if (*it == clause_ofs)
      {
        *it = watch_list_.back();
        watch_list_.pop_back();
        return;
      }
  }

  void replaceWatchLinkTo(ClauseOfs clause_ofs, ClauseOfs replace_ofs)
  {
    for (auto it = watch_list_.begin(); it != watch_list_.end(); it++)
      if (*it == clause_ofs)
      {
        *it = replace_ofs;
        return;
      }
  }

  void addWatchLinkTo(ClauseIndex clause_ofs)
  {
    watch_list_.push_back(clause_ofs);
  }

  void addBinLinkTo(LiteralID lit)
  {
    binary_links_.back() = lit;
    binary_links_.push_back(SENTINEL_LIT);
  }

  void resetWatchList()
  {
    watch_list_.clear();
    watch_list_.push_back(SENTINEL_CL);
  }

  bool hasBinaryLinkTo(LiteralID lit)
  {
    for (auto l : binary_links_)
    {
      if (l == lit)
        return true;
    }
    return false;
  }

  bool hasBinaryLinks()
  {
    return !binary_links_.empty();
  }
};

class Antecedent
{
  unsigned int val_;

public:
  Antecedent()
  {
    val_ = 1;
  }

  Antecedent(const ClauseOfs cl_ofs)
  {
    val_ = (cl_ofs << 1) | 1;
  }
  Antecedent(const LiteralID idLit)
  {
    val_ = (idLit.raw() << 1);
  }

  bool isAClause() const
  {
    return val_ & 0x01;
  }

  ClauseOfs asCl() const
  {
    return val_ >> 1;
  }

  LiteralID asLit()
  {
    LiteralID idLit;
    idLit.copyRaw(val_ >> 1);
    return idLit;
  }
  // A NON-Antecedent will only be A NOT_A_CLAUSE Clause Id
  bool isAnt()
  {
    return val_ != 1; // i.e. NOT a NOT_A_CLAUSE;
  }
};

struct Variable
{
  Antecedent ante;
  int decision_level = INVALID_DL;
};

// for now Clause Header is just a dummy
// we keep it for possible later changes
class ClauseHeader
{
  unsigned creation_time_; // number of conflicts seen at creation time
  unsigned glue_;
  unsigned length_;

public:
  void setGlue(unsigned glue)
  {
    glue_ = glue;
  }
  void see()
  {
    glue_ += 3;
  }
  unsigned glue() const
  {
    return glue_;
  }
  bool isLearned()
  {
    return glue_ >= 1;
  }
  unsigned creation_time()
  {
    return creation_time_;
  }
  unsigned length() { return length_; }
  void set_length(unsigned length) { length_ = length; }

  void set_creation_time(unsigned time)
  {
    creation_time_ = time;
  }
  static unsigned overheadInLits() { return sizeof(ClauseHeader) / sizeof(LiteralID); }
};

#endif /* STRUCTURES_H_ */
