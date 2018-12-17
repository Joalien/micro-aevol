//
// Created by arrouan on 01/10/18.
//

#pragma once

#include <cstdio>
#include <cstring>
#include <cassert>
#include <cstdint>
#include <vector>
#include <zlib.h>
#include <bitset>

#include "Threefry.h"

constexpr int8_t CODON_SIZE = 3;

const std::bitset<22> PROM_SEQ (std::string("0101011001110010010110"));
const std::bitset<9> SHINE_DAL_SEQ (std::string("0110110000000"));
const std::bitset<3> PROTEIN_END (std::string("001")); // CODON_STOP

class ExpManager;

class Dna {

 public:
  Dna() = default;

  Dna(const Dna& clone);

  Dna(int length, Threefry::Gen& rng);

  Dna(char* genome, int length);

  Dna(int length);

  ~Dna() = default;

  int length() const;

  void save(gzFile backup_file);
  void load(gzFile backup_file);

  void set(int pos, char c);

  void do_switch(int pos);

  int8_t promoter_at(int pos);

  int8_t terminator_at(int pos);

  bool shine_dal_start(int pos);

  bool protein_stop(int pos);

  int codon_at(int pos);

  std::bitset<5000> seq_;
};
