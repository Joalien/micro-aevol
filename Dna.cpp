//
// Created by arrouan on 01/10/18.
//

#include "Dna.h"
#include "ExpManager.h"

Dna::Dna(const Dna& clone) : seq_(clone.seq_) {
}

Dna::Dna(int length, Threefry::Gen& rng) {
    assert(length==5000);
    // Generate a random genome
    for (int32_t i = 0; i < length; i++) {
        seq_[i] = rng.random(NB_BASE);
    }
}

Dna::Dna(char* genome, int length) : seq_(genome) { //Assuming fixed length
    assert(length==5000);
}

Dna::Dna(int length) {
    assert(length==5000);
}

int Dna::length() const {
    return seq_.size();
}

void Dna::save(gzFile backup_file) {
    int dna_length = length();
    gzwrite(backup_file, &dna_length, sizeof(dna_length));
    gzwrite(backup_file, & seq_, seq_.size() * sizeof(seq_[0]));
}

void Dna::load(gzFile backup_file) {
    int dna_length;
    gzread(backup_file,&dna_length,sizeof(dna_length));
    gzread(backup_file, & seq_, sizeof(seq_));
}

void Dna::set(int pos, char c) {
    seq_[pos] = c;
}

void Dna::do_switch(int pos) {
    if (seq_[pos] == '0') seq_[pos] = '1';
    else seq_[pos] = '0';
}


int8_t Dna::promoter_at(int pos) {
    //We get a 22 length bitset from seq_
    std::bitset<PROM_SEQ.size()> prom_dist;
    std::bitset<PROM_SEQ.size()> seq_at_pos;
    for(int8_t motif_id=0; motif_id<PROM_SEQ.size(); motif_id++) {
        seq_at_pos[motif_id] = seq_[
                pos + motif_id >= seq_.size() ?// == (pos + motif_id) % seq_.size()
                pos + motif_id - seq_.size() :
                pos + motif_id];
    }

    //XOR operation
    prom_dist = (seq_at_pos^PROM_SEQ);

    //Debug line
//    std::cout << pos << "\n"
//    << seq_at_pos << " XOR\n"
//    << PROM_SEQ << " =\n"
//    << prom_dist << std::endl;

    return prom_dist[0] +
           prom_dist[1] +
           prom_dist[2] +
           prom_dist[3] +
           prom_dist[4] +
           prom_dist[5] +
           prom_dist[6] +
           prom_dist[7] +
           prom_dist[8] +
           prom_dist[9] +
           prom_dist[10] +
           prom_dist[11] +
           prom_dist[12] +
           prom_dist[13] +
           prom_dist[14] +
           prom_dist[15] +
           prom_dist[16] +
           prom_dist[17] +
           prom_dist[18] +
           prom_dist[19] +
           prom_dist[20] +
           prom_dist[21]; // return distance between seq_ and PROM_SEQ. 0 means equals
}

int8_t Dna::terminator_at(int pos) {
    const int TERMINATOR_SIZE = 4;
    std::bitset<TERMINATOR_SIZE> prom_dist;
    std::bitset<TERMINATOR_SIZE> seq_at_pos;
    std::bitset<TERMINATOR_SIZE> seq_at_pos_plus10;
    for(int8_t motif_id=0; motif_id<TERMINATOR_SIZE; motif_id++) {
        seq_at_pos[motif_id] = seq_[
                pos + motif_id >= seq_.size() ?// == (pos + motif_id) % seq_.size()
                pos + motif_id - seq_.size() :
                pos + motif_id];
        seq_at_pos_plus10[motif_id] = seq_[
                pos + motif_id + 10 >= seq_.size() ?// == (pos + motif_id) % seq_.size()
                pos + motif_id + 10 - seq_.size() :
                pos + motif_id + 10];
    }

    prom_dist = (seq_at_pos^seq_at_pos_plus10);

    return prom_dist[0]+
           prom_dist[1]+
           prom_dist[2]+
           prom_dist[3];
}

bool Dna::shine_dal_start(int pos) {
    std::bitset<SHINE_DAL_SEQ.size()> prom_dist;
    std::bitset<SHINE_DAL_SEQ.size()> seq_at_pos;
    int k_t;
    for(int8_t motif_id=0; motif_id<SHINE_DAL_SEQ.size(); motif_id++) {
        seq_at_pos[motif_id] = seq_[
                pos + motif_id >= seq_.size() ?// 101001****010
                pos + motif_id - seq_.size() :
                pos + motif_id];
    }

    prom_dist = (seq_at_pos^SHINE_DAL_SEQ);

    //Debug line
//    std::cout << pos << "\n"
//    << seq_at_pos << " XOR\n"
//    << SHINE_DAL_SEQ << " =\n"
//    << prom_dist << std::endl;

    return prom_dist[0] +
           prom_dist[1] +
           prom_dist[2] +
           prom_dist[3] +
           prom_dist[4] +
           prom_dist[5] +
           prom_dist[10]+
           prom_dist[11]+
           prom_dist[12]==0;
}

bool Dna::protein_stop(int pos) {
    std::bitset<PROTEIN_END.size()> prom_dist;
    std::bitset<PROTEIN_END.size()> seq_at_pos;
    for(int8_t motif_id=0; motif_id<PROTEIN_END.size(); motif_id++) {
        seq_at_pos[motif_id] = seq_[
                pos + motif_id >= seq_.size() ?// == (pos + motif_id) % seq_.size()
                pos + motif_id - seq_.size() :
                pos + motif_id];
    }

    prom_dist = (seq_at_pos^PROTEIN_END);

    return prom_dist[0] +
           prom_dist[1] +
           prom_dist[2]==0;
}

int Dna::codon_at(int pos) {
    int value = 0;

    int t_pos;

    for (int i = 0; i < 3; i++) {
        t_pos =
                static_cast<int>(pos + i >= seq_.size() ? pos + i -
                                                          seq_.size()
                                                        : pos + i);
        value += seq_[t_pos] << (CODON_SIZE - i - 1);
    }

//    std::cout << seq_[pos+0] << seq_[pos+1] << seq_[pos+2] << "\n"
//    << value << "\n";
    return value;
}