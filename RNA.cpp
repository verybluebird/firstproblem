#include <cstddef>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "RNA.h"
#include <chrono>

using namespace std;





RNA::reference:: reference(RNA *const rna, size_t ind) : rna(rna), ind(ind) {};
size_t RNA::reference:: change_block(size_t index_of_block, Nucleotid n, size_t index) {
    size_t shift_L = sizeof(size_t) * 8 - (index) * 2;
    size_t shift_R = (index + 1) * 2;
    size_t ed = O;
    size_t ed_R = B;
    if (shift_R != 64) {
        ed_R = ed >> shift_R;
    }
    //cout << "ed_r" << "\n";
    //out(ed_R);
    size_t ed_L = B;
    if (shift_L != 64) {
        ed_L = ed << shift_L;
    }
    //cout << "ed_l" << "\n";
    //out(ed_L);
    // cout << "mask1" << "\n";
    size_t mask1 = B | ed_L | ed_R;
    // out(mask1);
    size_t block = rna->rnaarr[index_of_block];
    block = rna->rnaarr[index_of_block] & mask1;
    size_t shift = 2 * (sizeof(size_t) * 4 - index - 1);
    size_t new_nucl = n;
    //cout<<"oldnucl\n";
    // out(n);
    //cout<<"newnucl\n";
    new_nucl = new_nucl << shift;
    //out(new_nucl);
    block = block | new_nucl;
    // cout<<"block\n";
    // out (block);
    return block;
}



RNA::reference & RNA::reference::operator=(Nucleotid nucl) {//reference
    size_t index_of_last_old_nucl_in_block = (rna->NuclNum - 1) % (sizeof(size_t) * 4);
    size_t index_of_last_old_nucl=(rna->NuclNum - 1);
    if (rna->NuclNum == 0) {
        index_of_last_old_nucl_in_block = 0;
        index_of_last_old_nucl=0;
    }
    if (rna->rnaarr == nullptr) {
        ++rna->NuclNum;
    }
    size_t index_of_changing_nucl_in_block = ind % (sizeof(size_t) * 4);
    size_t index_of_changing_block = (ind) / (sizeof(size_t) * 4);
    //size_t old_nucl=Nucleotid(reference(rna,index_of_changing_nucl_in_block));

    if ((rna->length < (index_of_changing_block + 1) && index_of_changing_nucl_in_block == 0) or
        rna->length == 0) {
        size_t *new_rnaarr;
        if (index_of_changing_block != 0) {
            new_rnaarr = new size_t[index_of_changing_block * 2];
        } else {
            new_rnaarr = new size_t[1];
        }
        for (size_t i = 0; i < rna->length; i++) {
            new_rnaarr[i] = rna->rnaarr[i];
        }
        delete[] rna->rnaarr;
        if (index_of_changing_block != 0) {
            rna->length = 2 * index_of_changing_block;
        } else {
            rna->length = 1;
        }
        rna->rnaarr = new_rnaarr;
    } else {
        assert(index_of_changing_nucl_in_block <= index_of_last_old_nucl_in_block + 1);
    }
    if (rna->rnaarr == nullptr) {
        rna->rnaarr = new size_t[1];
        ++rna->length;
        *rna->rnaarr = 0;
    }

    rna->rnaarr[index_of_changing_block] = change_block(index_of_changing_block, nucl,
                                                        index_of_changing_nucl_in_block);
    size_t index_of_changing_nucl=(index_of_changing_block)*(sizeof(size_t) * 4)+index_of_changing_nucl_in_block;
    if (index_of_last_old_nucl < index_of_changing_nucl) {
        rna->NuclNum++;
    }
    return *this;
}

static size_t get_mask(size_t index) {
    size_t mask = T;
    size_t shift = 2 * (sizeof(size_t) * 4 - index - 1);
    mask = mask << shift;
    return mask;
}


RNA::reference:: operator Nucleotid() {
    size_t k = sizeof(size_t);
    size_t index_of_block = ind / (k * 4);
    size_t index_of_nucl_in_block = ind % (sizeof(size_t) * 4);


    size_t shift = 2 * (sizeof(size_t) * 4 - index_of_nucl_in_block - 1);

    size_t block = *(rna->rnaarr + index_of_block);
    size_t new_block = block & get_mask(index_of_nucl_in_block);
    new_block = new_block >> shift;

    if (new_block == A) {
        return A;
    }
    if (new_block == G) {
        return G;
    }
    if (new_block == C) {
        return C;
    }
    if (new_block == T) {
        return T;
    }


}

void RNA::reference::out(size_t n) {
        for (int i = sizeof(size_t) * 8 - 1; i > 0; --i) // short 16 бит (00000000 00000000)
            cout << (n >> i & 1); // проходимся по битам и выводим через пробел
        cout << (n & 1) << " \n";
    }

size_t RNA::reference::get_mask(size_t index) {
    size_t mask = T;
    size_t shift = 2 * (sizeof(size_t) * 4 - index - 1);
    mask = mask << shift;
    return mask;
}


RNA::reference RNA::operator[](size_t ind) {
    assert(ind >= 0);
    return reference{this, ind};
};


RNA::RNA() : rnaarr(nullptr), NuclNum(0), length(0) {}

RNA::RNA(size_t NuclNum, size_t length) : NuclNum(NuclNum), length(length) {
    rnaarr = new size_t[length];
}

RNA::RNA(size_t num, Nucleotid n) {
    size_t nucl = n;
    assert(num >= 0);
    if (num > 0) {
        this->NuclNum = num;
        this->length = (num + 4 * sizeof(size_t) - 1) / (4 * sizeof(size_t));
        this->rnaarr = new size_t[length]; //memory allocation

        for (size_t i = 0; i < num; i++) {
            nucl = n;
            size_t index_of_changing_block = (i) / (sizeof(size_t) * 4);
            size_t index_of_nucl_in_block = i % (sizeof(size_t) * 4);
            size_t shift = 2 * (sizeof(size_t) * 4 - index_of_nucl_in_block - 1);
            nucl = nucl << shift;
            *(this->rnaarr + index_of_changing_block) = *(this->rnaarr + index_of_changing_block) | nucl;

        }
    } else {
        this->NuclNum = num;
        this->length = num;
        this->rnaarr = nullptr;
    }
}

RNA::~RNA() {
    delete[] rnaarr;
}

std::ostream &operator<<(ostream &out, RNA &other) {//вывод рнк на экран
    for (size_t i = 0; i < other.NuclNum; i++) {

        Nucleotid k = Nucleotid(other[i]);
        switch (k) {
            case A:
                cout << "A";
                break;
            case C:
                cout << "C";
                break;
            case G:
                cout << "G";
                break;
            case T:
                cout << "T";
                break;
            default:
                break;
        }
    }

    cout << "\n";
    return cout;
}

RNA::RNA(const RNA &rna) : NuclNum(rna.NuclNum), rnaarr(new size_t[rna.length]),
                           length(rna.length) {//copy-constr memory allocation
    std::copy(rna.rnaarr, rna.rnaarr + rna.length, rnaarr);
}

RNA &RNA::operator=(const RNA &r) { //к одной рнк приравниваем другую
    if (this->length == r.length) {
        for (size_t i = 0; i < r.length; i++) {
            this->rnaarr[i] = r.rnaarr[i];
        }
        this->NuclNum = r.NuclNum;
    } else {
        delete[] this->rnaarr;
        this->rnaarr = new size_t[r.length]; //memory allocation
        for (size_t i = 0; i < r.length; i++) { //copy data
            this->rnaarr[i] = r.rnaarr[i];
        }
        this->NuclNum = r.NuclNum;
        this->length = r.length;
    }
    return *this;
}


size_t RNA:: getsize(const RNA &r) { //возвращает размер цепочки в количествах нуклеотидов
    return r.NuclNum;
}

void RNA::add(Nucleotid k) { //доавляет нуклеотид в конец цепочки
    size_t number_nucl = this->NuclNum;
    reference r(this, number_nucl);
    r.operator=(k);
}

Nucleotid RNA::Nuclget(size_t index) {//по индексу возвращает нуклеотид
    reference r(this, index);
    return Nucleotid(r);
}

RNA RNA::operator+(RNA &r2) {

    size_t new_number_of_nucles = r2.NuclNum + this->NuclNum;
    size_t new_length = (r2.NuclNum + this->NuclNum + 4 * sizeof(size_t) - 1) / (4 * sizeof(size_t));
    RNA new_rna(new_number_of_nucles, new_length);
    for (size_t i = 0; i < this->length; i++) { //copy data first rna
        new_rna.rnaarr[i] = this->rnaarr[i];
    }
    for (size_t i = this->NuclNum; i < new_number_of_nucles; i++) { //copy data
        reference(&new_rna, i) = Nucleotid(reference(&r2, i - this->NuclNum));
    }

    return new_rna;
}

bool RNA::operator==(RNA &r2) {
    return this->NuclNum == r2.NuclNum && this->length == r2.length && *(this->rnaarr) == *(r2.rnaarr);
};


bool RNA:: operator!=(RNA &r2) {
    return !(this->NuclNum == r2.NuclNum && this->length == r2.length && *(this->rnaarr) == *(r2.rnaarr));

};

bool RNA:: isComplementary(RNA &r1) {
    size_t sum = 0;
    for (size_t i = 0; i < (r1.NuclNum); i++) {
        if (Nucleotid(reference(this, i)) == Nucleotid(reference(&r1, i))) {
            sum++;
            break;
        }
        sum += Nucleotid(reference(this, i)) & Nucleotid(reference(&r1, i));
    }
    return sum == 0;
}

RNA RNA:: trim(size_t ind) {
    if (ind == this->NuclNum + 1) {
        return *this;
    }
    size_t new_length = (ind + 4 * sizeof(size_t)) / (4 * sizeof(size_t));
    //auto *new_rnaarr = new size_t[new_length];
    RNA new_rna(ind + 1, new_length);
    for (size_t i = 0; i < ind + 1; i++) { //copy data
        reference(&new_rna, i) = Nucleotid(reference(this, i));
    }

    return new_rna;
}

RNA RNA:: split(size_t index) {

    assert (index >= 0);
    if (index == 0) {
        RNA rna_sp{this->NuclNum, this->length};
        rna_sp.rnaarr = this->rnaarr;
        this->rnaarr = nullptr;
        this->length = 0;
        this->NuclNum = 0;
        return rna_sp;
    }
    RNA rna_sp{};
    if (index >= this->NuclNum) {
        return rna_sp;
    }
    for (size_t i = index; i < (this->NuclNum); i++) { //copy data
        size_t b = i - index;

        rna_sp[b] = Nucleotid(this->operator[](i));
    }

    *this = this->trim(index);
    return rna_sp;
}

size_t RNA::  cardinality(Nucleotid value) {
    size_t returnNum = 0;
    for (size_t i = 0; i < this->NuclNum; i++) {
        if (Nucleotid(this->operator[](i)) == value) {
            returnNum++;
        }
    }
    return returnNum;
}

unordered_map<Nucleotid, size_t, hash<size_t> > RNA:: cardinality() {
    unordered_map<Nucleotid, size_t, hash<size_t> > returnNums;
    for (size_t i = 0; i < this->NuclNum; i++) {
        switch (Nucleotid(this->operator[](i))) {
            case A:
                returnNums[A]++;
                break;
            case C:
                returnNums[C]++;
                break;
            case G:
                returnNums[G]++;
                break;
            case T:
                returnNums[T]++;
                break;
            default:
                break;
        }
    }
    return returnNums;
}

size_t RNA::getlength() const {
    return this->length;
}
size_t RNA::capacity() const {
    return this->NuclNum;
}







// Created by Юлия Затолоцкая on 17.11.2019.
//

