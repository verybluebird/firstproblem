//
// Created by Юлия Затолоцкая on 17.11.2019.
//

#ifndef RNA1_RNA_H
#define RNA1_RNA_H



#include <cstddef>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <chrono>

using namespace std;
enum Nucleotid {
    A, G, C, T
};
enum mask {
    O, B
};

class RNA {
    friend class reference;

    class reference {
        RNA* rna;
        size_t ind;
    public:
        reference(RNA* const rna, size_t ind);

        size_t change_block(size_t index_of_block, Nucleotid n, size_t index);

        static void out(size_t n);

        reference& operator=(Nucleotid nucl);

        static size_t get_mask(size_t index);

        explicit operator Nucleotid();
    };


private:
    size_t NuclNum;
    size_t* rnaarr;
    size_t length;
public:
    reference operator[](size_t ind);

    RNA();

    RNA(size_t NuclNum, size_t length);

    RNA(size_t num, Nucleotid n);

    ~RNA();

    friend std::ostream& operator<<(ostream& out, RNA& other);

    RNA(const RNA& rna);

    RNA& operator=(const RNA& r);

    static size_t getsize(const RNA& r);

    void add(Nucleotid k);

    Nucleotid Nuclget(size_t index);

    RNA operator+(RNA& r2);

    bool operator==(RNA& r2);

    bool operator!=(RNA& r2);

    bool isComplementary(RNA& r1);

    RNA trim(size_t ind);

    RNA split(size_t index);

    size_t cardinality(Nucleotid value);
    size_t getlength() const;
    size_t capacity() const;

    unordered_map<Nucleotid, size_t, hash<size_t> > cardinality();


};






//
// Created by Юлия Затолоцкая on 16.11.2019.
//

#ifndef UNTITLED1_RNA_H
#define UNTITLED1_RNA_H

#endif //UNTITLED1_RNA_H
//
// Created by Юлия Затолоцкая on 17.11.2019.
//

#ifndef RNA1_RNA1_H
#define RNA1_RNA1_H

#endif //RNA1_RNA1_H




#endif //RNA1_RNA_H
