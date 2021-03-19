#ifndef ACOMS_PROBLEM
#define ACOMS_PROBLEM

#include "../standard.cpp"


const string DNA_ALPHABET = "ACGT";

const vector<double> DEFAULT_DNA_BACKGROUND = {0.25, 0.25, 0.25, 0.25};


struct Problem {
    int n;                              // number of sequences
    vector<string> sequence_names;
    vector<string> sequences;

    int w;                              // width of the motif

    int alp_size;                       // size of the alphabet
    string alphabet;
    vector<double> background;          // background probability of characters
    int total;                          // total number of characters

    Problem(const string &dataset_filename,
            const int &w,
            const string &alphabet = DNA_ALPHABET,
            const vector<double> &background = DEFAULT_DNA_BACKGROUND) :
            w(w),
            alphabet(alphabet),
            background(background) {
        alp_size = alphabet.length();

        read_dataset_from_file(dataset_filename);

        if (background.empty()) calculate_background();
    }

    void read_dataset_from_file(const string &dataset_filename) {
        n = 0;

        ifstream dataset_file(dataset_filename);
        string line;

        total = 0;

        while (dataset_file.good()) {
            getline(dataset_file, line);
            if (line.empty()) continue;
            if (line[0] == '>') {
                n++;
                sequence_names.emplace_back(line);
                sequences.emplace_back("");
            }
            else {
                refine_sequence(line);
                sequences.back().append(line);
                total += line.length();
            }
        }
    }

    void refine_sequence(string &s) {
        transform(s.begin(), s.end(), s.begin(), ::toupper);
        if (s.back() == '\r') s.pop_back();
        for (char &character : s)
            encode(character);
    }

    int encode(const char &character) {
        for (int i = 0; i < alp_size; ++i)
            if (alphabet[i] == character) return i;
        cerr << "Unidentified character " << character << endl;
        exit(1);
    }

    void calculate_background() {
        vector<int> occurrences(alp_size, 0);
        for (string &sequence : sequences)
            for (char &character : sequence)
                occurrences[encode(character)]++;

        for (int &occurrence : occurrences)
            background.emplace_back((double) occurrence / total);
    }

    vector<int> sequence_lengths() {
        vector<int> lengths;
        for (string &sequence : sequences) lengths.emplace_back(sequence.length());
        return lengths;
    }
};


#endif // ACOMS_PROBLEM