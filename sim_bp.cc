#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "sim_bp.h"
#include<vector>
#include<math.h>

using namespace std;
/*  argc holds the number of command line arguments
    argv[] holds the commands themselves

    Example:-
    sim bimodal 6 gcc_trace.txt
    argc = 4
    argv[0] = "sim"
    argv[1] = "bimodal"
    argv[2] = "6"
    ... and so on
*/



int main (int argc, char* argv[]) {
    FILE *FP;               // File handler
    char *trace_file;       // Variable that holds trace file name;
    bp_params params;       // look at sim_bp.h header file for the the definition of struct bp_params
    char outcome;           // Variable holds branch outcome
    unsigned long int addr; // Variable holds the address read from input file
    int sizeofbh=0;
    uint32_t new_add = 0;
    uint32_t new_add1=0;
    int mask = 0;
    int prediction_hits = 0;
    int prediction_misses = 0;
    bool already_found = false;
    bool taken = false; //basically prediction for regular bimodal and gshare
    bool g_taken=false;
    bool b_taken=false;
    int n=0;
    int det=0;
    uint32_t new_addg=0;
    uint32_t new_addb=0;
    int h=0;
    bool gflag=0;
    bool bflag=0;
    uint32_t new_add_c=0;
    string type="";

    if (!(argc == 4 || argc == 5 || argc == 7)) {
        printf("Error: Wrong number of inputs:%d\n", argc - 1);
        exit(EXIT_FAILURE);
    }

    params.bp_name = argv[1];

    // strtoul() converts char* to unsigned long. It is included in <stdlib.h>
    if (strcmp(params.bp_name, "bimodal") == 0)              // Bimodal
    {
        if (argc != 4) {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc - 1);
            exit(EXIT_FAILURE);
        }
        params.M2 = strtoul(argv[2], NULL, 10);
        trace_file = argv[3];
        printf("COMMAND\n%s %s %lu %s\n", argv[0], params.bp_name, params.M2, trace_file);
    } else if (strcmp(params.bp_name, "gshare") == 0)          // Gshare
    {
        if (argc != 5) {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc - 1);
            exit(EXIT_FAILURE);
        }
        params.M1 = strtoul(argv[2], NULL, 10);
        params.N = strtoul(argv[3], NULL, 10);
        trace_file = argv[4];
        printf("COMMAND\n%s %s %lu %lu %s\n", argv[0], params.bp_name, params.M1, params.N, trace_file);

    } else if (strcmp(params.bp_name, "hybrid") == 0)          // Hybrid
    {
        if (argc != 7) {
            printf("Error: %s wrong number of inputs:%d\n", params.bp_name, argc - 1);
            exit(EXIT_FAILURE);
        }
        params.K = strtoul(argv[2], NULL, 10);
        params.M1 = strtoul(argv[3], NULL, 10);
        params.N = strtoul(argv[4], NULL, 10);
        params.M2 = strtoul(argv[5], NULL, 10);
        trace_file = argv[6];
        printf("COMMAND\n%s %s %lu %lu %lu %lu %s\n", argv[0], params.bp_name, params.K, params.M1, params.N, params.M2,
               trace_file);
        //create table/vector of size 2^N and initialize all counters to 2
    } else {
        printf("Error: Wrong branch predictor name:%s\n", params.bp_name);
        exit(EXIT_FAILURE);
    }

    type=params.bp_name;
    if(type=="hybrid"){
        h=params.K;
    }
    vector<int> prediction_table(pow(2, det), 2); //bimodal prediction table
    vector<int> prediction_tableg(pow(2, det), 2); //gshare prediction table
    vector<int> chooser_table(pow(2, h), 1);
    if(type=="bimodal"){
        det=params.M2;
        prediction_table.resize(pow(2, det), 2);
    }
    if(type=="gshare"){
        det=params.M1;
        prediction_table.resize(pow(2, det), 2);
    }
    if(type=="hybrid"){
        det=params.M2;
        prediction_table.resize(pow(2, det), 2);
        det=params.M1;
        prediction_tableg.resize(pow(2, det), 2);
    }


    //creates 2^m2 sized prediction table



    //create hybrid prediction tables

    //create global branch history
    uint32_t branch_history=0; //initialize to 0

    // Open trace_file in read mode
    FP = fopen(trace_file, "r");
    if (FP == NULL) {
        // Throw error and exit if fopen() failed
        printf("Error: Unable to open file %s\n", trace_file);
        exit(EXIT_FAILURE);
    }

    char str[2];
    while (fscanf(FP, "%lx %s", &addr, str) != EOF) {

        outcome = str[0];
        //if (outcome == 't')
        //printf("%lx %s\n", addr, "t");           // Print and test if file is read correctly
        //else if (outcome == 'n')
        //printf("%lx %s\n", addr, "n");          // Print and test if file is read correctly
        /*************************************
            Add branch predictor code here
        **************************************/
        n++; //num of predictions
        //------------------bimodal---------------------
        if (type=="bimodal") {
            new_add = (addr >> 2); //yeets the last 2 bits
            mask = (1 << params.M2) - 1;
            new_add = new_add & mask;
            //make a prediction
            if (prediction_table[new_add] >= 2) {
                taken = true;
            } else {
                taken = false;
            }

            //if actually taken
            if (outcome == 't') {
                prediction_table[new_add]++;
                if (prediction_table[new_add] > 3) {
                    prediction_table[new_add] = 3;
                }
                if (taken == true) {
                    prediction_hits++;
                } else {
                    prediction_misses++;
                }
            }
            //if not taken
            if (outcome == 'n') {
                prediction_table[new_add]--;
                if (prediction_table[new_add] < 0) {
                    prediction_table[new_add] = 0;
                }
                if (taken == false) {
                    prediction_hits++;
                } else {
                    prediction_misses++;
                }
            }
        }

            //--------------------gshare-------------------------

        else if(type=="gshare") {
            //cout<<"entered gshare"<<endl;
            //the current n-bit global branch history register is
            // XORed with the uppermost n bits of the m PC bits
            new_add = (addr >> 2); //yeets the last 2 bits
            mask = (1 << params.M1) - 1;
            new_add = new_add & mask; //m PC bits
            if(n<=100){
                //printf("Address: %x\n", addr);
            }
            if(n<=100){
                //printf("New add: %x\n", new_add);
            }
            //just want first n bits of m
            mask = ~((1 << ((params.M1-params.N))) - 1);
            new_add1=new_add & mask;

            //xor with branch history
            new_add1 = new_add1 ^ branch_history;
            //clear top 3 bits of new_add- clear m to m-n
            mask = (1 << (params.M1-params.N)) - 1;
            new_add = new_add & mask;
            //or the top n bits with last m-n of addy
            new_add=new_add | new_add1;
            //make a prediction
            if (prediction_table[new_add] >= 2) {
                taken = true;
            } else {
                taken = false;
            }

            //if actually taken
            if (outcome == 't') {
                prediction_table[new_add]++;
                if (prediction_table[new_add] > 3) {
                    prediction_table[new_add] = 3;
                }
                if (taken == true) {
                    prediction_hits++;
                } else {
                    prediction_misses++;
                }

            }

            //if not taken
            if (outcome == 'n') {
                prediction_table[new_add]--;
                if (prediction_table[new_add] < 0) {
                    prediction_table[new_add] = 0;
                }
                if (taken == false) {
                    prediction_hits++;
                } else {
                    prediction_misses++;
                }
            }

            //update branch predictor
            //shift
            branch_history=branch_history>>1;

            if(outcome=='t'){
                //clear last m-n bits
                mask = ~((1 << ((params.M1-params.N))) - 1);
                //set first bit
                branch_history |= (1<<(params.M1-1));
                branch_history=branch_history & mask;

            }
            else{
                //clear last m-n bits
                mask = ~((1 << ((params.M1-params.N))) - 1);
                //clear first bit
                branch_history &= ~(1<<(params.M1-1));
                branch_history=branch_history & mask;
                //printf("%x\n", branch_history & ((1 << (params.M1-params.N)) - 1));

            }
        }

        //------------hybrid predictor---------------------
        else if(type=="hybrid"){
            //Obtain two predictions, one from the gshare predictor
            // and one from the bimodal predictor
            //-------OBTAIN GSHARE PREDICTION--------
                new_add = (addr >> 2); //yeets the last 2 bits
                mask = (1 << params.M1) - 1;
                new_add = new_add & mask; //m PC bits
                if(n<=100){
                    //printf("Address: %x\n", addr);
                }
                if(n<=100){
                    //printf("New add: %x\n", new_add);
                }
                //just want first n bits of m
                mask = ~((1 << ((params.M1-params.N))) - 1);
                new_add1=new_add & mask;
                if(n<=100){
                    //printf("New add 1 (first n bits of m): %x\n", new_add1);
                }

                //xor with branch history
                new_add1 = new_add1 ^ branch_history;
                if(n<=100){
                    //printf("New add 1 after xoring with bh %x : %x\n", branch_history, new_add1);
                }
                //clear top 3 bits of new_add- clear m to m-n
                mask = (1 << (params.M1-params.N)) - 1;
                new_add = new_add & mask;
                //or the top n bits with last m-n of addy
                new_add=new_add | new_add1;

                if (n<=100) {
                    //new_add = new_add >> (params.M1 - params.N); //gets uppermost n bits of the m PC bits
                    //printf("For gshare- Addr: %x, prediction index: %d, branch history: %x\n", addr, new_add, branch_history);
                }

                //make a prediction
                if (prediction_tableg[new_add] >= 2) {
                    //using gshare, resize prediction vector
                    g_taken = true;
                } else {
                    //using bimodal, resize prediction vector
                    g_taken = false;
                }
                new_addg=new_add;

                //-------GET BIMODAL PREDICTION--------
                new_add = (addr >> 2); //yeets the last 2 bits
                mask = (1 << params.M2) - 1;
                new_add = new_add & mask;
            if (n<=100) {
                //new_add = new_add >> (params.M1 - params.N); //gets uppermost n bits of the m PC bits
                //printf("For bimodal- Addr: %x, prediction index: %d, branch history: %x\n", addr, new_add, branch_history);
            }
                //make a prediction
                if (prediction_table[new_add] >= 2) {
                    //resize prediction table
                    b_taken = true;
                } else {
                    b_taken = false;
                }
                new_addb=new_add;

                //------------------MAKE OVERALL PREDICTION AND UPDATE---------
                new_add_c=(addr>>2);
                mask = (1 << params.K) - 1;
                new_add_c=new_add_c & mask;

                if(chooser_table[new_add_c]>=2){
                    taken=g_taken;
                if (n<=100) {
                    //new_add = new_add >> (params.M1 - params.N); //gets uppermost n bits of the m PC bits
                    //printf("Choosing gshare. Addr: %x, prediction index: %d, branch history: %x\n", addr, new_addg, branch_history);
                }
                //update branch predictor based on outcome
                if (outcome == 't') {
                    prediction_tableg[new_addg]++;
                    if (prediction_tableg[new_addg] > 3) {
                        prediction_tableg[new_addg] = 3;
                    }
                    if (taken == true) { //prediction=true
                        //edit chooser table
                        prediction_hits++;
                    } else { //prediction=false
                        prediction_misses++;
                    }

                }

                //if not taken
                if (outcome == 'n') {
                    prediction_tableg[new_addg]--;
                    if (prediction_tableg[new_addg] < 0) {
                        prediction_tableg[new_addg] = 0;
                    }

                    if (taken == false) { //true pred
                        prediction_hits++;
                    } else {            //false pred
                        prediction_misses++;
                    }
                }
            }

            else{
                taken=b_taken;
                //update branch predictor based on outcome
                if (outcome == 't') {
                    prediction_table[new_addb]++;
                    if (prediction_table[new_addb] > 3) {
                        prediction_table[new_addb] = 3;
                    }

                    if (taken == true) { //prediction=true
                        //edit chooser table
                        prediction_hits++;
                    } else { //prediction=false
                        prediction_misses++;
                    }

                }

                //if not taken
                if (outcome == 'n') {
                    prediction_table[new_addb]--;
                    if (prediction_table[new_addb] < 0) {
                        prediction_table[new_addb] = 0;
                    }

                    if (taken == false) { //true pred
                        prediction_hits++;
                    } else {            //false pred
                        prediction_misses++;
                    }
                }
            }

            //update branch predictor
            //shift
            branch_history=branch_history>>1;
            if(outcome=='t'){
                //clear last m-n bits
                mask = ~((1 << ((params.M1-params.N))) - 1);
                //set first bit
                branch_history |= (1<<(params.M1-1));
                branch_history=branch_history & mask;

            }
            else{
                //clear last m-n bits
                mask = ~((1 << ((params.M1-params.N))) - 1);
                //clear first bit
                branch_history &= ~(1<<(params.M1-1));
                branch_history=branch_history & mask;
            }

            //update chooser counter
            if(outcome=='t'){
                if((b_taken==true) && (g_taken==false)){ //g is wrong, b is right
                    chooser_table[new_add_c]--;
                    if(chooser_table[new_add_c] < 0){
                        chooser_table[new_add_c]=0;
                    }
                }
                if((b_taken==false) && (g_taken==true)){ //g is right, b is wrong
                    chooser_table[new_add_c]++;
                    if(chooser_table[new_add_c] > 3){
                        chooser_table[new_add_c]=3;
                    }
                }
            }

            if(outcome=='n'){
                if((b_taken==false) && (g_taken==true)){ //g is wrong, b is right
                    chooser_table[new_add_c]--;
                    if(chooser_table[new_add_c] < 0){
                        chooser_table[new_add_c]=0;
                    }
                }
                if((b_taken==true) && (g_taken==false)){ //g is right, b is wrong
                    chooser_table[new_add_c]++;
                    if(chooser_table[new_add_c] > 3){
                        chooser_table[new_add_c]=3;
                    }
                }
            }
        }
    }


    //-------------outputs----------------
    cout<<"OUTPUT"<<endl;
    cout<<"number of predictions:    "<<n<<endl;
    cout<<"number of mispredictions:    "<<prediction_misses<<endl;
    printf("misprediction rate:       %.2f", (double(prediction_misses) / double(n))*100);
    cout<<"%"<<endl;

    if(type=="gshare") {
        cout << "FINAL GSHARE CONTENTS" << endl;
    }
    else if(type=="bimodal"){
        cout << "FINAL BIMODAL CONTENTS" << endl;
    }

    if(type=="bimodal" || type=="gshare") {
        for (int i = 0; i < prediction_table.size(); i++) {
            cout << i << " " << prediction_table[i] << endl;
        }
    }

    if (type=="hybrid"){
        cout << "FINAL CHOOSER CONTENTS" << endl;
        for (int i = 0; i < chooser_table.size(); i++) {
            cout << i << " " << chooser_table[i] << endl;
        }

        cout << "FINAL GSHARE CONTENTS" << endl;
        for (int i = 0; i < prediction_tableg.size(); i++) {
            cout << i << " " << prediction_tableg[i] << endl;
        }

        cout << "FINAL BIMODAL CONTENTS" << endl;
        for (int i = 0; i < prediction_table.size(); i++) {
            cout << i << " " << prediction_table[i] << endl;
        }
    }



    return 0;

}