//
// Created by Martin Steinegger on 6/7/21.
//
#include "GemmiWrapper.h"
#include "mmread.hpp"
#ifdef HAVE_ZLIB
#include "gz.hpp"
#endif
#include "input.hpp"
#include "foldcomp.h"
#include "cif.hpp"

#include <algorithm>

GemmiWrapper::GemmiWrapper(){
    threeAA2oneAA = {{"ALA",'A'},  {"ARG",'R'},  {"ASN",'N'}, {"ASP",'D'},
                     {"CYS",'C'},  {"GLN",'Q'},  {"GLU",'E'}, {"GLY",'G'},
                     {"HIS",'H'},  {"ILE",'I'},  {"LEU",'L'}, {"LYS",'K'},
                     {"MET",'M'},  {"PHE",'F'},  {"PRO",'P'}, {"SER",'S'},
                     {"THR",'T'},  {"TRP",'W'},  {"TYR",'Y'}, {"VAL",'V'},
                     // additional standard res
                     {"SEC",'C'}, // U
                     // modified res
                     {"00C", 'C'}, {"02K", 'A'}, {"02Y", 'A'}, {"03Y", 'C'},
                     {"05N", 'P'}, {"05O", 'Y'}, {"07O", 'C'}, {"0A1", 'Y'},
                     {"0A8", 'C'}, {"0A9", 'F'}, {"0AF", 'W'}, {"0AH", 'S'},
                     {"0AK", 'D'}, {"0AR", 'R'}, {"0BN", 'F'}, {"0CS", 'A'},
                     {"0E5", 'T'}, {"0EA", 'Y'}, {"0FL", 'A'}, {"0LF", 'P'},
                     {"0QL", 'C'}, {"0TD", 'D'}, {"0UO", 'W'}, {"0WZ", 'Y'},
                     {"0Y8", 'P'}, {"11Q", 'P'}, {"143", 'C'}, {"1AC", 'A'},
                     {"1IP", 'N'}, {"1OP", 'Y'}, {"1PA", 'F'}, {"1TQ", 'W'},
                     {"1TY", 'Y'}, {"1X6", 'S'}, {"200", 'F'}, {"23F", 'F'},
                     {"23P", 'A'}, {"2AG", 'A'}, {"2CO", 'C'}, {"2GX", 'F'},
                     {"2HF", 'H'}, {"2JG", 'S'}, {"2KK", 'K'}, {"2KP", 'K'},
                     {"2LT", 'Y'}, {"2ML", 'L'}, {"2MR", 'R'}, {"2MT", 'P'},
                     {"2QZ", 'T'}, {"2R3", 'Y'}, {"2RA", 'A'}, {"2RX", 'S'},
                     {"2SO", 'H'}, {"2TL", 'T'}, {"2TY", 'Y'}, {"2ZC", 'S'},
                     {"30F", 'C'}, // U 
                     {"30V", 'C'}, {"31Q", 'C'}, {"33X", 'A'},
                     {"34E", 'V'}, {"3AH", 'H'}, {"3BY", 'P'}, {"3CF", 'F'},
                     {"3CT", 'Y'}, {"3GL", 'E'}, {"3MY", 'Y'}, {"3PX", 'P'},
                     {"3QN", 'K'}, {"3WX", 'P'}, {"3X9", 'C'}, {"3YM", 'Y'},
                     {"3ZH", 'H'}, {"41H", 'F'}, {"41Q", 'N'}, {"42Y", 'S'},
                     {"432", 'S'}, {"45F", 'P'}, {"4AF", 'F'}, {"4AK", 'K'},
                     {"4AR", 'R'}, {"4AW", 'W'}, {"4BF", 'Y'}, {"4CF", 'F'},
                     {"4CY", 'M'}, {"4D4", 'R'}, {"4DP", 'W'}, {"4FB", 'P'},
                     {"4FW", 'W'}, {"4GJ", 'C'}, {"4HH", 'S'}, {"4HL", 'Y'},
                     {"4HT", 'W'}, {"4II", 'F'}, {"4IN", 'W'}, {"4J4", 'C'},
                     {"4J5", 'R'}, {"4KY", 'P'}, {"4L0", 'P'}, {"4LZ", 'Y'},
                     {"4MM", 'M'}, {"4N7", 'P'}, {"4N8", 'P'}, {"4N9", 'P'},
                     {"4OG", 'W'}, {"4OU", 'F'}, {"4OV", 'S'}, {"4PH", 'F'},
                     {"4PQ", 'W'}, {"4WQ", 'A'}, {"51T", 'Y'}, {"54C", 'W'},
                     {"55I", 'F'}, {"56A", 'H'}, {"5CR", 'F'}, {"5CS", 'C'},
                     {"5CT", 'K'}, {"5CW", 'W'}, {"5GG", 'K'}, {"5GM", 'I'},
                     {"5JP", 'S'}, {"5MW", 'K'}, {"5OH", 'A'}, {"5OW", 'K'},
                     {"5PG", 'G'}, {"5R5", 'S'}, {"5T3", 'K'}, {"5VV", 'N'},
                     {"5XU", 'A'}, {"60F", 'C'}, {"66D", 'I'}, {"6BR", 'T'},
                     {"6CL", 'K'}, {"6CV", 'A'}, {"6CW", 'W'}, {"6DN", 'K'},
                     {"6G4", 'K'}, {"6M6", 'C'}, {"6V1", 'C'}, {"6WK", 'C'},
                     {"73C", 'S'}, {"73N", 'R'}, {"73O", 'Y'}, {"73P", 'K'},
                     {"74P", 'K'}, {"7ID", 'D'}, {"7N8", 'F'}, {"7O5", 'A'},
                     {"7OZ", 'A'}, {"7QK", 'K'}, {"7T2", 'F'}, {"7TK", 'D'},
                     {"7VN", 'A'}, {"7VU", 'F'}, {"7W2", 'F'}, {"7XC", 'F'},
                     {"823", 'N'}, {"85L", 'C'}, {"86N", 'E'}, {"8JB", 'C'},
                     {"8LJ", 'P'}, {"8RE", 'K'}, {"9E7", 'K'}, {"9IJ", 'F'},
                     {"9KP", 'K'}, {"9RI", 'K'}, {"9WV", 'A'}, {"A1ADO", 'L'},
                     {"A1ALE", 'A'}, {"A1AP1", 'W'}, {"A1AVU", 'V'}, {"A1AZ2", 'H'},
                     {"A1D5B", 'K'}, {"A1D5E", 'P'}, {"A1D5P", 'Y'}, {"A1D64", 'C'},
                     {"A1H2H", 'F'}, {"A1H2I", 'W'}, {"A1H45", 'W'}, {"A1H5W", 'A'},
                     {"A1L4B", 'Q'}, {"A1LTQ", 'R'}, {"A1LWV", 'K'}, {"A30", 'Y'},
                     {"A3U", 'F'}, {"A5N", 'N'}, {"A8E", 'V'}, {"AA4", 'A'},
                     {"AAR", 'R'}, {"ABA", 'A'}, {"ACB", 'D'}, {"AEA", 'C'},
                     {"AEI", 'D'}, {"AGM", 'R'}, {"AGQ", 'Y'}, {"AGT", 'C'},
                     {"AHB", 'N'}, {"AHO", 'A'}, {"AHP", 'A'}, {"AIB", 'A'},
                     {"AKZ", 'D'}, {"ALC", 'A'}, {"ALN", 'A'}, {"ALO", 'T'},
                     {"ALS", 'A'}, {"ALT", 'A'}, {"ALV", 'A'}, {"ALY", 'K'},
                     {"AME", 'M'}, {"API", 'K'}, {"APK", 'K'}, {"AR7", 'R'},
                     {"ARM", 'R'}, {"ARO", 'R'}, {"ASA", 'D'}, {"ASB", 'D'},
                     {"ASL", 'D'}, {"AYA", 'A'}, {"AZH", 'A'}, {"AZK", 'K'},
                     {"B27", 'T'}, {"B2A", 'A'}, {"B2F", 'F'}, {"B2I", 'I'},
                     {"B2V", 'V'}, {"B3A", 'A'}, {"B3D", 'D'}, {"B3E", 'E'},
                     {"B3K", 'K'}, {"B3S", 'S'}, {"B3X", 'N'}, {"B3Y", 'Y'},
                     {"BB6", 'C'}, {"BB7", 'C'}, {"BB8", 'F'}, {"BB9", 'C'},
                     {"BBC", 'C'}, {"BCS", 'C'}, {"BCX", 'C'}, {"BFD", 'D'},
                     {"BG1", 'S'}, {"BH2", 'D'}, {"BHD", 'D'}, {"BIF", 'F'},
                     {"BLE", 'L'}, {"BMT", 'T'}, {"BOR", 'R'}, {"BP5", 'A'},
                     {"BPE", 'C'}, {"BSE", 'S'}, {"BTK", 'K'}, {"BTR", 'W'},
                     {"BUC", 'C'}, {"BWV", 'R'}, {"BXT", 'S'}, {"BYR", 'Y'},
                     {"C1J", 'R'}, {"C1S", 'C'}, {"C1X", 'K'}, {"C3Y", 'C'},
                     {"C4G", 'R'}, {"C4R", 'C'}, {"C5C", 'C'}, {"C67", 'R'},
                     {"C6C", 'C'}, {"C6D", 'R'}, {"CAF", 'C'}, {"CAS", 'C'},
                     {"CCS", 'C'}, {"CE7", 'N'}, {"CG6", 'C'}, {"CGA", 'E'},
                     {"CGU", 'E'}, {"CGV", 'C'}, {"CHP", 'G'}, {"CIR", 'R'},
                     {"CME", 'C'}, {"CMH", 'C'}, {"CML", 'C'}, {"CMT", 'C'},
                     {"CR5", 'G'}, {"CS1", 'C'}, {"CS3", 'C'}, {"CS4", 'C'},
                     {"CSA", 'C'}, {"CSB", 'C'}, {"CSD", 'C'}, {"CSJ", 'C'},
                     {"CSO", 'C'}, {"CSP", 'C'}, {"CSR", 'C'}, {"CSS", 'C'},
                     {"CSU", 'C'}, {"CSX", 'C'}, {"CSZ", 'C'}, {"CTH", 'T'},
                     {"CWR", 'S'}, {"CXM", 'M'}, {"CY0", 'C'}, {"CY1", 'C'},
                     {"CY3", 'C'}, {"CY4", 'C'}, {"CYD", 'C'}, {"CYF", 'C'},
                     {"CYG", 'C'}, {"CYJ", 'K'}, {"CYQ", 'C'}, {"CYR", 'C'},
                     {"CYW", 'C'}, {"CZ2", 'C'}, {"CZS", 'A'}, {"CZZ", 'C'},
                     {"D11", 'T'}, {"D2T", 'D'}, {"D3P", 'G'}, {"D8R", 'K'},
                     {"DA2", 'R'}, {"DAB", 'A'}, {"DAH", 'F'},
                     {"DAL", 'A'}, {"DAR", 'R'}, {"DAS", 'D'}, {"DBU", 'T'},
                     {"DBY", 'Y'}, {"DBZ", 'A'}, {"DCY", 'C'},
                     {"DDE", 'H'}, {"DDZ", 'A'}, {"DGL", 'E'}, {"DGN", 'Q'},
                     {"DHA", 'S'}, {"DHI", 'H'}, {"DHV", 'V'}, {"DI7", 'Y'},
                     {"DIL", 'I'}, {"DIV", 'V'}, {"DJD", 'F'}, {"DLE", 'L'},
                     {"DLS", 'K'}, {"DLY", 'K'}, {"DM0", 'K'}, {"DMH", 'N'},
                     {"DMK", 'D'}, {"DNE", 'L'}, {"DNP", 'A'}, {"DNS", 'K'},
                     {"DNW", 'A'}, {"DPL", 'P'}, {"DPN", 'F'}, {"DPP", 'A'},
                     {"DPQ", 'Y'}, {"DPR", 'P'}, {"DSE", 'S'}, {"DSG", 'N'},
                     {"DSN", 'S'}, {"DTH", 'T'}, {"DTR", 'W'},
                     {"DTY", 'Y'}, {"DV9", 'E'}, {"DVA", 'V'}, {"DYA", 'D'},
                     {"DYJ", 'P'}, {"DYS", 'C'}, {"E9C", 'Y'}, {"E9M", 'W'},
                     {"E9V", 'H'}, {"ECC", 'Q'}, {"ECX", 'C'}, {"EFC", 'C'},
                     {"EI4", 'R'}, {"EJA", 'C'}, {"ELY", 'K'}, {"EME", 'E'},
                     {"ESB", 'Y'}, {"ESC", 'M'}, {"EU0", 'V'}, {"EUP", 'T'},
                     {"EW6", 'S'}, {"EXA", 'K'}, {"EXL", 'W'}, {"EZY", 'G'},
                     {"F2F", 'F'}, {"F2Y", 'Y'}, {"F7W", 'W'}, {"FAK", 'K'},
                     {"FC0", 'F'}, {"FCL", 'F'}, {"FDL", 'K'}, {"FF9", 'K'},
                     {"FFL", 'L'}, {"FGA", 'E'}, {"FGL", 'G'}, {"FGP", 'S'},
                     {"FH7", 'K'}, {"FHL", 'K'}, {"FHO", 'K'}, {"FIO", 'R'},
                     {"FL6", 'D'}, {"FLT", 'Y'}, {"FME", 'M'}, {"FOE", 'C'},
                     {"FP9", 'P'}, {"FQA", 'K'}, {"FT6", 'W'}, {"FTR", 'W'},
                     {"FTY", 'Y'}, {"FVA", 'V'}, {"FY2", 'Y'}, {"FY3", 'Y'},
                     {"FZN", 'K'}, {"G1X", 'Y'}, {"G5G", 'L'}, {"GAU", 'E'},
                     {"GFT", 'S'}, {"GGL", 'E'}, {"GHG", 'Q'}, {"GHP", 'G'},
                     {"GL3", 'G'}, {"GLJ", 'E'}, {"GLK", 'E'}, {"GLZ", 'G'},
                     {"GMA", 'E'}, {"GME", 'E'}, {"GNC", 'Q'}, {"GPL", 'K'},
                     {"GQI", 'C'}, {"GVL", 'S'}, {"H14", 'F'}, {"H5M", 'P'},
                     {"H7V", 'A'}, {"HAR", 'R'}, {"HIA", 'H'}, {"HIC", 'H'},
                     {"HIP", 'H'}, {"HIQ", 'H'}, {"HIX", 'A'}, {"HL2", 'L'},
                     {"HLU", 'L'}, {"HLY", 'K'}, {"HMR", 'R'}, {"HNC", 'C'},
                     {"HOO", 'H'}, {"HOX", 'F'}, {"HP9", 'F'}, {"HPE", 'F'},
                     {"HPH", 'F'}, {"HQA", 'A'}, {"HR7", 'R'}, {"HRG", 'R'},
                     {"HS8", 'H'}, {"HS9", 'H'}, {"HSE", 'S'}, {"HSK", 'H'},
                     {"HSL", 'S'}, {"HSV", 'H'}, {"HT7", 'W'}, {"HTI", 'C'},
                     {"HTN", 'N'}, {"HTR", 'W'}, {"HVA", 'V'}, {"HY3", 'P'},
                     {"HYP", 'P'}, {"HZP", 'P'}, {"I2F", 'K'}, {"I2M", 'I'},
                     {"I3D", 'W'}, {"I4G", 'G'}, {"I4O", 'H'}, {"I6O", 'P'},
                     {"I7F", 'S'}, {"IAE", 'F'}, {"IAM", 'A'}, {"IAR", 'R'},
                     {"IAS", 'D'}, {"IB9", 'Y'}, {"IGL", 'G'}, {"IIL", 'I'},
                     {"ILM", 'I'}, {"ILX", 'I'}, {"ILY", 'K'}, {"IML", 'I'},
                     {"IOY", 'F'}, {"IPG", 'G'}, {"IYR", 'Y'}, {"IZO", 'M'},
                     {"J2F", 'Y'}, {"J8W", 'S'}, {"JGO", 'K'}, {"JJJ", 'C'},
                     {"JJK", 'C'}, {"JJL", 'C'}, {"JKH", 'P'}, {"JLP", 'K'},
                     {"K5H", 'C'}, {"K5L", 'S'}, {"KBE", 'K'}, {"KCR", 'K'},
                     {"KCX", 'K'}, {"KEO", 'K'}, {"KFP", 'K'}, {"KGC", 'K'},
                     {"KHB", 'K'}, {"KKD", 'D'}, {"KOR", 'M'}, {"KPF", 'K'},
                     {"KPI", 'K'}, {"KPY", 'K'}, {"KR3", 'K'}, {"KST", 'K'},
                     {"KYN", 'W'}, {"KYQ", 'K'}, {"L3O", 'L'}, {"L5P", 'K'},
                     {"LA2", 'K'}, {"LAL", 'A'}, {"LAY", 'L'}, {"LBY", 'K'},
                     {"LBZ", 'K'}, {"LCK", 'K'}, {"LDH", 'K'}, {"LE1", 'V'},
                     {"LED", 'L'}, {"LEF", 'L'}, {"LEI", 'V'}, {"LEN", 'L'},
                     {"LET", 'K'}, {"LGY", 'K'}, {"LLO", 'K'}, {"LLP", 'K'},
                     {"LLY", 'K'}, {"LME", 'E'}, {"LMQ", 'Q'}, {"LP6", 'K'},
                     {"LPD", 'P'}, {"LRK", 'K'}, {"LSO", 'K'}, {"LTU", 'W'},
                     {"LVN", 'V'}, {"LWI", 'F'}, {"LYF", 'K'}, {"LYK", 'K'},
                     {"LYN", 'K'}, {"LYO", 'K'}, {"LYP", 'K'}, {"LYR", 'K'},
                     {"LYX", 'K'}, {"LYZ", 'K'}, {"M0H", 'C'}, {"M2L", 'K'},
                     {"M2S", 'M'}, {"M3L", 'K'}, {"MAA", 'A'}, {"MBQ", 'Y'},
                     {"MCS", 'C'}, {"MDF", 'Y'}, {"ME0", 'M'}, {"MEA", 'F'},
                     {"MED", 'M'}, {"MEN", 'N'}, {"MEQ", 'Q'}, {"MGG", 'R'},
                     {"MGN", 'Q'}, {"MH6", 'S'}, {"MHL", 'L'}, {"MHO", 'M'},
                     {"MHS", 'H'}, {"MHU", 'F'}, {"MIR", 'S'}, {"MIS", 'S'},
                     {"MK8", 'L'}, {"ML3", 'K'}, {"MLE", 'L'}, {"MLL", 'L'},
                     {"MLY", 'K'}, {"MLZ", 'K'}, {"MME", 'M'}, {"MMO", 'R'},
                     {"MND", 'N'}, {"MNL", 'L'}, {"MP8", 'P'}, {"MPQ", 'G'},
                     {"MSA", 'G'}, {"MSE", 'M'}, {"MSO", 'M'}, {"MTY", 'Y'},
                     {"MVA", 'V'}, {"MYK", 'K'}, {"MYN", 'R'}, {"N0A", 'F'},
                     {"N10", 'S'}, {"N65", 'K'}, {"N7P", 'P'}, {"N80", 'P'},
                     {"N9P", 'A'}, {"NA8", 'A'}, {"NAL", 'A'}, {"NBQ", 'Y'},
                     {"NC1", 'S'}, {"NCB", 'A'}, {"NEP", 'H'}, {"NFA", 'F'},
                     {"NIY", 'Y'}, {"NLB", 'L'}, {"NLE", 'L'}, {"NLF", 'W'},
                     {"NLN", 'L'}, {"NLO", 'L'}, {"NLW", 'L'}, {"NLY", 'G'},
                     {"NMC", 'G'}, {"NMM", 'R'}, {"NPH", 'C'}, {"NVA", 'V'},
                     {"NYB", 'C'}, {"NYS", 'C'}, {"NZC", 'T'}, {"NZH", 'H'},
                     {"O2E", 'S'}, {"O6H", 'W'}, {"O7D", 'W'}, {"O7G", 'V'},
                     {"OAR", 'R'}, {"OAS", 'S'}, {"OBS", 'K'}, {"OCS", 'C'},
                     {"OCY", 'C'}, {"OGC", 'H'}, {"OHI", 'H'}, {"OHS", 'D'},
                     {"OLD", 'H'}, {"OLT", 'T'}, {"OMH", 'S'}, {"OMT", 'M'},
                     {"OMX", 'Y'}, {"OMY", 'Y'}, {"OPR", 'R'}, {"ORN", 'A'},
                     {"ORQ", 'R'}, {"OSE", 'S'}, {"OTH", 'T'}, {"OYL", 'H'},
                     {"OZ3", 'G'}, {"OZW", 'F'}, {"P1L", 'C'}, {"P2Q", 'Y'},
                     {"P3Q", 'Y'}, {"P5U", 'S'}, {"P9S", 'C'}, {"PAQ", 'Y'},
                     {"PAT", 'W'}, {"PBF", 'F'}, {"PCA", 'Q'}, {"PCS", 'F'},
                     {"PEC", 'C'}, {"PF5", 'F'}, {"PFF", 'F'}, {"PG9", 'G'},
                     {"PH6", 'P'}, {"PHA", 'F'}, {"PHD", 'D'}, {"PHI", 'F'},
                     {"PHL", 'F'}, {"PJ3", 'H'}, {"PLJ", 'P'}, {"PM3", 'F'},
                     {"POK", 'R'}, {"POM", 'P'}, {"PPN", 'F'}, {"PR3", 'C'},
                     {"PR4", 'P'}, {"PR7", 'P'}, {"PR9", 'P'}, {"PRJ", 'P'},
                     {"PRK", 'K'}, {"PRS", 'P'}, {"PRV", 'G'}, {"PSA", 'F'},
                     {"PSH", 'H'}, 
                     {"PSW", 'C'}, // U 
                     {"PTH", 'Y'}, {"PTM", 'Y'},
                     {"PTR", 'Y'}, {"PXU", 'P'}, {"PYA", 'A'}, {"PYX", 'C'},
                     {"Q2E", 'W'}, {"Q3P", 'K'}, {"Q75", 'A'}, {"Q8X", 'C'},
                     {"QCI", 'Q'}, {"QCS", 'C'}, {"QIL", 'I'}, {"QM8", 'L'},
                     {"QMB", 'A'}, {"QMM", 'Q'}, {"QNQ", 'C'}, {"QNT", 'C'},
                     {"QNW", 'C'}, {"QO2", 'C'}, {"QO5", 'C'}, {"QO8", 'C'},
                     {"QPA", 'C'}, {"QPH", 'F'}, {"QQ8", 'Q'}, {"QVA", 'C'},
                     {"QX7", 'A'}, {"R0K", 'E'}, {"R1A", 'C'}, {"R4K", 'W'},
                     {"RE0", 'W'}, {"RE3", 'W'}, {"RGL", 'R'}, {"RPI", 'R'},
                     {"RVJ", 'A'}, {"RVX", 'S'}, {"RX9", 'I'}, {"RXL", 'V'},
                     {"S1H", 'S'}, {"SAC", 'S'}, {"SAR", 'G'}, {"SBG", 'S'},
                     {"SBL", 'S'}, {"SCH", 'C'}, {"SCS", 'C'}, {"SCY", 'C'},
                     {"SDP", 'S'}, 
                     {"SE7", 'C'}, // U 
                     {"SEB", 'S'}, {"SEE", 'S'},
                     {"SEL", 'S'}, {"SEM", 'S'}, {"SEN", 'S'}, {"SEP", 'S'},
                     {"SET", 'S'}, {"SGB", 'S'}, {"SLL", 'K'}, {"SLZ", 'K'},
                     {"SMC", 'C'}, {"SME", 'M'}, {"SMF", 'F'}, {"SNC", 'C'},
                     {"SNK", 'H'}, {"SNM", 'S'}, {"SNN", 'N'}, {"SRZ", 'S'},
                     {"SUN", 'S'}, {"SVA", 'S'}, {"SVV", 'S'}, {"SVX", 'S'},
                     {"SVY", 'S'}, {"SVZ", 'S'}, {"SWW", 'S'}, {"SXE", 'S'},
                     {"SYS", 'C'}, // U
                     {"T0I", 'Y'}, {"T3R", 'L'}, {"T8L", 'T'},
                     {"T9E", 'T'}, {"TBG", 'V'}, {"TCQ", 'Y'}, {"TDD", 'L'},
                     {"TFW", 'W'}, {"TGH", 'W'}, {"TH5", 'T'}, {"TH6", 'T'},
                     {"THC", 'T'}, {"THZ", 'R'}, {"TIH", 'A'}, {"TIS", 'S'},
                     {"TLY", 'K'}, {"TMD", 'T'}, {"TNQ", 'W'}, {"TOQ", 'W'},
                     {"TOX", 'W'}, {"TPL", 'W'}, {"TPO", 'T'}, {"TPQ", 'Y'},
                     {"TQQ", 'W'}, {"TQZ", 'C'}, {"TRF", 'W'}, {"TRN", 'W'},
                     {"TRO", 'W'}, {"TRQ", 'W'}, {"TRW", 'W'}, {"TRX", 'W'},
                     {"TS9", 'I'}, {"TSQ", 'F'}, {"TSY", 'C'}, {"TTQ", 'W'},
                     {"TTS", 'Y'}, {"TY2", 'Y'}, {"TY5", 'Y'}, {"TY8", 'Y'},
                     {"TYB", 'Y'}, {"TYC", 'Y'}, {"TYE", 'Y'}, {"TYI", 'Y'},
                     {"TYJ", 'Y'}, {"TYN", 'Y'}, {"TYO", 'Y'}, {"TYQ", 'Y'},
                     {"TYS", 'Y'}, {"TYT", 'Y'}, {"TYW", 'Y'}, {"TYY", 'Y'},
                     {"U2X", 'Y'}, {"U3X", 'F'}, {"UF0", 'S'}, {"UMA", 'A'},
                     {"UX8", 'W'}, {"UXQ", 'F'}, {"V1C", 'F'}, {"V3C", 'P'},
                     {"V44", 'C'}, {"V5N", 'H'}, {"V61", 'F'}, {"V6W", 'W'},
                     {"V7T", 'K'}, {"VAD", 'V'}, {"VAH", 'V'}, {"VAI", 'V'},
                     {"VHF", 'E'}, {"VI3", 'C'}, {"VPV", 'K'}, {"VR0", 'R'},
                     {"WFP", 'F'}, {"WLU", 'L'}, {"WPA", 'F'}, {"WRP", 'W'},
                     {"WVL", 'V'}, {"WWB", 'W'}, {"WYK", 'R'}, {"X5P", 'A'},
                     {"X60", 'V'}, {"XA6", 'F'}, {"XCN", 'C'}, // {"XPL", 'O'},
                     {"XPR", 'P'}, {"XRE", 'P'}, {"XSN", 'N'}, {"XW1", 'A'},
                     {"XX1", 'K'}, {"XYC", 'A'}, {"Y1V", 'L'}, {"Y57", 'K'},
                     {"YCM", 'C'}, {"YFQ", 'A'}, {"YHA", 'K'}, {"YOF", 'Y'},
                     {"YPZ", 'Y'}, {"YRV", 'C'}, {"YTF", 'Q'}, {"YTH", 'T'},
                     {"Z3E", 'T'}, {"ZAI", 'K'}, {"ZAL", 'A'}, {"ZBZ", 'C'},
                     {"ZCL", 'F'}, {"ZDJ", 'Y'}, {"ZJU", 'D'}, {"ZKO", 'Q'},
                     {"ZLF", 'C'}, {"ZPO", 'P'}, {"ZQN", 'V'}, {"ZSX", 'G'},
                     {"ZT1", 'K'}, {"ZT6", 'Y'}, {"ZTK", 'A'}, {"ZU0", 'T'},
                     {"ZV4", 'F'}, {"ZYJ", 'P'}, {"ZYK", 'P'}, {"ZZJ", 'A'},
                     // unknown
                     {"UNK",'X'}};
    fixupBuffer = NULL;
}

std::unordered_map<std::string, int> getEntityTaxIDMapping(gemmi::cif::Document& doc) {
    std::unordered_map<std::string, int> entity_to_taxid;
    static const std::vector<std::pair<std::string, std::string>> loops_with_taxids = {
        { "_entity_src_nat.", "?pdbx_ncbi_taxonomy_id"},
        { "_entity_src_gen.", "?pdbx_gene_src_ncbi_taxonomy_id"},
        { "_pdbx_entity_src_syn.", "?ncbi_taxonomy_id"}
    };
    for (gemmi::cif::Block& block : doc.blocks) {
        for (auto&& [loop, taxid] : loops_with_taxids) {
            for (auto row : block.find(loop, {"entity_id", taxid})) {
                if (row.has2(1) == false) {
                    continue;
                }
                std::string entity_id = gemmi::cif::as_string(row[0]);
                if (entity_to_taxid.find(entity_id) != entity_to_taxid.end()) {
                    continue;
                }
                const char* endptr = NULL;
                int taxId = gemmi::no_sign_atoi(row[1].c_str(), &endptr);
                if (endptr != NULL && *endptr == '\0') {
                    entity_to_taxid.emplace(entity_id, taxId);
                }
            }
        }
    }
    return entity_to_taxid;
}

GemmiWrapper::Format mapFormat(gemmi::CoorFormat format) {
    switch (format) {
        case gemmi::CoorFormat::Pdb:
            return GemmiWrapper::Format::Pdb;
        case gemmi::CoorFormat::Mmcif:
            return GemmiWrapper::Format::Mmcif;
        case gemmi::CoorFormat::Mmjson:
            return GemmiWrapper::Format::Mmjson;
        case gemmi::CoorFormat::ChemComp:
            return GemmiWrapper::Format::ChemComp;
        default:
            return GemmiWrapper::Format::Unknown;
    }
}

bool GemmiWrapper::load(const std::string& filename, Format format) {
    if ((format == Format::Foldcomp) || (format == Format::Detect && gemmi::iends_with(filename, ".fcz"))) {
        std::ifstream in(filename, std::ios::binary);
        if (!in) {
            return false;
        }
        return loadFoldcompStructure(in, filename);
    }
    try {
#ifdef HAVE_ZLIB
        gemmi::MaybeGzipped infile(filename);
#else
        gemmi::BasicInput infile(filename);
#endif
        if (format == Format::Detect) {
            format = mapFormat(gemmi::coor_format_from_ext(infile.basepath()));
        }
        gemmi::Structure st;
        std::unordered_map<std::string, int> entity_to_tax_id;
        switch (format) {
            case Format::Mmcif: {
                gemmi::CharArray mem = read_into_buffer(infile);
                char* data = mem.data();
                size_t dataSize = mem.size();

                // hack to fix broken _citation.title in AF3
                const char target0[] = "Accurate structure prediction of biomolecular interactions with AlphaFold 3\n";
                size_t target0Len = sizeof(target0) - 1;
                const char target1[] = "_citation.title";
                size_t target1Len = sizeof(target1) - 1;
                char* it = std::search(data, data + dataSize, target0, target0 + target0Len);
                if (it != data + dataSize) {
                    while (it > data && *(it - 1) != '\n') {
                        it--;
                    }
                    if (strncmp(it, target1, target1Len) == 0) {
                        it[0] = '#';
                        it[1] = ' ';
                    }
                }

                gemmi::cif::Document doc = gemmi::cif::read_memory(mem.data(), mem.size(), infile.path().c_str());
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                st = gemmi::make_structure(doc);
                break;
            }
            case Format::Mmjson: {
                gemmi::cif::Document doc = gemmi::cif::read_mmjson(infile);
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                st = gemmi::make_structure(doc);
                break;
            }
            case Format::ChemComp: {
                gemmi::cif::Document doc = gemmi::cif::read(infile);
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                st = gemmi::make_structure_from_chemcomp_doc(doc);
                break;
            }
            default:
                st = gemmi::read_pdb(infile);
        }
        updateStructure((void*) &st, filename, entity_to_tax_id);
    } catch (...) {
        return false;
    }
    return true;
}

// https://stackoverflow.com/questions/1448467/initializing-a-c-stdistringstream-from-an-in-memory-buffer/1449527
struct OneShotReadBuf : public std::streambuf
{
    OneShotReadBuf(char* s, std::size_t n)
    {
        setg(s, s, s + n);
    }
};

bool GemmiWrapper::loadFromBuffer(const char * buffer, size_t bufferSize, const std::string& name, GemmiWrapper::Format format) {
    if ((format == Format::Foldcomp) || (format == Format::Detect && (bufferSize > MAGICNUMBER_LENGTH && strncmp(buffer, MAGICNUMBER, MAGICNUMBER_LENGTH) == 0))) {
        OneShotReadBuf buf((char *) buffer, bufferSize);
        std::istream istr(&buf);
        if (!istr) {
            return false;
        }
        return loadFoldcompStructure(istr, name);
    }
    try {
#ifdef HAVE_ZLIB
        gemmi::MaybeGzipped infile(name);
#else
        gemmi::BasicInput infile(name);
#endif
        if (format == Format::Detect) {
            format = mapFormat(gemmi::coor_format_from_ext(infile.basepath()));
        }

        gemmi::Structure st;
        std::unordered_map<std::string, int> entity_to_tax_id;
        switch (format) {
            case Format::Pdb:
                st = gemmi::pdb_impl::read_pdb_from_stream(gemmi::MemoryStream(buffer, bufferSize), name, gemmi::PdbReadOptions());
                break;
            case Format::Mmcif: {
                const char* targetBuffer = buffer;
                // hack to fix broken _citation.title in AF3
                const char target0[] = "Accurate structure prediction of biomolecular interactions with AlphaFold 3\n";
                size_t target0Len = sizeof(target0) - 1;
                const char target1[] = "_citation.title";
                size_t target1Len = sizeof(target1) - 1;
                const char* it = std::search(targetBuffer, targetBuffer + bufferSize, target0, target0 + target1Len);
                if (it != targetBuffer + bufferSize) {
                    if (fixupBuffer == NULL) {
                        fixupBufferSize = bufferSize;
                        fixupBuffer = (char*)malloc(fixupBufferSize);
                    } else if (bufferSize > fixupBufferSize) {
                        fixupBufferSize = bufferSize * 1.5;
                        fixupBuffer = (char*)realloc(fixupBuffer, fixupBufferSize);
                    }
                    memcpy(fixupBuffer, targetBuffer, bufferSize);
                    while (it > targetBuffer && *(it - 1) != '\n') {
                        it--;
                    }
                    if (strncmp(it, target1, target1Len) == 0) {
                        *(fixupBuffer + (it - targetBuffer)) = '#';
                        *(fixupBuffer + (it - targetBuffer) + 1) = ' ';
                    }
                    targetBuffer = fixupBuffer;
                }
                gemmi::cif::Document doc = gemmi::cif::read_memory(targetBuffer, bufferSize, name.c_str());
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                st = gemmi::make_structure(doc);
                break;
            }
            case Format::Mmjson: {
                char* bufferCopy = (char*)malloc(bufferSize + 1 * sizeof(char));
                if (bufferCopy == NULL) {
                    return false;
                }
                if (memcpy(bufferCopy, buffer, bufferSize) == NULL) {
                    free(bufferCopy);
                    return false;
                }
                bufferCopy[bufferSize] = '\0';
                gemmi::cif::Document doc = gemmi::cif::read_mmjson_insitu(bufferCopy, bufferSize, name);
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                st = gemmi::make_structure(doc);
                free(bufferCopy);
                break;
            }
            case Format::ChemComp: {
                gemmi::cif::Document doc = gemmi::cif::read_memory(buffer, bufferSize, name.c_str());
                entity_to_tax_id = getEntityTaxIDMapping(doc);
                st = gemmi::make_structure_from_chemcomp_doc(doc);
                break;
            }
            default:
                return false;
        }
        updateStructure((void*) &st, name, entity_to_tax_id);
    } catch (...) {
        return false;
    }
    return true;
}

bool GemmiWrapper::loadFoldcompStructure(std::istream& stream, const std::string& filename) {
    std::cout.setstate(std::ios_base::failbit);
    Foldcomp fc;
    int res = fc.read(stream);
    if (res != 0) {
        return false;
    }
    std::vector<AtomCoordinate> coordinates;
    fc.useAltAtomOrder = false;
    res = fc.decompress(coordinates);
    if (res != 0) {
        return false;
    }
    std::cout.clear();
    if (coordinates.size() == 0) {
        return false;
    }

    title.clear();
    chain.clear();
    names.clear();
    chainNames.clear();
    modelIndices.clear();
    ca.clear();
    ca_bfactor.clear();
    c.clear();
    cb.clear();
    n.clear();
    ami.clear();
    title.append(fc.strTitle);
    names.push_back(filename);
    const AtomCoordinate& first = coordinates[0];
    chainNames.push_back(first.chain);
    modelCount = 1;
    modelIndices.push_back(modelCount);
    int residueIndex = INT_MAX;
    Vec3 ca_atom = {NAN, NAN, NAN};
    Vec3 cb_atom = {NAN, NAN, NAN};
    Vec3 n_atom  = {NAN, NAN, NAN};
    Vec3 c_atom  = {NAN, NAN, NAN};
    float ca_atom_bfactor = 0.0;
    for (std::vector<AtomCoordinate>::const_iterator it = coordinates.begin(); it != coordinates.end(); ++it) {
        const AtomCoordinate& atom = *it;
        if (atom.residue_index != residueIndex) {
            if (residueIndex != INT_MAX) {
                ca.push_back(ca_atom);
                cb.push_back(cb_atom);
                n.push_back(n_atom);
                c.push_back(c_atom);
                ca_bfactor.push_back(ca_atom_bfactor);
                ca_atom = {NAN, NAN, NAN};
                cb_atom = {NAN, NAN, NAN};
                n_atom  = {NAN, NAN, NAN};
                c_atom  = {NAN, NAN, NAN};
                ca_atom_bfactor = 0.0;
            }
            if (threeAA2oneAA.find(atom.residue) == threeAA2oneAA.end()) {
                ami.push_back('X');
            } else {
                ami.push_back(threeAA2oneAA[atom.residue]);
            }
            residueIndex = atom.residue_index;
        }

        if (atom.atom == "CA") {
            ca_atom = { atom.coordinate.x, atom.coordinate.y, atom.coordinate.z };
            ca_atom_bfactor = atom.tempFactor;
        } else if (atom.atom == "CB") {
            cb_atom = { atom.coordinate.x, atom.coordinate.y, atom.coordinate.z };
        } else if (atom.atom == "N") {
            n_atom = { atom.coordinate.x, atom.coordinate.y, atom.coordinate.z };
        } else if (atom.atom == "C") {
            c_atom = { atom.coordinate.x, atom.coordinate.y, atom.coordinate.z };
        }
    }
    ca.push_back(ca_atom);
    cb.push_back(cb_atom);
    n.push_back(n_atom);
    c.push_back(c_atom);
    ca_bfactor.push_back(ca_atom_bfactor);
    chain.emplace_back(0, ca.size());
    return true;
}

void GemmiWrapper::updateStructure(void * void_st, const std::string& filename, std::unordered_map<std::string, int>& entity_to_tax_id) {
    gemmi::Structure * st = (gemmi::Structure *) void_st;

    title.clear();
    chain.clear();
    names.clear();
    chainNames.clear();
    modelIndices.clear();
    modelCount = 0;
    ca.clear();
    ca_bfactor.clear();
    c.clear();
    cb.clear();
    n.clear();
    ami.clear();
    taxIds.clear();
    title.append(st->get_info("_struct.title"));
    size_t currPos = 0;
    for (gemmi::Model& model : st->models){
        modelCount++;
        for (gemmi::Chain& ch : model.chains) {
            size_t chainStartPos = currPos;
            size_t pos = filename.find_last_of("\\/");
            std::string name = (std::string::npos == pos)
                               ? filename
                               : filename.substr(pos+1, filename.length());
            //name.push_back('_');
            chainNames.push_back(ch.name);
            char* rest;
            errno = 0;
            unsigned int modelNumber = strtoul(model.name.c_str(), &rest, 10);
            if ((rest != model.name.c_str() && *rest != '\0') || errno == ERANGE) {
                modelIndices.push_back(modelCount);
            }else{
                modelIndices.push_back(modelNumber);
            }

            names.push_back(name);
            int taxId = -1;
            for (gemmi::Residue &res : ch.first_conformer()) {
                if (taxId == -1) {
                    auto it = entity_to_tax_id.find(res.entity_id);
                    if (it != entity_to_tax_id.end()) {
                        taxId = it->second;
                    }
                }
                bool notPolymer = res.entity_type != gemmi::EntityType::Polymer;
                if (notPolymer == true) {
                        continue;
                }
                Vec3 ca_atom = {NAN, NAN, NAN};
                Vec3 cb_atom = {NAN, NAN, NAN};
                Vec3 n_atom  = {NAN, NAN, NAN};
                Vec3 c_atom  = {NAN, NAN, NAN};
                float ca_atom_bfactor;
                bool hasCA = false;
                for (gemmi::Atom &atom : res.atoms) {
                    if (atom.name == "CA") {
                        ca_atom.x = atom.pos.x;
                        ca_atom.y = atom.pos.y;
                        ca_atom.z = atom.pos.z;
                        ca_atom_bfactor = atom.b_iso;
                        hasCA = true;
                    } else if (atom.name == "CB") {
                        cb_atom.x = atom.pos.x;
                        cb_atom.y = atom.pos.y;
                        cb_atom.z = atom.pos.z;
                    } else if (atom.name == "N") {
                        n_atom.x = atom.pos.x;
                        n_atom.y = atom.pos.y;
                        n_atom.z = atom.pos.z;
                    } else if (atom.name == "C") {
                        c_atom.x = atom.pos.x;
                        c_atom.y = atom.pos.y;
                        c_atom.z = atom.pos.z;
                    }
                }
                if(hasCA == false){
                    continue;
                }
                ca_bfactor.push_back(ca_atom_bfactor);
                ca.push_back(ca_atom);
                cb.push_back(cb_atom);
                n.push_back(n_atom);
                c.push_back(c_atom);
                currPos++;
                if (threeAA2oneAA.find(res.name) == threeAA2oneAA.end()) {
                    ami.push_back('X');
                } else {
                    ami.push_back(threeAA2oneAA[res.name]);
                }
            }
            taxIds.push_back(taxId == -1 ? 0 : taxId);
            chain.push_back(std::make_pair(chainStartPos, currPos));
        }
    }
}
