The code is used to verify the complexity evaluation of our algebraic attacks on full LowMC with 2 chosen plaintexts.

The code can be compiled using the command g++ LowMC.cpp LowMC.h main.cpp -03 -std=c++11

Then, when running the code, you are required to input a command (1 or 2) to perform experiments on different constructions.

1 -> test the construction with a partial S-box layer (1000 tests in a few minutes)

2 -> test the construction with a full S-box layer (10000 tests in about 2 minute)

Our paper is available at https://eprint.iacr.org/2020/1034
