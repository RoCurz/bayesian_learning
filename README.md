Usage is instructed below:
- Change permission of `compiler.sh` and `run.sh` if required
- Run `./compiler.sh`
- Run `./run.sh <filename>.bif <dataset>.bif`

To verify the score:
- Run `g++ -std=c++20 format_checker.cpp -o format_checker && ./format_checker`