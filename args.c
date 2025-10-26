#include <args.h>

#include <iostream>
#include <string>

double square_of_max(double a, double b) {
    double m = std::max(a, b);
    return m * m;
}

Args parse_args(int argc, char *argv[]) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "N") {
            args.N = std::stoi(argv[i+1]);
            ++i;
            continue;
        }
        if (std::string(argv[i]) == "M") {
            args.M = std::stoi(argv[i+1]);
            ++i;
            continue;
        }
        if (std::string(argv[i]) == "count_iter") {
            args.count_iter = std::stoi(argv[i+1]);
            ++i;
            continue;
        }
        if (std::string(argv[i]) == "delta") {
            args.delta = std::stof(argv[i+1]);
            ++i;
            continue;
        }
        if (std::string(argv[i]) == "--test") {
            args.with_test = true;
            continue;
        }
    }
    args.h1 = (args.B1 - args.A1) / (1.0 * args.M);
    args.h2 = (args.B2 - args.A2) / (1.0 * args.N);
    args.eps = square_of_max(args.h1, args.h2);
    std::cout << "args: N == " << args.N <<
                       " M == " << args.M << 
                       " count_iter == " << args.count_iter << 
                       " delta == " << args.delta << 
                       " test == " << args.with_test << std::endl;
    return args;
}