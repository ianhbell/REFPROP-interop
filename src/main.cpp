#include <iostream>
#include <initializer_list>
#include <args.hxx>

#include "RPinterop/interop.hpp"

int main(int argc, char **argv)
{
    args::ArgumentParser p("C++ interop converter");
    args::Group arguments(p, "arguments", args::Group::Validators::DontCare, args::Options::Global);
    args::ValueFlag<std::string> infile(arguments, "path", "The path to the input file", {"infile"});
    args::ValueFlag<std::string> outfile(arguments, "path", "The path to the output file", {"outfile"});
    args::HelpFlag h(arguments, "help", "help", {'h', "help"});

    try
    {
        p.ParseCLI(argc, argv);
        std::cout << infile.Get() << std::endl;
        std::cout << outfile.Get() << std::endl;
        RPinterop::FLDfile FLD(infile.Get());
        std::string name = "CO2";
        
        if (!outfile.Get().empty()){
            std::ofstream ofs(outfile.Get());
            ofs << FLD.make_json(name).dump(1);
        }
        else{
            std::cout << FLD.make_json(name).dump(1) << std::endl;
        }
    }
    catch (args::Help)
    {
        std::cout << p;
    }
    catch (args::Error& e)
    {
        std::cerr << e.what() << std::endl << p;
        return 1;
    }
    return 0;
}