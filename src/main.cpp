#include <iostream>
#include <initializer_list>
#include <args.hxx>

#include "RPinterop/interop.hpp"

int main(int argc, char **argv)
{
    args::ArgumentParser p("C++ interop converter");
    args::Group group(p, "Specification of input (pick one):", args::Group::Validators::Xor);
    args::ValueFlag<std::string> infile(group, "path", "The path to the input file", {"infile"});
    args::ValueFlag<std::string> infolder(group, "path", "The path to the input folder within which the .FLD files will be processed ", {"infolder"});
    
    args::Group outgroup(p, "Specification of output (pick one, or none to stdout):", args::Group::Validators::Xor);
    args::ValueFlag<std::string> outfile(outgroup, "path", "The path to the output file", {"outfile"});
    args::ValueFlag<std::string> outfolder(outgroup, "path", "The path to the output folder", {"outfolder"});
    args::ValueFlag<int> outstdout(outgroup, "", "Print to console on stdout", {"outstdout"});
    
    args::Group arguments(p, "Other arguments", args::Group::Validators::DontCare, args::Options::Global);
    args::HelpFlag h(arguments, "help", "help", {'h', "help"});
    args::ValueFlag<int> verbosity(arguments, "", "The verbosity of the output, 0 is no output", {"v"});
    

    try
    {
        p.ParseCLI(argc, argv);
       
        if (infolder){
            for (auto entry : std::filesystem::directory_iterator(infolder.Get())){
                if (entry.path().filename().extension() != ".FLD"){
                    continue;
                }
                std::string name = entry.path().filename().replace_extension("").string();
                
                try{
                    RPinterop::FLDfile FLD(entry.path());
                    std::string j = FLD.make_json(name).dump(1);
                    
                    auto warnings = FLD.get_warnings();
                    auto errors = FLD.get_errors();
                    
                    if (!warnings.empty()){
                        std::cout << name << std::endl;
                        for (auto &warn : warnings){
                            std::cout << "    [WARN]: " << warn << std::endl;
                        }
                    }
                    
                    if (!outfolder.Get().empty()){
                        if (!std::filesystem::exists(outfolder.Get())){
                            std::cerr << "Output folder of " << outfolder.Get() << " does not exist" << std::endl;
                        }
                        else{
                            std::ofstream ofs(outfolder.Get()+"/"+name+".json");
                            ofs << j << std::endl;
                        }
                    }
                    else{
                        std::cout << j << std::endl;
                    }
                }
                catch(std::exception &e){
                    //throw e;
                    std::cout << name << std::endl;
                    std::cout << "    [ERROR]: " << e.what() << std::endl;
                }
            }
        }
        else if (infile){
            RPinterop::FLDfile FLD(infile.Get());
            std::string name = std::filesystem::path(infile.Get()).filename().replace_extension("").string();
            if (!outfile.Get().empty()){
                std::ofstream ofs(outfile.Get());
                ofs << FLD.make_json(name).dump(1);
            }
            else{
                std::cout << FLD.make_json(name).dump(1) << std::endl;
            }
        }
        else{
            throw std::invalid_argument("Must provide either infile or infolder argument");
        }
    }
    catch (args::Help)
    {
        std::cout << p;
    }
    catch (args::Error& e)
    {
        std::cerr << "Error: " << e.what() << std::endl << p;
        return 1;
    }
    return 0;
}
