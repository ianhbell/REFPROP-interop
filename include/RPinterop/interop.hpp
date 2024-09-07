#pragma once

#include <filesystem>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <regex>
#include <optional>

#include "nlohmann/json.hpp"

namespace RPinterop{
namespace internal{

    // for string delimiter (https://stackoverflow.com/a/46931770)
    inline std::vector<std::string> strsplit (const std::string& s, const std::string& delimiter) {
        std::size_t pos_start = 0,
               pos_end = std::string::npos,
               delim_len = delimiter.length();
        std::string token;
        std::vector<std::string> result;
        
        while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
            token = s.substr (pos_start, pos_end - pos_start);
            pos_start = pos_end + delim_len;
            result.push_back (token);
        }
        
        result.push_back (s.substr (pos_start));
        return result;
    }
    inline std::string strjoin(const std::vector<std::string> &lines, const std::string& delim){
        if (lines.empty()){ return ""; }
        std::string result = lines[0];
        for (auto i = 1U; i < lines.size(); ++i){
            result += delim + lines[i];
        }
        return result;
    }

    /**
    For a line with "!", return all characters before (non-inclusive) the "!", otherwise return the whole string
     */
    inline std::string strip_line_comment(const std::string& line){
        auto ipos = line.find_first_of("!");
        if (ipos == std::string::npos){
            // No comment found, return whole string
            return line;
        }
        else{
            // Return the characters preceding the comment character of "!"
            return line.substr(0, ipos);
        }
    }
    /**
    For a line with a character, return all characters after the character, otherwise return the whole string
     */
    inline std::string strip_line_leading_comment(const std::string& line, const std::string& comment_char = "!"){
        auto ipos = line.find_first_of(comment_char);
        if (ipos == std::string::npos){
            // No comment found, return whole string
            return line;
        }
        else{
            // Return the characters after the comment character of "!"
            return line.substr(ipos+1);
        }
    }

    inline std::string strip_trailing_whitespace(const std::string& line){
        auto ipos = line.find_last_not_of(" \t\n\v\f\r");
        if (ipos == std::string::npos){
            // No whitespace found, return whole string
            return line;
        }
        else{
            return line.substr(0, ipos+1);
        }
    }
    inline std::string strip_leading_whitespace(const std::string& line){
        auto ipos = line.find_first_not_of(" \t\n\v\f\r");
        if (ipos == std::string::npos){
            // No whitespace found, return whole string
            return line;
        }
        else{
            return line.substr(ipos);
        }
    }

    /// Return a chunk of lines starting from a line that begins with starting_key and ends before the first empty line
    inline std::vector<std::string> get_line_chunk(
                                            const std::vector<std::string>& lines,
                                            const std::string& starting_key,
                                            const std::optional<std::size_t>& initial_index = std::nullopt
                                            )
    {
        // block starts with the key
        size_t istart = std::string::npos;
        for (auto i = initial_index.value_or(0); i < lines.size(); ++i){
            if (lines[i].find(starting_key) == 0){
                istart = i; break;
            }
        }
        if (istart == std::string::npos){
            throw std::invalid_argument("Block starting with "+starting_key+" was not found");
        }
        // block ends with an empty line
        size_t iend = std::string::npos;
        for (auto i = istart+1; i < lines.size(); ++i){
            if (lines[i].find_first_not_of(" \t\n\v\f\r") == std::string::npos){
                iend = i; break;
            }
        }
        if (iend == std::string::npos){
            throw std::invalid_argument("Could not find an empty line at end of block");
        }
        return std::vector<std::string>(lines.begin()+istart, lines.begin()+iend);
    }

    inline std::string to_upper(std::string s){
        std::locale locale;
        auto f = [&locale] (char ch) { return std::use_facet<std::ctype<char>>(locale).toupper(ch); };
        std::transform(s.begin(), s.end(), s.begin(), f);
        return s;
    }

    inline auto read_allnum_from_line(const std::string& line){
        auto vals = strsplit(strip_line_comment(line), " ");
        std::vector<double> o;
        for (auto v : vals){
            if (!v.empty()){
                // Replace the D or d in FORTRAN scientific notation with e
                // because C++ (and everything else) uses e for exponentiation in
                // scientific notation
                std::string v_e = std::regex_replace(v, std::regex("[Dd]"), "e");
                // Replace spaces with nothing
                std::string v_ew = std::regex_replace(v_e, std::regex("\\s+"), "");
                if (v_ew.empty()){
                    continue;
                }
                o.push_back(strtod(v_ew.c_str(), nullptr));
            }
        }
        return o;
    }

    inline std::string read_file(const std::filesystem::path &path, bool replace_tabs = true){
        
        std::ifstream ifs(path);
        if (!ifs){
            throw std::invalid_argument("Path to be loaded is not readable: " + path.string());
        }
        std::stringstream buffer;  buffer << ifs.rdbuf();
        
        std::string contents = buffer.str();
        
        if (contents.find("\t") != std::string::npos){
            if(replace_tabs){
                contents = std::regex_replace(contents, std::regex("\t"), "    ");
            }
            else{
                throw std::invalid_argument("Tabs found in the file "+path.string()+". Replace them with spaces.");
            }
        }
        return contents;
    }
}

/***
 A parsing class that takes in a set of lines and line-by-line reads in things
 There are methods that increment the line counter, using the next line
 */
class LineParser{
private:
    const std::vector<std::string>& lines_;
    std::size_t i;
public:
    LineParser(const std::vector<std::string>& lines) : lines_(lines), i(0){}
    
    /// Get a reference to the next line to be parsed
    const auto& get_next_line() const {
        if (i > lines_.size()-1){ throw std::invalid_argument("no next line to return"); }
        return lines_[i];
    }
    /// Set the index of the next line
    void set_i(std::size_t i_){ this->i = i_; }
    /// Get the index of the next line
    std::size_t get_i() const { return i; }
    
    auto read_allnum_from_line(const std::string& line) const{
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), " ");
        std::vector<double> o;
        for (auto v : vals){
            if (!v.empty()){
                // Replace the D or d in FORTRAN scientific notation with e
                // because C++ (and everything else) uses e for exponentiation in
                // scientific notation
                std::string v_e = std::regex_replace(v, std::regex("[Dd]"), "e");
                // Replace spaces with nothing
                std::string v_ew = std::regex_replace(v_e, std::regex("\\s+"), "");
                if (v_ew.empty()){
                    continue;
                }
                o.push_back(strtod(v_ew.c_str(), nullptr));
            }
        }
        return o;
    }
    /// Read in a fixed number of numbers (floats or integers) from the line
    auto read_Nnum_from_line(const std::string& line, std::size_t n) const{
        auto vals = read_allnum_from_line(line);
        if (n > 0){
            if (vals.size() != n){ throw std::invalid_argument("Unable to read "+std::to_string(n)+" numbers from this line:"+line); }
        }
        return vals;
    }
    /// Read in the strings from the line, assuming (possibly multiple) " " as the string delimiter
    auto read_strs_from_line(const std::string &line) const {
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), " ");
        return vals;
    }
    /// Read in one string from the line
    auto read_1str_from_line(const std::string &line) const {
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), " ");
        if (vals.empty()){ throw std::invalid_argument("Unable to read at least one string from this line:"+line); }
        return vals[0];
    }
    
    /// Read in a fixed number of numbers (floats or integers) from the next line and increment the counter
    auto read_Nnum_and_increment(std::size_t n){
        auto vals = read_Nnum_from_line(get_next_line(), n);
        i++;
        return vals;
    };
    /// Read in all the numbers (floats or integers) from the next line and increment the counter
    auto read_allnum_and_increment(){
        auto vals = read_allnum_from_line(get_next_line());
        i++;
        return vals;
    };
    /// Read in 1 space-delimited string from the next line and increment the counter
    auto read_1str_and_increment(){
        auto val = read_1str_from_line(get_next_line());
        i++;
        return val;
    }
    /// Read in one number from the next line and increment the counter
    auto read_1num_and_increment(){
        auto vals = read_Nnum_from_line(get_next_line(), 1);
        i++;
        return vals[0];
    };
    
};

/**
 * Things that are obtained from parsing a residual Helmholtz
 * energy EOS
 */
struct ResidualResult{
    double Tmin_K,Tmax_K,pmax_kPa,rhomax_molL;
    std::string cp0_pointer;
    double MM_kgkmol = -1,Ttriple_K,ptriple_kPa,rhotriple_molL,Tnbp_K,acentric,
    Tcrit_K,pcrit_kPa,rhocrit_molL,Tred_K,rhored_molL,R=-1;
    nlohmann::json alphar;
    std::string DOI_EOS = "";
};

/**
 * Method for converting modified-Benedict-Webb-Rubin EOS form to conventional 
 * Helmholtz energy form
 */ 
inline auto BWR2FEQ(const std::vector<std::string>& lines){
    ResidualResult res;
    LineParser parser(lines);
    
    auto model_key = parser.read_1str_from_line(lines[1]);
    
    // Find the DOI if present
    std::string DOI_EOS = "??";
    
    // Find the first non-header row;
    size_t i = std::string::npos;
    for (auto j = 2U; j < lines.size(); ++j){
        // First line not started by element in {:,?,!}, stop
        auto ind = lines[j].find_first_of(":?!");
        
        auto indEOS = lines[j].find(":DOI: ");
        if(!lines[j].empty() && (indEOS == 0)){
            std::string init = lines[j].substr(indEOS+6, lines[j].size()-6);
            DOI_EOS = internal::strip_trailing_whitespace(init);
        }
                                              
        if(!lines[j].empty() && (ind == 0)){
            continue;
        }
        i = j; break;
    }
    parser.set_i(i);
    
    // These lambdas are just for concision in the next block
    auto readallline = [&](const std::string& line){ return parser.read_allnum_from_line(line); };
    auto readn = [&](std::size_t n){ return parser.read_Nnum_and_increment(n); };
    auto read1 = [&](){ return parser.read_Nnum_and_increment(1)[0]; };
    auto read1str = [&](){ return parser.read_1str_and_increment(); };
    
    // Read in all the metadata parameters
    res.Tmin_K = read1();
    res.Tmax_K = read1();
    res.pmax_kPa = read1();
    res.rhomax_molL = read1();
    res.cp0_pointer = read1str();
    res.MM_kgkmol = read1();
    res.Ttriple_K = read1();
    res.ptriple_kPa = read1();
    res.rhotriple_molL = read1();
    res.Tnbp_K = read1();
    res.acentric = read1();
    auto c = readn(3);
    res.Tcrit_K = c[0];
    res.pcrit_kPa = c[1];
    res.rhocrit_molL = c[2];
    auto reds = readn(2);
    res.Tred_K = reds[0];
    res.rhored_molL = reds[1];
    double gamma_FLD = read1();
    double R_barL = read1();
    res.R = R_barL*100.0; // MBWR in FLD files uses L-bar for energy unit
    res.DOI_EOS = DOI_EOS;
    
//    res.R = 8.3144598; // TODO remove!!!!!
    
    double gamma_aspublished = -1/(gamma_FLD*gamma_FLD);
    // Leading coefficient in exp() in exponential term is then:
    double g = -gamma_aspublished*(res.rhocrit_molL*res.rhocrit_molL);
    
    double l = 2.0; // exponent on delta in exponential term

    // And now we read the EOS coefficients in
    auto termcounts = readn(2); // Should always be 32, 1
    // Collect all the coefficients
    std::vector<double> cc = {0.0}; // We are going to 1-index for simplicity in what follows
    for (auto j = parser.get_i(); j < lines.size(); ++j){
        for (auto& val : readallline(lines[j])){
            cc.push_back(val*100); // Pressures in the MBWR formulation in REFPROP are in bar, so convert to kPa
        }
    }
    
    // Collectors for generated terms
    std::vector<double> n_, t_, d_, l_, g_;
    auto addterm = [&](double n, double t, double d, double l, double g){
        n_.push_back(n);
        t_.push_back(t);
        d_.push_back(d);
        l_.push_back(l);
        g_.push_back(g);
//        std::cout << n << "    " << t << "    " << d << "    " << l << "    " << g << "    " << std::endl;
    };
    
    // Powers of gamma^k, to be used in the denominators of exponential terms
    double g1 = g, g2 = g1*g, g3 = g2*g, g4 = g3*g, g5 = g4*g, g6 = g5*g;
    
    std::vector<double> t2 = {0,1,0.5,0,-1,-2,1,0,-1,-2,1,0,-1,0,-1,-2,-1,-1,-2,-2}; // exponents on T in normal part
    std::vector<double> d2 = {0,2,2,2,2,2,3,3,3,3,4,4,4,5,6,6,7,8,8,9}; // exponents on density in normal part
    std::vector<double> s = {-2,-3,-2,-4,-2,-3,-2,-4,-2,-3,-2,-3,-4}; // exponents on T in exponential part
    std::vector<double> r = {3,3,5,5,7,7,9,9,11,11,13,13,13}; // exponents on density in exponential part

    // For exponential part of Z-1...
    // See Table 3.5 from Span book. But first convert leading coefficients to the form of Eq 3.26,
    // but here *without* the mysterious term in the denominator
    std::vector<double> n(20, 0.0);
    for (auto k = 20; k < 33; ++k)
        n.push_back(cc[k]*pow(res.rhored_molL, r[k-20]-1)*pow(res.Tred_K, s[k-20]-1)/res.R);
    
    addterm(n[20]/(2*g1) + n[22]/(2*g2) + n[24]/g3 + 3*n[26]/g4 + 12*n[28]/g5 + 60*n[30]/g6, 3, 0, 0, 0);
    addterm(n[21]/(2*g1) + n[25]/g3 + 12*n[29]/g5 + 60*n[31]/g6, 4, 0, 0, 0);
    addterm(n[23]/(2*g2) + 3*n[27]/g4 + 60*n[32]/g6, 5, 0, 0, 0);
    
    // For polynomial term, see page 26 of Span monograph, straightforward conversion to alphar contribution
    // Converting from Eq. 3.28 to 3.26 in the Span formulation, but with an additional factor of d2[i]-1 in the denominator
    // to match Younglove and McLinden, Eq. B6. This is also what is done in REFPROP
    for (int k = 1; k < 20; ++k){
        addterm(cc[k]*pow(res.rhored_molL, d2[k]-1)*pow(res.Tred_K, t2[k]-1)/res.R/(d2[k]-1), 1-t2[k], d2[k]-1, 0, 0);
    }
    
    addterm(-(n[20]/(2*g1) + n[22]/(2*g2) + n[24]/g3 + 3*n[26]/g4 + 12*n[28]/g5 + 60*n[30]/g6), 3, 0, l, g);
    addterm(-(n[21]/(2*g1) + n[25]/g3 + 12*n[29]/g5 + 60*n[31]/g6), 4, 0, l, g);
    addterm(-(n[23]/(2*g2) + 3*n[27]/g4 + 60*n[32]/g6), 5, 0, l, g);

    addterm(-n[22]/(2*g1) - n[24]/g2 -3*n[26]/g3 - 12*n[28]/g4 - 60*n[30]/g5, 3, 2, l, g);
    addterm(-n[25]/g2 - 12*n[29]/g4 -60*n[31]/g5, 4, 2, l, g);
    addterm(-n[23]/2/g1 - 3*n[27]/g3 -60*n[32]/g5, 5, 2, l, g);

    addterm(-n[24]/2/g1 - 3*n[26]/2/g2 - 6*n[28]/g3-30*n[30]/g4, 3, 4, l, g);
    addterm(-n[25]/2/g1 - 6*n[29]/g3 - 30*n[31]/g4, 4, 4, l, g);
    addterm(-3*n[27]/2/g2 - 30*n[32]/g4, 5, 4, l, g);

    addterm(-n[26]/2/g1 - 2*n[28]/g2 - 10*n[30]/g3, 3, 6, l, g);
    addterm(-2*n[29]/g2 - 10*n[31]/g3, 4, 6, l, g);
    addterm(-n[27]/2/g1 - 10*n[32]/g3, 5, 6, l, g);

    addterm(-n[28]/2/g1 - 5*n[30]/2/g2, 3, 8, l, g);
    addterm(-n[29]/2/g1 - 5*n[31]/2/g2, 4, 8, l, g);
    addterm(-5*n[32]/2/g2, 5, 8, l, g);

    addterm(-n[30]/2/g1, 3, 10, l, g);
    addterm(-n[31]/2/g1, 4, 10, l, g);
    addterm(-n[32]/2/g1, 5, 10, l, g);

    res.alphar = nlohmann::json::array();
    nlohmann::json term = {{"type", "ResidualHelmholtzExponential"}, {"n", n_}, {"t", t_}, {"d", d_}, {"l", l_}, {"g", g_}};
    res.alphar.push_back(term);
    return res;
}

inline std::optional<ResidualResult> parse_ECS(const std::vector<std::string>& lines){
    ResidualResult res;
    LineParser parser(lines);
    
    auto model_key = parser.read_1str_from_line(lines[1]);
    
    // Find the DOI if present
    std::string DOI_EOS = "??";
    
    // Find the first non-header row;
    size_t i = std::string::npos;
    for (auto j = 2U; j < lines.size(); ++j){
        // First line not started by element in {:,?,!}, stop
        auto ind = lines[j].find_first_of(":?!");
        
        auto indEOS = lines[j].find(":DOI: ");
        if(!lines[j].empty() && (indEOS == 0)){
            std::string init = lines[j].substr(indEOS+6, lines[j].size()-6);
            DOI_EOS = internal::strip_trailing_whitespace(init);
        }
                                              
        if(!lines[j].empty() && (ind == 0)){
            continue;
        }
        i = j; break;
    }
    parser.set_i(i);
    
    res.Tmin_K = parser.read_1num_and_increment();
    res.Tmax_K = parser.read_1num_and_increment();
    res.pmax_kPa = parser.read_1num_and_increment();
    res.rhomax_molL = parser.read_1num_and_increment();
    res.cp0_pointer = parser.read_1str_and_increment();
    std::string reference_fluid = internal::strip_trailing_whitespace(parser.read_1str_and_increment());
    std::string ref_model_pointer = parser.read_1str_and_increment();
    double acentric_reference = parser.read_1num_and_increment();
    double Z_crit_reference = parser.read_1num_and_increment();
    double acentric_fluid = parser.read_1num_and_increment();
    double T_crit_fluid_K = parser.read_1num_and_increment();
    double p_crit_fluid_kPa = parser.read_1num_and_increment();
    double rho_crit_fluid_molL = parser.read_1num_and_increment();
    
    auto read_coeffs_exps = [&](int N){
        std::vector<double> coeffs, exps;
        for (auto i = 0; i < N; ++i){
            auto coefexp = parser.read_Nnum_and_increment(2);
            coeffs.push_back(coefexp[0]);
            exps.push_back(coefexp[1]);
        }
        return std::make_tuple(coeffs, exps);
    };
    int NTF = static_cast<int>(parser.read_1num_and_increment());
    auto [fT_coeffs, fT_exps] = read_coeffs_exps(NTF);
    int NDF = static_cast<int>(parser.read_1num_and_increment());
    auto [fD_coeffs, fD_exps] = read_coeffs_exps(NDF);
    int NTH = static_cast<int>(parser.read_1num_and_increment());
    auto [hT_coeffs, hT_exps] = read_coeffs_exps(NTH);
    int NDH = static_cast<int>(parser.read_1num_and_increment());
    auto [hD_coeffs, hD_exps] = read_coeffs_exps(NDH);
    
    if (fT_coeffs.size() == 2 && fD_coeffs.empty() && hT_coeffs.size() == 2 && hD_coeffs.empty()){
        nlohmann::json jref = {
            {"name", reference_fluid},
            {"acentric", acentric_reference},
            {"Z_crit", Z_crit_reference}
        };
        nlohmann::json jfluid = {
            {"acentric", acentric_fluid},
            {"f_T_coeffs", fT_coeffs},
            {"h_T_coeffs", hT_coeffs},
            {"rhomolar_crit / mol/m^3", rho_crit_fluid_molL*1000},
            {"T_crit / K", T_crit_fluid_K},
            {"Z_crit", p_crit_fluid_kPa*1000/(8.31446261815324*T_crit_fluid_K*rho_crit_fluid_molL*1000)}
        };
        
        nlohmann::json j = {
            {"kind", "multifluid-ECS-HuberEly1994"},
            {"model", {
                {"reference_fluid", jref},
                {"fluid", jfluid}
            }}
        };
        // std::cout << j.dump(1) << std::endl;
        res.alphar = j;
    }
    else{
        throw std::invalid_argument("I don't know yet how to parse this kind of ECS term");
    }
    return res;
}

/**
Given a string that contains an EOS block, return its
JSON representation
*/
inline std::optional<ResidualResult> convert_FEQ(const std::vector<std::string>& lines){
    ResidualResult res;
    LineParser parser(lines);
    
    // These lambdas are just for concision in the next block
    auto readn = [&](std::size_t n){ return parser.read_Nnum_and_increment(n); };
    auto read1 = [&](){ return parser.read_Nnum_and_increment(1)[0]; };
    auto read1str = [&](){ return parser.read_1str_and_increment(); };
    
    auto model_key = parser.read_1str_from_line(lines[1]);
    
    if (model_key == "BWR"){
        return BWR2FEQ(lines);
    }
    else if (model_key == "ECS"){
        return parse_ECS(lines);
    }
    else if (model_key != "FEQ"){
        throw std::invalid_argument("Cannot parse this EOS type:" + model_key);
    }
    
    // Find the DOI if present
    std::string DOI_EOS = "??";
    
    // Find the first non-header row;
    size_t i = std::string::npos; 
    for (auto j = 2U; j < lines.size(); ++j){
        // First line not started by element in {:,?,!}, stop
        auto ind = lines[j].find_first_of(":?!");
        
        auto indEOS = lines[j].find(":DOI: ");
        if(!lines[j].empty() && (indEOS == 0)){
            std::string init = lines[j].substr(indEOS+6, lines[j].size()-6);
            DOI_EOS = internal::strip_trailing_whitespace(init);
        }
                                              
        if(!lines[j].empty() && (ind == 0)){
            continue;
        } 
        i = j; break;
    }
    parser.set_i(i);
    
    // Read in all the metadata parameters
    res.Tmin_K = read1();
    res.Tmax_K = read1();
    res.pmax_kPa = read1();
    res.rhomax_molL = read1();
    res.cp0_pointer = read1str();
    res.MM_kgkmol = read1();
    res.Ttriple_K = read1();
    res.ptriple_kPa = read1();
    res.rhotriple_molL = read1();
    res.Tnbp_K = read1();
    res.acentric = read1();
    auto c = readn(3);
    res.Tcrit_K = c[0];
    res.pcrit_kPa = c[1];
    res.rhocrit_molL = c[2];
    auto r = readn(2);
    res.Tred_K = r[0];
    res.rhored_molL = r[1];
    res.R = read1();
    res.DOI_EOS = DOI_EOS;

    // And now we read the EOS
    auto termcounts = readn(12);
    std::size_t Nnormal = static_cast<std::size_t>(termcounts[0]);
    std::size_t Nnormalcount = static_cast<std::size_t>(termcounts[1]); // Parameters per term
    std::size_t NGaussian = static_cast<std::size_t>(termcounts[2]);
    std::size_t NGaussiancount = static_cast<std::size_t>(termcounts[3]); // Parameters per term
    std::size_t NGao = static_cast<std::size_t>(termcounts[4]);

    auto read_normal = [&readn](const std::vector<std::string> &lines, std::size_t Nnormalcount) -> nlohmann::json {
        std::vector<double> n, t, d, l, g;
        for (auto line : lines){
            auto z = readn(Nnormalcount);
            n.push_back(z[0]);
            t.push_back(z[1]);
            d.push_back(z[2]);
            if (Nnormalcount == 4){
                l.push_back(z[3]);
            }
            if (Nnormalcount == 5){
                double l_ = z[3];
                l.push_back(l_);
                if (l_ == 0){
                    g.push_back(0);
                }
                else{
                    g.push_back(z[4]);
                }
            }
        }
        if (Nnormalcount <= 4){
            return {{"type", "ResidualHelmholtzPower"}, {"n", n}, {"t", t}, {"d", d}, {"l", l}};
        }
        else if (Nnormalcount == 5){
            return {{"type", "ResidualHelmholtzExponential"}, {"n", n}, {"t", t}, {"d", d}, {"l", l}, {"g", g}};
        }
        else{
            throw std::invalid_argument("Normalcould is not correct");
        }
    };

    auto read_Gaussian = [&readn](const std::vector<std::string> &lines, size_t NGaussiancount) -> nlohmann::json {
        // These terms can be of a number of different kinds
        // and they are disambiguated based upon the parameters in the term
        // (this is not documented anywhere in REFPROP)

        nlohmann::json normal, nonanalyt, R125, doubleexp;
        auto init_normal = [&normal](){
            std::vector<double> _;
            normal = {
                {"n", _}, {"t", _}, {"d", _}, 
                {"eta", _}, {"beta", _}, {"gamma", _}, {"epsilon", _},
                {"type", "ResidualHelmholtzGaussian"}
            };
        };
        auto init_doubleexp = [&doubleexp](){
            std::vector<double> _;
            doubleexp = {
                {"n", _}, {"t", _}, {"d", _},
                {"ld", _}, {"gd", _}, {"lt", _}, {"gt", _},
                {"type", "ResidualHelmholtzDoubleExponential"}
            };
        };
        auto init_nonanalyt = [&nonanalyt](){
            std::vector<double> _;
            nonanalyt = {
                {"n", _}, {"a", _}, {"b", _}, {"beta", _}, 
                {"A", _}, {"B", _}, {"C", _}, {"D", _}, 
                {"type", "ResidualHelmholtzNonAnalytic"}
            };
        };
        auto init_R125= [&R125](){
            std::vector<double> _;
            R125 = {
                {"n", _}, {"t", _}, {"d", _}, {"l", _}, {"m", _},
                {"type", "ResidualHelmholtzLemmon2005"}
            };
        };
        
        for (auto i = 0U; i < lines.size(); ++i){
            auto z = readn(NGaussiancount);
            assert(NGaussiancount == 12);
            // Determine what kind of term it is and store accordingly
            if (z[3] == 2 && z[4] == 2 && z[9] == 0 && z[10] == 0 && z[11] == 0){
                // It is a normal Gaussian term
                if (normal.empty()){ init_normal(); }
                normal["n"].push_back(z[0]);
                normal["t"].push_back(z[1]);
                normal["d"].push_back(z[2]);
                normal["eta"].push_back(-z[5]);
                normal["beta"].push_back(-z[6]);
                normal["gamma"].push_back(z[7]);
                normal["epsilon"].push_back(z[8]);
            }
            else if (z[5] == -1 && z[6] == -1 && z[7] == 0 && z[8] == 0 && z[9] == 0 && z[10] == 0 && z[11] == 0){
                // It is a R125 term from Lemmon, 2005
                if (R125.empty()){ init_R125(); }
                R125["n"].push_back(z[0]);
                R125["t"].push_back(z[1]);
                R125["d"].push_back(z[2]);
                R125["l"].push_back(z[3]);
                R125["m"].push_back(z[4]);
            }
            else if (z[4] == 1 && z[1] == 0 && z[8] == 0 && z[9] == 0 && z[10] == 0 && z[11] == 0){
                // It is a methanol term from de Reuck IUPAC monograph
                // Conversions from de Reuck were done to get in REFPROP-compatible format,
                // and then the format is again converted to go from REFPROP format to double-exponential
                // form accepted in teqp.
                // The term started life in the form N_i*delta^(r_i)*exp(c_i*C_0*tau - b_i - (j_i*C_1*delta)^(k_i))
                // and then converted to the form
                // N_i*delta^d_i*EXP[eta_i*(delta)^de_i+beta_i*(tau_i-gamma_i)]
                // and is once again converted, where the leading coefficient becomes
                // N_i*exp(-beta_i*gamma_i)
                // td are 0
                // lt are 1
                // The term in teqp is defined as:
                //
                //\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\gamma_{d,i}\delta^{l_{d,i}}-\gamma_{t,i}\tau^{l_{t,i}})
                // where gt and
                double N=z[0], d=z[2], de=z[3], eta=z[5], beta=z[6], gamma=z[7];
                if (doubleexp.empty()){ init_doubleexp(); }
                doubleexp["n"].push_back(N*exp(-beta*gamma));
                doubleexp["t"].push_back(0);
                doubleexp["d"].push_back(d);
                doubleexp["gd"].push_back(-eta);
                doubleexp["ld"].push_back(de);
                doubleexp["gt"].push_back(-beta);
                doubleexp["lt"].push_back(1);
            }
            else if (z[3] == 2 && z[4] == 2 && z[9] != 0 && z[10] != 0 && z[11] != 0){
                // It is a non-analytic term
                if (nonanalyt.empty()){ init_nonanalyt(); }
                nonanalyt["n"].push_back(z[0]);
                nonanalyt["b"].push_back(z[5]);
                nonanalyt["beta"].push_back(z[6]);
                nonanalyt["A"].push_back(z[7]);
                nonanalyt["C"].push_back(z[8]);
                nonanalyt["D"].push_back(z[9]);
                nonanalyt["B"].push_back(z[10]);
                nonanalyt["a"].push_back(z[11]);
            }
            else{
                throw std::invalid_argument("Bad Gaussian term");
            }
        }
        nlohmann::json o = nlohmann::json::array();
        if(!normal.empty()){ o.push_back(normal); }
        if(!nonanalyt.empty()){ o.push_back(nonanalyt); }
        if(!R125.empty()){ o.push_back(R125); }
        if(!doubleexp.empty()){ o.push_back(doubleexp); }
        return o;
    };
    auto read_Gao = [&readn](const std::vector<std::string> &lines, size_t NGaocount) -> nlohmann::json{
        std::vector<double> n,t,d,eta,beta,gamma,epsilon,b;
        for (auto i = 0U; i < lines.size(); ++i){
            auto z = readn(NGaocount);
            assert(NGaocount == 12);
            n.push_back(z[0]);
            t.push_back(z[1]);
            d.push_back(z[2]);
            eta.push_back(z[5]);
            beta.push_back(z[6]);
            gamma.push_back(z[7]);
            epsilon.push_back(z[8]);
            b.push_back(z[9]);
        }
        return {{"n", n}, {"t", t}, {"d", d},
            {"eta", eta}, {"beta", beta}, {"gamma", gamma}, {"epsilon", epsilon}, {"b", b},
            {"type", "ResidualHelmholtzGaoB"}};
    };

    res.alphar = nlohmann::json::array();
    if (Nnormal > 0){
        auto l = std::vector<std::string>(lines.begin()+i, lines.begin()+i+Nnormal);
        res.alphar.push_back(read_normal(l, Nnormalcount));
    }
    if (NGaussian > 0){
        auto l = std::vector<std::string>(lines.begin()+i+Nnormal, lines.begin()+i+NGaussian+Nnormal);
        for (auto thing : read_Gaussian(l, NGaussiancount)){
            res.alphar.push_back(thing);
        }
    }
    if (NGao > 0){
        auto l = std::vector<std::string>(lines.begin()+i+NGaussian+Nnormal, lines.begin()+i+Nnormal+NGaussian+NGao);
        res.alphar.push_back(read_Gao(l, 12));
    }
    return res;
}

/**
 * Things that are obtained from parsing an ideal-gas Helmholtz
 * energy EOS
 */
struct Alpha0Result{
    double Tmin_K,Tmax_K,pmax_kPa,rhomax_molL;
    std::string cp0_pointer;
    double Tred_K, cp0red_JmolK, R;
    nlohmann::json alpha0;
};

/**
 \param lines The set of lines used for the ideal-gas part
 \param Tri The reducing temperature used in the contribution, in K
 \param R_alphar_JmolK The gas constant coming from the default residual EOS contribution, in J/mol/K. This may (or may not) be used depending on whether a scaling parameter for R is provided
 */
inline Alpha0Result convert_CP0(const std::vector<std::string>& lines, double Tri, double R_alphar_JmolK){
    Alpha0Result a;
    
    LineParser parser(lines);
    
    // These lambdas are just for concision in the next block
    auto readnline = [&](const std::string& line, std::size_t n){ return parser.read_Nnum_from_line(line, n); };
    auto readn = [&](std::size_t n){ return parser.read_Nnum_and_increment(n); };
    auto read1 = [&](){ return parser.read_Nnum_and_increment(1)[0]; };
    
    a.cp0_pointer = parser.read_1str_from_line(lines[1]);
    
    // Find the first non-header row;
    size_t i = std::string::npos;
    for (auto j = 2U; j < lines.size(); ++j){
        // First line not started by element in {:,?,!}, stop
        auto ind = lines[j].find_first_of(":?!");
        if(!lines[j].empty() && (ind == 0)){
            continue;
        }
        i = j; break;
    }
    parser.set_i(i);
    
    // Read in all the metadata parameters
    a.Tmin_K = read1();
    a.Tmax_K = read1();
    a.pmax_kPa = read1();
    a.rhomax_molL = read1();
    auto r = readn(2);
    a.Tred_K = r[0];
    a.cp0red_JmolK = r[1];
    double R_used_JmolK = R_alphar_JmolK;
    
    auto N = readn(7);
    size_t Npoly = static_cast<std::size_t>(N[0]);
    size_t NPlanck = static_cast<std::size_t>(N[1]);
    for (auto k = 2U; k < N.size(); ++k){
        if (N[k] != 0){
            throw std::invalid_argument("Found a cp0 term that cannot currently be parsed");
        }
    }

    // cp0/R = sum_i (c_i*T^t_i)
    // the scaling factor on cp0 might not be the gas constant, so you need to correct because
    // REFPROP files use cp0/K = sum_i (c_i*(T/Q)^t_i) = sum_i (c_i/(Q^t_i) * T^t_i)
    // To get to cp0/R use (cp0/K)*(K/R)
    auto read_polynomial = [&readnline, &a, &Tri, &R_used_JmolK](const std::vector<std::string> &lines) -> nlohmann::json {
        std::vector<double> c, t;
        for (auto line : lines){
            auto z = readnline(line, 2);
            c.push_back(z[0]*a.cp0red_JmolK/R_used_JmolK/pow(a.Tred_K, z[1]));
            t.push_back(z[1]);
        }
        return {{"type", "IdealGasHelmholtzCP0PolyT"}, {"c", c}, {"t", t}, {"T0", 298.15}, {"Tc", Tri}, {"R", R_used_JmolK}};
    };
    
    auto read_Planck = [&readnline, &a, &Tri, &R_used_JmolK](const std::vector<std::string> &lines) -> nlohmann::json {
        std::vector<double> n, v;
        for (auto line : lines){
            auto z = readnline(line, 2);
            n.push_back(z[0]*a.cp0red_JmolK/R_used_JmolK);
            v.push_back(z[1]*a.Tred_K);
        }
        return {{"type", "IdealGasHelmholtzPlanckEinsteinFunctionT"}, {"n", n}, {"v", v}, {"T0", 298.15}, {"Tcrit", Tri}, {"R", R_used_JmolK} };
    };

    a.alpha0 = nlohmann::json::array();
    // The log(rho) term that is always included
    a.alpha0.push_back({{"a1", 0.0}, {"a2", 0.0}, {"type", "IdealGasHelmholtzLead"}});
    a.alpha0.push_back({{"a", -1}, {"type", "IdealGasHelmholtzLogTau"}});
    i = parser.get_i();
    if (Npoly > 0){
        auto l = std::vector<std::string>(lines.begin()+i, lines.begin()+i+Npoly);
        a.alpha0.push_back(read_polynomial(l));
    }
    if (NPlanck > 0){
        auto l = std::vector<std::string>(lines.begin()+i+Npoly, lines.begin()+i+NPlanck+Npoly);
        a.alpha0.push_back(read_Planck(l));
    }
    
    return a;
}

/**
 * Things that are obtained from parsing a residual Helmholtz
 * energy EOS
 */
struct HeaderResult{
    std::string short_name, CASnum, full_name, chemical_formula, synonym;
    double MM_kgkmol = -1, Ttriple_K, Tnbp_K, Tcrit_K, pcrit_kPa, rhocrit_molL, acentric, dipole_D;
    std::string refstate, RPversion;
    std::string UNnumber = "?";
    std::string family = "?";
    double HVupper_kJmol = 99999999, GWP100=99999999, RCL_vv = 99999999, ODP=99999999;
    std::string ASHRAE34 = "?", StdInChIstr = "?", StdInChIKey = "?", altid = "?", hash = "?";
    auto to_int(const std::string& s) -> int {
        return static_cast<int>(strtol(s.c_str(), nullptr, 10));
    }
    auto to_double(const std::string& s) -> double {
        return static_cast<int>(strtod(s.c_str(), nullptr));
    }
    void set(const std::string &k, const std::string &val){
        if (k == "FAMILY"){
            family = val;
        }
        else if (k == "UN"){
            UNnumber = val;
        }
        else if (k == "INCHI"){
            StdInChIstr = val;
            // The string must start with InChI=, otherwise it is not a standard InChI string
            // so prepend the necessary prefix
            if (StdInChIstr.find("InChI=") != 0){
                StdInChIstr = "InChI=" + StdInChIstr;
            }
        }
        else if (k == "INCHIKEY"){
            StdInChIKey = val;
        }
        else if (k == "SAFETY"){
            ASHRAE34 = val;
        }
        else if (k == "ALTID"){
            altid = val;
        }
        else if (k == "HASH"){
            hash = val;
        }
        else if (k == "HEAT"){
            HVupper_kJmol = to_double(val);
        }
        else if (k == "GWP"){
            GWP100 = to_double(val);
        }
        else if (k == "ODP"){
            ODP = to_double(val);
        }
        else if (k == "RCL"){
            RCL_vv = to_double(val);
        }
        else{
            throw std::invalid_argument(k);
        }
    }
};

inline HeaderResult convert_header(const std::vector<std::string>& lines){
    HeaderResult h;
    std::size_t i = 0;
    auto _read1str = [](const std::string &line){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), "  ");
        if (vals.empty()){throw std::invalid_argument("Unable to read one string from this line:"+line); }
        return vals[0];
    };
    auto read1str = [&lines, &i, &_read1str](){
        std::string val = _read1str(lines[i]); i++; return val;
    };
    auto _read1 = [](const std::string &line){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), " ");
        if (vals.empty()){throw std::invalid_argument("Unable to read one number from this line:"+line); }
        return strtod(vals[0].c_str(), nullptr);
    };
    auto read1 = [&lines, &i, &_read1](){
        double val = _read1(lines[i]); i++; return val;
    };
    h.short_name = read1str();
    h.CASnum = read1str();
    h.full_name = read1str();
    h.chemical_formula = read1str();
    h.synonym = read1str();
    h.MM_kgkmol = read1();
    h.Ttriple_K = read1();
    h.Tnbp_K = read1();
    h.Tcrit_K = read1();
    h.pcrit_kPa = read1();
    h.rhocrit_molL = read1();
    h.acentric = read1();
    h.dipole_D = read1();
    h.refstate = read1str();
    if (h.refstate.find("OT") == 0){ read1(); }
    h.RPversion = read1str();
    
    if (h.RPversion == "10.0"){
        for (; i < lines.size(); ++i){
            std::smatch match;
            // If an empty string
            if (lines[i].find_first_not_of(" \t\n\v\f\r") == std::string::npos){
                break;
            }
            if (!std::regex_search(lines[i], match, std::regex(R"(:(\w+):)"))){
                throw std::invalid_argument("This line doesn't contain a key contained in colons: "+lines[i]+"; it has length of "+std::to_string(lines[i].size()));
            }
            using namespace internal;
            h.set(internal::to_upper(match[1]), strip_trailing_whitespace(strip_line_comment(lines[i])));
        }
    }
    else if (h.RPversion == "9.0" || h.RPversion == "8.0"){
        h.UNnumber = read1str();
        h.family = read1str();
        if (i < lines.size()-1){
            h.HVupper_kJmol = read1();
        }
    }
    
    return h;
}

inline auto get_ancillary_description(const std::string& key){
    const std::map<std::string, std::string> keys = {
        // Most are in this form:
        {"PS5",  "P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T]"},
        {"DL1",  "D=Dc*[1+SUM(Ni*Theta^ti)]"},
        {"DV3",  "D=Dc*EXP[SUM(Ni*Theta^ti)]"},
        
        // Other ones might be of this form:
        {"DL2",  "D=Dc*[1+SUM(Ni*Theta^(ti/3))]"},
        {"DL3",  "D=Dc*EXP[SUM(Ni*Theta^ti)]"},
        {"DL4",  "D=Dc*EXP[SUM(Ni*Theta^(ti/3))]"},
        {"DL5",  "D=Dc*EXP[SUM(Ni*Theta^ti)*Tc/T]"},
        {"DL6",  "D=Dc*EXP[SUM(Ni*Theta^(ti/3))*Tc/T]"},
        {"DV4",  "D=Dc*EXP[SUM(Ni*Theta^(ti/3))]"},
        {"DV5",  "D=Dc*EXP[SUM(Ni*Theta^ti)*Tc/T]"},
        {"DV6",  "D=Dc*EXP[SUM(Ni*Theta^(ti/3))*Tc/T]"},
        {"PS6",  "P=Pc*EXP[SUM(Ni*Theta^(ti/3))*Tc/T]"},
    };
    if (keys.count(key) > 0){
        return keys.at(key);
    }
    else{
        throw std::invalid_argument("Bad ancillary model for getting description: " + key);
    }
};

inline nlohmann::json get_ancillary(const std::vector<std::string>& lines){
    
    LineParser parser(lines);
    
    // These lambdas are just for concision in the next block
    auto readn = [&](std::size_t n){ return parser.read_Nnum_and_increment(n); };
    auto read1 = [&](){ return parser.read_Nnum_and_increment(1)[0]; };
    
    auto modelname = lines[1].substr(0, 3);
    
    // Find the first non-header row;
    size_t i = std::string::npos;
    for (auto j = 2U; j < lines.size(); ++j){
        // First line not started by element in {:,?,!}, stop
        auto ind = lines[j].find_first_of(":?!");
        if(!lines[j].empty() && (ind == 0)){
            continue;
        }
        i = j; break;
    }
    parser.set_i(i);
    
    double Tmin_K = read1();
    double Tmax_K = read1();
    read1(); // placeholder
    read1(); // placeholder
    auto reducing = readn(2);
    reducing[1] *= 1000; // REFPROP uses mol/L & kPa, CoolProp uses mol/m^3 & Pa, multiply by 1000 to convert
    auto Ncoeffs = readn(6);
    
    std::vector<double> n, t;
    for (auto k = 0; k < Ncoeffs[0]; ++k){
        auto z = readn(2);
        n.push_back(z[0]);
        t.push_back(z[1]);
    }
    //assert(len(n) == Ncoeffs[0] and len(t) == Ncoeffs[0]);
    
    auto desc = get_ancillary_description(modelname);
    
    using seti = std::set<int>;
    using sets = std::set<std::string>;
    std::string model_key = modelname.substr(0,2);
    int key_index = modelname[2] - '0'; // See https://stackoverflow.com/a/628766
    
    std::string type_key;
    if (model_key == "PS"){
        type_key = "pL";
    }
    else if (sets{"DL1", "DL2"}.count(modelname) > 0){
        type_key = "rhoLnoexp";
        // Convert DL2 to DL1
        if (seti{2,4,6}.count(key_index) > 0){
            for (auto&t_ : t){ t_ /= 3; }
            desc = get_ancillary_description("DL1");
        }
    }
    else if (model_key == "DL"){
        type_key = "rhoL";
        if (seti{2,4,6}.count(key_index) > 0){
            for (auto&t_ : t){ t_ /= 3; }
            desc = get_ancillary_description("DL" + std::to_string(key_index-1));
        }
    }
    else if (model_key == "DV"){
        type_key = "rhoV";
        if (seti{2,4,6}.count(key_index) > 0){
            for (auto&t_ : t){ t_ /= 3; }
            desc = get_ancillary_description("DV" + std::to_string(key_index-1));
        }
    }
    else{
        throw std::invalid_argument(model_key);
    }
    
    return {
        {"T_r", reducing[0]},
        {"Tmax", std::min(Tmax_K, reducing[0])},
        {"Tmin", Tmin_K},
        {"description", desc},
        {"n", n},
        {"reducing_value", reducing[1]},
        {"t", t},
        {"type", type_key},
        {"using_tau_r", (desc.find("Tc/T") != std::string::npos)}
    };
}

inline nlohmann::json get_all_ancillaries(const std::vector<std::string>& lines){
    return {
        {"pS", get_ancillary(internal::get_line_chunk(lines, "#PS"))},
        {"rhoV", get_ancillary(internal::get_line_chunk(lines, "#DV"))},
        {"rhoL", get_ancillary(internal::get_line_chunk(lines, "#DL"))}
    };
}

inline auto get_rhoLV_ancillaries(double T, const nlohmann::json& jancillaries){
    auto eval_ancillary = [T](const nlohmann::json &j){
        std::valarray<double> n = j.at("n");
        std::valarray<double> t = j.at("t");
        double reducing_value = j.at("reducing_value");
        double T_r = j.at("T_r");
        auto Theta = 1.0 - T/T_r;
        auto RHS = (pow(Theta, t)*n).sum();
        if (j.at("using_tau_r")){
            RHS *= T_r/T;
        }
        if (j.at("type") == "rhoLnoexp"){
            return reducing_value*(1+RHS);
        }
        else{
            return exp(RHS)*reducing_value;
        }
    };
    return std::make_tuple(eval_ancillary(jancillaries.at("rhoL")), eval_ancillary(jancillaries.at("rhoV")));
}

class FLDfile{
private:
    const std::string contents;
    const std::vector<std::string> lines;
    std::vector<std::string> warnings, errors;
public:
    FLDfile(const std::filesystem::path &path) : contents(internal::read_file(path)), lines(internal::strsplit(contents, "\n")){ };
    
    auto get_warnings(){ return warnings; }
    auto get_errors(){ return errors; }

    nlohmann::json convert_EOS(const RPinterop::ResidualResult& feq, const RPinterop::Alpha0Result& alpha0, const HeaderResult& head, const nlohmann::json& ancillaries){
        auto [rhoLtriple, rhoVtriple] = get_rhoLV_ancillaries(feq.Tmin_K, ancillaries);
        nlohmann::json STATES = {
            {"reducing", {
              {"T", feq.Tred_K},
              {"T_units", "K"},
              {"hmolar", -99999999999.0},
              {"hmolar_units", "J/mol"},
              {"p", feq.pcrit_kPa*1e3},
              {"p_units", "Pa"},
              {"rhomolar", feq.rhored_molL*1e3},
              {"rhomolar_units", "mol/m^3"},
              {"smolar", 999999999999999.0},
              {"smolar_units", "J/mol/K"}
            }},
            {"sat_min_liquid", {
              {"T", feq.Tmin_K},
              {"T_units", "K"},
              {"hmolar", -999999999999.0},
              {"hmolar_units", "J/mol"},
              {"p", feq.ptriple_kPa*1e3},
              {"p_units", "Pa"},
              {"rhomolar", rhoLtriple},
              {"rhomolar_units", "mol/m^3"},
              {"smolar", 99999999999999999999.0},
              {"smolar_units", "J/mol/K"}
            }},
            {"sat_min_vapor", {
              {"T", feq.Tmin_K},
              {"T_units", "K"},
              {"hmolar", 9999999999999.0},
              {"hmolar_units", "J/mol"},
              {"p", feq.ptriple_kPa*1e3},
              {"p_units", "Pa"},
              {"rhomolar", rhoVtriple},
              {"rhomolar_units", "mol/m^3"},
              {"smolar", 9999999999999999999.0},
              {"smolar_units", "J/mol/K"}
            }}
        };

        nlohmann::json EOS = {
          {"BibTeX_CP0", ""},
          {"BibTeX_EOS", feq.DOI_EOS},
          {"STATES", STATES},
          {"T_max", feq.Tmax_K},
          {"T_max_units", "K"},
          {"Ttriple", feq.Ttriple_K},
          {"Ttriple_units", "K"},
          {"acentric", feq.acentric},
          {"acentric_units", "-"},
          {"alpha0", alpha0.alpha0},
          {"alphar", feq.alphar},
          {"gas_constant", feq.R},
          {"gas_constant_units", "J/mol/K"},
          {"molar_mass", feq.MM_kgkmol/1e3},
          {"molar_mass_units", "kg/mol"},
          {"p_max", feq.pmax_kPa*1e3},
          {"p_max_units", "Pa"},
          {"pseudo_pure", false} /// TODO
        };
        if (feq.MM_kgkmol < 0){
            EOS["molar_mass"] = head.MM_kgkmol/1e3;
        }
        if (feq.R < 0){
            EOS["gas_constant"] = 8.31446261815324;
        }
        return EOS;
    }

    nlohmann::json get_critical_state(const RPinterop::HeaderResult &head){
        return {
            {"T", head.Tcrit_K},
            {"T_units", "K"},
            {"hmolar", -99999999999.0},
            {"hmolar_units", "J/mol"},
            {"p", head.pcrit_kPa*1e3},
            {"p_units", "Pa"},
            {"rhomolar", head.rhocrit_molL*1e3},
            {"rhomolar_units", "mol/m^3"},
            {"smolar", 999999999999999.0},
            {"smolar_units", "J/mol/K"}
        };
    };

    /** 
     * Parse the standard InChI string to extract the chemical formula in Hill order 
    */
    // auto formula_from_inchi(const string& stdinchistring){
        
    //     if '/' not in stdinchistring:
    //         raise ValueError(f'{stdinchistring} is not a valid standard InChI key')
    //     formula = stdinchistring.split('/')[1]
    //     matches = re.findall(r'[A-Z][a-z]*[0-9]*', formula)
    //     assert(sum([len(m) for m in matches]) == len(formula)) # make sure all parts are matched
    //     o = ''
    //     for match in matches:
    //         def replacer(m):
    //             N = len(m.groups())
    //             # print(m, N)
    //             if N == 0:
    //                 return m.group(0)+'_{1}'
    //             else:
    //                 i = m.group(1)
    //                 return f'_{{{i}}}'
    //         newmatch = re.sub(r'([0-9]+)', replacer, match)
    //         if '_{' not in newmatch:
    //             newmatch += '_{1}'
    //         o += newmatch
    //     return o
    // }

    nlohmann::json get_info(const std::string& name, const RPinterop::HeaderResult &head){
        int CHEMSPIDER_ID = -1;
        return {
            {"2DPNG_URL", "http://www.chemspider.com/ImagesHandler.ashx?id=" + std::to_string(CHEMSPIDER_ID)},
            {"ALIASES", nlohmann::json::array()},
            {"CAS", head.CASnum},
            {"CHEMSPIDER_ID", CHEMSPIDER_ID},
            {"ENVIRONMENTAL", {
              {"ASHRAE34", head.ASHRAE34},
              {"FH", 99999999},
              {"GWP100", head.GWP100},
              {"GWP20", -99999999999},
              {"GWP500", -9999999999},
              {"HH", -999999999},
              {"Name", name},
              {"ODP", head.ODP},
              {"PH", -999999}
            }},
            // {"FORMULA", formula_from_inchi(head.StdInChIstr)},
            {"FORMULA", head.chemical_formula},
            {"INCHI_KEY", head.StdInChIKey},
            {"INCHI_STRING", head.StdInChIstr},
            {"NAME", name},
            {"REFPROP_NAME", name},
            {"HASH", head.hash},
            {"SMILES", "?"},
            {"DIPOLE", head.dipole_D*3.33564e-30},
            {"DIPOLE_units", "C*m"}
        };
    }

    auto make_json(const std::string& name){
        auto head = convert_header(lines);
        
        nlohmann::json ancillaries = nlohmann::json::object();
        try{
            ancillaries = get_all_ancillaries(lines);
        }
        catch(std::exception &e){
            warnings.push_back(std::string("Could not load ancillaries; message:") + e.what());
//            std::cerr << e.what() << std::endl;
        }
        
        auto feq = convert_FEQ(internal::get_line_chunk(lines, "#EOS"));
        auto alpha0 = convert_CP0(internal::get_line_chunk(lines, "#AUX"), feq.value().Tred_K, feq.value().R);
        auto EOS = convert_EOS(feq.value(), alpha0, head, ancillaries);

        nlohmann::json f = {
            {"EOS", {EOS}},
            {"ANCILLARIES", ancillaries},
            {"INFO", get_info(name, head)}
        };
        f["STATES"] = {
            {"critical", get_critical_state(head)},
            {"triple_liquid", f["EOS"][0]["STATES"]["sat_min_liquid"]},
            {"triple_vapor", f["EOS"][0]["STATES"]["sat_min_vapor"]}
        };
        return f;
    }
};


class HMXBNCfile{
private:
    const std::string contents;
    const std::vector<std::string> lines;
    std::vector<std::string> warnings, errors;
public:
    HMXBNCfile(const std::filesystem::path &path) : contents(internal::read_file(path)), lines(internal::strsplit(contents, "\n")){ };
    
    auto get_warnings(){ return warnings; }
    auto get_errors(){ return errors; }
    
    auto convert_BNC(){
        auto BNCchunk = internal::get_line_chunk(lines, "#BNC", 0);
        
        int j = static_cast<int>(BNCchunk.size()-1);
        int jmax = j;
        nlohmann::json entries = nlohmann::json::array();
        while (j > 0){
            // Locate the next pair
            for (; j > 0; --j){
                if(BNCchunk[j].find("!") == 0 && j != jmax)
                    break;
            }
            std::vector<std::string> header, hashes;
            if (j == 0){
                break;
            }
            
            for (auto k = j+1; k < jmax; ++k){
                auto line = internal::strip_trailing_whitespace(internal::strip_leading_whitespace(BNCchunk[k]));
                if (line.empty()){ continue; }
                if (line.find("?") == 0){
                    header.push_back(line.substr(1));
                }
                else if (line.find("/") != std::string::npos){
                    hashes = internal::strsplit(line, "/");
                }
                else if (line.find("TC") == 0 || line.find("VC") == 0){
                }
                else if (line.find("PR1") == 0){
                }
                else if (line.find("TRN") == 0){
                }
                else if (line.find("ST1") == 0){
                }
                else{
                    auto elements = internal::strsplit(internal::strip_line_comment(line), " ");
                    
                    std::vector<double> nums;
                    std::string function_name;
                    for (const auto& el : elements){
                        // First element is the name of the function
                        if (function_name.empty()){
                            function_name = el;
                            continue;
                        }
                        if (el.size() > 0){
                            nums.push_back(strtod(el.c_str(), nullptr));
                        }
                    }
                    
                    entries.push_back({
                        {"info", internal::strjoin(header,"\n")},
                        {"BiBTeX", "?"},
                        {"hash1", hashes[0]},
                        {"hash2", hashes[1]},
                        {"betaT", nums[0]},
                        {"gammaT", nums[1]},
                        {"betaV", nums[2]},
                        {"gammaV", nums[3]},
                        {"F", nums[4]},
                        {"function", function_name}
                    });
                }
            }
            jmax = j;
        }
        return entries;
    }
    
    auto convert_MXM(){
        std::vector<std::size_t> starting_indices;
        for (auto i = 0U; i < lines.size(); ++i){
            if (lines[i].find("#MXM") == 0){
                starting_indices.push_back(i);
            }
        }
        
        nlohmann::json terms = nlohmann::json::array();
        
        for (auto istart: starting_indices){
            auto chunk = internal::get_line_chunk(lines, "#MXM", istart);
            
            if (chunk[1].find("XR0 ") == 0){
                continue;
            }
            std::string function_name = chunk[1].substr(0, 3);
            
            std::vector<std::string> header_lines;
            // Get all the header information
            std::size_t i = 2;
            for (; i < chunk.size(); ++i){
                if (chunk[i].find("!") == 0){
                    break;
                }
                auto line = internal::strip_line_leading_comment(chunk[i], "?");
                if (line.find("`") != 0){
                    header_lines.push_back(line);
                }
            }
            if (i == chunk.size()){
                continue;
            }
            std::string header = internal::strip_trailing_whitespace(internal::strjoin(header_lines, "\n"));
            
            // Figure out the number and types of terms
            auto elements = internal::strsplit(internal::strip_line_comment(chunk[i+3]), " ");
            
            std::vector<int> numbers;
            for (const auto& el : elements){
                if (el.size() > 0){
                    numbers.push_back(static_cast<int>(strtod(el.c_str(), nullptr)));
                }
            }
            if (numbers.size() != 11){
                throw std::invalid_argument("??");
            }
            std::size_t Npower = numbers[0], Npower_coeffs = numbers[1], N_GERG=numbers[3], N_GERG_coeffs=numbers[4], N_exp=numbers[5], N_exp_coeffs=numbers[6];
            
            std::vector<double> n,t,d,c,eta,beta,gamma,epsilon;
            
            if (Npower > 0 && N_GERG == 0 && N_exp == 0){
                for (auto k = 1U; k <= Npower; ++k){
                    auto nums = internal::read_allnum_from_line(chunk[i+3+k]);
                    if (nums.size() != Npower_coeffs){
                        throw std::invalid_argument("wrong length");
                    }
                    n.push_back(nums[0]);
                    t.push_back(nums[1]);
                    d.push_back(nums[2]);
                    c.push_back(nums[3]);
                }
                if (n.size() != Npower){
                    throw std::invalid_argument("Didn't collect the right number of terms");
                }
                terms.push_back({
                    {"Name", function_name},
                    {"type", "Exponential"},
                    {"kind", "alpharij = sum_i n_i*tau^t_i*delta^d_i*exp(-delta^c_i)"},
                    {"description", header},
                    {"n",n}, {"t",t}, {"d",d}, {"l", c}
                });
            }
            else if (Npower > 0 && N_GERG > 0 && N_exp == 0){
                for (auto k = 1U; k <= Npower; ++k){
                    auto nums = internal::read_allnum_from_line(chunk[i+3+k]);
                    if (nums.size() != Npower_coeffs){
                        throw std::invalid_argument("wrong length");
                    }
                    n.push_back(nums[0]);
                    t.push_back(nums[1]);
                    d.push_back(nums[2]);
                    c.push_back(nums[3]);
                    eta.push_back(0);
                    beta.push_back(0);
                    gamma.push_back(0);
                    epsilon.push_back(0);
                }
                for (auto k = Npower+1; k <= Npower+N_GERG; ++k){
                    auto nums = internal::read_allnum_from_line(chunk[i+3+k]);
                    if (nums.size() != N_GERG_coeffs){
                        throw std::invalid_argument("wrong length");
                    }
                    n.push_back(nums[0]);
                    t.push_back(nums[1]);
                    d.push_back(nums[2]);
                    eta.push_back(nums[3]);
                    epsilon.push_back(nums[4]);
                    beta.push_back(nums[5]);
                    gamma.push_back(nums[6]);
                }
                if (n.size() != Npower + N_GERG){
                    throw std::invalid_argument("Didn't collect the right number of terms");
                }
                terms.push_back({
                    {"Name", function_name},
                    {"type", "GERG-2008"},
                    {"Npower", Npower},
                    {"kind", "alpharij = sum_i n_i*tau^t_i*delta^d_i*exp(-beta_i*(tau-gamma_i)^2-eta_i*(delta-epsilon_i)^2)"},
                    {"description", header},
                    {"n",n}, {"t",t}, {"d",d}, {"eta",eta}, {"beta",beta}, {"gamma",gamma}, {"epsilon",epsilon},
                });
            }
            else if (Npower > 0 && N_GERG == 0 && N_exp > 0){
                for (auto k = 1U; k <= Npower; ++k){
                    auto nums = internal::read_allnum_from_line(chunk[i+3+k]);
                    if (nums.size() != Npower_coeffs){
                        throw std::invalid_argument("wrong length");
                    }
                    n.push_back(nums[0]);
                    t.push_back(nums[1]);
                    d.push_back(nums[2]);
                    c.push_back(nums[3]);
                    eta.push_back(0);
                    beta.push_back(0);
                    gamma.push_back(0);
                    epsilon.push_back(0);
                }
                for (auto k = Npower+1; k <= Npower+N_exp; ++k){
                    auto nums = internal::read_allnum_from_line(chunk[i+3+k]);
                    if (nums.size() != N_exp_coeffs){
                        throw std::invalid_argument("wrong length");
                    }
                    n.push_back(nums[0]);
                    t.push_back(nums[1]);
                    d.push_back(nums[2]);
                    c.push_back(0);
                    eta.push_back(nums[3]);
                    epsilon.push_back(nums[4]);
                    beta.push_back(nums[5]);
                    gamma.push_back(nums[6]);
                    
                }
                if (n.size() != Npower + N_exp){
                    throw std::invalid_argument("Didn't collect the right number of terms");
                }
                
                terms.push_back({
                    {"Name", function_name},
                    {"type", "Gaussian+Exponential"},
                    {"Npower", Npower},
                    {"kind", "alpharij = sum_i n_i*tau^t_i*delta^d_i*exp(-beta_i*(tau-gamma_i)^2-eta_i*(delta-epsilon_i)^2)"},
                    {"description", header},
                    {"n",n}, {"t",t}, {"d",d}, {"eta",eta}, {"beta",beta}, {"gamma",gamma}, {"epsilon",epsilon},
                });
            }
            else{
                throw std::invalid_argument("Unable to make sense of how to parse this set of terms");
            }
        }
        return terms;
    }
    
    auto make_jsons(){
        auto BIP = convert_BNC();
        auto deps = convert_MXM();
        return std::make_tuple(BIP, deps);
    }
};


}
