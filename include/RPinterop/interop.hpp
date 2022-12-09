#pragma once

#include <filesystem>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <regex>

#include "nlohmann/json.hpp"

//using namespace std;
using std::vector;
using std::string;

namespace RPinterop{

namespace internal{

// for string delimiter (https://stackoverflow.com/a/46931770)
vector<string> strsplit (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;
    
    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }
    
    res.push_back (s.substr (pos_start));
    return res;
}
string strjoin(const vector<string> &lines, const std::string& delim){
    if (lines.empty()){ return ""; }
    else{
        std::string o = lines[0];
        for (auto i = 1; i < lines.size(); ++i){
            o += delim + lines[i];
        }
        return o;
    }
}

string strip_line_comment(const string& line){
    auto ipos = line.find_first_of("!");
    if (ipos == std::string::npos){
        // No comment found, return whole string
        return line;
    }
    else{
        return line.substr(0, ipos);
    }
}
string strip_trailing_whitespace(const string& line){
    auto ipos = line.find_first_of(" \t\n\v\f\r");
    if (ipos == std::string::npos){
        // No whitespace found, return whole string
        return line;
    }
    else{
        return line.substr(0, ipos);
    }
}

std::vector<std::string> get_line_chunk(
                                        const std::vector<std::string>& lines,
                                        const string& starting_key
                                        )
{
    // block starts with the key
    size_t istart = std::string::npos;
    for (auto i = 0; i < lines.size(); ++i){
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
    return std::vector<std::string>(lines.begin()+istart, lines.begin()+iend);
}

string to_upper(const string& str){
    auto s = str; // copy :(
    std::locale locale;
    auto f = [&locale] (char ch) { return std::use_facet<std::ctype<char>>(locale).toupper(ch); };
    std::transform(s.begin(), s.end(), s.begin(), f);
    return s;
}
}

/**
 * Things that are obtained from parsing a residual Helmholtz
 * energy EOS
 */
struct ResidualResult{
    double Tmin_K,Tmax_K,pmax_kPa,rhomax_molL;
    std::string cp0_pointer;
    double MM_kgkmol,Ttriple_K,ptriple_kPa,rhotriple_molL,Tnbp_K,acentric,
    Tcrit_K,pcrit_kPa,rhocrit_molL,Tred_K,rhored_molL,R;
    nlohmann::json alphar;
    std::string DOI_EOS = "";
};

/**
Given a string that contains an EOS block, return its
JSON representation
*/
ResidualResult convert_FEQ(const vector<string>& lines){
    ResidualResult res;
    
    auto read1strline = [](const string &line){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), " ");
        if (vals.empty()){throw std::invalid_argument("Unable to read one string from this line:"+line); }
        return vals[0];
    };
    auto model_key = read1strline(lines[1]);
    if (model_key != "FEQ"){
        throw std::invalid_argument("Cannot parse this EOS type:" + model_key);
    }
    
    // Find the DOI if present
    std::string DOI_EOS = "??";
    
    // Find the first non-header row;
    size_t i = std::string::npos; 
    for (auto j = 2; j < lines.size(); ++j){
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
    auto readnline = [](const string& line, int n){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), " ");
        std::vector<double> o;
        for (auto v : vals){
            if (!v.empty()){
                std::string v_e = std::regex_replace(v, std::regex("[Dd]"), "e");
                std::string v_ew = std::regex_replace(v_e, std::regex("\\s+"), "");
                if (v_ew.empty()){
                    continue;
                }
                o.push_back(strtod(v_ew.c_str(), nullptr));
            }
        }
        if (n > 0){
            if (o.size() != n){ throw std::invalid_argument("Unable to read "+std::to_string(n)+" numbers from this line:"+line); }
        }
        return o;
    };
    auto readn = [&lines, &i, &readnline](int n){
        auto vals = readnline(lines[i], n);
        i++;
        return vals;
    };
    auto read1str = [&lines, &i](){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(lines[i]), " ");
        i++;
        if (vals.empty()){throw std::invalid_argument("Unable to read one string from this line:"+lines[i]); }
        return vals[0];
    };
    auto read1 = [&lines, &i](){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(lines[i]), " ");
        i++;
        if (vals.empty()){throw std::invalid_argument("Unable to read one number from this line:"+lines[i]); }
        return strtod(vals[0].c_str(), nullptr);
    };
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
    auto termcounts = readn(-1);
    std::size_t Nnormal = termcounts[0];
    std::size_t Nnormalcount = termcounts[1]; // Parameters per term
    std::size_t NGaussian = termcounts[2];
    std::size_t NGaussiancount = termcounts[3]; // Parameters per term
    std::size_t NGao = termcounts[4];

    auto read_normal = [&readnline](const vector<string> &lines, size_t Nnormalcount) -> nlohmann::json {
        std::vector<double> n, t, d, l, g;
        for (auto line : lines){
            auto z = readnline(line, Nnormalcount);
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
    };

    auto read_Gaussian = [&readnline](const vector<string> &lines, size_t NGaussiancount) -> nlohmann::json {
        // These terms can be of a number of different kinds
        // and they are disambiguated based upon the parameters in the term
        // (this is not documented anywhere in REFPROP)

        nlohmann::json normal, nonanalyt, R125;
        auto init_normal = [&normal](){
            std::vector<double> _;
            normal = {
                {"n", _}, {"t", _}, {"d", _}, 
                {"eta", _}, {"beta", _}, {"gamma", _}, {"epsilon", _},
                {"type", "ResidualHelmholtzGaussian"}
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
        
        for (auto& line : lines){
            auto z = readnline(line, NGaussiancount);
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
                normal["n"].push_back(z[0]);
                normal["t"].push_back(z[1]);
                normal["d"].push_back(z[2]);
                normal["l"].push_back(z[3]);
                normal["m"].push_back(z[4]);
            }
            else if (z[3] == 2 && z[4] == 2 && z[9] != 0 && z[10] != 0 && z[11] != 0){
                // It is a non-analytic term
                if (nonanalyt.empty()){ init_nonanalyt(); }
                nonanalyt["n"].push_back(z[0]);
                nonanalyt["b"].push_back(z[5]);
                nonanalyt["beta"].push_back(z[6]);
                nonanalyt["A"].push_back(z[7]);
                nonanalyt["B"].push_back(z[8]);
                nonanalyt["C"].push_back(z[9]);
                nonanalyt["D"].push_back(z[10]);
                nonanalyt["a"].push_back(z[11]);
            }
            else{
                throw std::invalid_argument("Bad Gaussian term");
            }
        }
        nlohmann::json o = nlohmann::json::array();
        if(!normal.empty()){ o.push_back(normal); }
        if(!nonanalyt.empty()){ o.push_back(nonanalyt); }
        return o;
    };
    auto read_Gao = [&readnline](const vector<string> &lines, size_t NGaocount) -> nlohmann::json{
        std::vector<double> n,t,d,eta,beta,gamma,epsilon,b;
        for (auto& line : lines){
            auto z = readnline(line, NGaocount);
            assert(NGaussiancount == 12);
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
 * Things that are obtained from parsing a residual Helmholtz
 * energy EOS
 */
struct Alpha0Result{
    double Tmin_K,Tmax_K,pmax_kPa,rhomax_molL;
    std::string cp0_pointer;
    double Tred_K, cp0red_JmolK, R;
    nlohmann::json alpha0;
};

Alpha0Result convert_CP0(const vector<string>& lines, double Tri){
    Alpha0Result a;
    
    auto read1strline = [](const string &line){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), " ");
        if (vals.empty()){throw std::invalid_argument("Unable to read one string from this line:"+line); }
        return vals[0];
    };
    a.cp0_pointer = read1strline(lines[1]);
    
    // Find the first non-header row;
    size_t i = std::string::npos;
    for (auto j = 2; j < lines.size(); ++j){
        // First line not started by element in {:,?,!}, stop
        auto ind = lines[j].find_first_of(":?!");
        if(!lines[j].empty() && (ind == 0)){
            continue;
        }
        i = j; break;
    }
    auto readnline = [](const string& line, int n){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), " ");
        std::vector<double> o;
        for (auto v : vals){
            if (!v.empty()){
                std::string v_e = std::regex_replace(v, std::regex("[Dd]"), "e");
                std::string v_ew = std::regex_replace(v_e, std::regex("\\s+"), "");
                if (v_ew.empty()){
                    continue;
                }
                o.push_back(strtod(v_ew.c_str(), nullptr));
            }
        }
        if (n > 0){
            if (o.size() != n){ throw std::invalid_argument("Unable to read "+std::to_string(n)+" numbers from this line:"+line); }
        }
        return o;
    };
    auto readn = [&lines, &i, &readnline](int n){
        auto vals = readnline(lines[i], n);
        i++;
        return vals;
    };
    auto read1str = [&lines, &i](){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(lines[i]), " ");
        i++;
        if (vals.empty()){throw std::invalid_argument("Unable to read one string from this line:"+lines[i]); }
        return vals[0];
    };
    auto read1 = [&lines, &i](){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(lines[i]), " ");
        i++;
        if (vals.empty()){throw std::invalid_argument("Unable to read one number from this line:"+lines[i]); }
        return strtod(vals[0].c_str(), nullptr);
    };
    
    // Read in all the metadata parameters
    a.Tmin_K = read1();
    a.Tmax_K = read1();
    a.pmax_kPa = read1();
    a.rhomax_molL = read1();
    auto r = readn(2);
    a.Tred_K = r[0];
    a.cp0red_JmolK = r[1];
    
    auto N = readn(7);
    auto Npoly = N[0];
    auto NPlanck = N[1];

    // cp0/R = c_i*T^t_i
    auto read_polynomial = [&readnline, &a, &Tri](const vector<string> &lines) -> nlohmann::json {
        std::vector<double> c, t;
        for (auto line : lines){
            auto z = readnline(line, 2);
            c.push_back(z[0]);
            t.push_back(z[1]);
        }
        return {{"type", "IdealGasHelmholtzCP0PolyT"}, {"c", c}, {"t", t}, {"T0", 298.15}, {"Tc", Tri}, {"R", a.cp0red_JmolK}};
    };
    
    auto read_Planck = [&readnline, &a, &Tri](const vector<string> &lines) -> nlohmann::json {
        std::vector<double> n, v;
        for (auto line : lines){
            auto z = readnline(line, 2);
            n.push_back(z[0]);
            v.push_back(z[1]);
        }
        return {{"type", "IdealGasHelmholtzPlanckEinsteinFunctionT"}, {"n", n}, {"v", v}, {"T0", 298.15}, {"Tcrit", Tri}, {"R",a.cp0red_JmolK} };
    };

    a.alpha0 = nlohmann::json::array();
    // The log(rho) term that is always included
    a.alpha0.push_back({{"a1", 0.0}, {"a2", 0.0}, {"type", "IdealGasHelmholtzLead"}});
    a.alpha0.push_back({{"a", -1}, {"type", "IdealGasHelmholtzLogTau"}});
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
    double MM_kgkmol, Ttriple_K, Tnbp_K, Tcrit_K, pcrit_kPa, rhocrit_molL, acentric, dipole_D;
    std::string refstate, RPversion;
    std::string UNnumber = "?";
    std::string family = "?";
    double HVupper_kJmol = 99999999, GWP100=99999999, RCL_vv = 99999999, ODP=99999999;
    std::string ASHRAE34 = "?", StdInChIstr = "?", StdInChIKey = "?", altid = "?", hash = "?";
    auto to_int(const std::string& s) -> int {
        return strtol(s.c_str(), nullptr, 10);
    }
    auto to_double(const std::string& s) -> int {
        return strtod(s.c_str(), nullptr);
    }
    void set(const std::string &k, const std::string &val){
        if (k == "FAMILY"){
            family = val;
        }
        else if (k == "UN"){
            UNnumber = to_int(val);
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

HeaderResult convert_header(const vector<string>& lines){
    HeaderResult h;
    std::size_t i = 0;
    auto _read1str = [](const std::string &line){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), "  ");
        if (vals.empty()){throw std::invalid_argument("Unable to read one string from this line:"+line); }
        return vals[0];
    };
    auto read1str = [&lines, &i, &_read1str](){
        string val = _read1str(lines[i]); i++; return val;
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
        h.UNnumber = read1();
        h.family = read1str();
        if (i < lines.size()-1){
            h.HVupper_kJmol = read1();
        }
    }
    
    return h;
}

auto get_ancillary_description(const string& key){
    const std::map<std::string, std::string> keys = {
        // Most are in this form:
        {"PS5",  "P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T]"},
        {"DL1",  "D=Dc*[1+SUM(Ni*Theta^ti)]"},
        {"DV3",  "D=Dc*EXP[SUM(Ni*Theta^ti)]"},
        
        // Other ones might be of this form:
        {"DL2",  "D=Dc*[1+SUM(Ni*Theta^(ti/3))]"},
        {"DL3",  "D=Dc*EXP[SUM(Ni*Theta^ti)]"},
        {"DL4",  "D=Dc*EXP[SUM(Ni*Theta^(ti/3))]"},
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

nlohmann::json get_ancillary(const vector<string>& lines){
    
    auto modelname = lines[1].substr(0, 3);
    
    // Find the first non-header row;
    size_t i = std::string::npos;
    for (auto j = 2; j < lines.size(); ++j){
        // First line not started by element in {:,?,!}, stop
        auto ind = lines[j].find_first_of(":?!");
        if(!lines[j].empty() && (ind == 0)){
            continue;
        }
        i = j; break;
    }
    
    auto readnline = [](const string line, int n){
        using namespace internal;
        auto vals = strsplit(strip_line_comment(line), " ");
        std::vector<double> o;
        for (auto v : vals){
            if (!v.empty()){
                o.push_back(strtod(v.c_str(), nullptr));
            }
        }
        if (o.size() != n){ throw std::invalid_argument("Unable to read "+std::to_string(n)+" numbers from this line:"+line);  }
        return o;
    };
    auto readn = [&lines, &i, &readnline](int n){
        auto vals = readnline(lines[i], n);
        i++;
        return vals;
    };
    auto read1 = [&lines, &i, &readn](){ return readn(1)[0]; };
    
    double Tmin_K = read1();
    double Tmax = read1();
    double ph1 = read1();
    double ph2 = read1();
    auto reducing = readn(2);
    reducing[1] *= 1000; // REFPROP uses mol/L & kPa, CoolProp uses mol/m^3 & Pa, multiply by 1000 to convert
    auto Ncoeffs = readn(6);
    
    std::vector<double> n, t;
    for (auto i = 0; i < Ncoeffs[0]; ++i){
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
    }
    else if (model_key == "DL"){
        type_key = "rhoL";
        if (seti{2,4,6}.count(key_index) > 0){
            for (auto&t_ : t){ t_ /= 3; }
        }
    }
    else if (model_key == "DV"){
        type_key = "rhoV";
        if (seti{2,4,6}.count(key_index) > 0){
            for (auto&t_ : t){ t_ /= 3; }
        }
    }
    else{
        throw std::invalid_argument(model_key);
    }
    
    return {
        {"T_r", reducing[0]},
        {"Tmax", reducing[0]},
        {"Tmin", Tmin_K},
        {"description", desc},
        {"n", n},
        {"reducing_value", reducing[1]},
        {"t", t},
        {"type", type_key},
        {"using_tau_r", (desc.find("Tc/T") != std::string::npos)}
    };
}

nlohmann::json get_all_ancillaries(const vector<string>& lines){
    return {
        {"pS", get_ancillary(internal::get_line_chunk(lines, "#PS"))},
        {"rhoV", get_ancillary(internal::get_line_chunk(lines, "#DV"))},
        {"rhoL", get_ancillary(internal::get_line_chunk(lines, "#DL"))}
    };
}

class FLDfile{
private:
    const std::string contents;
    const std::vector<std::string> lines;
    std::string read_file(const std::filesystem::path &path, bool replace_tabs = true){
        
        std::ifstream ifs(path);
        if (!ifs){
            throw std::invalid_argument(path);
        }
        std::stringstream buffer;  buffer << ifs.rdbuf();
        
        std::string contents = buffer.str();
        
        if (contents.find("\t") != std::string::npos){
            if(replace_tabs){
                std::regex_replace(contents, std::regex("\t"), "    ");
            }
            else{
                throw std::invalid_argument("Tabs found in the file "+path.string()+". Replace them with spaces.");
            }
        }
        return buffer.str();
    }
    std::vector<std::string> warnings, errors;
public:
    FLDfile(const std::filesystem::path &path) : contents(read_file(path)), lines(internal::strsplit(contents, "\n")){ };
    
    auto get_warnings(){ return warnings; }
    auto get_errors(){ return errors; }

    nlohmann::json convert_EOS(const RPinterop::ResidualResult& feq, const RPinterop::Alpha0Result& alpha0){

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
              {"T", feq.Ttriple_K},
              {"T_units", "K"},
              {"hmolar", -999999999999.0},
              {"hmolar_units", "J/mol"},
              {"p", feq.ptriple_kPa*1e3},
              {"p_units", "Pa"},
              {"rhomolar", 9999999999999999.0},
              {"rhomolar_units", "mol/m^3"},
              {"smolar", 99999999999999999999.0},
              {"smolar_units", "J/mol/K"}
            }},
            {"sat_min_vapor", {
              {"T", feq.Ttriple_K},
              {"T_units", "K"},
              {"hmolar", 9999999999999.0},
              {"hmolar_units", "J/mol"},
              {"p", feq.ptriple_kPa*1e3},
              {"p_units", "Pa"},
              {"rhomolar", 99999999999999.0},
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
            {"SMILES", "?"}
        };
    }

    auto make_json(const string& name){
        auto head = convert_header(lines);
        auto feq = convert_FEQ(internal::get_line_chunk(lines, "#EOS"));
        auto alpha0 = convert_CP0(internal::get_line_chunk(lines, "#AUX"), feq.Tred_K);
        auto EOS = convert_EOS(feq, alpha0);
        
        nlohmann::json ancillaries = nlohmann::json::object();
        try{
            ancillaries = get_all_ancillaries(lines);
        }
        catch(std::exception &e){
            warnings.push_back("Could not load ancillaries");
//            std::cerr << e.what() << std::endl;
        }

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

}
