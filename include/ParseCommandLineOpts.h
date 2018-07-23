#ifndef INCLUDE_PARSECOMMANDLINE_OPTS_H_
#define INCLUDE_PARSECOMMANDLINE_OPTS_H_

#include <string>
#include <iterator>
#include <vector>

class ParseCommandLineOpts{
    public:
        ParseCommandLineOpts (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                tokens.push_back(std::string(argv[i]));
        }
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(tokens.begin(), tokens.end(), option);
            if (itr != tokens.end() && ++itr != tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }
        bool cmdOptionExists(const std::string &option) const{
            return std::find(tokens.begin(), tokens.end(), option)
                   != tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};

#endif
