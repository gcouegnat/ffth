#include "meta.h"

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>

class Option {
  public:
    Option(const std::string& longopt, std::string shortopt, const std::string type,
           const std::string desc, bool required, const std::string& def)
        : longopt_(longopt), shortopt_(shortopt), desc_(desc), required_(required), value_(def),
          set_(false), type_(type) 
    {}

    bool is_required() {
        return required_;
    }
    bool is_set() {
        return set_;
    }
    
    std::string value() {
        return value_;
    }
    
    std::string description() {
        return desc_;
    }
    
    std::string name() {
        return longopt_;
    }
    
    std::string shortname() {
        return shortopt_;
    }
    
    void set_value(std::string value) {
        value_ = value;
        set_ = true;
    }

    std::string type() {
        return type_;
    }

  private:
    std::string longopt_;
    std::string shortopt_;
    std::string desc_;
    bool required_;
    std::string value_;
    bool set_;
    std::string type_;

};

class Argparser {

  public:
    Argparser() {};

    template<class T>
    void add_option(const std::string& name, const std::string& desc) {
        add_option(name, "", desc, true, T());
    }

    template<class T>
    void add_option(const std::string& name, const std::string& desc, const T value) {
        add_option(name, "", desc, false, value);
    }

    void parse(int argc, const char** argv) {
        progname_ = argv[0];
        for(int i=1; i<argc; ++i) {
            std::string argn = argv[i];
            std::string argname = argn.substr(2);
            if (argname == "help") {
                print_usage();
                std::exit(0);
            }
            
            opt_t::const_iterator it = options_.find(argname);
            if (it != options_.end()) {
                it->second->set_value( argv[++i] );
            } else {
                error(argn + " is not a valid option");
            }
        }
        check();
    }
    
    template<class T>
    T get(const std::string& name) {
        T value;
        opt_t::const_iterator it = options_.find(name);
        value = meta::lexical_cast<T>(it->second->value());
        return value;
    }

    void print_options() {
        for(opt_t::iterator it = options_.begin(); it != options_.end(); ++it) {
            Option* opt = it->second;
            std::cout << "--" << it->first << ", -" << opt->shortname() << std::endl << "       " << opt->description();
            if (opt->is_required()) {
                std::cout << " [required]";
            } else {
                std::cout << " [optional]" << opt->value() << "]";
            }
            std::cout << std::endl;
        }
    }

    void print_values() {
        for(std::vector<Option*>::iterator it = ordered_.begin(); it != ordered_.end(); ++it) {
            std::cout << (*it)->name() << ": " << (*it)->value() << std::endl;
        }
    }

    void print_usage() {
        std::cout << "\nUSAGE:  " << progname_ ;
        for(opt_t::iterator it = options_.begin(); it != options_.end(); ++it) {
            Option* opt = it->second;
            if (opt->is_required()) {
                std::cout << " --" << opt->name() << " <" << opt->type() << ">";
            }
        }
        std::cout << " [options]" << std::endl;
        std::cout << "\nOPTIONS:" << std::endl;
        for(std::vector<Option*>::iterator it = ordered_.begin(); it != ordered_.end(); ++it) {
            Option* opt = *it;
            std::cout << "  --" << opt->name() << " <" << opt->type() << ">     " << opt->description();
            if (!opt->is_required()) {
                std::cout << " (default value: " << opt->value() << ")";
            }
            std::cout << std::endl;
        }
    }

  private:
    void error(std::string msg) {
        std::cerr << "ERROR: ";
        std::cerr << msg << std::endl;
        print_usage();
        std::exit(1);
    }

    void check() {
        for(opt_t::iterator it = options_.begin(); it != options_.end(); ++it) {
            Option* opt = it->second;
            if (opt->is_required() and !opt->is_set()) {
                error("Mandatory argument --" + opt->name() + " is not set.");
            }
        }
    }

    template<class T>
    void add_option(const std::string& longopt, std::string shortopt, const std::string& desc, bool required, const T value) {
        std::string def = meta::lexical_cast<std::string>(value);
        std::string type = meta::readable_typename<T>();
        Option* opt = new Option(longopt, shortopt, type, desc, required, def);
        options_[longopt] = opt;
        ordered_.push_back(opt);
    }

  private:
    typedef std::map<std::string, Option*> opt_t;
    opt_t options_;
    std::vector<Option*> ordered_;
    std::string progname_;

};
