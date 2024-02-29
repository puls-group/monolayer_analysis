/* This is a small tool for easy command line option handling.
 * For a quick start, copy the following code into your
 * 	int main (int argc, char *argv[]) :

op::OptionHandler OH ("This is the helptext");
int istore = 0;
OH.addOption((new op::SingleValueOption<int>("i", istore))->description("Integer
option")); op::pRes res = OH.procOptions(argc, argv); if (res != op::ok) return
res;

 * This registers an option -i, which takes an integer argument and stores it in
the variable istore.
 * Furthermore, the options -h[elp] and -? are always available; they display
the help text.
 */
#ifndef __OPTION_HANDLER_H__
#define __OPTION_HANDLER_H__

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace op {

typedef enum { ok = 1, ret = 0, ret_err = -1 } pRes;
typedef enum { opt_ok = 1, opt_error = -1 } optRes;

class PositionalOption {
  public:
    PositionalOption(const std::string& name, bool required = true,
                     int count = 0)
        : argc_(count), name_(name), desc_(""), required_(required) {}

    virtual ~PositionalOption() {}

    std::string name() const { return name_; }

    int argc() const { return argc_; }

    bool is_required() const { return required_; }
    virtual bool is_satisfied() const {
        if (!is_required()) {
            return true;
        }
        return false;
    }
    virtual std::string get_default_text() const = 0;

    std::string helpline() const {
        std::string res = std::string("");

        if (!is_required()) {
            res.append("[");
        }

        res.append(name_);

        if (!is_required()) {
            res.append("]");
        }
        res.append("\t");
        if (desc_.length() > 0) {
            res.append(desc_);
        } else {
            res.append("No known description");
        }
        if (!is_required()) {
            res.append(get_default_text());
        }
        res.append("\n\n");
        return res;
    }

    virtual PositionalOption* description(std::string desc) {
        desc_ = desc;
        return this;
    }

    virtual pRes proc(char* argv[]) = 0;

  protected:
    int argc_;

  private:
    const std::string name_;
    std::string desc_;
    bool required_;
}; // class Option

class Option {
  public:
    Option(const std::string& key, int count = 0)
        : argc_(count), name_(key), desc_("") {}

    virtual ~Option() {}

    std::string name() const { return name_; }

    int argc() const { return argc_; }

    virtual std::string get_default_text() const = 0;

    std::string helpline() const {
        if (desc_.length() > 0) {
            std::string res = std::string("  -").append(name_).append("\t");
            res.append(desc_);
            res.append(get_default_text());
            res.append("\n\n");
            return res;
        }
        return std::string("");
    }

    virtual Option* description(std::string desc) {
        desc_ = desc;
        return this;
    }

    virtual pRes proc(char* argv[]) = 0;

  protected:
    int argc_;

  private:
    const std::string name_;
    std::string desc_;
}; // class Option

class OptionHandler {
  public:
    OptionHandler(std::string helptext) : desc_(helptext) {}
    ~OptionHandler() {
        /*for (auto& x : opMap_)
            delete x.second;*/
    }
    optRes addOption(Option& opt) {
        opMap_[opt.name()] = &opt;
        return optRes::opt_ok;
    }
    optRes addOption(PositionalOption& opt) {
        if (positional_required && !opt.is_required()) {
            positional_required = false;
        } else if (!positional_required && opt.is_required()) {
            std::cerr << "Error: cannot have required positional argument "
                         "after optional one"
                      << std::endl;
            std::cerr << "Optional: " << pos_opt_[pos_opt_.size() - 1]->name()
                      << "\t"
                      << "Required" << opt.name() << std::endl;
            return optRes::opt_error;
        }
        pos_opt_.push_back(&opt);
        return optRes::opt_ok;
    }

    pRes procOptions(int argc, char* argv[]);
    std::string getHelp() {
        std::string ret = "\n";
        ret += desc_;
        ret += "\n\nAvailable options:\n\n";
        for (auto& x : pos_opt_)
            ret += x->helpline();
        for (auto& x : opMap_)
            ret += x.second->helpline();
        ret += "  -h[elp]\n  -?\tDisplays this help\n\n";
        return ret;
    }

  private:
    std::map<std::string, Option*> opMap_;
    std::vector<PositionalOption*> pos_opt_;
    const std::string desc_;
    bool positional_required = true;
    unsigned int curr_positional = 0;
}; // class OptionHandler

template <class value_type> class SingleValueOption : public Option {
  public:
    SingleValueOption(const std::string& key, value_type& storage)
        : Option(key, 1), store_(storage) {}

    virtual ~SingleValueOption() {}

    virtual pRes proc(char* argv[]);

    virtual std::string get_default_text() const;

  private:
    value_type& store_;
}; // class SingleValueOption

template <class value_type>
pRes SingleValueOption<value_type>::proc(char* argv[]) {
    std::istringstream oss(argv[1]);
    oss >> store_;
    return ok;
}

template <class value_type>
std::string SingleValueOption<value_type>::get_default_text() const {
    return std::string("(Default: ").append(std::to_string(store_)).append(")");
}

template <>
std::string SingleValueOption<std::string>::get_default_text() const {
    return std::string("(Default: ").append(store_).append(")");
}

template <class value_type> class AppendableValueOption : public Option {
  public:
    AppendableValueOption(const std::string& key,
                          std::vector<value_type>& storage)
        : Option(key, 1), store_(storage) {}

    virtual ~AppendableValueOption() {}

    virtual pRes proc(char* argv[]);

    virtual std::string get_default_text() const;

  private:
    std::vector<value_type>& store_;
}; // class SingleValueOption

template <class value_type>
pRes AppendableValueOption<value_type>::proc(char* argv[]) {
    value_type dummy;
    std::istringstream oss(argv[1]);
    oss >> dummy;
    store_.push_back(dummy);
    return ok;
}

template <class value_type>
std::string AppendableValueOption<value_type>::get_default_text() const {
    return "";
}

template <>
std::string AppendableValueOption<std::string>::get_default_text() const {
    return "";
}

template <class value_type>
class SingleValuePositionalOption : public PositionalOption {
  public:
    SingleValuePositionalOption(const std::string& name, value_type& storage,
                                bool required = true)
        : PositionalOption(name, 1, required), store_(storage),
          satisfied_(false) {}

    virtual ~SingleValuePositionalOption() {}

    virtual pRes proc(char* argv[]);
    virtual bool is_satisfied() const;

    virtual std::string get_default_text() const;

  private:
    value_type& store_;
    bool satisfied_;
}; // class SingleValuePositionalOption

template <class value_type>
pRes SingleValuePositionalOption<value_type>::proc(char* argv[]) {
    std::istringstream oss(argv[0]);
    oss >> store_;
    satisfied_ = true;
    return ok;
}

template <class value_type>
bool SingleValuePositionalOption<value_type>::is_satisfied() const {
    return satisfied_;
}

template <class value_type>
std::string SingleValuePositionalOption<value_type>::get_default_text() const {
    return std::string("(Default: ").append(std::to_string(store_)).append(")");
}

template <>
std::string SingleValuePositionalOption<std::string>::get_default_text() const {
    return std::string("(Default: ").append(store_).append(")");
}

pRes OptionHandler::procOptions(int argc, char* argv[]) {
    bool error = false;
    int errorAt = 0;
    for (int i = 1; i < argc; i++) {
        if (!(argv[i][0] == '-')) {
            // Try treating as positional argument
            if (curr_positional >= pos_opt_.size()) {
                // No positionals left
                error = true;
                errorAt = i;
                break;
            } else {
                // Feed to next positional
                PositionalOption* posopt = pos_opt_[curr_positional];
                pRes posres = posopt->proc((char**)(argv + i));
                if (posres != ok)
                    return posres;

                if (posopt->is_satisfied()) {
                    curr_positional++;
                }
                continue;
            }
        }

        char* curName = argv[i] + 1;
        if (std::string(curName) == "h" || std::string(curName) == "help" ||
            std::string(curName) == "?") {
            std::cout << getHelp();
            return ret;
        }

        try {
            Option* opt = opMap_.at(argv[i] + 1);
            if (i + opt->argc() >= argc) {
                std::cout << "Too few arguments specified for option "
                          << argv[i] << "; expected " << opt->argc()
                          << std::endl;
                return ret_err;
            }
            pRes res = opt->proc((char**)(argv + i));
            if (res != ok)
                return res;

            i += opt->argc();
        } catch (const std::out_of_range& oor) {
            error = true;
            errorAt = i;
            break;
        }
    }

    if (error) {
        std::cout << "Invalid option: " << argv[errorAt] << std::endl;
        return ret_err;
    } else if (curr_positional < pos_opt_.size()) {
        if (pos_opt_[curr_positional]->is_required()) {
            std::cout << "Missing positional option: "
                      << pos_opt_[curr_positional]->name() << std::endl;
            return ret_err;
        }
    }

    return ok;
} // OptionHandler::procOptions

} // namespace op

#endif