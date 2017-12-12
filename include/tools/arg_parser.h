#ifndef _ARG_PARSER_H_
#define _ARG_PARSER_H_

#include <cassert>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

/// Parse command-line arguments
/**
 * A simple command-line parser.
 *
 * Supports named arguments and switches as well as unnamed data arguments
 *
 * Example: ./foo -v -o outfolder in1.xml in2.xml
 */
class ArgParser {
 public:
  /// Parse command-line arguments
  ArgParser(int argc, char **argv) : named_args_(), data_args_() {
    int pos = 1;
    while (pos < argc) {
      if (argv[pos][0] == '-') {
        // Read argument and advance
        std::string arg(&argv[pos++][1]);
        // check whether argument has a value (i.e. not just a flag)
        if (pos < argc && argv[pos][0] != '-') {
          named_args_[arg] =
              std::string{argv[pos++]};  // assign and advance to next
        } else {
          named_args_[arg] = "";
        }
      } else {
        data_args_.emplace_back(argv[pos++]);
      }
    }
  }

  /// Get a named argument's value
  /// \param key the argument name
  /// \param default_value the value to return if the argument wasn't set
  template <typename T>
  T Get(const std::string &key, const T default_value = T()) const {
    T retval;
    auto it = named_args_.find(key);
    if (it != named_args_.end()) {
      std::istringstream s(it->second);
      s >> retval;
    } else {
      // do this in the else case, otherwise empty string arguments
      // would return the default value instead of ""
      retval = default_value;
    }
    return retval;
  }

  /// check whether an argument was set
  bool IsSet(const std::string &arg) const {
    return named_args_.find(arg) != named_args_.end();
  }

  /// the number of unnamed data arguments
  size_t NumDataArgs() const { return data_args_.size(); }

  /// get a data argument by its index (among the data arguments)
  std::string DataArg(const size_t index) const {
    assert(index < NumDataArgs());
    return data_args_[index];
  }

 protected:
  std::unordered_map<std::string, std::string> named_args_;
  std::vector<std::string> data_args_;
};

#endif
