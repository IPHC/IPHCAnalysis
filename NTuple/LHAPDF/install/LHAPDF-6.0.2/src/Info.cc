// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/Info.h"
#include "LHAPDF/PDFIndex.h"
#include "yaml-cpp/yaml.h"

namespace LHAPDF {


  void Info::load(const path& filepath) {
    // Silently do nothing if the provided path is empty
    if (filepath.empty()) return;

    // But complain if a non-empty path is provided, but it's invalid
    if (!exists(filepath)) throw ReadError("PDF data file '" + filepath.string() + "' not found");

    // Read the YAML part of the file into the metadata map
    try {
      // Do the parsing "manually" up to the first doc delimiter
      std::ifstream file(filepath.c_str());
      YAML::Node doc;

      #if YAMLCPP_API == 3

      YAML::Parser parser(file);
      parser.GetNextDocument(doc);
      for (YAML::Iterator it = doc.begin(); it != doc.end(); ++it) {
        string key, val;
        it.first() >> key;
        try {
          // Assume the value is a scalar type -- it'll throw an exception if not
          it.second() >> val;
        } catch (const YAML::InvalidScalar& ex) {
          // It's a list: process the entries individually into a comma-separated string
          string subval;
          for (size_t i = 0; i < it.second().size(); ++i) {
            it.second()[i] >> subval;
            val += subval + ((i < it.second().size()-1) ? "," : "");
          }
        }
        //cout << key << ": " << val << endl;
        _metadict[key] = val;
      }

      #elif YAMLCPP_API == 5

      string docstr, line;
      while (getline(file, line)) {
        //cout << "@ " << line << endl;
        if (line == "---") break;
        docstr += line + "\n";
      }
      doc = YAML::Load(docstr);
      for (YAML::const_iterator it = doc.begin(); it != doc.end(); ++it) {
        const string key = it->first.as<string>();
        const YAML::Node& val = it->second;
        if (val.IsScalar()) {
          // Scalar value
          _metadict[key] = val.as<string>();
        } else {
          // Process the sequence entries into a comma-separated string
          string seqstr = "";
          for (size_t i = 0; i < val.size(); ++i)
            seqstr += val[i].as<string>() + ((i < val.size()-1) ? "," : "");
          _metadict[key] = seqstr;
        }
      }

      #endif
    } catch (const YAML::ParserException& ex) {
      throw ReadError("YAML parse error in " + filepath.string() + " :" + ex.what());
    } catch (const LHAPDF::Exception& ex) {
      throw;
    } catch (const std::exception& ex) {
      throw ReadError("Trouble when reading " + filepath.string() + " :" + ex.what());
    }

  }



  // /// @todo Only support loading via PDF set name and member ID, not explicit paths
  // /// @todo Replace the loading of the set metadata into the member info with set-level Info singletons
  // void Info::loadFull(const path& mempath) { //< @todo Need a better method name!
  //   // Extract the set name from the member data file path
  //   const path memberdata = findFile(mempath);
  //   if (memberdata.empty() || !exists(memberdata)) throw ReadError("Could not find PDF data file '" + mempath.string() + "'");
  //   const string memname = memberdata.filename().string(); //< Can use this to alternatively work out the set name...
  //   const path setdir = memberdata.parent_path();
  //   const string setname = setdir.filename().string();
  //   path setinfo = findpdfsetinfopath(setname);
  //   // Load the set info
  //   if (exists(setinfo)) load(setinfo.string());
  //   // Load the member info (possibly overriding the set-level metadata)
  //   load(memberdata.string());
  // }


}
