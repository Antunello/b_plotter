#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <sstream>
#include <vector>
#include <iostream>

#include <TFile.h>
#include <TObject.h>

namespace Utils {

  std::vector<std::string> tokenize(std::string str, std::string delim);

  bool pathExists(std::string pathName);
  bool isFile(std::string pathName);
  bool isDir(std::string pathName);

  template <class T>
  T *readObject(TDirectory *d, std::string name);
  template <class T>
  T *readObject(std::string fname, std::string name);

}

#define UTILS__ERROR(messageExp)			      \
  do {							      \
    std::ostringstream oss;				      \
    oss << "Error in " << __FILE__ << ":" << __LINE__ << ": " \
	<< messageExp;					      \
    std::cout << oss.str().c_str() << std::endl;	      \
    throw oss.str();					      \
  }while(false)


template <class T>
T *Utils::readObject(TDirectory *d, std::string name) {
  TObject *obj = d->Get(name.c_str());
  if(!obj) {
    UTILS__ERROR("Object " << name << " not found in file " << d->GetName());
  }
  return (T*)obj->Clone();
}

template <class T>
T *Utils::readObject(std::string fname, std::string name) {
  TFile *f = TFile::Open(fname.c_str(), "READ");
  if(!f || f->IsZombie()) {
    UTILS__ERROR("Could not open file " << fname << " for reading");
  }
  T *obj = readObject<T>(f, name);
  f->Clear();
  f->Close();
  delete f;
  return obj;
}

#endif
