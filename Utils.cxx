#include "Utils.h"

#include <sys/stat.h>

std::vector<std::string> Utils::tokenize(std::string str, std::string delim){

  std::vector<std::string> tokens;
  std::string::size_type sPos, sEnd, sLen;
  sPos = 0;
  while ( sPos != std::string::npos ) {
    sEnd = str.find_first_of(delim, sPos);
    if(sEnd==std::string::npos) sEnd = str.length();
    sLen = sEnd - sPos;
    std::string word = str.substr(sPos,sLen);
    tokens.push_back(word);
    sPos = str.find_first_not_of(delim, sEnd);
  }
  return tokens;
}

bool Utils::pathExists(std::string pathName) {

  struct stat fileAtt; 

  return stat(pathName.c_str(), &fileAtt) == 0; 
}

bool Utils::isFile(std::string pathName) {

  struct stat fileAtt; 

  int exists = stat(pathName.c_str(), &fileAtt);

  if(exists != 0) return false;

  return S_ISREG(fileAtt.st_mode);
}

bool Utils::isDir(std::string pathName) {

  struct stat fileAtt; 

  int exists = stat(pathName.c_str(), &fileAtt);

  if(exists != 0) return false;

  return S_ISDIR(fileAtt.st_mode);
}
