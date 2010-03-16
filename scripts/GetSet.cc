////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Mon Mar  1 09:08:37 CET 2010
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <string>
#include <fstream>
#include <vector>


int GetSet (std::string const Class, std::string const FileName)
{
  std::ifstream InFile(FileName.c_str());

  std::string Type;
  std::string Name;

  std::vector<std::string> FuncDef;
  char Def[500];
  while ( InFile >> Type >> Name ) {
    Name = std::string(Name.begin(), Name.end()-1);

    sprintf(Def, "%s Get%s ();", Type.c_str(), Name.c_str());
    FuncDef.push_back(Def);
    printf("%s %s::Get%s ()\n", Type.c_str(), Class.c_str(), Name.c_str());
    printf("{\n");
    printf("  return %s;\n", Name.c_str());
    printf("}\n\n\n");

    sprintf(Def, "void Set%s (%s const);", Name.c_str(), Type.c_str());
    FuncDef.push_back(Def);
    printf("void %s::Set%s (%s const in)\n", Class.c_str(), Name.c_str(), Type.c_str());
    printf("{\n");
    printf("  %s = in;\n", Name.c_str());
    printf("  return;\n");
    printf("}\n\n\n");
  }

  for (size_t i = 0; i != FuncDef.size(); ++i) {
    std::cout << FuncDef[i] << std::endl;
  }
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [ClassName] [FileName]" << std::endl;
    return 1;
  }

  GetSet(argv[1], argv[2]);

  return 0;
}
