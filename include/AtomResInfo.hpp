#ifndef ATOMRESINFO_HPP
#define ATOMRESINFO_HPP

#include <string>
#include <iostream>
#include <map>

/* SCREAM AtomRes info: fields: Resname (3 letter), AtomLabel */

class AtomResInfo {
public:
  AtomResInfo();
  AtomResInfo(std::string, std::string);
  ~AtomResInfo();

  AtomResInfo& operator=(const AtomResInfo&);
  bool operator==(const AtomResInfo&) const;
  bool operator<(const AtomResInfo&) const;
  friend std::ostream& operator<<(std::ostream&, const AtomResInfo&);

  std::string resName; // 3 letter resname
  std::string atomLabel;

private:
  static std::map<std::string, double> atomLabel_value_map;
  //std::map<std::string, double> atomLabel_value_map;
};


#endif /* ATOMRESINFO_HPP */
