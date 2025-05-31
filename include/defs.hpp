#ifndef DEFS_HPP
#define DEFS_HPP

#include "scream_atom.hpp"
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

typedef vector<SCREAM_ATOM*> ScreamAtomV;
typedef vector<SCREAM_ATOM*>::const_iterator ScreamAtomVConstItr;
typedef vector<SCREAM_ATOM*>::iterator ScreamAtomVItr;

typedef vector<string> stringV;
typedef vector<string>::const_iterator stringVConstItr;
typedef vector<string>::iterator stringVItr;

/* Functions for easy deletion of SCREAM_ATOM*'s in a ScreamAtomV. */

void deleteScreamAtom(SCREAM_ATOM*);

/* String Manipulation functions. */

void split (const std::string& text, const std::string& separators,
            std::vector<std::string, std::allocator<std::string> >& words);
void split (const std::string& text, 
	    std::vector<std::string, std::allocator<std::string> >& words); // splits to " " and "\t"

string itoa(int, unsigned int = 10);

/* Below: http://www.parashift.com/c++-faq-lite/misc-technical-issues.html#faq-39.2 */

class BadConversion : public std::runtime_error {
 public:
   BadConversion(const std::string& s)
     : std::runtime_error(s)
     { }
 };

inline std::string stringify(double x)
{
  std::ostringstream o;
  if (!(o << x))
    throw BadConversion("stringify(double)");
  return o.str();
} 

/* String to int or float conversion */
template <class T> 
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&)) 
{
  // the third parameter of from_string() should be 
  // one of std::hex, std::dec or std::oct
  std::istringstream iss(s);
  return !(iss>> f >> t).fail();

}


/* a Pair of atoms */


class AtomPair {
public:
  AtomPair();
  AtomPair(SCREAM_ATOM*, SCREAM_ATOM*);
  ~AtomPair();
  
  AtomPair& operator=(const AtomPair&);

  SCREAM_ATOM* a1;
  SCREAM_ATOM* a2;
};

void print(string str);
void dump(SCREAM_ATOM*);

/* DEBUG class.  Borrowed from: http://www.parasoft.com/jsp/products/article.jsp?articleId=489
 * changed const char * to string.
 */

#ifndef DEBUG

class Debug
{
public:
  Debug(const string) {}
  static bool active()               { return false; }
  static void out(string, string = "stdout") { return; }
};

#endif /* ifndef DEBUG */

#ifdef DEBUG

class Debug
{
public:
  Debug(string category);
  bool active() const { return _active; }
  //void out(const char *file, long line, const char * format, ...) const;
  void out(string line, string outFilename = "stdout") const;
  ~Debug();

public:
  static void trap();
  
private:
  char * _category; // Not a string, so no overhead until used
  bool _active;
  
  static bool _showAll;
};

#endif /* ifdef DEBUG */


#endif /* DEFS_HPP */
