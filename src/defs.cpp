#include "defs.hpp"
#include "scream_atom.hpp"

void deleteScreamAtom(SCREAM_ATOM* atom) {
  /* do i need to do type checking? */
  delete atom;
}

/* below: from www.roguewave.com/support/docs/sourcepro/stdlibug/12-3.html */
void split (const std::string& text, const std::string& separators,
            std::vector<std::string,
	    std::allocator<std::string> >& words) {

    size_t n     = text.length ();
    size_t start = text.find_first_not_of (separators);

    while (start < n) {
        size_t stop = text.find_first_of (separators, start);
        if (stop > n) stop = n;
        words.push_back (text.substr (start, stop-start));
        start = text.find_first_not_of (separators, stop+1);
    }
}


/**
 * C++ version std::string style "itoa":
 */
	
std::string itoa(int value, unsigned int base) {
	
	const char digitMap[] = "0123456789abcdef";
	std::string buf;
	// Guard:
	if (base == 0 || base > 16) {
		// Error: may add more trace/log output here
		return buf;
	}
	// Take care negative int:
	std::string sign;
	int _value = value;
	if (value < 0) {
		_value = -value;
		sign = "-";
	}
	// Translating number to string with base:
	for (int i = 30; _value && i ; --i) {
		buf = digitMap[ _value % base ] + buf;
		_value /= base;
	}
	return sign.append(buf);
}


/* two helper functions */

size_t get_first_of_either(string text, vector<string> white_spaces, size_t pos) {

  size_t n     = text.length ();
  size_t start_1 = text.find_first_of (white_spaces[0], pos);
  size_t start_2 = text.find_first_of (white_spaces[1], pos);
  size_t start = (start_1 < start_2) ? start_1 : start_2;

  return start;

}

size_t get_first_not_of_either(string text, vector<string> white_spaces, size_t pos) {
  size_t n     = text.length ();
  size_t start_1 = text.find_first_not_of (white_spaces[0], pos);
  size_t start_2 = text.find_first_not_of (white_spaces[1], pos);
  size_t start = (start_1 < start_2) ? start_1 : start_2;
  size_t start_tmp = start;

  for (int i = 0; i < n; i++) {
    start_tmp = text.find_first_not_of (white_spaces[0], start_tmp);
    start_tmp = text.find_first_not_of (white_spaces[1], start_tmp);
  }

  start = start_tmp;
  return start;

}

/* end of the two helper functions */

void split(const std::string& text, 
	   std::vector<std::string, std::allocator<std::string> >& words) {

  /* This splits according to white spaces, including tabs. */
  size_t n     = text.length ();
  std::vector<std::string> white_spaces;
  white_spaces.push_back(" ");
  white_spaces.push_back("\t");
  
  size_t start = get_first_not_of_either(text, white_spaces, 0);  // start_1 == start_2 at this point.

  while (start < n) {
    //    size_t stop = text.find_first_of (separators, start);
    size_t stop = get_first_of_either(text, white_spaces, start);
    if (stop > n) stop = n;
    words.push_back (text.substr (start, stop-start));
    //    start = text.find_first_not_of (separators, stop+1);
    start = get_first_not_of_either(text, white_spaces, stop+1);
  }
  
  
}

AtomPair::AtomPair() {

}

AtomPair::AtomPair(SCREAM_ATOM* a1, SCREAM_ATOM* a2) {

  this->a1 = a1;
  this->a2 = a2;

}

AtomPair::~AtomPair() {

}

AtomPair& AtomPair::operator=(const AtomPair& ap) {
  if (this == &ap)
    return *this;
  this->a1 = ap.a1;
  this->a2 = ap.a2;
  return *this;
}

void print(string str) {
  cout << str << endl;
}

void dump(SCREAM_ATOM* atom) {
  atom->dump();
}


/* DEBUG class stuff. */

#include <string.h>

#ifdef DEBUG

Debug::Debug(string category) {
  this->_category = new char[400];
  strcpy(this->_category , category.c_str());
  this->_active = true;
}

Debug::~Debug() {
  delete this->_category;
}

void Debug::out(string line, string outFilename) const {
  string outstream(outFilename);
  ostream* OUTSTREAM = new ofstream;
  if (outstream == "stdout" or outstream == "STDOUT") {
    OUTSTREAM = &cout;
  } else {
    ((ofstream*)OUTSTREAM)->open(outFilename.c_str());
  }
  *OUTSTREAM << string(this->_category) << "   " << line << endl;

  if (outstream != "stdout") {
    ((ofstream*)OUTSTREAM)->close();
    delete OUTSTREAM;
  }
}

void Debug::trap() {
  // do nothing for now
}

#endif /* DEBUG */
