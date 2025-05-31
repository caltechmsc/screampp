%module py_scream_ee

%{

#include "scream_atom.hpp"

#include "scream_model.hpp"

#include "scream_ctl_reader.hpp"
#include "Rotamer.hpp"
#include "AARotamer.hpp"
#include "sc_Protein.hpp"
#include "Rotlib.hpp"
#include "RotlibCollection.hpp"
#include "ClashCollection.hpp"
#include "RotamerNeighborList.hpp"
#include "RotamerCluster.hpp"
#include "RotConnInfo.hpp"

#include "sc_bgf_handler.hpp"
#include "scream_rtf.hpp"
#include "scream_EE.hpp"
#include "MutInfo.hpp"
%}

// Parse the original header file

%include "stl.i"
%include std_vector.i
%include std_string.i

class SCREAM_ATOM;

%template() std::vector<SCREAM_ATOM*>;
typedef std::vector<SCREAM_ATOM*> ScreamAtomV;

%include "scream_model.hpp"

%include "scream_ctl_reader.hpp"
%include "scream_atom.hpp"
%include "Rotamer.hpp"
%include "AARotamer.hpp"
%include "sc_Protein.hpp"
%include "Rotlib.hpp"
%include "RotlibCollection.hpp"
%include "ClashCollection.hpp"
%include "RotamerNeighborList.hpp"
%include "RotamerCluster.hpp"
%include "RotConnInfo.hpp"

%include "sc_bgf_handler.hpp"
%include "scream_rtf.hpp"
%include "scream_EE.hpp"
%include "MutInfo.hpp"





// Instantiate some templates
%template() std::vector<int>; 

%template() std::pair<SCREAM_ATOM*, int>; // for scream_atom connectivity_m
%template(ConnectivityMap) std::map<SCREAM_ATOM*, int>;

//%apply std::string& {std::string *};

//%template(ScreamAtomV) std::vector< SCREAM_ATOM* >;

%template(stringV) std::vector<std::string>;

//%template() std::pair<std::string, int>;
//%template(ExcitationEnumeration) std::map<std::string, int>;
%template() std::pair<std::string, unsigned short>;
%template(ExcitationEnumeration) std::map<std::string, unsigned short>;
%template() std::pair<std::string, Rotamer*>;
%template(ExcitedRotamers) std::map<std::string, Rotamer*>;

%template(RotamerV) std::vector<Rotamer*>;

%template(pairds) std::pair<double, std::string>;

%template(MutInfoListPy) std::vector<MutInfo>;
%template(MutInfoPairListPy) std::vector<MutInfoPair>;



%inline %{
    string derefString(string* str) {
	return *str;	
    }
%}


%inline %{
    Rotamer* derefRotamer(Rotamer** x) {
       return *x;
    }
%}

%inline %{
    AARotamer* derefAARotamer(AARotamer** x) {
	return *x;
    }

%}

%inline %{
    AARotamer *castRotamerToAARotamer(Rotamer* x) {
	return (AARotamer*)x;
    }

%}


