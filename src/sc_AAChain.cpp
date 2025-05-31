#include "defs.hpp"
#include "scream_atom.hpp"
#include "sc_ProteinComponent.hpp"
#include "sc_AminoAcid.hpp"
#include "scream_tools.hpp"
//#include "GenericRotamer_defs.hpp"

#include <cassert>
#include <map>
#include <iostream>
using namespace std;
#include "sc_AAChain.hpp"
//vcvicek
#include <cstdlib>

AAChain::AAChain() {

}

AAChain::AAChain(const ScreamAtomV& atom_v) {

  nterm_mm_on_free_store = false;
  cterm_mm_on_free_store = false;

  this->InitFromAtomList(atom_v);

}

AAChain::AAChain(const AAChain& ch) {

}


AAChain::~AAChain() {
  /* Cleaning up AminoAcid data structures. */
  map<int, AminoAcid*>::const_iterator itr;
  for (itr = aa_m.begin(); itr != aa_m.end(); ++itr) {
    delete itr->second;
  }
    
  multimap<string, SCREAM_ATOM*>::const_iterator m_itr;
  if (nterm_mm_on_free_store) {
    for (m_itr = nterm_mm.begin(); m_itr != nterm_mm.end(); ++m_itr) {
      delete m_itr->second;
    }
  }

  if (cterm_mm_on_free_store) {
    for (m_itr = cterm_mm.begin(); m_itr != cterm_mm.end(); ++m_itr) {
      delete m_itr->second;
    }
  }

}

ScreamAtomV AAChain::getAtomList() const {

  ScreamAtomV returnList;
  ScreamAtomV tempList;
  
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = nterm_mm.begin(); itr != nterm_mm.end(); ++itr) {
    tempList.push_back(itr->second);
  }
  tempList.clear();
  
  returnList.insert(returnList.begin(), tempList.begin(), tempList.end());
  
  /* AminoAcid traversal */
  for (map<int, AminoAcid*>::const_iterator itr = aa_m.begin(); itr != aa_m.end(); ++itr) {
    ScreamAtomV aaAtomList = itr->second->getAtomList();
    tempList.insert(tempList.end(), aaAtomList.begin(), aaAtomList.end());
  }
  returnList.insert(returnList.begin(), tempList.begin(), tempList.end());
  tempList.clear();

  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = cterm_mm.begin(); itr != cterm_mm.end(); ++itr) {
    tempList.push_back(itr->second);
  }
  
  returnList.insert(returnList.begin(), tempList.begin(), tempList.end());

  return returnList;

}


int AAChain::length() const {

  return this->aa_m.size();

}

int AAChain::number_of_atoms() const {
  
  int nterm_n = nterm_mm.size();
  int cterm_n = cterm_mm.size();
  int main_chain_atoms = 0;
  for (map<int, AminoAcid*>::const_iterator itr = aa_m.begin(); itr != aa_m.end(); ++itr) {
    main_chain_atoms += itr->second->number_of_atoms();
  }

  return (nterm_n + cterm_n + main_chain_atoms);

}

double AAChain::total_charge() const {

  
  double total_charge = 0;
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = nterm_mm.begin(); itr != nterm_mm.end(); ++itr) {
    total_charge += itr->second->q[0];
  }
  for (multimap<string, SCREAM_ATOM*>::const_iterator itr = cterm_mm.begin(); itr != cterm_mm.end(); ++itr) {
    total_charge += itr->second->q[0];
  }

  for (map<int, AminoAcid*>::const_iterator itr = aa_m.begin(); itr != aa_m.end(); ++itr) {
    total_charge += itr->second->total_charge();
  }
  
  return total_charge;

}

vector<int> AAChain::getResidueNumbers() const {
  /* Returns a list of residue numbers that are well-defined on this list. */

  vector<int> residueNumbers;
  residueNumbers.clear();

  map<int, AminoAcid*>::const_iterator itr = this->aa_m.begin();
  for (; itr != aa_m.end(); ++itr) {
    int residueN = itr->first;
    residueNumbers.push_back(residueN);
  }
  return residueNumbers;

}

AminoAcid* AAChain::operator[](int res_n) const {

  //  cout << "in operator[] " << endl; flush(cout);
  //  cout << "aa_m.size is " << aa_m.size() << endl; flush(cout);
 
  map<int, AminoAcid*>::const_iterator itr = this->aa_m.find(res_n);
  
  if (itr != aa_m.end()) {
    return itr->second;
  } else {
    cerr << "Error: Residue " << res_n << " not found in this chain." << endl;
    exit(8);
    return NULL;
  }
}

ScreamAtomV AAChain::get_sc_atoms(int pstn) const {

  

  return this->operator[](pstn)->get_sc_atoms();

}

double AAChain::PHI(int res_n) const {

  AminoAcid* aa = this->aa_m.find(res_n)->second;
  if (aa == this->aa_m.begin()->second) {
    return 0.0001;

  } else {

    AminoAcid* prev_aa = aa_m.find(res_n-1)->second;
    return aa->PHI(prev_aa);
  }

}

double AAChain::PSI(int res_n) const {

  AminoAcid* aa = aa_m.find(res_n)->second;

  if (aa == this->aa_m.end()->second) {
    return 0.0001;

  } else {
    AminoAcid* next_aa = aa_m.find(res_n+1)->second;
    return aa->PSI(next_aa);
  }

}

ScreamMatrix AAChain::set_PHI(int res_n, double angle) {

  ScreamVector origin_offset = ScreamVector( ( (AABackBone*) (this->aa_m.find(res_n)->second->get_rot()->get_bb()) )->N());
  ScreamMatrix R = aa_m.find(res_n)->second->set_PHI(angle, this->aa_m.find(res_n-1)->second);
  

  map<int, AminoAcid*>::const_iterator ptr;

  for (ptr = aa_m.begin(); ptr!= aa_m.end(); ++ptr) {
    if (ptr->first > res_n) {
      ptr->second->translate(ScreamVector(0,0,0) - origin_offset);
      ptr->second->transform(R);  // ptr->second is a AminoAcid*
      ptr->second->translate(origin_offset);
    }
  }

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = this->cterm_mm.begin(); itr!=cterm_mm.end(); ++itr) {
    ScreamVector this_atom = ScreamVector(itr->second);
    this_atom = this_atom - origin_offset;
    this_atom = R * this_atom;
    this_atom = this_atom + origin_offset;
    for (int i = 0; i<=2; ++i) {
      itr->second->x[i] = this_atom[i];
    }
  }
  
  return R;

}

ScreamMatrix AAChain::set_PSI(int res_n, double angle) {

  ScreamVector origin_offset = ScreamVector( ( (AABackBone*) (this->aa_m.find(res_n)->second->get_rot()->get_bb()) )->CA());
  ScreamMatrix R = aa_m.find(res_n)->second->set_PSI(angle, this->aa_m.find(res_n+1)->second);

  map<int, AminoAcid*>::const_iterator ptr;

  for (ptr = aa_m.begin(); ptr != aa_m.end(); ++ptr) {
    if (ptr->first > res_n) {
      ptr->second->translate(ScreamVector(0,0,0) - origin_offset);
      ptr->second->transform(R);  // ptr->second is a AminoAcid*
      ptr->second->translate(origin_offset);
    }
  }

  multimap<string, SCREAM_ATOM*>::const_iterator itr;
  for (itr = this->cterm_mm.begin(); itr!=cterm_mm.end(); ++itr) {
    ScreamVector this_atom = ScreamVector(itr->second);
    this_atom = this_atom - origin_offset;
    this_atom = R * this_atom;
    this_atom = this_atom + origin_offset;
    for (int i = 0; i<=2; ++i) {
      itr->second->x[i] = this_atom[i];
    }
  }
  
  return R;

}

void AAChain::replace_AminoAcid(int pstn, AminoAcid* new_AA) {
  // Replaces an AminoAcid at one position with another AminoAcid, deletes old AMinoAcid from memory..

  AminoAcid *old_AA = (this->aa_m)[pstn];
  (this->aa_m)[pstn] = new_AA;
  delete old_AA;

}

void AAChain::fix_toggle(bool value) {

  multimap<string, SCREAM_ATOM*>::iterator itr;
  for (itr = nterm_mm.begin(); itr != nterm_mm.end(); ++itr) {
    itr->second->fix_atom(value);
  }
  for (itr = cterm_mm.begin(); itr != cterm_mm.end(); ++itr) {
    itr->second->fix_atom(value);
  }

  map<int, AminoAcid*>::iterator itr2;
  for (itr2 = aa_m.begin(); itr2 != aa_m.end(); ++itr2) {
    itr2->second->fix_toggle(value);
  }

}

void AAChain::fix_toggle_sc_pstn(bool value, int res_n) {

  this->aa_m.find(res_n)->second->fix_sc_toggle(value);

}

void AAChain::SC_replacement(int res_n, const AARotamer* const aa_rot, string placementMethod, vector<double>& CreateCBParameters) {
  
  string res_name = scream_tools::three_letter_AA(aa_rot->get_resName());
  if (aa_m.find(res_n)->second->get_rot()->get_resName() == res_name) {
    //    cout << "AAChain SC_replacement? " << endl;
    if (res_name == "GLY") 
      this->aa_m.find(res_n)->second->SC_replacement(aa_rot, "Default", CreateCBParameters); // GLY case.
    else
      this->aa_m.find(res_n)->second->SC_replacement(aa_rot, placementMethod, CreateCBParameters);              // this would be SC_replacment without mutation.

  } else {
    cout << "Mutation??? This should never happen.  If mutation, should have been taken care of in Protein." << endl;
    exit(2);
    mutation_replacement(res_n, aa_rot);
  }

}

void AAChain::print_Me() const {

  //  cout << "in AAChain::print_Me() " << endl;

  map<int, AminoAcid*>::const_iterator itr;
  for (itr = this->aa_m.begin(); itr != aa_m.end(); ++itr) {
    itr->second->print_Me();
  }
  multimap<string, SCREAM_ATOM*>::const_iterator c_itr;
  for (c_itr = this->nterm_mm.begin(); c_itr != nterm_mm.end(); ++c_itr) {
    c_itr->second->dump();
  }
  for (c_itr = this->cterm_mm.begin(); c_itr != cterm_mm.end(); ++c_itr) {
    c_itr->second->dump();
  }

}

void AAChain::print_ordered_by_n() const {

  map<int, AminoAcid*>::const_iterator itr;
  for (itr = this->aa_m.begin(); itr != aa_m.end(); ++itr) {
    itr->second->print_ordered_by_n();

  }
  
  multimap<string, SCREAM_ATOM*>::const_iterator c_itr;
  map<int, SCREAM_ATOM*> ordered_m;

  // printing nterm_mm

  for (c_itr = this->nterm_mm.begin(); c_itr != nterm_mm.end(); ++c_itr) {
    ordered_m.insert(make_pair(c_itr->second->n, c_itr->second));
  }


  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->dump();
  }

  // printing c_term_mm
  
  ordered_m.clear();

  for (c_itr = this->cterm_mm.begin(); c_itr != cterm_mm.end(); ++c_itr) {
    ordered_m.insert(make_pair(c_itr->second->n, c_itr->second));
  }

  for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->dump();
  }

}

void AAChain::append_to_filehandle(ostream* ofstream_p) const {


  map<int, AminoAcid*>::const_iterator itr;
  for (itr = this->aa_m.begin(); itr != aa_m.end(); ++itr) {
    itr->second->append_to_filehandle(ofstream_p);

  }
  
  multimap<string, SCREAM_ATOM*>::const_iterator c_itr;
  map<int, SCREAM_ATOM*> ordered_m;

  // printing nterm_mm

  for (c_itr = this->nterm_mm.begin(); c_itr != nterm_mm.end(); ++c_itr) {
    ordered_m.insert(make_pair(c_itr->second->n, c_itr->second));
  }


  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->append_to_filehandle(ofstream_p);
  }

  // printing c_term_mm
  
  ordered_m.clear();

  for (c_itr = this->cterm_mm.begin(); c_itr != cterm_mm.end(); ++c_itr) {
    ordered_m.insert(make_pair(c_itr->second->n, c_itr->second));
  }

  for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->append_to_filehandle(ofstream_p);
  }
  

}


void AAChain::pdb_append_to_filehandle(ostream* ofstream_p) const {


  map<int, AminoAcid*>::const_iterator itr;
  for (itr = this->aa_m.begin(); itr != aa_m.end(); ++itr) {
    itr->second->pdb_append_to_filehandle(ofstream_p);

  }
  
  multimap<string, SCREAM_ATOM*>::const_iterator c_itr;
  map<int, SCREAM_ATOM*> ordered_m;

  // printing nterm_mm

  for (c_itr = this->nterm_mm.begin(); c_itr != nterm_mm.end(); ++c_itr) {
    ordered_m.insert(make_pair(c_itr->second->n, c_itr->second));
  }


  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->pdb_append_to_filehandle(ofstream_p);
  }

  // printing c_term_mm
  
  ordered_m.clear();

  for (c_itr = this->cterm_mm.begin(); c_itr != cterm_mm.end(); ++c_itr) {
    ordered_m.insert(make_pair(c_itr->second->n, c_itr->second));
  }

  for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->pdb_append_to_filehandle(ofstream_p);
  }
 
}


void AAChain::append_to_ostream_connect_info(ostream* ofstream_p) const {


  map<int, AminoAcid*>::const_iterator itr;
  for (itr = this->aa_m.begin(); itr != aa_m.end(); ++itr) {
    itr->second->append_to_ostream_connect_info(ofstream_p);

  }
  
  multimap<string, SCREAM_ATOM*>::const_iterator c_itr;
  map<int, SCREAM_ATOM*> ordered_m;

  // printing nterm_mm

  for (c_itr = this->nterm_mm.begin(); c_itr != nterm_mm.end(); ++c_itr) {
    ordered_m.insert(make_pair(c_itr->second->n, c_itr->second));
  }


  map<int, SCREAM_ATOM*>::const_iterator itr_m;
  for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->append_to_ostream_connect_info(ofstream_p);
  }

  // printing c_term_mm
  
  ordered_m.clear();

  for (c_itr = this->cterm_mm.begin(); c_itr != cterm_mm.end(); ++c_itr) {
    ordered_m.insert(make_pair(c_itr->second->n, c_itr->second));
  }

  for (itr_m = ordered_m.begin(); itr_m != ordered_m.end(); ++itr_m) {
    itr_m->second->append_to_ostream_connect_info(ofstream_p);
  }
 
}


void AAChain::print_chi_angle_spread(ostream* ofstream_p) const {

  double chi1, chi2, chi3, chi4, chi5;

  setiosflags(std::ios::fixed);
  setprecision(5);

  for (int i = aa_m.begin()->first; i < (aa_m.size() + aa_m.begin()->first); ++i) {
    *ofstream_p << setw(5) << i << ' ' << aa_m.find(i)->second->residue_type();

    chi1 = this->operator[](i)->chi1();
    if (chi1 == 1000) {
      *ofstream_p << endl;
      continue;
    } else {
      *ofstream_p << ' ' << setprecision(5) << setiosflags(std::ios::fixed) << setw(10) << chi1;
    }
    chi2 = this->operator[](i)->chi2();
    if (chi2 == 1000) {
      *ofstream_p << endl;
      continue;
    } else {
      *ofstream_p << ' ' << setprecision(5) << setiosflags(std::ios::fixed) << setw(10) << chi2;
    }
    chi3 = this->operator[](i)->chi3();
    if (chi3 == 1000) {
      *ofstream_p << endl;
      continue;
    } else {
      *ofstream_p << ' ' << setprecision(5) << setiosflags(std::ios::fixed) << setw(10) << chi3;
    }
    chi4 = this->operator[](i)->chi4();
    if (chi4 == 1000) {
      *ofstream_p << endl;
      continue;
    } else {
      *ofstream_p << ' ' << setprecision(5) << setiosflags(std::ios::fixed) << setw(10) << chi4;
    }  
    chi5 = this->operator[](i)->chi5();
    if (chi5 == 1000) {
      *ofstream_p << endl;
      continue;
    } else {
      *ofstream_p << ' ' << setprecision(5) << setiosflags(std::ios::fixed) << setw(10) << chi5 << endl;;
    }  
  } // end for loop

}

void AAChain::InitFromAtomList(const ScreamAtomV& atom_list) {
  /* Check input data. */
  Debug debugInfo("AAChain::InitFromAtomList(const ScreamAtomV& atom_list)");

  debugInfo.out("InitFromAtomList() in AAChain... ");

  if (atom_list.size() == 0) {
    cout << "Warning: AAChain to be initiated has zero atoms. " << endl;
  }

  /* Counting number of residues in this AAChain for information. */

  ScreamAtomVConstItr itr;
  int last_res_n = -1;
  int res_c = 0;
  
  for (itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    SCREAM_ATOM* this_atom = *itr;

    /* Tripeptide cases.  Highly esoteric--should encapsulate this somehow. */
    if ((*itr)->resName == string("ACE") || (*itr)->resName == string("NME") ) {
      continue;
    }

    if (this_atom->resNum != last_res_n) {
      ++res_c;                            // counting res_c
      last_res_n = this_atom->resNum;
    }
  }
  /* This is a horrible idea. */
  debugInfo.out("Total number of residues in this chain: ");
  debugInfo.out(itoa(res_c));
  
  last_res_n = -1;                            // reset parameters
  res_c = 0;
  
  // start populating map<int, AminoAcid*> aa_m.
  
  ScreamAtomV one_aa_v;
  
  for (itr = atom_list.begin(); itr != atom_list.end(); ++itr) {
    
    if ((*itr)->resName == string("ACE")) {
      nterm_mm.insert(make_pair((*itr)->atomLabel, *itr));
      continue;
    } else if ((*itr)->resName == string("NME")) {
      cterm_mm.insert(make_pair((*itr)->atomLabel, *itr));
      continue;
    }
    
    int cur_res_n = (*itr)->resNum;

    debugInfo.out(" resNum resName");
    debugInfo.out(" resNum:  " + itoa((*itr)->resNum) + " resName: " + (*itr)->resName );
    debugInfo.out(" curResN: " + itoa(cur_res_n) );

    // initial condition
    if (last_res_n == -1) {
      one_aa_v.push_back(*itr);
      last_res_n = cur_res_n;
      continue;
    } else if (last_res_n == cur_res_n) { // loop body
      one_aa_v.push_back(*itr);
      continue;
    } else if (last_res_n != cur_res_n) {
      debugInfo.out(" amino acid # about to be added: " + itoa(last_res_n));
      AminoAcid* new_AA = new AminoAcid(one_aa_v);

      debugInfo.out(" AminoAcid successfully added. ");
      this->aa_m.insert(make_pair(last_res_n, new_AA));
      
      ++res_c;
      one_aa_v.clear();
      one_aa_v.push_back(*itr);
      last_res_n = cur_res_n;
      continue;
    }
  }

  // end condition
  AminoAcid* new_AA = new AminoAcid(one_aa_v);
  debugInfo.out(" AminoAcid successfully added. ");
  this->aa_m.insert(make_pair(last_res_n, new_AA));
  
  /* Assigning chain designator.
   */

  debugInfo.out("Getting chain desination...");
  AminoAcid* aa = this->aa_m.begin()->second;
  AARotamer* aaRot = aa->get_rot();
  AABackBone* aaBB = aaRot->get_bb();
  SCREAM_ATOM* N = aaBB->N();

  if ( N == NULL ) {
    debugInfo.out("Using alternative atom (CA) for chain designation ");
    this->chain_desig = string(aaBB->CA()->chain);
  } else {
      this->chain_desig = string(N->chain);
  }

  debugInfo.out("Number of amino acids in this chain: " + itoa(aa_m.size()) );
  debugInfo.out("Chain name: " + this->chain_desig );
  
}

void AAChain::mutation_replacement(int res_n, const AARotamer* const in_rot) {

  /* first copy backbone info because the amino acid will be deleted.
   */

  assert(this->aa_m.find(res_n)->second->get_rot()->get_resName() != in_rot->get_resName() );

  //  AminoAcid* old_AA = aa_m.find(res_n)->second;
  //  BackBone* old_bb  = old_AA->get_rot()->get_bb();
  //  SideChain* old_sc = old_AA->get_rot()->get_sc();

  AARotamer in_aarot_deepcopy;
  in_aarot_deepcopy.deepcopy(*((AARotamer*)in_rot));

  // working on this!!!!

}
