use phf::{Map, phf_map};

#[rustfmt::skip]
pub static ATOM_ORDER_WEIGHTS: Map<&'static str, i32> = phf_map! {
    // --- N-Terminus Variants (0-99) ---
    "H1" => 0, "H2" => 1, "H3" => 2,
    "NT" => 10, "HT1" => 11, "HT2" => 12, "HT3" => 13,

    // --- Main Backbone (100-199) ---
    "N"   => 100,
    "HN"  => 110,
    "CA"  => 120,
    "HA"  => 130,
    "HA1" => 131, "HA2" => 132, // Glycine

    // --- Side Chain Atoms (ordered by topological distance from backbone) ---
    // -- Beta Position (200-299) --
    "CB"  => 200, "HB" => 210, "HB1" => 211, "HB2" => 212, "HB3" => 213,

    // -- Gamma Position (300-399) --
    "OG"  => 300, "HG" => 301, // SER, CYS (SG/HSG are aliased to OG/HG)
    "CG"  => 310, "HG1" => 311, "HG2" => 312, // PRO, ARG, GLN, GLU, LEU, LYS
    "OG1" => 320, // THR
    // HG1 on THR's OG1 is aliased to HG1
    "CG1" => 330, "HG11" => 331, "HG12" => 332, "HG13" => 333, // VAL, ILE
    "CG2" => 340, "HG21" => 341, "HG22" => 342, "HG23" => 343, // THR, ILE, LEU

    // -- Delta Position (400-499) --
    "SD"  => 400, // MET
    "CD"  => 410, "HD1" => 411, "HD2" => 412, // PRO, LYS, ARG
    "ND1" => 420, // HIS, TRP
    "OD1" => 425, // ASP, ASN
    "CD1" => 430, "HD11" => 431, "HD12" => 432, "HD13" => 433, // ILE, LEU
    "CD2" => 440, // HIS, PHE, TYR, TRP
    "ND2" => 450, "HD21" => 451, "HD22" => 452, // ASN
    "OD2" => 455, // ASP

    // -- Epsilon Position (500-599) --
    "NE"  => 500, "HE" => 501, // ARG
    "CE"  => 510, "HE1" => 511, "HE2" => 512, "HE3" => 513, // LYS, MET
    "OE1" => 515, // GLN, GLU
    "CE1" => 520, // HIS
    "NE1" => 530, // TRP
    "CE2" => 540, // HIS, PHE, TYR, TRP
    "NE2" => 550, "HE21" => 551, "HE22" => 552, // GLN
    "OE2" => 555, // GLU
    "CE3" => 560, // TRP

    // -- Zeta Position (600-699) --
    "NZ"  => 600, "HZ1" => 601, "HZ2" => 602, "HZ3" => 603, // LYS
    "CZ"  => 610, "HZ" => 611, // PHE, TYR, ARG
    "CZ2" => 620, // TRP
    "CZ3" => 630, // TRP

    // -- Eta Position (700-799) --
    "OH"  => 700, "HH" => 701, // TYR
    "NH1" => 710, "HH11" => 711, "HH12" => 712, // ARG
    "NH2" => 720, "HH21" => 722, "HH22" => 723, // ARG
    "CH2" => 730, "HH2" => 731, // TRP

    // --- Main Backbone Carbonyl & C-Terminus (900+) ---
    "C"   => 900,
    "O"   => 910,
    "OXT" => 920, "HXT" => 921,
};

pub static ATOM_NAME_ALIASES: Map<&'static str, &'static str> = phf_map! {
    // --- General ---
    "H" => "HN", // The most common backbone H alias

    // --- Backbone and Alpha Hydrogens ---
    "HCA" => "HA", "1HA" => "HA1", "2HA" => "HA2",

    // --- Beta Hydrogens ---
    "HCB" => "HB", "1HB" => "HB1", "2HB" => "HB2", "3HB" => "HB3",

    // --- Gamma Hydrogens ---
    "HCG" => "HG", "1HG" => "HG1", "2HG" => "HG2",
    "HCG1" => "HG11", "HCG2" => "HG23", // For VAL/ILE vs THR
    "SG" => "OG", "HSG" => "HG",     // CYS aliased to SER for topology ordering
    "HG1" => "HG", // THR's OG1-H aliased to generic gamma H

    // --- Delta Hydrogens ---
    "HCD" => "HD", "1HD" => "HD1", "2HD" => "HD2",
    "HCD1" => "HD11", "HCD2" => "HD21",
    "HND1" => "HD1",
    "HND2" => "HD21",

    // --- Epsilon Hydrogens ---
    "HCE"  => "HE", "1HE" => "HE1", "2HE" => "HE2",
    "HCE1" => "HE1", "HCE2" => "HE2", "HCE3" => "HE3",
    "HNE"  => "HE",
    "HNE1" => "HE1",
    "HNE2" => "HE21",

    // --- Zeta Hydrogens ---
    "HCZ"  => "HZ", "HCZ2" => "HZ2", "HCZ3" => "HZ3",
    "HNZ"  => "HZ1",

    // --- Eta Hydrogens ---
    "HOH"  => "HH",
    "HNH1" => "HH11",
    "HNH2" => "HH21",
    "HCH2" => "HH2",

    // --- Protonated State Aliases (from delta file analysis) ---
    "HSD" => "HD1",
    "HSE" => "HE2",
    "HOD1" => "OD1",
    "HOD2" => "OD2",
    "HOE1" => "OE1",
    "HOE2" => "OE2",

    // --- C-Terminus Aliases ---
    "HOXT" => "HXT",
    "OT1"  => "OXT",
    "OT2"  => "O",
};
