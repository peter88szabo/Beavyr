// Auto-defined fragment geometries (XYZ coords only, Angstrom).
// The first atom is the connector atom for fragment attachment.

pub struct FragmentDef {
    pub name: &'static str,
    pub xyz: &'static str,
}

/// Source: -CCH.xyz
pub const FRAG_CCH: &str = r#"
 C     0.000000     0.000000     0.000000
 C     1.189000     0.000000     0.000000
 H     2.240000     0.000000     0.000000
"#;

/// Source: -CH=CH2.xyz
pub const FRAG_CH_CH2: &str = r#"
 C     0.000000     0.000000     0.000000
 C     0.000000     0.000000     1.335000
 H     0.943102     0.000000    -0.544500
 H     0.943102     0.000000     1.879500
 H    -0.943102     0.000000     1.879500
"#;

/// Source: -CH=O.xyz
pub const FRAG_CH_O: &str = r#"
 C     0.000000     0.000000     0.000000
 O     0.000000     0.000000     1.220000
 H     0.943102     0.000000    -0.544500
"#;

/// Source: -COOH.xyz
pub const FRAG_COOH: &str = r#"
 C     0.000000     0.000000     0.000000
 O     0.000000     0.000000     1.400000
 O     1.056551     0.000000    -0.610000
 H     0.895669     0.000000     1.716667
"#;

/// Source: -CycloHexane.xyz
pub const FRAG_CYCLOHEXANE: &str = r#"
 C     0.000000     0.000000     0.000000
 C     0.000000     0.000000     1.500000
 H     1.008807     0.000000    -0.356663
 C    -0.707107     1.224745     2.000000
 H    -0.504403    -0.873651     1.856667
 H     1.008806     0.000000     1.856667
 C    -2.121320     1.224745     1.500000
 H    -0.707107     1.224745     3.070000
 H    -0.202704     2.098396     1.643333
 C    -2.121320     1.224745     0.000000
 H    -2.625723     0.351093     1.856667
 H    -2.625723     2.098396     1.856667
 C    -1.414214     0.000000    -0.500000
 H    -3.130126     1.224745    -0.356667
 H    -1.616917     2.098396    -0.356667
 H    -1.918616    -0.873651    -0.143333
 H    -1.414214     0.000000    -1.570000
"#;

/// Source: -CycloPentane.xyz
pub const FRAG_CYCLOPENTANE: &str = r#"
 C     0.000000     0.000000     0.000000
 C     0.000000     0.000000     1.546000
 C     1.509354     0.000000     1.880616
 C     2.065338     1.075069     0.918730
 C     1.202221     0.869825    -0.347374
 H    -0.459782     0.929692     1.921343
 H    -0.476802    -0.885923     1.998115
 H     1.940819    -0.984445     1.633000
 H     1.648971     0.233932     2.949443
 H     3.151686     0.930424     0.794090
 H     1.887893     2.076822     1.344891
 H     1.809833     0.389425    -1.132633
 H     0.840971     1.847422    -0.708499
 H     0.164473    -1.026837    -0.367645
"#;

/// Source: -NH2.xyz
pub const FRAG_NH2: &str = r#"
 N     0.000000     0.000000     0.000000
 H     0.000000     0.000000     1.030000
 H     0.971095     0.000000    -0.343330
"#;

/// Source: -NO2.xyz
pub const FRAG_NO2: &str = r#"
 N     0.000000     0.000000     0.000000
 O     0.000000     0.000000     1.210000
 O     0.988136     0.000000    -0.698346
"#;

/// Source: -OCH3.xyz
pub const FRAG_OCH3: &str = r#"
 O     0.000000     0.000000     0.000000
 C     0.000000     0.000000     1.400000
 H     1.026720     0.000000     1.762996
 H    -0.513360     0.889165     1.763000
 H    -0.513360    -0.889165     1.763000
"#;

/// Source: -OH.xyz
pub const FRAG_OH: &str = r#"
 O     0.000000     0.000000     0.000000
 H     0.000000     0.000000     0.947000
"#;

/// Source: -OOH.xyz
pub const FRAG_OOH: &str = r#"
 O     0.000000     0.000000     0.000000
 O     0.000000     0.000000     1.260000
 H     0.892841     0.000000     1.575663
"#;

/// Source: -Phenyl.xyz
pub const FRAG_PHENYL: &str = r#"
 C     0.000000     0.000000     0.000000
 C     0.000000     0.000000     1.400000
 C     1.212436     0.000000     2.100000
 C     2.424871     0.000000     1.400000
 C     2.424871     0.000000     0.000000
 C     1.212436     0.000000    -0.700000
 H    -0.943102     0.000000     1.944500
 H     1.212436     0.000000     3.189000
 H     3.367973     0.000000     1.944500
 H     3.367973     0.000000    -0.544500
 H     1.212436     0.000000    -1.789000
"#;

/// Source: -Pyrrole.xyz
pub const FRAG_PYRROLE: &str = r#"
 N     0.000000     0.000000     0.000000
 C     0.000000     0.000000     1.301000
 C     1.272510     0.000000     1.704664
 C     2.074733     0.000000     0.474055
 C     1.189444     0.000000    -0.527089
 H     1.429077     0.000000    -1.590422
 H     3.161806     0.000000     0.394225
 H     1.639943     0.000000     2.730867
 H    -0.875071     0.000000     1.950885
"#;

/// Source: CH3.xyz
pub const FRAG_CH3: &str = r#"
 C     0.000000     0.000000     0.000000
 H     0.000000     0.000000     1.070000
 H     1.008807     0.000000    -0.356663
 H    -0.504403    -0.873651    -0.356667
"#;

pub const FRAGMENTS: &[FragmentDef] = &[
    FragmentDef { name: "-CH3", xyz: FRAG_CH3 },
    FragmentDef { name: "-CH=CH2", xyz: FRAG_CH_CH2 },
    FragmentDef { name: "-OH", xyz: FRAG_OH },
    FragmentDef { name: "-NH2", xyz: FRAG_NH2 },
    FragmentDef { name: "-Phenyl", xyz: FRAG_PHENYL },
    FragmentDef { name: "-Pyrrole", xyz: FRAG_PYRROLE },
    FragmentDef { name: "-OCH3", xyz: FRAG_OCH3 },
    FragmentDef { name: "-CH=O", xyz: FRAG_CH_O },
    FragmentDef { name: "-COOH", xyz: FRAG_COOH },
    FragmentDef { name: "-NO2", xyz: FRAG_NO2 },
    FragmentDef { name: "-OOH", xyz: FRAG_OOH },
    FragmentDef { name: "-CCH", xyz: FRAG_CCH },
    FragmentDef { name: "-CycloPentane", xyz: FRAG_CYCLOPENTANE },
    FragmentDef { name: "-CycloHexane", xyz: FRAG_CYCLOHEXANE },
];
