# Ancient DNA authentication: per-sample verdicts

| Sample | Species | Protocol | Verdict | 5' C>T excess | 3' G>A excess | Control | Endogenous | Merge rate | Median frag |
|---|---|---|---|---|---|---|---|---|---|
| TRP002.A | ypestis | partial_udg | **ANCIENT** | 0.1918 | 0.1811 | 0.0032 | 0.2329 | n/a | 42.0 |
| LBG002.A | ypestis | full_udg | **UDG_TREATED_RESIDUAL_SIGNAL** | 0.0007 | 0.0101 | 0.0028 | 0.3198 | n/a | 56.0 |
| LVC005.A | ypestis | undeclared | **WEAK_ANCIENT_SIGNAL** | 0.0054 | 0.0181 | 0.0072 | 0.1834 | n/a | 53 |
| LVC001.C | ypestis | undeclared | **NO_DAMAGE_DETECTED** | 0.0080 | 0.0083 | 0.0094 | 0.0699 | 0.916 | 43 |
| LP31b | ypestis | undeclared | **UNRELIABLE** | 0.0150 | -0.0027 | 0.0112 | 0.0246 | 0.915 | 49 |
| MOD1953 | ypestis | modern | **MODERN_AS_EXPECTED** | 0.0000 | -0.0000 | 0.0015 | 0.9472 | 0.858 | 305 |
| Motala1 | human | undeclared | **NO_DAMAGE_DETECTED** | 0.0031 | 0.0034 | 0.0035 | 0.0011 | n/a | 54.0 |

## Reasoning

**TRP002.A** (partial_udg) - ANCIENT  
terminal C>T excess +0.1918 and G>A excess +0.1811 both clear 0.02, they agree within 1.1x, and the control stays flat at +0.0032

**LBG002.A** (full_udg) - UDG_TREATED_RESIDUAL_SIGNAL  
3' G>A excess +0.0101 decays monotonically from the terminus (5' C>T 0.0012, 0.0011, 0.0010; 3' G>A 0.0106, 0.0049, 0.0030) with the control flat at +0.0028. A small residual after full UDG treatment is expected, since the enzyme does not excise every uracil, so this is not evidence against the treatment having worked

**LVC005.A** (undeclared) - WEAK_ANCIENT_SIGNAL  
3' G>A excess +0.0181 decays monotonically from the terminus (5' C>T 0.0062, 0.0031, 0.0020; 3' G>A 0.0188, 0.0099, 0.0057) with the control flat at +0.0072; below the strong threshold (0.02) but far above the noise floor and damage-shaped. Consistent with a partly UDG-treated or well-preserved library, or with a modest ancient fraction diluted by modern DNA; damage alone cannot say which

**LVC001.C** (undeclared) - NO_DAMAGE_DETECTED  
no terminal damage (C>T +0.0080); the library protocol is not declared in the archive, so this is consistent with a UDG-treated library, with modern contamination throughout, or with a signal too weak for this depth

**LP31b** (undeclared) - UNRELIABLE  
the A>G control rises by +0.0112 at the terminus, above the 0.01 tolerance, so terminal mismatches are not driven by deamination alone and the C>T estimate (+0.0150) cannot be separated from that artefact

**MOD1953** (modern) - MODERN_AS_EXPECTED  
no terminal damage (C>T +0.0000, control +0.0015); this is the negative control and the estimator correctly stayed silent

**Motala1** (undeclared) - NO_DAMAGE_DETECTED  
no terminal damage (C>T +0.0031); the library protocol is not declared in the archive, so this is consistent with a UDG-treated library, with modern contamination throughout, or with a signal too weak for this depth

## Did the estimator behave?

- True negative control (modern, must stay silent):
  - MOD1953: MODERN_AS_EXPECTED, C>T excess 0.0000 -> silent, as required
- Positive control (declared damaged, must fire):
  - TRP002.A: ANCIENT, C>T excess 0.1918 -> detected
- Protocol undeclared (no expectation, reported as found):
  - LVC005.A: WEAK_ANCIENT_SIGNAL
  - LVC001.C: NO_DAMAGE_DETECTED
  - LP31b: UNRELIABLE
  - Motala1: NO_DAMAGE_DETECTED
