#!/usr/bin/python3
from collections import defaultdict
##################################################################################################
class Atoms(object):
    def __init__(self):
        self.atom_dict = {1:'H', 3:'Li', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 11:'Na', 12:'Mg', 14:'Si', 15:'P', 16:'S', 17:'Cl', 19:'K',  20:'Ca', 35:'Br', 53:'I'} 
        self.SUPPORTED_1 = {"H","B","C","N","O","F","P","S","I","K"}  # one letter only
        self.SUPPORTED_2 = {"Li","Na","Mg","Si","Cl","Ca","Br"}      # two letter only
        self.DUMMY_PREFIXES = {"DUM":"DUM", "VS":"VS", "EP":"EP", "LP":"LP", "MW":"MW"}        # extend as needed
        self.MASS_REF = {
                        "H": 1.008, "Li": 6.941, "B": 10.811, "C": 12.011, "N": 14.007, "O": 15.999,
                        "F": 18.998, "Na": 22.990, "Mg": 24.305, "Si": 28.086, "P": 30.974, "S": 32.065,
                        "Cl": 35.453,  "K": 39.098, "Ca": 40.078, "Br": 79.904, "I": 126.904
                    } # Probably need to handle deviations from this
        self.AMBIGUOUS_PREFIXES = self._build_ambiguous_prefixes()
        self.ORGANIC_DEFAULTS = {
                        'NA': 'N',    # alpha-N in protein convention more common than sodium in organics
                        'CA': 'C',    # alpha-carbon
                        'CL': 'Cl',   # chlorine
                        'BR': 'Br',   # bromine
                        'SI': 'Si',   # silicon
                    }
        pass
##################################################################################################  
    def _build_ambiguous_prefixes(self):
        ambiguous = {}
        for two_letter in self.SUPPORTED_2:
            if two_letter[0] in self.SUPPORTED_1:
                ambiguous[two_letter.upper()] = [two_letter[0], two_letter]
        return ambiguous
##################################################################################################  
    def _closest_by_mass(self, mass: float, candidates: list[str], tol: float) -> str | None:
        best = None
        best_err = None
        for sym in candidates:
            err = abs(mass - self.MASS_REF[sym])
            if best_err is None or err < best_err:
                best, best_err = sym, err
        if best_err is not None and best_err <= tol:
            return best
        return None
##################################################################################################         
    def gaussian_symbol_from_atomname(self, atomname: str, mass: float | None = None) -> str:
        name = atomname
        if not name:
            raise ValueError("Empty atomname cannot be resolved")
    
        up = name.upper()

        # 1) explicit dummy sites -> Gaussian Bq
        for pref in self.DUMMY_PREFIXES:
            if up.startswith(pref):
                return "Bq"

        # 2) Ambiguous two-letter prefixes (e.g. CA, NA, CL, SI, BR). When the
        #    first two letters form a valid two-letter element AND the first
        #    letter alone is also a valid one-letter element, mass disambiguates.
        for prefix, candidates in self.AMBIGUOUS_PREFIXES.items():
            if up.startswith(prefix):
                if mass is None:
                    raise ValueError(
                        f"Ambiguous atomname '{atomname}' (could be one of "
                        f"{candidates}) requires a mass for disambiguation."
                    )
                pick = self._closest_by_mass(mass, candidates, tol=2.0)
                if pick is None:
                    raise ValueError(
                        f"Ambiguous atomname '{atomname}' with mass {mass:.3f} "
                        f"matches none of {candidates} within tolerance."
                    )
                return pick

        # 3) try two-letter element from first two characters, case insensitive
        if len(name) >= 2 and name[0].isalpha() and name[1].isalpha():
            cand2 = name[0].upper() + name[1].lower()
            if cand2 in self.SUPPORTED_2:
                return cand2

        # 4) try one-letter element from first character
        cand1 = name[0].upper()
        if cand1 in self.SUPPORTED_1:
            return cand1

        # 5) mass fallback across supported elements, if available
        if mass is not None:
            all_supported = list(self.SUPPORTED_1 | self.SUPPORTED_2)
            pick = self._closest_by_mass(mass, all_supported, tol=1.0)
            if pick is not None:
                return pick

        raise ValueError(f"Could not resolve element from atomname '{atomname}'")
##################################################################################################    
    def gaussian_symbol_from_gro_name(self, gro_name):
        if not gro_name:
            raise ValueError("Empty gro atom name cannot be resolved")
    
        up = gro_name.upper()
    
        # 1) Dummy/virtual sites → Bq
        for pref in self.DUMMY_PREFIXES:
            if up.startswith(pref):
                return "Bq"
    
        # 2) Ambiguous two-letter prefix → use organic-chemistry default
        for prefix, default in self.ORGANIC_DEFAULTS.items():
            if up.startswith(prefix):
                return default
    
        # 3) Two-letter element match
        if len(gro_name) >= 2 and gro_name[0].isalpha() and gro_name[1].isalpha():
            cand2 = gro_name[0].upper() + gro_name[1].lower()
            if cand2 in self.SUPPORTED_2:
                return cand2
    
        # 4) One-letter element match
        cand1 = gro_name[0].upper()
        if cand1 in self.SUPPORTED_1:
            return cand1
    
        raise ValueError(f"Could not resolve element from gro atom name '{gro_name}'")
##################################################################################################
    def Atom_Types(self,Gros,Masses):
        Gaus_List=[]
        for i,m in zip(Gros,Masses): # If the atomic number is present we could just use atom_dict...
            Gaus = self.gaussian_symbol_from_atomname(i, m)
            Gaus_List.append(Gaus)

        Gro_List=[]
        element_counters = defaultdict(int)
        Gro_List = []
        for elem in Gaus_List:
            element_counters[elem] += 1
            Gro_List.append(f"{elem}{element_counters[elem]}")
        
        Total_Atoms = len(Gaus_List)
        Dummy_List=[]
        for i in Gaus_List:
            Dummy_List.append(1 if i == 'Bq' else 0)
            
        return  Gro_List, Gaus_List, Dummy_List, Total_Atoms
##################################################################################################
    def Atom_Symbol(self,Atom):
        atom=self.atom_dict[int(Atom)]
        return atom
##################################################################################################
