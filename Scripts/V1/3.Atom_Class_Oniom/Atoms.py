#!/usr/bin/python3
import re
##################################################################################################
class Atoms(object):
    def __init__(self):
        self.atom_dict = {1:'H', 3:'Li', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 11:'Na', 12:'Mg', 14:'Si', 15:'P', 16:'S', 17:'Cl', 19:'K',  20:'Ca', 35:'Br', 53:'I'} 
        self.SUPPORTED_1 = {"H","B","C","N","O","F","P","S","I","K"}  # one letter only
        self.SUPPORTED_2 = {"Li","Na","Mg","Si","Cl","Ca","Br"}      # two letter only
        self.DUMMY_PREFIXES = {"DUM":"DUM", "VS":"VS", "EP":"EP", "LP":"LP", "MW":"MW"}        # extend as needed
        self.MASS_REF = {
    "H": 1.008, "B": 10.811, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "P": 30.974, "S": 32.065, "Cl": 35.453, "Br": 79.904, "I": 126.904,
    "Li": 6.941, "Na": 22.990, "Mg": 24.305, "Si": 28.086, "K": 39.098, "Ca": 40.078
        }

        pass
################################################################################################## 
#What I don't know how to do is handling finding sigma and epsilon...
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
        """
        Returns Gaussian symbol: element like 'C'/'Cl' or 'Bq' for recognised dummy sites.
        Raises if cannot be resolved confidently.
        """
        name = atomname.strip()
        if not name:
            raise ValueError("Empty atomname cannot be resolved")

        up = name.upper()

        # 1) explicit dummy sites -> Gaussian Bq
        # eg MW, EP, VS, LP, DUM (with optional digits: MW1, EP2 etc)
        for pref in self.DUMMY_PREFIXES:
            if up.startswith(pref):
                return "Bq"

        # 2) handle the NA ambiguity using mass if present
        if up.startswith("NA") and mass is not None:
            pick = self._closest_by_mass(mass, ["N", "Na"], tol=1.5)
            if pick is None:
                raise ValueError(f"Ambiguous atomname '{atomname}' with mass {mass:.3f}")
            return pick

        # 3) try two-letter element from first two characters, case insensitive
        if len(name) >= 2 and name[0].isalpha() and name[1].isalpha():
            cand2 = name[0].upper() + name[1].lower()  # CL -> Cl
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
    def Atom_Types(self,Gros,Masses):
        Gaus_List=[]
        for i,m in zip(Gros,Masses): # If the atomic number is present we could just use atom_dict...
            Gaus = self.gaussian_symbol_from_atomname(i, m)
            Gaus_List.append(Gaus)
        print(Gaus_List)

        Gro_List=[]
        for n, i in enumerate(Gaus_List, start=1):
            Gro = f"{i}{n}"
            Gro_List.append(Gro)
        
        Heavy_Atoms=0
        Total_Atoms=0
        Dummy_List=[]
        for i in Gaus_List:
            print(f'{i}')
            if i =='H':
                Dummy_List.append(0)
            elif i=='Bq':
                Dummy_List.append(1)
            else:
                Heavy_Atoms+=1
                Dummy_List.append(0)
            Total_Atoms+=1
            
        return Heavy_Atoms, Total_Atoms, Gro_List, Gaus_List, Dummy_List
##################################################################################################
    def Atom_Symbol(self,Atom):
        atom=self.atom_dict[int(Atom)]
        return atom
##################################################################################################
string='1,6,6,7,8,1,53,35,16,11,9'
Strings=[s.strip() for s in string.split(',')]
Atoms=Atoms()
for S in Strings:
    atom=Atoms.Atom_Symbol(Atom=S)

Gros=['HNC','CN','ClCH3','MW','NCH','NAC3','NACl','NaCl']
Masses=[1.008, 12.011, 35.453, 0.000, 14.007, 14.007, 22.990, 22.990]
Heavy_Atoms,Gro_List,Gaus_List,Dummy_List=Atoms.Atom_Types(Gros,Masses)
for Gro, gaus,gro,dummy in zip(Gros,Gaus_List, Gro_List,Dummy_List):
    print(f'Gros:{Gro} gaus:{gaus} gro:{gro} dummy:{dummy}')
