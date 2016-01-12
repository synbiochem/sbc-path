# This tool imports a pathway into a chassis organisms
# - It has to deal with missing metabolites
# - It can be extended in order to perform comparative analysis, simulations, etc.
from os import path, getenv
import sys
SBCPATH = path.join(getenv('HOME'), 'mibsspc2', 'synbiochemdb_scripts')
sys.path.insert(0, SBCPATH)
import sbcpsql
import cobra

def parse_equation(eq):
    d = []
    m = eq.split()
    i = 0
    c = []
    while i < len(m):
        if m[i] == '=':
            d.append(c)
            c = []
            i += 1
        elif m[i] == '+':
            i += 1
        else:
            try:
                c.append((m[i], m[i+1]))
                i += 2
            except:
                break
    if len(c) > 0:
        d.append(c)
    return d
        


class Path:
    
    def __init__(self):
        self.rlist = {}
        self.precursors = set()
        self.products = set()
        self.sbcdb = sbcpsql.sbcdb()
        self.sbcdb.chem_xref()
        self.sbcdb.reac_xref()
        self.sbcdb.dict_reac()


    # reaction id (it needs to be in the database)
    # info: additional info (dictionary)
    def add_reaction(self, r, info):
        if r in self.sbcdb.rxrefreac:
            self.rlist[r] = info
            try:
                eq = self.sbcdb.dictreac[r]
                req = parse_equation(eq[1])
                self.rlist[r]['subs'] = req[0]
                self.rlist[r]['prods'] = req[1]
            except:
                pass

    # Compute overall path stoichiometry
    def path_stoichiometry(self):
        stoi = {}
        for r in self.rlist:
            factor = 1.0
            for side in ['subs', 'prods']:
                factor *= -1 
                for c in self.rlist[r]['subs']:
                    try:
                        n = float(c[0])
                        cmp = c[1]
                    except:
                        # Take into account only well-defined stoichiometry
                        continue
                    if cmp not in stoi:
                        stoi[cmp] = 0.0
                    stoi[cmp] += factor*n
        import pdb
        pdb.set_trace()
        bal = set()
        for c in stoi:
            if stoi[c] == 0:
                bal.add(c)
        for c in bal:
            del stoi[c]
        self.stoi = stoi

    def map_stoi(self, db='bigg'):
        mstoi = {}
        for c in stoi:
            mfound = False
            if c in self.sbcdb.xrefchem:
                for xref in self.sbcdb.xrefchem[c]:
                    if xref.startswith('bigg'):
                        mstoi[self.sbcdb.xrefchem[c]] = stoi[c]
                        mfound = Trud
            if not mfound:
                mstoi[c] = stoi[c]
        return mstoi

            

# Import the model
def import_model(modelfile, format='JSON'):
    if format == 'JSON':
        return cobra.io.load_json_model(modelfile)
    elif format == 'SBML':
        return cobra.io.load_sbml_model(modelfile)        
    elif format == 'MATLAB':
        return cobra.io.load_matlab_model(modelfile)

# Import the pathway
def import_pathway(pathfile):
    p = Path()
    # Basic format: reaction_id, enzyme_id
    for l in open(pathfile):
        m = l.rstrip().split()
        p.add_reaction(m[0], {'sequence': m[1]})
        
    return p
        
