# This tool imports a pathway into a chassis organisms
# - It has to deal with missing metabolites
# - It can be extended in order to perform comparative analysis, simulations, etc.
from os import path, getenv
import sys
# This path would need to be updated to a more stable location
SBCPATH = path.join(getenv('HOME'), 'mibsspc2', 'synbiochemdb_scripts')
sys.path.insert(0, SBCPATH)
import sbcpsql
import cobra
import copy
import re

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
        
# Find bigg metabolite in list
def find_cid(cid, cx):
    if cid in cx:
        return cx[cid]
    else:
        cid = re.sub('__', '_', cid)
        m = cid.split('_')
        cid1 = '_'.join(m[0:(len(m)-1)])
        cid2 = '-'.join(m[0:(len(m)-1)])
        if cid1 == 'sertrna_sec':
            cid1 = 'sertrna(sec)'
            return cx[cid1]
        if cid1 in cx:
            return cx[cid1]
        elif cid2 in cx:
            return cx[cid2]
    return None

def remove_compartment(cid, compartment='c'):
    return re.sub('_'+compartment, '', cid)


def new_metabolite(mid, formula, name, compartment):
    # Create the metabolites
    met = cobra.Metabolite(mid, formula=formula,
                     name=name, compartment=compartment)
    return(met)

def new_reaction(rid, name, react, subsys='SYNBIOCHEM', rev=0, obj=0.0, lb=0, ub=9999):
    metab = {}
    for x in react['subs']:
        if x[2] not in metab:
            metab[x[2]] = 0
        try:
            metab[x[2]] -= float(x[0])
        except:
            pass
    for x in react['prods']:
        if x[2] not in metab:
            metab[x[2]] = 0
        try:
            metab[x[2]] += float(x[0])
        except:
            pass
    rm = set()
    for c in metab:
        if metab[c] == 0:
            rm.add(c)
    for c in rm:
        del metab[c]
    if rev != 0:
        lb = -99999
    reaction = cobra.Reaction(rid) # 'my_new_reaction'
    reaction.name = name #'3 oxoacyl acyl carrier protein synthase n C140'
    reaction.subsystem = subsys #'Cell Envelope Biosynthesis'
    reaction.lower_bound = lb  # This is the default
    reaction.upper_bound = ub  # This is the default
#    reaction.reversibility = rev  # This is the default
    reaction.objective_coefficient = obj  # This is the default
    reaction.add_metabolites(metab)
#    print rid,reaction.reaction
    return reaction

def add_transport(model, met):
    mid = met.id
    midp = met.id+'_p'
    mide = met.id+'_e'
    formula = met.formula
    name = met.name          
    met_p = new_metabolite(midp, formula, name, compartment='Periplasm')
    model.add_metabolites(met_p)
    met_e = new_metabolite(mide, formula, name, compartment='Extra_organism')
    model.add_metabolites(met_e)
    stoit = {'subs': {met: 1}, 'prods': {met_p: 1}}
    r1 = new_reaction(mid+'Ppp',name='Transport',
                 subsys='SYNBIOCHEM Transport', react=stoit, rev=1)
    stoit = {'subs': {met_p: -1}, 'prods': {met_e: 1}}
    r2 = new_reaction(mid+'tex',name='Transport',
                 subsys='SYNBIOCHEM Transport', react=stoit, rev=1)
    stoit = {'subs': {met_e: -1}, 'prods':{}}
    r3 = new_reaction('EX_'+mide,name='Transport',
                 subsys='SYNBIOCEM Transport', react=stoit, rev=1)
    model.add_reactions( [r1, r2, r3])


class Path:
    

    def __init__(self):
        self.rlist = {}
        self.rlistchassis = {}
        self.precursors = set()
        self.products = set()
        self.sbcdb = sbcpsql.sbcdb()
        self.sbcdb.chem_xref()
        self.sbcdb.reac_xref()
        self.sbcdb.dict_reac()
        self.chassis = None
        self.metdb = {}
        self.dbmet = {}
        self.pathstoi = None

    def __repr__(self):
        return str(self.rlist)


        
   # Map metabolites in the pathway into chassis
    def plug_path_chassis(self, db='bigg', map_reaction=False):
        stats = {'ml': set(), 'mdbl': set()}
        for r in self.rlist:
            rx = r
            # If we map the reactions into the ones in bigg,
            # then do not add them to the model if they already exist
            if map_reaction:
                for ri in self.sbcdb.rxrefreac[r]:
                    if ri.startswith(db):
                        rx = re.sub(db+':', '', ri)
                        break
            self.rlistchassis[rx] = copy.deepcopy(self.rlist[r])
            for cc in ['subs', 'prods']:
                self.rlistchassis[rx][cc] = []
                for c in self.rlist[r][cc]:
                    cx = c[1]
                    # Check if the metabolite is already in the E. coli model
#                    if cx == 'MNXM1':
#                        import pdb
#                        pdb.set_trace()
                    if cx in self.dbmet:
                        for ci in self.dbmet[cx]:
                            if ci.endswith('_'+self.rlistchassis[rx]['compartment']):
                                cx = ci
                                stats['mdbl'].add(cx)
                                break
                    stats['ml'].add(cx)
                    if cx not in self.chassis.metabolites:
                        met = new_metabolite(cx, formula=None, name=cx, compartment=self.rlistchassis[rx]['compartment'])
                        self.chassis.add_metabolites([met])
                    self.rlistchassis[rx][cc].append((c[0], cx, getattr(self.chassis.metabolites, cx)))
            nrx = new_reaction(rx, rx, self.rlistchassis[rx])
            self.chassis.add_reaction(nrx)
        self.message(['Reactions:', len(self.rlist), 'Metabolites:', len(stats['ml']),
                      'New metabolites:', len(stats['ml'] - stats['mdbl'])])
            


    def add_path(self, pathfile):
        # Basic format: reaction_id, direction, reversibility, enzyme_id (optional)
        for l in open(pathfile):
            if l.startswith('#'):
                continue
            m = l.rstrip().split('\t')
            reaction = m[0]
            direction = 1
            reversibility = 0
            compartment = 'c'
            enzyme_id = 'e_'+reaction
            comment = None
            try:
                direction = m[1]
                reversibility = m[2]
            except:
                pass
            try:
                compartment = m[3]
            except:
                pass
            try:
                enzyme_id = m[4]
            except:
                pass
            try:
                comment = m[5]
            except:
                pass
            self.add_reaction(reaction, {'sequence': enzyme_id, 'compartment': compartment, 'comment': comment})


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
        for r in self.rlistchassis:
            factor = 1.0
            for side in ['subs', 'prods']:
                factor *= -1
                for c in self.rlistchassis[r][side]:
                    try:
                        n = float(c[0])
                        cmp = c[1]
                    except:
                        # Take into account only well-defined stoichiometry
                        continue
                    if cmp not in stoi:
                        stoi[cmp] = 0.0
                    stoi[cmp] += factor*n
        bal = set()
        for c in stoi:
            if stoi[c] == 0:
                bal.add(c)
        for c in bal:
            del stoi[c]
        stats = set()
        for c in stoi:
            if c not in self.metdb:
                add_transport(self.chassis, getattr(self.chassis.metabolites, c))
                stats.add(c)
        self.message(['Added transport for metabolites:', stats])
        self.pathstoi = stoi

            

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



    def map_chassis_met_sbcdb(self, db='bigg:'):
        cnx = {}
        for x in self.sbcdb.xrefchem:
            if x.startswith(db):
                cnx[re.sub(db, '', x)] = self.sbcdb.xrefchem[x]
        for m in self.chassis.metabolites:
            cid = remove_compartment(m.id, m.compartment)
            msbc = find_cid(cid, cnx)
            if msbc is not None:
                self.metdb[m.id] = msbc
                if msbc not in self.dbmet:
                    self.dbmet[msbc] = set()
                self.dbmet[msbc].add(m.id)
        self.message(['Species in chassis:', len(self.chassis.metabolites)])
        self.message(['Mapped metabolites:',  len(self.dbmet)])
        self.message(['Missing species', len(self.chassis.metabolites) -len(self.metdb)])


    # Before mapping the pathway to the database, we need to declare our chassis
    def add_chassis(self, modelfile, format='JSON'):
        if format == 'JSON':
            self.chassis = cobra.io.load_json_model(modelfile)
        elif format == 'SBML':
            self.chassis = cobra.io.load_sbml_model(modelfile)        
        elif format == 'MATLAB':
            self.chassis = cobra.io.load_matlab_model(modelfile)
        self.message(["Chassis", modelfile, "added"])
        self.map_chassis_met_sbcdb()

        
    def message(self, l):
        print ' '.join(map(str, l))

def ptest():
    p = Path()
    p.add_chassis('../../data/strains/iAF1260.json')
    p.add_path('limonene.path')
    p.plug_path_chassis()
    p.path_stoichiometry()
    # TO DO: routine that verifies that a pathway is viable
    # for r in p.chassis.reactions:
    #    if r.objective_coefficient != 0:
    #        x.append(r)
    # biom = x[0]
    # biom.objective_coefficient = 0
    # limon = p.chassis.reactions.get_by_id('MNXR70692')
    # limon.objective_coefficient = 1
    # p.chassis.optimize()
    return p
