import codecs
import re
import json
import glob
import os
import math
import hashlib

import numpy as np
import CoolProp.CoolProp as CP

def lines_starts_with(s, lines, *, start_index=0):
    matches = []
    for i in range(start_index, len(lines)):
        if lines[i].startswith(s):
            matches.append(i)
    return matches

def lines_contains(s, lines, *, start_index=0):
    matches = []
    for i in range(start_index, len(lines)):
        if s in lines[i]:
            matches.append(i)
    return matches

def lines_equals(s, lines, *, start_index=0):
    matches = []
    for i in range(start_index, len(lines)):
        if re.fullmatch(s, lines[i]):
            matches.append(i)
    return matches


class FLDDeconstructor:
    """
    Read a REFPROP fluid file, and convert to the 
    JSON format of CoolProp
    """
    def __init__(self, path, encoding='cp1252'):
        """
        Load in the FLD file, do some pre-processing
        """
        self.path = path
        with codecs.open(path, encoding=encoding) as fp:
            self.contents = fp.read()
            self.lines = self.contents.split('\n')

    def get_keyed_line(self, block, key, conv, *, alias=None):
        keys = [key]
        if alias:
            keys += [alias]
        for k in keys:
            if not k in ''.join(block): continue
            lines = list(filter(lambda l: k in l, block))
            if len(lines) == 1:
                els = lines[0].split('!')[0]
                return [conv(el) for el in els.split(' ') if el]
        raise KeyError(key)

    
    def get_block(self, startreg, endreg, start_index=0, indices=False):
        matches = lines_starts_with(startreg, self.lines, start_index = start_index)
        istart = matches[0]
        matches = lines_equals(endreg, self.lines, start_index = istart+1)
        iend = matches[0]
        if indices:
            return self.lines[istart:iend], istart, iend
        else:
            return self.lines[istart:iend]

    def strip_commented(self, block):
        return list(filter(lambda l: not l.startswith('?') and not l.startswith('!'), block))

    def strip_comments(self, block):
        block = list(filter(lambda l: not l.startswith('?') and not l.startswith('!'), block))
        for iline, line in enumerate(block):
            if '!' in line:
                block[iline] = line.split('!')[0].rstrip() + '\n'
        return block

    def get_ancillary_description(self, key):
        return {
        # Most are in this form:
        'PS5':  'P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T]',
        'DL1':  'D=Dc*[1+SUM(Ni*Theta^ti)]',
        'DV3':  'D=Dc*EXP[SUM(Ni*Theta^ti)]',

        # Other ones might be of this form:
        'DL2':  'D=Dc*[1+SUM(Ni*Theta^(ti/3))]',
        'DL4':  'D=Dc*EXP[SUM(Ni*Theta^(ti/3))]',
        'DV6':  'D=Dc*EXP[SUM(Ni*Theta^(ti/3))*Tc/T]',
        'DV4':  'D=Dc*EXP[SUM(Ni*Theta^(ti/3))]',
        'DL6':  'D=Dc*EXP[SUM(Ni*Theta^(ti/3))*Tc/T]',
        }[key]

    def get_ancillary(self, key):
        block = self.get_block('#'+key, r'\s*')
        block = self.strip_comments(block)
        model_key = block[1].strip()
        limits = [float(block[i]) for i in range(2, 6)]
        reducing = [float(v) for v in block[6].strip().split(' ') if v]
        Ncoeffs = [int(v) for v in block[7].split(' ') if v]
        n, t = [], []
        reducing[1] *= 1000 # REFPROP uses mol/L & kPa, CoolProp uses mol/m^3 & Pa, multiply by 1000 to convert
        for row in block[8::]:
            _n, _t = [_ for _ in row.lstrip().split(' ') if _][0:2]
            n.append(float(_n.replace('D+','e+')))
            t.append(float(_t))
        assert(len(n) == Ncoeffs[0] and len(t) == Ncoeffs[0])
        desc = self.get_ancillary_description(model_key.strip())
        if model_key[0:2] == 'PS':
            type_key = 'pL'
        elif model_key in ['DL1','DL2']:
            type_key = 'rhoLnoexp'
        elif model_key[0:2] == 'DL':
            type_key = 'rhoL'
            t = t if model_key[2] in ['1','3','5'] else [t_/3 for t_ in t]
        elif model_key[0:2] == 'DV':
            type_key = 'rhoV'
            t = t if model_key[2] in ['1','3','5'] else [t_/3 for t_ in t]
        else:
            raise KeyError(model_key)

        return {
            "T_r": reducing[0],
            "Tmax": reducing[0],
            "Tmin": limits[0],
            'description': desc,
            'n': n,
            'reducing_value': reducing[1],
            't': t,
            'type': type_key,
            'using_tau_r': 'Tc/T' in desc
        }

    def get_all_ancillaries(self):
        return{
            'PS': self.get_ancillary('PS'),
            'DV': self.get_ancillary('DV'),
            'DL': self.get_ancillary('DL')
        }

    def get_ancillaries(self):
        anc = self.get_all_ancillaries()
        return {'pS': anc['PS'], 'rhoL': anc['DL'], 'rhoV': anc["DV"]}

    def formula_from_inchi(self):
        """ Parse the standard inchi string to extract the chemical formula in Hill order """
        try:
            stdinchistring = self.get_keyed_line(self.lines, '!Standard InChI String', lambda x: x)[0]
        except:
            return 'COULD NOT PARSE InChI STRING'
        if '/' not in stdinchistring:
            raise ValueError(f'{stdinchistring} is not a valid standard InChI key')
        formula = stdinchistring.split('/')[1]
        matches = re.findall(r'[A-Z][a-z]*[0-9]*', formula)
        assert(sum([len(m) for m in matches]) == len(formula)) # make sure all parts are matched
        o = ''
        for match in matches:
            def replacer(m):
                N = len(m.groups())
                # print(m, N)
                if N == 0:
                    return m.group(0)+'_{1}'
                else:
                    i = m.group(1)
                    return f'_{{{i}}}'
            newmatch = re.sub(r'([0-9]+)', replacer, match)
            if '_{' not in newmatch:
                newmatch += '_{1}'
            o += newmatch
        return o

    def get_info(self):
        CHEMSPIDER_ID = -1
        i = lines_contains('!Short name', self.lines)
        short_name = self.lines[i[0]].split('!')[0].strip()
        name = os.path.split(self.path)[1].split('.')[0]

        i = lines_contains('!Full name', self.lines)
        full_name = self.lines[i[0]].split('!')[0].strip()

        name = os.path.split(self.path)[1].split('.')[0]
        aliases = [short_name, short_name.lower(), full_name, full_name.lower(), full_name.replace(' ','').upper(), short_name.replace(' ','').upper()]
        aliases = [alias for alias in aliases if alias != name]
        
        StdInChIstr = 'UNKNOWN'
        try:
            StdInChIstr = self.get_keyed_line(self.lines, '!Standard InChI String', lambda x: x)[0]
            if not StdInChIstr.startswith('InChI='):
                StdInChIstr = "InChI=" + StdInChIstr
        except KeyError:
            print('Unable to find InChI key')

        InChIKey = 'UNKNOWN'
        try:
            InChIKey = self.get_keyed_line(self.lines, '!Standard InChI Key', lambda x: x)[0]
        except KeyError:
            print('Unable to find InChI key')

        GWP100 = 'UNKNOWN'
        try:
            GWP100 = float(self.get_keyed_line(self.lines, '!GWP', lambda x: x)[0])
        except KeyError:
            print('Unable to find GWP')

        ODP = 'UNKNOWN'
        try:
            ODP = float(self.get_keyed_line(self.lines, '!ODP', lambda x: x)[0])
        except KeyError:
            print('Unable to find ODP')

        ASHRAE34 = 'UNKNOWN'
        try:
            ASHRAE34 = self.get_keyed_line(self.lines, '!Safety Group', lambda x: x)[0]
        except KeyError:
            print('Unable to find Safety group')
        
        return {
            "2DPNG_URL": f"http://www.chemspider.com/ImagesHandler.ashx?id={CHEMSPIDER_ID}",
            "ALIASES": list(set(aliases)),
            "CAS": self.get_keyed_line(self.lines, '!CAS number', lambda x: x)[0],
            "CHEMSPIDER_ID": CHEMSPIDER_ID,
            "ENVIRONMENTAL": {
              "ASHRAE34": ASHRAE34,
              "FH": 0,
              "GWP100": GWP100,
              "GWP20": 0.0,
              "GWP500": 0.0,
              "HH": 0,
              "Name": name,
              "ODP": ODP,
              "PH": 1e30
            },
            "FORMULA": self.formula_from_inchi(),
            "INCHI_KEY": InChIKey,
            "INCHI_STRING": StdInChIstr,
            "NAME": name,
            "REFPROP_NAME": name,
            "SMILES": "?"
        }

    def get_critical_state(self):

        i = lines_contains('!Critical temperature [K]', self.lines)
        Tc = float(self.lines[i[0]].split('!')[0].strip())

        i = lines_contains('!Critical pressure [kPa]', self.lines)
        pc_kPa = float(self.lines[i[0]].split('!')[0].strip())

        i = lines_contains('!Critical density [mol/L]', self.lines)
        rhoc_moldm3 = float(self.lines[i[0]].split('!')[0].strip())
        
        return {
              "T": Tc,
              "T_units": "K",
              "hmolar": -99999999999,
              "hmolar_units": "J/mol",
              "p": pc_kPa*1e3,
              "p_units": "Pa",
              "rhomolar": rhoc_moldm3*1e3,
              "rhomolar_units": "mol/m^3",
              "smolar": 999999999999999,
              "smolar_units": "J/mol/K"
        }

    def ancillary_evaluator(self, anc, T):
        T_r = anc['T_r']
        theta = 1-T/T_r
        RHS = sum([_n*theta**_t for _n, _t in zip(anc['n'], anc['t'])])
        if anc['using_tau_r']:
            RHS *= T_r/T
        if anc['type'] == 'rhoLnoexp':
            return anc['reducing_value']*(1+RHS)
        elif anc['type'] in ['pV','pL','rhoV','rhoL']:
            return math.exp(RHS)*anc['reducing_value']
        else:
            raise KeyError(anc['type'])

    def check_ancillaries(self, anc, name, T, backend = 'REFPROP'):
        AS = CP.AbstractState(backend, name)
        AS.update(CP.QT_INPUTS, 0, T)
        rhoL = AS.rhomolar()
        rhoV = AS.saturated_vapor_keyed_output(CP.iDmolar)
        p = AS.p()

        try:
            rhoL_anc = self.ancillary_evaluator(anc['DL'], T)
            rhoV_anc = self.ancillary_evaluator(anc['DV'], T)
            p_anc = self.ancillary_evaluator(anc['PS'], T)

            def check_ok(k, err, thresh):
                if err > thresh:
                    raise ValueError(name, k, err, thresh)
            check_ok('DL',abs(rhoL_anc-rhoL)/rhoL, 1e-3)
            check_ok('DV',abs(rhoV_anc-rhoV)/rhoV, 1e-3)
            check_ok('PS',abs(p_anc-p)/p, 1e-3)

        except BaseException as BE:
            print(BE)
            # raise

    def get_alpha0(self, Tc):
        def get_PX0_block():
            imin, imax = 0, len(self.lines)
            for i in range(100):
                try:
                    block, imin, imaxnew = self.get_block('#AUX', r'^\s*[\n\r]', imin, indices=True) # end of block is either empty line or three or more spaces
                    imin = imaxnew
                    kind = block[1].split(' ')[0]
                    if kind == 'PX0':
                        return block
                except IndexError:
                    pass
            raise IndexError("Cannot find PX0 block")
        PX0_block = get_PX0_block()
        # print(''.join(PX0_block))
        term_spec = self.get_keyed_line(PX0_block,' !Nterms',int)
        names = ['ai*log(tau**ti)','ai*tau**ti','ai*log(1-exp(bi*tau))']

        alpha0 = []

        # To get the ln(delta) term that doesn't appear in the block
        # anywhere
        alpha0.append({
            "a1": 0.0,
            "a2": 0.0,
            "type": "IdealGasHelmholtzLead"
        })

        i = lines_contains('!Nterms', PX0_block)[0]+1
        for Nterms, name in zip(term_spec, names):
            lines = PX0_block[i:i+Nterms]
            if name == 'ai*log(tau**ti)':
                if len(lines) == 1:
                    els = [el.strip() for el in lines[0].split(' ') if el]
                    if float(els[1]) != 1.0: 
                        raise ValueError("ti must be 1.0 in ai*log(tau**ti)")
                    alpha0.append({
                        "a": float(els[0]),
                        'type': 'IdealGasHelmholtzLogTau'
                    })
                else:
                    raise IndexError()
            elif name == 'ai*tau**ti':
                # This is a combination of the IdealGasHelmholtzLogTau and IdealGasHelmholtzLogTau terms in CoolProp
                n,t = [],[]
                for line in lines:
                    els = [el.strip() for el in line.split(' ') if el]
                    ai = float(els[0])
                    ti = float(els[1])
                    n.append(ai)
                    t.append(ti)
                alpha0.append({
                    "n": n,
                    "t": t,
                    "type": "IdealGasHelmholtzPower"
                })
            elif name == 'ai*log(1-exp(bi*tau))':
                n, t = [],[]
                for line in lines:
                    els = [el.strip() for el in line.split(' ') if el]
                    ni,ti = [float(el) for el in els[0:2]]
                    n.append(ni)
                    t.append(ti/Tc)
                alpha0.append({
                    "n": n, "t": t,
                    "type": "IdealGasHelmholtzPlanckEinstein"
                })
            elif name == 'Hmm':
                pass
            else:
                raise ValueError("Bad alpha0 contribution" + name)

            i += Nterms

        return alpha0

    def get_EOS(self, ideal=True):
        block = self.get_block('#EOS', r'^\s*[\n\r]') # end of block is a line with only whitespace and terminated with a newline, or a carriage return
        block = self.strip_commented(block)
        indices = lines_contains('eta      beta    gamma   epsilon', block)
        if indices:
            block = block[0:indices[0]]
        # print(''.join(block))
        kind = block[1].split(' ')[0]
        if kind != 'FEQ':
            raise ValueError(f'Kind of EOS [{kind}] is not supported; only valid options are: {{FEQ}}')

        def get_keyed_line(key,conv,*,alias=None):
            return self.get_keyed_line(block, key, conv, alias=alias)

        # Various constants. Anything without a unit suffix is in base SI units
        Tmax = get_keyed_line('!Upper temperature limit [K]',float)[0]
        Tt = get_keyed_line('!Triple point temperature',float)[0]
        acentric = get_keyed_line('!Acentric factor',float)[0]
        R = get_keyed_line('!Gas constant',float,alias='!gas constant')[0]
        molemass_kgkmol = get_keyed_line('!Molar mass',float)[0]
        pt_kPa = get_keyed_line('!Pressure at triple point',float)[0]
        pmax_kPa = get_keyed_line('!Upper pressure limit',float)[0]
        # Reducing parameters
        _Tr, pc_kPa, _rhoc_moldm3 = get_keyed_line('!Tc [K], pc [kPa], rhoc [mol/L]',float)
        Tr, rhor_moldm3 = get_keyed_line('!Reducing parameters [K, mol/L]',float)
        # Term specification
        try:
            terms = get_keyed_line(' !# terms and # coefs/term',int)
            i0 = lines_contains('!# terms and # coefs/term', block)[0]
        except KeyError:
            terms = get_keyed_line(' !First digit: # of normal terms',int)
            i0 = lines_contains(' !First digit: # of normal terms', block)[0]
            
        Npoly, Ncoefpoly, NGaussian, NcoefGaussian, NGao = terms[0:5]
        assert(sum(terms[5::]) == 0)

        pairs = []
        names = []
        if Npoly > 0:
            names.append('Polynomial')
            pairs.append((Npoly, Ncoefpoly))
        if NGaussian > 0:
            names.append('Gaussian')
            pairs.append((NGaussian, NcoefGaussian))
        if NGao > 0:
            names.append('Gao')
            pairs.append((NGao, 12))

        alphar = []
        i = i0 + 1
        for name, pair in zip(names, pairs):

            Nterms, Ncoefperterm = pair
            if Ncoefperterm == 0 and name in ['Polynomial','Gaussian']:
                continue
            lines = block[i:i+Nterms]
            if name == 'Polynomial':
                if Ncoefperterm == 4:
                    n,t,d,l = [],[],[],[]
                    for line in lines:
                        els = [el.strip() for el in line.split(' ') if el]
                        ni,ti,di,li = [float(el) for el in els[0:4]]
                        n.append(ni)
                        t.append(ti)
                        d.append(di)
                        l.append(li)
                    alphar.append({"n": n, "t": t, "d": d, "l": l,
                                   "type": "ResidualHelmholtzPower"})
                elif Ncoefperterm == 5:
                    n,t,d,l,g = [],[],[],[],[]
                    for line in lines:
                        els = [el.strip() for el in line.split(' ') if el]
                        ni,ti,di,li,gi = [float(el) for el in els[0:5]]
                        n.append(ni)
                        t.append(ti)
                        d.append(di)
                        l.append(li)
                        g.append(gi)
                    alphar.append({"n": n, "t": t, "d": d, "l": l, "g": g,
                                   "type": "ResidualHelmholtzExponential"})
                else:
                    raise ValueError()

            elif name == 'Gaussian':

                def is_normal_gaussian(lines):
                    for line in lines:
                        els = [el.strip() for el in line.split(' ') if el]
                        if not (all([float(el)==0 for el in els[9:12]]) and all([float(el)==2.0 for el in els[3:5]])):
                            return False
                    return True

                def is_R125_gaussian(lines):
                    for line in lines:
                        els = [el.strip() for el in line.split(' ') if el]
                        if not (all([float(el)==0 for el in els[9:12]]) and all([float(el)==-1.0 for el in els[5:7]])):
                            return False
                    return True

                def is_nonanalytic_gaussian(lines):
                    numgaussian = 0; numNA = 0
                    for line in lines:
                        els = [el.strip() for el in line.split(' ') if el]
                        if all([float(el)==0 for el in els[9:12]]) and all([float(el)==2.0 for el in els[3:5]]):
                            numgaussian += 1
                        if all([float(el)!=0 for el in els[9:12]]) and all([float(el)==2.0 for el in els[3:5]]):
                            numNA += 1
                    return (numNA > 0, numgaussian, numNA)

                def add_normal_gaussian(lines):
                    n,t,d,eta,beta,gamma,epsilon = [],[],[],[],[],[],[]
                    for line in lines:
                        els = [el.strip() for el in line.split(' ') if el]
                        ni,ti,di,ph1,ph2,negetai,negbetai,gammai,epsiloni,ph3,ph4,ph5 = [float(el) for el in els[0:12]]
                        n.append(ni)
                        t.append(ti)
                        d.append(di)
                        eta.append(-negetai)
                        beta.append(-negbetai)
                        gamma.append(gammai)
                        epsilon.append(epsiloni)
                    alphar.append({"n": n, "t": t, "d": d, 
                               "eta": eta, "beta": beta, "gamma": gamma, "epsilon": epsilon,
                               "type": "ResidualHelmholtzGaussian"})

                def add_nonanalytic(lines):
                    n,a,b,beta,A,B,C,D = [],[],[],[],[],[],[],[]
                    for line in lines:
                        els = [el.strip() for el in line.split(' ') if el]
                        ni,ph1,ph2,ph3,ph4,bi,betai,Ai,Ci,Di,Bi,ai = [float(el) for el in els[0:12]]
                        n.append(ni)
                        a.append(ai)
                        b.append(bi)
                        beta.append(betai)
                        A.append(Ai)
                        B.append(Bi)
                        C.append(Ci)
                        D.append(Di)
                    alphar.append({"n": n, "b": b, "a": a, "b": b, "A": A, "B": B, "C": C, "D": D, "beta": beta,
                               "type": "ResidualHelmholtzNonAnalytic"})

                if is_normal_gaussian(lines):
                    add_normal_gaussian(lines)
                elif is_R125_gaussian(lines):
                    n,t,d,l,m = [],[],[],[],[]
                    for line in lines:
                        els = [el.strip() for el in line.split(' ') if el]
                        # Normal gaussian term if last three placeholders are zero
                        ni,ti,di,li,mi,ph1,ph2,ph3,ph4,ph5,ph6,ph7 = [float(el) for el in els[0:12]]
                        n.append(ni)
                        t.append(ti)
                        d.append(di)
                        l.append(li)
                        m.append(mi)
                    alphar.append({"n": n, "t": t, "d": d, "l": l, "m": m, "type": "ResidualHelmholtzLemmon2005"})
                else:
                    isnonanalytic, NGaussian, NNA = is_nonanalytic_gaussian(lines)
                    if isnonanalytic:
                        add_normal_gaussian(lines[0:NGaussian])
                        add_nonanalytic(lines[NGaussian::])
                    else:
                        raise ValueError("I don't yet understand this kind of Gaussian:"+str(lines))
            elif name == 'Gao':
                n,t,d,eta,beta,gamma,epsilon,b = [],[],[],[],[],[],[],[]
                for line in lines:
                    els = [el.strip() for el in line.split(' ') if el]
                    ni,ti,di,ph1,ph2,etai,betai,gammai,epsiloni,bi,ph4,ph5 = [float(el) for el in els[0:12]]
                    n.append(ni)
                    t.append(ti)
                    d.append(di)
                    eta.append(etai)
                    beta.append(betai)
                    gamma.append(gammai)
                    epsilon.append(epsiloni)
                    b.append(bi)
                alphar.append({"n": n, "t": t, "d": d, 
                           "eta": eta, "beta": beta, "gamma": gamma, "epsilon": epsilon, "b": b,
                           "type": "ResidualHelmholtzGaoB"})
            else:
                raise ValueError("I don't yet understand:"+name)
                
            i += Nterms

        pr = pc_kPa*1e3

        # Get the terms
        EOS = {
          "BibTeX_CP0": "",
          "BibTeX_EOS": "XXX-X-XXXX",
          "STATES": {
            "reducing": {
              "T": Tr,
              "T_units": "K",
              "hmolar": -99999999999,
              "hmolar_units": "J/mol",
              "p": pr,
              "p_units": "Pa",
              "rhomolar": rhor_moldm3*1e3,
              "rhomolar_units": "mol/m^3",
              "smolar": 999999999999999,
              "smolar_units": "J/mol/K"
            },
            "sat_min_liquid": {
              "T": Tt,
              "T_units": "K",
              "hmolar": -999999999999,
              "hmolar_units": "J/mol",
              "p": pt_kPa*1e3,
              "p_units": "Pa",
              "rhomolar": 9999999999999999,
              "rhomolar_units": "mol/m^3",
              "smolar": 99999999999999999999,
              "smolar_units": "J/mol/K"
            },
            "sat_min_vapor": {
              "T": Tt,
              "T_units": "K",
              "hmolar": 9999999999999,
              "hmolar_units": "J/mol",
              "p": pt_kPa*1e3,
              "p_units": "Pa",
              "rhomolar": 99999999999999,
              "rhomolar_units": "mol/m^3",
              "smolar": 9999999999999999999,
              "smolar_units": "J/mol/K"
            }
          },
          "T_max": Tmax,
          "T_max_units": "K",
          "Ttriple": Tt,
          "Ttriple_units": "K",
          "acentric": acentric,
          "acentric_units": "-",
          "alpha0": self.get_alpha0(Tc=Tr) if ideal else [],
          "alphar": alphar,
          "gas_constant": R,
          "gas_constant_units": "J/mol/K",
          "molar_mass": molemass_kgkmol/1e3,
          "molar_mass_units": "kg/mol",
          "p_max": pmax_kPa*1e3,
          "p_max_units": "Pa",
          "pseudo_pure": False
        }
        return EOS

    def write_JSON(self, jsonpath, *, ideal=True, ancillaries=False):
        """ """
        f = {
            'EOS': [self.get_EOS(ideal=ideal)],
            'ANCILLARIES': self.get_ancillaries() if ancillaries else {},
            'INFO': self.get_info(),
        }
        f['STATES'] = {
            'critical': self.get_critical_state(),
            'triple_liquid': f["EOS"][0]['STATES']['sat_min_liquid'],
            'triple_vapor': f["EOS"][0]['STATES']['sat_min_vapor']
        }
        with open(jsonpath, 'w') as fp:
            fp.write(json.dumps(f, indent=2, sort_keys=True))

class HMXDeconstructor:
    """ 
    Take a REFPROP 10 compatible HMX.BNC, convert to 
    CoolProp format departure and binary interaction files
    """
    def __init__(self, path, RProot):
        self.path = path
        self.RProot = RProot
        self.lines = [line.strip() for line in codecs.open(path, encoding='cp1252').readlines()]
        self.CASlookup = self.get_CASlookup()
        self.RPlookup = self.get_RPlookup()
        self.dep = self.get_departure_functions()
        self.BIP = self.get_binary_interaction()

    def get_CASlookup(self):
        o = {}
        for FLDpath in glob.glob(os.path.join(self.RProot, 'FLUIDS', '*.FLD')):
            FLD = os.path.split(FLDpath)[1].split('.')[0]
            hash_ = None
            for line in open(FLDpath).readlines():
                if 'Hash' in line:
                    hash_ = line.split('!')[0].strip()
                    break
            for line in open(FLDpath).readlines():
                if 'CAS' in line:
                    o[hash_] = line.split('!')[0].strip()
                    break
        return o

    def get_RPlookup(self):
        o = {}
        for FLDpath in glob.glob(os.path.join(self.RProot, 'FLUIDS', '*.FLD')):
            FLD = os.path.split(FLDpath)[1].split('.')[0]
            hash_ = None
            for line in open(FLDpath).readlines():
                if 'Hash' in line:
                    hash_ = line.split('!')[0].strip()
                    break
            o[hash_] = FLD
        return o

    def get_binary_interaction(self):
        BIP = []
        iBNC = lines_equals("BNC", self.lines)[0]
        iend = iBNC
        while iend < iBNC + 10000:
            if self.lines[iend].strip() == '':
                break
            else:
                iend += 1
        for chunk in '\n'.join(self.lines[iBNC+1:iend]).split('!\n'):
            if all([line.strip().startswith('?') or not line.strip() for line in chunk.split('\n')]):
                continue
            else:
                chunklines = chunk.split('\n')
                Name1, Name2 = chunklines[0][1::].split('[')[0].strip().split('/')
                ilastcomment = lines_starts_with('?', chunklines)[-1]
                refrow = chunklines[1:ilastcomment+1]
                hash1, hash2 = [el.strip() for el in chunklines[ilastcomment+1].split('/')]
                def to_float(s):
                    try:
                        return float(s)
                    except:
                        return s
                # start at the bottom and work backwards to find the last mixture model definition
                for irow in reversed(range(len(chunklines)-1)):
                    model, betaT, gammaT, betaV, gammaV, Fij = [to_float(el) for el in chunklines[irow].split(' ') if el][0:6]
                    if model in (['PR1','ST1','TRN'] + [f'TC{i}' for i in range(1,8)] + [f'VC{i}' for i in range(1,8)]):
                        # print('XX', chunklines[irow])
                        pass
                    else:
                        # print('  ', chunklines[irow])
                        break

                try:
                    BIP.append(dict(
                        betaT=betaT,
                        gammaT=gammaT,
                        betaV=betaV,
                        gammaV=gammaV,
                        F=Fij,
                        BibTeX='?',
                        Name1=self.RPlookup[hash1],
                        Name2=self.RPlookup[hash2],
                        CAS1=self.CASlookup[hash1],
                        CAS2=self.CASlookup[hash2],
                        function=model
                    ))
                except KeyError as ke:
                    print(ke)
        return BIP

    def get_departure_functions(self):
        iMXM = lines_starts_with("#MXM", self.lines)
        departure_functions = []

        # This is a placeholder for XR0 models
        departure_functions.append(dict(
            type='GERG-2008',
            n = [],
            t = [],
            d = [],
            l = [],
            eta = [],
            beta = [],
            gamma = [],
            epsilon = [],
            Npower=0,
            Name='XR0',
            aliases=[],
            BibTeX='??'
        ))

        for istart in iMXM:
            iend = istart
            while iend < istart + 30:
                if self.lines[iend].strip() == '':
                    break
                else:
                    iend += 1
            blocklines = self.lines[istart:iend]
            code = blocklines[1].split(' ')[0]
            if code in ['XR0','LIN','TR1']:
                continue
            istart = lines_contains('Descriptors for ', blocklines)[0] + 2 # index of row with spec
            specnumbers = [int(el) for el in blocklines[istart].split("!")[0].split(' ') if el]
            Npoly, Mpoly, dummy, Nspecial, Mspecial, NGaussian, MGaussian = specnumbers[0:7]
            def split_line(line):
                return [float(el.replace('d','e')) for el in line.split('!')[0].split()]

            if Nspecial > 0:
                # These GERG terms are combinations of polynomial and special terms
                Nlines = Npoly + Nspecial
                if Mpoly == 4:
                    blockpoly = blocklines[istart+1:istart+1+Npoly]
                    npoly, tpoly, dpoly, lpoly = zip(*[split_line(line) for line in blockpoly])
                    assert(len(npoly)==Npoly)
                else:
                    if Mpoly != 0:
                        raise ValueError(Mpoly)

                if Mspecial == 7:
                    blockspe = blocklines[istart+1+Npoly: istart+1+Npoly+Nspecial]
                    nspe, tspe, dspe, eta, epsilon, beta, gamma  = zip(*[split_line(line) for line in blockspe])
                else:
                    nspe = [0.0]*Nspecial
                    tspe = [0.0]*Nspecial
                    dspe = [0.0]*Nspecial
                    eta = [0.0]*Nspecial
                    beta = [0.0]*Nspecial
                    gamma = [0.0]*Nspecial
                    epsilon = [0.0]*Nspecial
                    if Mspecial != 0:
                        raise ValueError(Mspecial)

                departure_functions.append(dict(
                    type='GERG-2008',
                    n = list(npoly) + list(nspe),
                    t = list(tpoly) + list(tspe),
                    d = list(dpoly) + list(dspe),
                    l = list(lpoly) + [0.0]*Nspecial,
                    eta = [0.0]*Npoly + list(eta),
                    beta = [0.0]*Npoly + list(beta),
                    gamma = [0.0]*Npoly + list(gamma),
                    epsilon = [0.0]*Npoly + list(epsilon),
                    Npower=Npoly,
                    Name=code,
                    aliases=[],
                    BibTeX='??'
                ))

            else:
                # Polynomial and Gaussian terms
                Nlines = Npoly + NGaussian
                blockpoly = blocklines[istart+1 : istart+1+Npoly]
                if Mpoly == 4:
                    npoly, tpoly, dpoly, lpoly = zip(*[split_line(line) for line in blockpoly])
                elif Mpoly == 3:
                    npoly, tpoly, dpoly = zip(*[split_line(line) for line in blockpoly])
                    lpoly = [0.0]*len(dpoly)
                else:
                    if Mpoly != 0:
                        raise ValueError(Mpoly)

                if MGaussian == 12:
                    blockgau = blocklines[istart+1+Npoly:istart+1+Npoly+NGaussian+1]
                    ngau, tgau, dgau, dum1gau, dum2gau, eta, beta, gamma, epsilon, dum3, dum4, dum5  = zip(*[split_line(line) for line in blockgau])
                else:
                    ngau = [0.0]*NGaussian
                    tgau = [0.0]*NGaussian
                    dgau = [0.0]*NGaussian
                    eta = [0.0]*NGaussian
                    beta = [0.0]*NGaussian
                    gamma = [0.0]*NGaussian
                    epsilon = [0.0]*NGaussian
                    if MGaussian != 0:
                        raise ValueError(MGaussian)

                departure_functions.append(dict(
                    type='Gaussian+Exponential',
                    n=list(npoly) + list(ngau),
                    t=list(tpoly) + list(tgau),
                    d=list(dpoly) + list(dgau),
                    l=list(lpoly) + [0.0]*NGaussian,
                    eta=[0.0]*Npoly + (-np.array(eta)).tolist(),
                    beta=[0.0]*Npoly + (-np.array(beta)).tolist(),
                    gamma=[0.0]*Npoly + list(gamma),
                    epsilon=[0.0]*Npoly + list(epsilon),
                    Npower=Npoly,
                    Name=code,
                    aliases=[],
                    BibTeX='??'
                ))

        return departure_functions

    def write(self, *, base):
        with open(base + '/mixture_departure_functions.json','w') as fp:
            fp.write(json.dumps(self.dep, indent=2))
        with open(base + '/mixture_binary_pairs.json','w') as fp:
            fp.write(json.dumps(self.BIP, indent=2))

HMX_header = """HMX               !Mnemonic for mixture model, must match hfmix on call to SETUP.
4                 !Version number

! Changelog:
! ---------

#BNC              !Binary mixing coefficients
BNC
? Binary mixing coefficients for the various mixing rules used with the HMX model:
?
? KWi:  (i = 1,2,3,...,A,B,...)  --->  Kunz-Wagner mixing rules
?   model     BetaT     GammaT    BetaV     GammaV    Fij      not used
?
!"""

HMX_footer = """
@END
c        1         2         3         4         5         6         7         8
c2345678901234567890123456789012345678901234567890123456789012345678901234567890
"""

BIN_template = """
?{name1}/{name2}                        [{name1}/{name2}]
?FITTED {annotation}
  {hash1}/{hash2}
    {model}     {betaT:14.11f} {gammaT:14.11f} {betaV:14.11f} {gammaV:14.11f} {Fij:14.11f}  0.             0. 0. 0. 0. 0. 0.
    TC5    0.0      0.0     0.0       0.             0.             0.             0. 0. 0. 0. 0. 0.
    VC5    0.0      0.0     0.0       0.             0.             0.             0. 0. 0. 0. 0. 0.
!"""

MODEL_template = """#MXM              !Mixture model specification
{model} {annotation}
?
!```````````````````````````````````````````````````````````````````````````````
 BetaT    GammaT   BetaV    GammaV    Fij    not used      !Descriptors for binary-specific parameters
  1.0      1.0      1.0      1.0      0.0      0.0         !Default values (i.e. ideal-solution)
  {Nexp} {Ntermsexp}      0        {NKW} {NtermsKW}      {NGaussian} {NtermsGaussian}      0 0      0 0         !# terms and # coefs/term for normal terms, Kunz-Wagner terms, and Gaussian terms.  3rd column is not used.
  """

class HMXBuilder:
    """
    Read in CoolProp-format departure term JSON structure and build a REFPROP-format HMX.BNC
    """
    def __init__(self, jbin, jdep):

        """ Determine acceptable names for functions """
        def get_function_names():
            function_names = set([el['Name'] for el in jdep])
            return function_names
        function_names = get_function_names()
        if len(function_names) > 100:
            raise ValueError("Max of 100 functions allowed")
        RP_function_names = {f:f'B{i:2d}'.replace(' ','0') for i,f in enumerate(function_names)}
        def get_hash(CAS, name):
            try:
                # Use CoolProp to obtain the standard inchi key
                AS = CP.AbstractState('HEOS', CAS)
                StdInChIkey = CP.get_fluid_param_string(CAS, 'INCHIKEY')
                uid = hashlib.sha256(StdInChIkey.encode('UTF-8')).hexdigest()[2:9] + '0'
            except ValueError:
                # Use REFPROP to do so
                import ctREFPROP.ctREFPROP as c
                RP = c.REFPROPFunctionLibrary(os.getenv('RPPREFIX'))
                r = RP.REFPROPdll(name,'','HASH',RP.MOLAR_BASE_SI,0,0,0,0,[1.0])
                assert(r.ierr == 0)
                uid = r.hUnits
            return uid

        out = HMX_header

        # The top part with the definitions of what model to use for each binary pair and the 
        # model name
        for el in jbin:
            if 'function' not in el:
                RP_function_number = 'BA0'
            else:
                if el['function'] == 'XR0':
                    RP_function_number = 'XR0'
                else:
                    RP_function_number = RP_function_names[el['function']]
            if 'xi' in el and 'zeta' in el:
                # convert values to gammaT and gammaV
                # yadda...
                print("TODO: Need to convert values to gammaT and gammaV......................................")
                # values = 1.0, el['gammaT'], 1.0, el['gammaV'], el['F']
                continue
            else:
                template_values = {
                    'betaT': el['betaT'], 
                    'gammaT': el['gammaT'], 
                    'betaV': el['betaV'], 
                    'gammaV': el['gammaV'], 
                    'Fij': el['F'],
                    'name1': el['Name1'],
                    'name2': el['Name2'],
                    'annotation': 'time/date',
                    'hash1': get_hash(el['CAS1'], el['Name1']),
                    'hash2': get_hash(el['CAS2'], el['Name2']),
                    'model': RP_function_number
                }
                entry = BIN_template.format(**template_values)
                # print(entry)
                out += entry

        out += '\n\n'

        for el in jdep:
            # Build the departure term
            # print(el)
            model = RP_function_names[el["Name"]]
            annotation = f'{el["Name"]} '

            Nexp = Ntermsexp = NKW = NtermsKW = NGaussian = NtermsGaussian = 0
            if el['type'] == 'Exponential':
                Nexp = len(el['n'])
                Ntermsexp = 4 if 'l' in el and isinstance(el['l'], list) else 3
                
                # print(Nexp, Ntermsexp)
                if Ntermsexp == 4:
                    n, t, d, l = el['n'], el['t'], el['d'], el['l']
                    rows = []
                    for i in range(len(t)):
                        rows.append(f'{n[i]} {t[i]} {d[i]:0.10f} {l[i]:0.10f}')
                        if i == 0:
                            rows[-1] += f' ! n(i),t(i),d(i),l(i) in term n_i*tau^t_i*delta^d_i*exp(-delta^l_i)'
                elif Ntermsexp == 3:
                    n, t, d = el['n'], el['t'], el['d']
                    first_row = f'{n[0]} {t[0]} {d[0]:0.1f} ! n(i),t(i),d(i) in term n_i*tau^t_i*delta^d_i'
                    rows = []
                    for i in range(len(t)):
                        rows.append(f'{n[i]} {t[i]} {d[i]:0.10f}')
                        if i == 0:
                            rows[-1] += f' ! n(i),t(i),d(i) in term n_i*tau^t_i*delta^d_i'
                else:
                    raise ValueError()

            elif el['type'] == 'GERG-2008':
                # print(el)
                Nexp = el['Npower']
                Ntermsexp = 3
                assert('l' not in el)
                NKW = len(el['n'])-Nexp
                NtermsKW = 7
                
                n, t, d, eta, epsilon, beta, gamma = el['n'], el['t'], el['d'], el['eta'], el['epsilon'], el['beta'], el['gamma']
                l = None
                rows = []
                if Nexp > 0:
                    for i in range(Nexp):
                        rows.append(f'{n[i]} {t[i]} {d[i]:0.10f} ')
                        if i == 0:
                            rows[-1] += '! n(i),t(i),d(i) in term n_i*tau^t_i*delta^d_i'
                if NKW > 0:
                    for i in range(Nexp,len(t)):
                        rows.append(f'{n[i]} {t[i]} {d[i]:0.10f} {eta[i]} {epsilon[i]} {beta[i]} {gamma[i]} ')
                        if i == Nexp:
                            rows[-1] += '! n(i),t(i),d(i),eta(i),epsilon(i),beta(i),gamma(i) in term n_i*tau^t_i*delta^d_i*exp(-eta*(delta-epsilon)^2-beta*(delta-gamma))'

            elif el['type'] == 'Gaussian+Exponential':
                # print(el)
                Nexp = el['Npower']
                NGaussian = len(el['n'])-Nexp
                
                n, t, d, l, eta, epsilon, beta, gamma = el['n'], el['t'], el['d'], el['l'], el['eta'], el['epsilon'], el['beta'], el['gamma']
                rows = []
                if Nexp > 0:
                    Ntermsexp = 4
                    for i in range(Nexp):
                        rows.append(f'{n[i]} {t[i]} {d[i]:0.16f} {l[i]:0.16f} ')
                        if i == 0:
                            rows[-1] += '! n(i),t(i),d(i),l(i) in term n_i*tau^t_i*delta^d_i*exp(-sgn(l_i)*delta^l_i)'
                if NGaussian > 0:
                    NtermsGaussian = 12
                    for i in range(Nexp,len(t)):
                        negetai = -eta[i]
                        negbetai = -beta[i]
                        rows.append(f'{n[i]} {t[i]} {d[i]:0.16f} 2.0 2.0 {negetai} {negbetai} {gamma[i]} {epsilon[i]} 0.0 0.0 0.0 0.0')
                        if i == Nexp:
                            rows[-1] += '! n(i),t(i),d(i),_,_,eta(i),beta(i),gamma(i),epsilon(i),_,_,_,_ in term n_i*tau^t_i*delta^d_i*exp(eta*(delta-epsilon)^2+beta*(tau-gamma)^2)'
            else:
                raise KeyError(el['type'])
                
            out += MODEL_template.format(**locals())
            out += '\n'.join(rows) + '\n\n'

        out += '\n\n' + HMX_footer

        self.out = out

    def get(self):
        return self.out

def to_teqp(RProot, outroot):
    os.makedirs(os.path.join(outroot,'dev','mixtures'), exist_ok=True)
    os.makedirs(os.path.join(outroot,'dev','fluids'), exist_ok=True)
    hmx = HMXDeconstructor(RProot+'/FLUIDS/HMX.BNC', RProot)
    hmx.write(base=outroot+'/dev/mixtures')
    for fld in glob.glob(RProot+'/FLUIDS/*.FLD'):
        name = os.path.split(fld)[1].split('.')[0]
        try:
            FLD = FLDDeconstructor(fld).write_JSON(outroot+f'/dev/fluids/{name}.json')
        except BaseException as BE:
            print(fld, BE)
            # raise
            pass

def test_CoolProp(root):
    CP.set_config_bool(CP.OVERWRITE_FLUIDS, True)
    for path in glob.glob(os.path.join(root,'dev','fluids','*.json')):
        FLD = os.path.split(path)[1].split('.')[0]
        CP.add_fluids_as_JSON('HEOS', open(path).read())
        AS = CP.AbstractState('HEOS', FLD)
    CP.set_departure_functions(open(root+'/dev/mixtures/mixture_departure_functions.json').read())
    CP.set_interaction_parameters(open(root+'/dev/mixtures/mixture_binary_pairs.json').read())

def get_all_ancillaries():
    ancillaries = {}
    for fld in glob.glob('C:/Program Files (x86)/REFPROP/FLUIDS/*.FLD'):
        name = os.path.split(fld)[1].split('.')[0]
        FLD = FLDDeconstructor(fld)
        print(fld)
        FLD.get_EOS()
        anc = FLD.get_all_ancillaries()
        T = anc['DV']['T_r']*0.8
        FLD.check_ancillaries(anc, name, T)
        ancillaries[name] = anc

    with open('all_ancillaries.json','w') as fp:
        fp.write(json.dumps(ancillaries, indent=2))

if __name__ == '__main__':
    RProot = os.getenv('RPPREFIX')
    to_teqp(RProot, 'teqp')
    test_CoolProp('teqp')

    import ctREFPROP.ctREFPROP as c, os
    RP = c.REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    print(RP.RPVersion())