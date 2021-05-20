import codecs
import re
import json
import glob
import os
import math

import CoolProp.CoolProp as CP

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

    def lines_starts_with(self, s, lines, *, start_index=0):
        matches = []
        for i in range(start_index, len(lines)):
            if lines[i].startswith(s):
                matches.append(i)
        return matches


    def lines_contains(self, s, lines, *, start_index=0):
        matches = []
        for i in range(start_index, len(lines)):
            if s in lines[i]:
                matches.append(i)
        return matches

    def lines_equals(self, s, lines, *, start_index=0):
        matches = []
        for i in range(start_index, len(lines)):
            if re.fullmatch(s, lines[i]):
                matches.append(i)
        return matches

    def get_block(self, startreg, endreg, start_index=0, indices=False):
        matches = self.lines_starts_with(startreg, self.lines, start_index = start_index)
        istart = matches[0]
        matches = self.lines_equals(endreg, self.lines, start_index = istart+1)
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
            n.append(float(_n))
            t.append(float(_t))
        assert(len(n) == Ncoeffs[0] and len(t) == Ncoeffs[0])
        desc = self.get_ancillary_description(model_key.strip())
        if model_key[0:2] == 'PS':
            type_key = 'pL'
        elif model_key in ['DL1','DL2']:
            type_key = 'rhoLnoexp'
        elif model_key[0:2] == 'DL':
            type_key = 'rhoL'
        elif model_key[0:2] == 'DV':
            type_key = 'rhoV'
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

    def get_info(self):
        CHEMSPIDER_ID = -1
        name = os.path.split(self.path)[1].split('.')[0]
        return {
            "2DPNG_URL": f"http://www.chemspider.com/ImagesHandler.ashx?id={CHEMSPIDER_ID}",
            "ALIASES": [],
            "CAS": self.get_keyed_line(self.lines, '!CAS number', lambda x: x)[0],
            "CHEMSPIDER_ID": CHEMSPIDER_ID,
            "ENVIRONMENTAL": {
              "ASHRAE34": "?",
              "FH": 0,
              "GWP100": 0.0,
              "GWP20": 0.0,
              "GWP500": 0.0,
              "HH": 0,
              "Name": "",
              "ODP": -1e30,
              "PH": 0
            },
            "FORMULA": self.get_keyed_line(self.lines, '!Chemical formula', lambda x: x)[0],
            "INCHI_KEY": self.get_keyed_line(self.lines, '!Standard InChI Key', lambda x: x)[0],
            "INCHI_STRING": self.get_keyed_line(self.lines, '!Standard InChI String', lambda x: x)[0],
            "NAME": name,
            "REFPROP_NAME": name,
            "SMILES": "?"
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

    def check_ancillaries(self, ancillaries, name, T, backend = 'REFPROP'):
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
        print(''.join(PX0_block))
        term_spec = self.get_keyed_line(PX0_block,' !Nterms',int)
        names = ['ai*log(tau**ti)','ai*tau**ti','ai*log(1-exp(bi*tau))']

        alpha0 = []
        i = self.lines_contains('!Nterms', PX0_block)[0]+1
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
                a1 = None # constant term
                a2 = None # term multiplying tau
                if len(lines) != 2:
                    raise ValueError("Don't understand what to do when more than two terms")
                for line in lines:
                    els = [el.strip() for el in line.split(' ') if el]
                    ai = float(els[0])
                    ti = float(els[1])
                    if ti == 0.0: 
                        a1 = ai
                    elif ti == 1.0: 
                        a2 = ai
                    else:
                        raise ValueError(els, ai, ti)
                assert(a1 is not None)
                assert(a2 is not None)
                alpha0.append({
                    "a1": a1,
                    "a2": a2,
                    "type": "IdealGasHelmholtzLead"
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

    def get_EOS(self):
        block = self.get_block('#EOS', r'^\s*[\n\r]') # end of block is a line with only whitespace and terminated with a newline, or a carraige return
        block = self.strip_commented(block)
        indices = self.lines_contains('eta      beta    gamma   epsilon', block)
        if indices:
            block = block[0:indices[0]]
        # print(''.join(block))
        kind = block[1].split(' ')[0]
        if kind != 'FEQ':
            return None

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
        Tr, rhor_moldm3 = get_keyed_line('!Reducing parameters',float)
        # Term specification
        terms = get_keyed_line(' !# terms and # coefs/term',int)
        def chunkify(y, *, n):
            return [y[i*n:(i+1)*n] for i in range((len(y)+ n - 1)//n)]
        pairs = chunkify(terms, n=2)
        names = ['Polynomial','Gaussian','Hmm','Hmm','Hmm','Hmm']
        alphar = []
        i = self.lines_contains('!# terms and # coefs/term', block)[0]+1
        for name, pair in zip(names, pairs):
            Nterms, Ncoefperterm = pair
            lines = block[i:i+Nterms]
            if name == 'Polynomial':
                n,t,d,l = [],[],[],[]
                assert(Ncoefperterm==4)
                for line in lines:
                    els = [el.strip() for el in line.split(' ') if el]
                    ni,ti,di,li = [float(el) for el in els[0:4]]
                    n.append(ni)
                    t.append(ti)
                    d.append(di)
                    l.append(li)
                alphar.append({"n": n, "t": t, "d": d, "l": l,
                               "type": "ResidualHelmholtzPower"})
            elif name == 'Gaussian':
                n,t,d,eta,beta,gamma,epsilon = [],[],[],[],[],[],[]
                assert(Ncoefperterm==12)
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
                
            i += Nterms

        pr = 1e30

        # Get the terms
        EOS = {
          "BibTeX_CP0": "",
          "BibTeX_EOS": "XXX-X-XXXX",
          "STATES": {
            "reducing": {
              "T": Tr,
              "T_units": "K",
              "hmolar": -172.39792903434846,
              "hmolar_units": "J/mol",
              "p": pr,
              "p_units": "Pa",
              "rhomolar": rhor_moldm3*1e3,
              "rhomolar_units": "mol/m^3",
              "smolar": 89.79140692970108,
              "smolar_units": "J/mol/K"
            },
            "sat_min_liquid": {
              "T": Tt,
              "T_units": "K",
              "hmolar": -4851.143195734251,
              "hmolar_units": "J/mol",
              "p": pt_kPa*1e3,
              "p_units": "Pa",
              "rhomolar": 35465.24383230989,
              "rhomolar_units": "mol/m^3",
              "smolar": 53.11037416000078,
              "smolar_units": "J/mol/K"
            },
            "sat_min_vapor": {
              "T": Tt,
              "T_units": "K",
              "hmolar": 1689.054333680143,
              "hmolar_units": "J/mol",
              "p": pt_kPa*1e3,
              "p_units": "Pa",
              "rhomolar": 101.4989524639359,
              "rhomolar_units": "mol/m^3",
              "smolar": 131.15006869148877,
              "smolar_units": "J/mol/K"
            }
          },
          "T_max": Tmax,
          "T_max_units": "K",
          "Ttriple": Tt,
          "Ttriple_units": "K",
          "acentric": acentric,
          "acentric_units": "-",
          "alpha0": self.get_alpha0(Tc=Tr),
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

    def write_JSON(self, jsonpath):
        """ """
        f = {
            'EOS': [self.get_EOS()],
            'ANCILLARIES': self.get_ancillaries(),
            'INFO': self.get_info(),
        }
        f['STATES'] = {
            'critical': f["EOS"][0]['STATES']['reducing'],
            'triple_liquid': f["EOS"][0]['STATES']['sat_min_liquid'],
            'triple_vapor': f["EOS"][0]['STATES']['sat_min_vapor']
        }
        with open(jsonpath, 'w') as fp:
            fp.write(json.dumps(f, indent=2))

if __name__ == '__main__':

    for FLD in '3METHYLPENTANE','23DIMETHYLBUTANE','22DIMETHYLBUTANE':
        FLDDeconstructor(f'{FLD}.FLD').write_JSON(f'{FLD}.json')
    quit()

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