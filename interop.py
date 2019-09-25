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

    def lines_starts_with(self, s, lines, *, start_index=0):
        matches = []
        for i in range(start_index, len(lines)):
            if lines[i].startswith(s):
                matches.append(i)
        return matches

    def lines_equals(self, s, lines, *, start_index=0):
        matches = []
        for i in range(start_index, len(lines)):
            if re.fullmatch(s, lines[i]):
                matches.append(i)
        return matches

    def get_block(self, startrow, endrow):
        matches = self.lines_starts_with(startrow, self.lines, start_index = 0)
        istart = matches[0]
        matches = self.lines_equals(endrow, self.lines, start_index = istart+1)
        iend = matches[0]
        return self.lines[istart:iend]

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

if __name__ == '__main__':
    ancillaries = {}
    for fld in glob.glob('D:/Program Files (x86)/REFPROP/FLUIDS/*.FLD'):
        name = os.path.split(fld)[1].split('.')[0]
        FLD = FLDDeconstructor(fld)
        anc = FLD.get_all_ancillaries()
        T = anc['DV']['T_r']*0.8
        FLD.check_ancillaries(anc, name, T)
        ancillaries[name] = anc

    with open('all_ancillaries.json','w') as fp:
        fp.write(json.dumps(ancillaries, indent=2))