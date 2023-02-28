import os
import glob
import subprocess
import json

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

class AncillaryFitter:
    
    def __init__(self, *, FLD, prop, RP):
        """
        FLD: the FLD name for REFPROP
        RP: the REFPROP instance
        prop: What to fit, options are {'rhoL','rhoV','pS'}
        """
        self.FLD = FLD
        Tt, Tc, Tr, Tmin, rho_c, p_c, rhor = RP.REFPROPdll('','','TTRP;TC;TRED;TMIN;DC;PC;DRED', RP.MOLAR_BASE_SI, 0,0,0,0,[1.0]).Output[0:7]
        self.Ts = np.linspace(Tmin, 0.99*Tc)
        self.pS = np.array([RP.REFPROPdll('','TQ','P',RP.MOLAR_BASE_SI,0,0,T,0,[1.0]).Output[0] for T in self.Ts])
        self.rhoL = np.array([RP.REFPROPdll('','TQ','D',RP.MOLAR_BASE_SI,0,0,T,0,[1.0]).Output[0] for T in self.Ts])
        self.rhoV = np.array([RP.REFPROPdll('','TQ','D',RP.MOLAR_BASE_SI,0,0,T,1,[1.0]).Output[0] for T in self.Ts])
        T_r = Tc if Tr==Tc else Tr
        rho_r = rho_c if Tr==Tc else rhor
        self.theta = 1-self.Ts/T_r

        using_tau_r = False
        if prop in ['pS','pL','pV']:
            self.LHS = np.log(self.pS/p_c)*self.Ts/T_r
            using_tau_r = True
            reducing_value = p_c
            description = 'P=Pc*exp[sum(Ni*Theta^ti)*Tred/T]; Theta=1-T/Tred'
        elif prop in ['rhoL']:
            self.LHS = (self.rhoL/rho_r)-1
            reducing_value = rho_r
            description = 'D=Dred*[1+sum(Ni*Theta^ti)]; Theta=1-T/Tred'
        elif prop in ['rhoV']:
            self.LHS = np.log(self.rhoV/rho_r)
            reducing_value = rho_r
            description = 'D=Dred*exp[sum(Ni*Theta^ti)]; Theta=1-T/Tred'
        else:
            raise KeyError(prop)

        self.base_anc = {
            'T_r': T_r,
            'type': prop + ('noexp' if prop =='rhoL' else ''),
            'description': description,
            'using_tau_r': using_tau_r,
            'reducing_value': reducing_value,
            'Tmin': Tt,
            'Tmax': T_r
        }
        # print(self.LHS)

    def plot_values(self):
        plt.plot(self.theta, self.LHS)
        plt.gca().set(xlabel=r'$\theta=1-T/T_{\rm red}$', ylabel='LHS')
        plt.show()

    def evaluator(self, *, anc, T):
        """
        # Most are in this form:
        'PS5':  'P=Pc*EXP[SUM(Ni*Theta^ti)*Tc/T]',
        'DL1':  'D=Dc*[1+SUM(Ni*Theta^ti)]',
        'DV3':  'D=Dc*EXP[SUM(Ni*Theta^ti)]',
        """
        T_r = anc['T_r']
        theta = 1-T/T_r
        RHS = sum([_n*theta**_t for _n, _t in zip(anc['n'], anc['t'])])
        if anc['using_tau_r']:
            RHS *= T_r/T
        
        if anc['type'] == 'rhoLnoexp':
            return anc['reducing_value']*(1+RHS)
        elif anc['type'] in ['pV','pL','pS','rhoV','rhoL']:
            return np.exp(RHS)*anc['reducing_value']
        else:
            raise KeyError(anc['type'])

    def get_ancillary(self, *, n, t):
        anc = self.base_anc.copy()
        anc['n'] = np.array(n).tolist()
        anc['t'] = np.array(t).tolist()
        return anc
    
    def get_n(self, t):
        """
        This function takes the exponents t_i and does the 
        least squares solve for the coefficients n_i
        """
        A = np.zeros((len(self.theta), len(t)))
        for i in range(len(t)):
            A[:, i] = self.theta**t[i]
        n = np.linalg.lstsq(A, self.LHS, rcond=None)[0]
        return n, A

    def split_objective(self, params):
        n, A = self.get_n(t=params)
        predicted = np.dot(A, n)
        err = (predicted - self.LHS)/self.LHS
        val = np.sum(err**2)**0.5
        return val
    
    def optimize(self, bounds):
        return scipy.optimize.differential_evolution(self.split_objective, bounds=bounds)

def ancillary_injecter(*, RP, FLD, teqp, plot=True):
    ierr = RP.SETFLUIDSdll(FLD)
    if ierr > 0:
        print(RP.ERRMSGdll(ierr))
        quit()

    with open(teqp) as fp:
        j = json.load(fp)

    if 'rhoL' in j['ANCILLARIES']:
        print('Ancillaries already exist for FLD')
        return

    for propkey in ['pS', 'rhoL', 'rhoV']:
        anc = AncillaryFitter(FLD=FLD, prop=propkey, RP=RP)
        # anc.plot_values()
        opt = scipy.optimize.differential_evolution(anc.split_objective, bounds=[(0.1, 5) for i in range(5)])
        n, A = anc.get_n(t=opt.x)
        err = np.dot(A, n)-anc.LHS
        a = anc.get_ancillary(n=n, t=opt.x)
        y = np.array([anc.evaluator(anc=a, T=T_) for T_ in anc.Ts])
        if plot:
            plt.plot(anc.Ts, getattr(anc, propkey)/y-1)
            plt.title(f'{propkey} / {os.path.basename(FLD)}')
            plt.show()
        j['ANCILLARIES'][propkey] = a
        
    with open(teqp, 'w') as fp:
        fp.write(json.dumps(j,indent=2))