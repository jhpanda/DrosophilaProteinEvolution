""" 
    Asymptotic MK:
    alpha(x) = 1 - Ds/Dn * Pn(x)/Ps(x)
    Fit alpha(x) to alpha(x) = a(1-exp(-b*x))+c; constraints: b>0,alpha(x)<=1
    Reference: 
        1. Uricchio, Petrov and Enard, 2019, Nat Ecol Evol
        2. Messer and Petrov, 2013, PNAS
        3. Haller and Messer, 2017, G3
    Added by Junhui Peng, 31/01/2020
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from lmfit import Minimizer, Parameters

### fitting functions ###
def lmfunc(pars, x, y=None):
    a, b, c = pars['a'], pars['b'], pars['c']
    #model = a + b*np.exp(c*x)
    model = a*np.exp(-b*x)+c
    if y is None:
        return model
    #return np.sqrt(abs(model-y))
    return model-y

def lmdfunc(pars, x, y=None):
    a,b = pars['a'], pars['b']
    exp_v = np.exp(-b*x)
    da = exp_v
    db = -a*x*exp_v
    dc = np.ones(len(x))
    return np.array([da,db,dc])

def expfit(x_tofit,y_tofit,c_mka):
    ### c_mka is used as the initial value of alpha ###
    pars = Parameters()
    pars.add('alpha',value=0.,min=0.,max=1.,vary=True)
    #pars.add('c', value=0.,min=0.,max=1.)
    pars.add('a', value=0.,min=-10,max=10)
    pars.add('b', value=.1,min=0.,max=10)
    pars.add('c', expr='alpha-a*exp(-b)')
    #pars.add('a', expr='(alpha-c)/exp(-b)')
    mfit = Minimizer(lmfunc,pars,fcn_args=(x_tofit,),fcn_kws={'y':y_tofit},reduce_fcn='neglogcauchy')
    nfit = mfit.minimize(method='nedeler')
    dfit = mfit.minimize(method='leastsq',params=nfit.params)

    #dfit = mfit.leastsq(Dfun=lmdfunc, col_deriv=1)
    #dfit = mfit.leastsq()

    a = dfit.params['a']
    b = dfit.params['b']
    c = dfit.params['c']

    alpha = lmfunc(dfit.params,1.)
    alpha1 = dfit.params['alpha']
    std = alpha1.stderr
    if std==None:
        std = 0.05*alpha
    #print(alpha,alpha1,std)
    expr  = 'Fitted line:y = %f*exp(-%f*x)+%f'%(a,b,c)
    #chisqr = dfit.redchi
    chisqr = dfit.chisqr
    ndf = len(x_tofit-2)
    pvalue = stats.chi2.cdf(chisqr,ndf)

    return alpha,chisqr,pvalue,expr,dfit,std


def fit(Pn,Ps,polymseq,nallele,lcut=0.1,hcut=0.9,debug=False,plot=False):
        ### to fit from Pn_sfs and Ps_sfs ###
        #Dn = Pn[-1]
        #Ds = Ps[-1]
        Dn = polymseq[1]
        Ds = polymseq[3]
        n_allele = nallele
 
        x_tofit = []
        y_tofit = []
        xp      = [s/n_allele for s in range(n_allele+1)]
        if Dn and np.sum(Ps):
            ds_over_dn  = Ds/Dn
            #print(Dn,Ds,ds_over_dn)
            #if debug:
            #    print('Dn:',Dn,'Ds:',Ds,'Ds/Dn:',ds_over_dn)
            lcut = int(lcut*n_allele)
            hcut = int(hcut*n_allele)
 
            N = 0
            totPn = np.sum(Pn)
            totPs = np.sum(Ps)
            scale = totPn/totPs
            c_mka = 1-ds_over_dn*scale
            for i in range(lcut,hcut):
                pn_x = np.sum(Pn[i:-1])
                ps_x = np.sum(Ps[i:-1])
                #pn_x = Pn[i]
                #ps_x = Ps[i]
                #if pn_x>0 and ps_x>0:
                if pn_x>0 and ps_x>0:
                    if Pn[i]:
                    #print(pn_x,ps_x)
                        y = 1 - ds_over_dn*pn_x/ps_x
                        x_tofit += [xp[i]]
                        y_tofit += [y]
                        N += 1
                #else:
                #    #if Pn[i] or Ps[i]:
                #    y = c_mka
                #    x_tofit += [xp[i]]
                #    y_tofit += [y]
                #    N += 1
 
            x_tofit = np.array(x_tofit)
            y_tofit = np.array(y_tofit)
            if debug:
                #print(Pn)
                #print(Ps)
                ndata = len(x_tofit)
                print('Ds:',Ds,'Dn:',Dn,'Pn:',totPn,'Ps:',totPs)
                print('Ds/Dn:',ds_over_dn,'Pn/Ps:',scale,'ndata',ndata)
                print('classical MK alpha:',1-ds_over_dn*scale)

            if N>=3:
                alpha,chisqr,pvalue,expr,dfit,std = expfit(x_tofit,y_tofit,c_mka)
            elif N<3:
                alpha = np.mean(y_tofit)
                chisqr = 0
                pvalue = 0
                expr   = 'mean'
                dfit   = None
                std    = np.std(alpha)
            if debug:
                print(expr)
                print('Alpha:%5.3f Chi:%5.3f P:%5.3f'%(alpha,chisqr,pvalue))

            if plot:
                xfit = np.arange(0.1,0.95,0.05)
                if N>=3:
                    yfit = lmfunc(dfit.params, xfit)
                else:
                    yfit = np.array([alpha for s in xfit])
                fig,ax = plt.subplots()
                ax.scatter(x_tofit,y_tofit)
                #print(np.max(x_tofit),np.min(x_tofit))
                #ax.plot(x_tofit,expfunc(x_tofit,*popt),'r-')
                ax.plot(xfit,yfit,'b-')
                #ax.plot(x_tofit,fit1,'k-')
                plt.show()
                plt.close()
        else:
            alpha  = -99.99
            pvalue = -99.99 
            c_mka  = -99.99

        return alpha,std,c_mka

if __name__ == '__main__':
    Pn = [0, 0, 191710, 0, 72130, 0, 39715, 0, 27375, 0, 19415, 0, 14710, 0, 11990, 0, 9760, 0, 8360, 0, 7290, 0, 6785, 0, 5955, 0, 5475, 0, 4600, 0, 4370, 0, 4280, 0, 3665, 0, 3470, 0, 3155, 0, 3020, 0, 2865, 0, 2930, 0, 2465, 0, 2315, 0, 2400, 0, 2015, 0, 2045, 0, 1890, 0, 1770, 0, 1990, 0, 1605, 0, 1760, 0, 1610, 0, 1750, 0, 1495, 0, 1425, 0, 1420, 0, 1505, 0, 1460, 0, 1305, 0, 1205, 0, 1210, 0, 1255, 0, 1180, 0, 1175, 0, 1105, 0, 1110, 0, 1015, 0, 1110, 0, 1005, 0, 1115, 0, 1050, 0, 965, 0, 880, 0, 1015, 0, 915, 0, 865, 0, 865, 0, 845, 0, 895, 0, 890, 0, 865, 0, 835, 0, 895, 0, 840, 0, 855, 0, 775, 0, 755, 0, 765, 0, 695, 0, 725, 0, 775, 0, 735, 0, 765, 0, 710, 0, 730, 0, 865, 0, 720, 0, 605, 0, 720, 0, 670, 0, 700, 0, 600, 0, 545, 0, 660, 0, 570, 0, 625, 0, 625, 0, 610, 0, 660, 0, 695, 0, 525, 0, 535, 0, 615, 0, 630, 0, 515, 0, 505, 0, 615, 0, 560, 0, 570, 0, 545, 0, 590, 0, 610, 0, 500, 0, 450, 0, 505, 0, 645, 0, 530, 0, 530, 0, 515, 0, 465, 0, 560, 0, 510, 0, 405, 0, 540, 0, 555, 0, 615, 0, 520, 0, 565, 0, 460, 0, 460, 0, 465, 0, 520, 0, 530, 0, 545, 0, 440, 0, 505, 0, 465, 0, 540, 0, 540, 0, 490, 0, 500, 0, 540, 0, 520, 0, 480, 0, 485, 0, 490, 0, 470, 0, 470, 0, 490, 0, 510, 0, 430, 0, 500, 0, 500, 0, 510, 0, 395, 0, 440, 0, 505, 0, 450, 0, 455, 0, 495, 0, 390, 0, 390, 0, 555, 0, 475, 0, 510, 0, 460, 0, 470, 0, 445, 0, 490, 0, 480, 0, 505, 0, 445, 0, 425, 0, 485, 0, 495, 0, 505, 0, 515, 0, 565, 0, 620, 0, 480, 0, 470, 0, 570, 0, 505, 0, 530, 0, 535, 0, 525, 0, 695, 0, 565, 0, 515, 0, 570, 0, 630, 0, 525, 0, 580, 0, 580, 0, 630, 0, 555, 0, 550, 0, 575, 0, 625, 0, 665, 0, 680, 0, 800, 0, 865, 0, 815, 0, 750, 0, 985, 0, 930, 0, 1020, 0, 1285, 0, 1440, 0, 1890, 0, 2460, 0, 4080, 0, 1244780]
    Ps = [0, 0, 200120, 0, 97160, 0, 63740, 0, 46995, 0, 36900, 0, 30770, 0, 25920, 0, 22970, 0, 20455, 0, 18385, 0, 17155, 0, 16335, 0, 14460, 0, 13305, 0, 13175, 0, 11665, 0, 10905, 0, 10985, 0, 10220, 0, 9645, 0, 9125, 0, 8920, 0, 8845, 0, 8685, 0, 8170, 0, 7890, 0, 7285, 0, 7400, 0, 6940, 0, 6970, 0, 6760, 0, 6500, 0, 6220, 0, 6075, 0, 6550, 0, 6105, 0, 5685, 0, 5290, 0, 5690, 0, 5180, 0, 5260, 0, 5400, 0, 5020, 0, 4830, 0, 5060, 0, 4785, 0, 4965, 0, 4845, 0, 4695, 0, 4490, 0, 4390, 0, 4545, 0, 4300, 0, 4800, 0, 4360, 0, 4225, 0, 4295, 0, 3820, 0, 4000, 0, 4060, 0, 3760, 0, 3845, 0, 3790, 0, 3725, 0, 3775, 0, 3810, 0, 3720, 0, 3570, 0, 3900, 0, 3740, 0, 3660, 0, 3425, 0, 3420, 0, 3770, 0, 3455, 0, 3360, 0, 3545, 0, 3295, 0, 3365, 0, 3425, 0, 3670, 0, 3375, 0, 3595, 0, 3540, 0, 3240, 0, 3635, 0, 3130, 0, 3085, 0, 2955, 0, 3215, 0, 3140, 0, 3050, 0, 3020, 0, 2975, 0, 3200, 0, 3125, 0, 3035, 0, 2705, 0, 3000, 0, 2910, 0, 2860, 0, 2850, 0, 2790, 0, 2760, 0, 3005, 0, 2870, 0, 2725, 0, 2735, 0, 2780, 0, 2815, 0, 2780, 0, 2975, 0, 2650, 0, 2700, 0, 2645, 0, 2730, 0, 2845, 0, 3015, 0, 2735, 0, 2795, 0, 2710, 0, 2695, 0, 2725, 0, 2835, 0, 2505, 0, 2645, 0, 2985, 0, 2465, 0, 2850, 0, 2805, 0, 2695, 0, 2535, 0, 2680, 0, 2735, 0, 2655, 0, 2430, 0, 2630, 0, 2695, 0, 2610, 0, 2615, 0, 2900, 0, 2880, 0, 2575, 0, 2670, 0, 2600, 0, 2440, 0, 2755, 0, 2735, 0, 2725, 0, 2810, 0, 2715, 0, 2585, 0, 2820, 0, 2735, 0, 2885, 0, 2920, 0, 2840, 0, 2520, 0, 3010, 0, 2755, 0, 2730, 0, 2730, 0, 2580, 0, 2865, 0, 2680, 0, 2690, 0, 2810, 0, 2880, 0, 3230, 0, 2990, 0, 2775, 0, 2925, 0, 3270, 0, 2940, 0, 2850, 0, 2935, 0, 3190, 0, 2950, 0, 3080, 0, 3130, 0, 3255, 0, 3295, 0, 3255, 0, 3365, 0, 3250, 0, 3485, 0, 3655, 0, 3565, 0, 3655, 0, 3770, 0, 4020, 0, 3925, 0, 4020, 0, 4420, 0, 4580, 0, 4645, 0, 4930, 0, 5385, 0, 6220, 0, 6475, 0, 7680, 0, 9290, 0, 11365, 0, 18225, 0, 2389605]
    nallele = 410
    fit(Pn,Ps,nallele,lcut=0.1,hcut=0.9,debug=True,plot=True)
