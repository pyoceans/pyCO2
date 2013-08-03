#!/usr/bin/env python
#
#
# alkalinity.py
#
# purpose:  Compute Alkalinity for CO2SYS
# author:   Filipe P. A. Fernandes
# e-mail:   ocefpaf@gmail
# web:      http://ocefpaf.tiddlyspot.com/
# created:  28-Jan-2012
# modified: Sat 28 Jan 2012 07:37:12 PM EST
#
# obs: Enter a vector of pH, Acid Vol, temp, and cond
#      Enter a scalar of sal
#      Outputs a scalar for alkalinity


K = 273

def conts_G(vol, pot, t):
    r"""
    Constant "G"
    Parameters
    ----------
    vol : array_like
        Volume of acid added to the sample [ml]
    pot : array_like
        Potential measured during the titration [mV]
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]
    """

def at(S=35, T=25, C=0.1, d=1, pHTris=NULL, ETris=NULL, weight, E, volume):
    r"""
    S salinity    (constant)
    T temperature (degC)  (vector or constant)
    C molarity of acid   (constant)
    d density of acid    (constant)
    pHTris, pH used for the calibration of the electrode  (constant)
    ETris, voltage used for the calibration of the electrode (constant)
    E  voltage recorded during the titration in milliVolt(vector)
    volume : volume of acid (vector)
    weight : mass of the sample  (constant)
    """

    # Constants.
    R, F = -8.31447215, -96485.339924

    Tk = T + 273.15

    # Creation of a table p.
    p = data.frame(E=E, volume=volume, Tk=Tk)
    z = p

    # Transform mV in pH (total scale).
    if (!is.null(pHTris) & !is.null(ETris)){
        pH <- pHTris + (ETris/1000-E/1000)/(R*Tk*log(10)/F)

    # creation of a table p
    p <- data.frame(p, pH=pH)
    #to use only the value between pH 3.5 and 3 (according to Dickson, 2007):
    iii <- which((3<= p$pH)&(p$pH<=3.5))
    z<- p[iii,]
    }

    options(digits=9)
    m <- z$volume*d # Mass of acid
    m0 <- weight    # Mass of the sample

    #linear estimation of the total alkalinity (gran function):
    F1<-(m0+m)*exp((z$E/1000)/(R*(z$Tk)/F))
    f<-lm(m~F1)
    TA<-f$coefficients[1]*C/m0[1]

    #non linear estimation:
    E0 <- z$E/1000-(R*z$Tk/F)*log((-m0*TA+m*C)/(m0+m))
    Hprime <- exp((z$E/1000-E0)/(R*z$Tk/F))

    St <- (0.14/96.062)*(S/1.80655)
    Ks <- Ks(S,T[1], 0)
    Z <- 1+ St/Ks
    Ft <- (0.000067/18.998)*(S/1.80655)
    Kf <- exp(874/z$Tk-9.68+0.111*S**(1/2))
    y <- (m/m0)
    regr <- nls(y ~ ((At + (St/(1 + Ks*Z/(f*Hprime)))+(Ft/ (1+Kf/(f*Hprime)))+(f*Hprime/Z))/(C-f*Hprime/Z)),start=list(At=TA, f=1))

    ALK <- summary(regr)$parameters[1]
    attr(ALK,"unit") <- "mol/kg-soln"
    attr(ALK,"name") <- "Total Alkalinity"

    # At= Total Alkalinity (mol/kg)
    ALK = (100 * vol) * (10 ** ((pot + 96.4853 * 10 ** 3) / (np.log(10) * 8.31451 * 1000 * (t + K)))) * 0.001
    return ALK

# Using UNESCO 1983 (EOS 1980) polynomial.
def d_smow(t):
    r"""
    Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.

    Parameters
    ----------
    t : array_like
        in situ temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    Dsmow(t) : array_like
        density [kg m^-3]

    Notes
    -----
    Standard Mean Ocean Water (SMOW) is the water collected in the deep ocean
    used as a reference.

    Examples
    --------
    Data from UNESCO Tech. Paper in Marine Sci. No. 44, p22

    >>> import seawater.csiro as sw
    >>> t = T90conv([0, 0, 30, 30, 0, 0, 30, 30])
    >>> sw.smow(t)
    array([ 999.842594  ,  999.842594  ,  995.65113374,  995.65113374,
            999.842594  ,  999.842594  ,  995.65113374,  995.65113374])

    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
    of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
    doi:10.1016/0198-0149(81)90122-9
    """

    t = np.asanyarray(t)

    # Constants
    a0 = 999.842594
    a1 = 6.793952e-2
    a2 = -9.095290e-3
    a3 = 1.001685e-4
    a4 = -1.120083e-6
    a5 = 6.536332e-9

    T68 = t * 1.00024
    Dsmow = a0 + (a1 + (a2 + (a3 + (a4 + a5 * T68) * T68) * T68) * T68) * T68

    return Dsmow


def dens0(s, t):
    r"""Density of Sea Water at atmospheric pressure.

    Parameters
    ----------
    s(p=0) : array_like
        salinity [psu (PSS-78)]
    t(p=0) : array_like
        temperature [:math:`^\circ` C (ITS-90)]

    Returns
    -------
    dens0(s, t) : array_like
        density  [kg m :sup:`3`] of salt water with properties
        (s, t, p=0) 0 db gauge pressure


    References
    ----------
    .. [1] Fofonoff, P. and Millard, R.C. Jr UNESCO 1983. Algorithms for
    computation of fundamental properties of seawater. UNESCO Tech. Pap. in
    Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39.
    http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    .. [2] Millero, F.J. and  Poisson, A. International one-atmosphere equation
    of state of seawater. Deep-Sea Res. 1981. Vol28A(6) pp625-629.
    doi:10.1016/0198-0149(81)90122-9
    """

    s, t = np.asanyarray(s), np.asanyarray(t)

    T68 = T * 1.00024

    # UNESCO 1983 eqn(13) p17.
    b0 = 8.24493e-1
    b1 = -4.0899e-3
    b2 = 7.6438e-5
    b3 = -8.2467e-7
    b4 = 5.3875e-9

    A = (b0 + (b1 + (b2 + ( b3 + b4 * T68) * T68) * T68) * T68)

    c0 = -5.72466e-3
    c1 = 1.0227e-4
    c2 = -1.6546e-6

    B = (c0 + (c1 + c2 * T68) * T68)

    C = 4.8314e-4

    dens0 = d_smow(t) + A * s + B * s ** 1.5 + C * s ** 2

    return dens0 / 1000.