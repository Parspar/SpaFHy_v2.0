# -*- coding: utf-8 -*-
"""
Muunnosfunktiot koordinaattiprojektioille
ETRS89-TM35FIN, geodeettisista tasokoordinaateiksi ja takaisin

Lähde: JHS 154, 6.6.2008
http://www.jhs-suositukset.fi/web/guest/jhs/recommendations/154

2013-04-29/JeH, loukko (at) loukko (dot) net
http://www.loukko.net/koord_proj/
Vapaasti käytettävissä ilman toimintatakuuta.

25.1.2019 khaahti to python
13.5.2019 khaahti: muokattu muuttamaan KKJ-grid, zone 3/Uniform (YKJ) -> KKJ geografical
"""

# Muuntaa desimaalimuotoiset leveys- ja pituusasteet YKJ tasokoordinaateiksi
# koordGT(60.565894, 24.822422) -->  (6719258, 3380581)

import numpy as np

def koordGT(lev_aste, pit_aste, desimals=0):

    # Vakiot
    f = 1 / 297.0  # Ellipsoidin litistyssuhde
    a = 6378388  # Isoakselin puolikas
    lmbda_nolla = 0.471238898  # Keskimeridiaani (rad), 27 astetta
    k_nolla = 1.0  # Mittakaavakerroin
    E_nolla = 3500000  # Itäkoordinaatti

    # Kaavat
    # Muunnetaan astemuotoisesta radiaaneiksi
    fii = np.pi / 180.0  * lev_aste
    lmbda = np.pi / 180.0  * pit_aste

    n = f / (2-f)
    A1 = (a/(1+n)) * (1 + (pow(n, 2)/4) + (pow(n, 4)/64))
    e_toiseen = (2 * f) - pow(f, 2)
    e_pilkku_toiseen = e_toiseen / (1 - e_toiseen)
    h1_pilkku = (1/2)*n - (2/3)*pow(n, 2) + (5/16)*pow(n, 3) + (41/180)*pow(n, 4)
    h2_pilkku = (13/48)*pow(n, 2) - (3/5)*pow(n, 3) + (557/1440)*pow(n, 4)
    h3_pilkku =(61/240)*pow(n, 3) - (103/140)*pow(n, 4)
    h4_pilkku = (49561/161280)*pow(n, 4)
    Q_pilkku = np.arcsinh( np.tan(fii))
    Q_2pilkku = np.arctanh(np.sqrt(e_toiseen) * np.sin(fii))
    Q = Q_pilkku - np.sqrt(e_toiseen) * Q_2pilkku
    l = lmbda - lmbda_nolla
    beeta = np.arctan(np.sinh(Q))
    eeta_pilkku = np.arctanh(np.cos(beeta) * np.sin(l))
    zeeta_pilkku = np.arcsin(np.sin(beeta)/(1/np.cosh(eeta_pilkku)))
    zeeta1 = h1_pilkku * np.sin( 2 * zeeta_pilkku) * np.cosh( 2 * eeta_pilkku)
    zeeta2 = h2_pilkku * np.sin( 4 * zeeta_pilkku) * np.cosh( 4 * eeta_pilkku)
    zeeta3 = h3_pilkku * np.sin( 6 * zeeta_pilkku) * np.cosh( 6 * eeta_pilkku)
    zeeta4 = h4_pilkku * np.sin( 8 * zeeta_pilkku) * np.cosh( 8 * eeta_pilkku)
    eeta1 = h1_pilkku * np.cos( 2 * zeeta_pilkku) * np.sinh( 2 * eeta_pilkku)
    eeta2 = h2_pilkku * np.cos( 4 * zeeta_pilkku) * np.sinh( 4 * eeta_pilkku)
    eeta3 = h3_pilkku * np.cos( 6 * zeeta_pilkku) * np.sinh( 6 * eeta_pilkku)
    eeta4 = h4_pilkku * np.cos( 8 * zeeta_pilkku) * np.sinh( 8 * eeta_pilkku)
    zeeta = zeeta_pilkku + zeeta1 + zeeta2 + zeeta3 + zeeta4
    eeta = eeta_pilkku + eeta1 + eeta2 + eeta3 + eeta4

    # Tulos tasokoordinaatteina
    N = A1 * zeeta * k_nolla
    E = A1 * eeta * k_nolla + E_nolla

    return np.round(N, desimals), np.round(E, desimals)


# koordTG
# Muuntaa YKJ tasokoordinaatit desimaalimuotoisiksi leveys- ja pituusasteiksi
# koordTG(6719258, 3380581) -> (60.565894, 24.822422)

def  koordTG(N, E, desimals=2):

    # Vakiot
    f = 1 / 297.0  # Ellipsoidin litistyssuhde
    a = 6378388  # Isoakselin puolikas
    lmbda_nolla = 0.471238898  # Keskimeridiaani (rad), 27 astetta
    k_nolla = 1.0  # Mittakaavakerroin
    E_nolla = 3500000  # Itäkoordinaatti

    # Kaavat
    n = f / (2-f)
    A1 = (a/(1+n)) * (1 + (pow(n, 2)/4) + (pow(n, 4)/64))
    e_toiseen = (2 * f) - pow(f, 2)
    h1 = (1/2)*n - (2/3)*pow(n, 2) + (37/96)*pow(n, 3) - (1/360)*pow(n, 4)
    h2 = (1/48)*pow(n, 2) + (1/15)*pow(n, 3) - (437/1440)*pow(n, 4)
    h3 =(17/480)*pow(n, 3) - (37/840)*pow(n, 4)
    h4 = (4397/161280)*pow(n, 4)
    zeeta = N / (A1 * k_nolla)
    eeta = (E - E_nolla) / (A1 * k_nolla)
    zeeta1_pilkku = h1 * np.sin( 2 * zeeta) * np.cosh( 2 * eeta)
    zeeta2_pilkku = h2 * np.sin( 4 * zeeta) * np.cosh( 4 * eeta)
    zeeta3_pilkku = h3 * np.sin( 6 * zeeta) * np.cosh( 6 * eeta)
    zeeta4_pilkku = h4 * np.sin( 8 * zeeta) * np.cosh( 8 * eeta)
    eeta1_pilkku = h1 * np.cos( 2 * zeeta) * np.sinh( 2 * eeta)
    eeta2_pilkku = h2 * np.cos( 4 * zeeta) * np.sinh( 4 * eeta)
    eeta3_pilkku = h3 * np.cos( 6 * zeeta) * np.sinh( 6 * eeta)
    eeta4_pilkku = h4 * np.cos( 8 * zeeta) * np.sinh( 8 * eeta)
    zeeta_pilkku = zeeta - (zeeta1_pilkku + zeeta2_pilkku + zeeta3_pilkku + zeeta4_pilkku)
    eeta_pilkku = eeta - (eeta1_pilkku + eeta2_pilkku + eeta3_pilkku + eeta4_pilkku)
    beeta = np.arcsin((1/np.cosh(eeta_pilkku)*np.sin(zeeta_pilkku)))
    l = np.arcsin(np.tanh(eeta_pilkku)/(np.cos(beeta)))
    Q = np.arcsinh(np.tan(beeta))
    Q_pilkku = Q + np.sqrt(e_toiseen) * np.arctanh(np.sqrt(e_toiseen) * np.tanh(Q))

    for kierros in range(2):
        Q_pilkku = Q + np.sqrt(e_toiseen) * np.arctanh(np.sqrt(e_toiseen) * np.tanh(Q_pilkku))

    # Tulos radiaaneina
    fii = np.arctan(np.sinh(Q_pilkku))
    lmbda = lmbda_nolla + l

    # Tulos asteina
    fii = fii / np.pi * 180.0
    lmbda = lmbda / np.pi * 180.0

    return np.round(fii, desimals), np.round(lmbda, desimals)